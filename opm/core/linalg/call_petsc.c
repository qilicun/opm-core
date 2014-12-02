/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.
  Copyright 2014 STATOIL ASA.

  This file is part of the Open Porous Media Project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "config.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <opm/core/linalg/call_petsc.h>

static PetscErrorCode code;

typedef struct {
    Vec     b;      /*b = rhs*/
    Vec     u;      /*u = solution*/
    Mat     A;      /*A = matrix*/
    KSP     ksp;    /*ksp = solver*/
    PC      pc;     /*pc = preconditioner*/
} OEM_DATA;

typedef struct {
    KSPType     method;  /*ksp method*/
    PCType      pcname;  /*pc method*/
    int         view_ksp;    /*whether view ksp detail information*/
    double      rtol;
    double      atol;
    double      dtol;
    int         maxits;
} KSP_OPT;

static int
init(OEM_DATA* t, KSP_OPT* opts)
{
    t->b = PETSC_NULL;
    t->u = PETSC_NULL;
    t->A = PETSC_NULL;
    t->ksp = PETSC_NULL;
    t->pc = PETSC_NULL;

    /*set default options for ksp solvers*/
    opts->method = KSPGMRES;
    opts->pcname = PCJACOBI;
    opts->view_ksp = 0;
    opts->rtol = PETSC_DEFAULT;
    opts->atol = PETSC_DEFAULT;
    opts->dtol = PETSC_DEFAULT;
    opts->maxits = PETSC_DEFAULT;

    return 0;
}

static int
create_mat_vec(const int size, OEM_DATA* t)
{
    code = VecCreate(PETSC_COMM_WORLD,&t->u);CHKERRQ(code);
    code = VecSetSizes(t->u,PETSC_DECIDE, size);CHKERRQ(code);
    code = VecSetFromOptions(t->u);CHKERRQ(code);
    code = VecDuplicate(t->u, &t->b);CHKERRQ(code);

    code = MatCreate(PETSC_COMM_WORLD, &t->A);CHKERRQ(code);
    code = MatSetSizes(t->A,PETSC_DECIDE,PETSC_DECIDE, size, size);CHKERRQ(code);
    code = MatSetFromOptions(t->A);CHKERRQ(code);
    code = MatSetUp(t->A);CHKERRQ(code);
    
    return 0;
}

static int
to_petsc_vec(const double *x, Vec v)
{
    PetscScalar *vec;
    PetscInt size;
    int i;

    if (v == PETSC_NULL) {
        PetscPrintf(PETSC_COMM_WORLD,"PETSc CopySolution: invalid PETSc vector.\n");
        CHKERRQ(code);
    }
    code = VecGetLocalSize(v, &size); CHKERRQ(code);
    code = VecGetArray(v, &vec); CHKERRQ(code);
    for (i = 0; i < size; ++i) {
        vec[i] = x[i];
    }
    code = VecRestoreArray(v, &vec); CHKERRQ(code);

    return 0;
}

static int
from_petsc_vec(double *x, Vec v)
{
    PetscScalar *vec;
    PetscInt size;
    int i;
    
    if (v == PETSC_NULL) {
        PetscPrintf(PETSC_COMM_WORLD,"PETSc CopySolution: invalid PETSc vector.\n");
        CHKERRQ(code);
    }
    code = VecGetLocalSize(v, &size); CHKERRQ(code);
    code = VecGetArray(v, &vec); CHKERRQ(code);
    for (i = 0; i < size; ++i) {
        x[i] = vec[i];
    }
    code = VecRestoreArray(v, &vec); CHKERRQ(code);

    return 0;
}

static int
to_petsc_mat(const int size, const int nonzeros, const int* ia, const int* ja, const double* sa, Mat A)
{
    for (int i = 0; i < size; ++i) {
        int nrows = ia[i+1] - ia[i];
        if (nrows == 0) {
            continue;
        }
        for (int j = ia[i]; j < ia[i+1]; ++j) {
            code = MatSetValues(A,1,&i,1,&ja[j],&sa[j],INSERT_VALUES);
            CHKERRQ(code);
        }
    } 
    code = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(code);
    code = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(code);

    return 0;
}

static int
solve_system(OEM_DATA* t, KSP_OPT* opts, PetscInt* its, PetscReal* res, KSPConvergedReason* reason)
{
    PetscInt it;
    PetscReal residual;
    code = KSPCreate(PETSC_COMM_WORLD,&t->ksp);CHKERRQ(code);
    code = KSPSetOperators(t->ksp,t->A,t->A,DIFFERENT_NONZERO_PATTERN);CHKERRQ(code);
    code = KSPGetPC(t->ksp,&t->pc);CHKERRQ(code);
    code = KSPSetType(t->ksp, opts->method);
    code = PCSetType(t->pc, opts->pcname);CHKERRQ(code);
    code = KSPSetTolerances(t->ksp,opts->rtol,opts->atol,opts->dtol,opts->maxits);CHKERRQ(code);
    code = KSPSetFromOptions(t->ksp);CHKERRQ(code);
    code = KSPSetInitialGuessNonzero(t->ksp,PETSC_TRUE);CHKERRQ(code);
    code = KSPSolve(t->ksp,t->b,t->u);CHKERRQ(code);
    code = KSPGetConvergedReason(t->ksp, reason); CHKERRQ(code);
    code = KSPGetIterationNumber(t->ksp, &it); CHKERRQ(code);
    code = KSPGetResidualNorm(t->ksp, &residual); CHKERRQ(code);
    if (opts->view_ksp) {
        code = KSPView(t->ksp,PETSC_VIEWER_STDOUT_WORLD);
        CHKERRQ(code);
    }

    code = PetscPrintf(PETSC_COMM_WORLD,"KSP Iterations %D, Final Residual %G\n",it,residual);CHKERRQ(code);
    *its = it;
    *res = residual;

    return 0;
}

static int
destory(OEM_DATA* t, KSP_OPT* opts)
{
    if(t == NULL && opts == NULL) {
       return 0;
    }
    if(t->u != PETSC_NULL) {
        code = VecDestroy(&t->u);
        CHKERRQ(code); 
    }
    if(t->b != PETSC_NULL) {
        code = VecDestroy(&t->b);
        CHKERRQ(code);
    }
    if(t->A != PETSC_NULL) {
        code = MatDestroy(&t->A);
        CHKERRQ(code); 
    }
    if(t->ksp != PETSC_NULL) {
        code = KSPDestroy(&t->ksp);
        CHKERRQ(code);
    }

    free(t);
    free(opts);

    return 0; 
}

static int 
set_ksp_opts(const KSPType ksp_type, const PCType pc_type, const double rtol, const double atol, const double dtol, const int maxits, const int view_ksp, KSP_OPT* opts)
{
    opts->method = ksp_type;
    opts->pcname = pc_type;
    opts->view_ksp = view_ksp;
    opts->rtol = rtol;
    opts->atol = atol;
    opts->dtol = dtol;
    opts->maxits = maxits;

    return 0;
}

KSPConvergedReason
call_Petsc(const int size, const int nonzeros, const int* ia, const int* ja, const double* sa, const double* b, double* x, const KSPType ksp_type, const PCType pc_type, const double rtol, const double atol, const double dtol, const int maxits, const int view_ksp, int* it, double* res)
{
    OEM_DATA* t;
    KSP_OPT* opts;
    KSPConvergedReason reason;
    t = malloc(sizeof(OEM_DATA));
    opts = malloc(sizeof(KSP_OPT));
    assert(t);
    assert(opts);
    init(t, opts);
    create_mat_vec(size, t);
    to_petsc_mat(size, nonzeros, ia, ja, sa, t->A);
    to_petsc_vec(b, t->b);
    set_ksp_opts(ksp_type, pc_type, rtol, atol, dtol, maxits, view_ksp, opts);
    solve_system(t, opts, it, res, &reason);  
    from_petsc_vec(x, t->u);
    destory(t, opts);

    return reason;
}

