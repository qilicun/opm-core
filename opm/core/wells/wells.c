/*===========================================================================
//
// File: newwells.c
//
// Created: 2012-02-03 11:28:40+0100
//
// Authors: Knut-Andreas Lie      <Knut-Andreas.Lie@sintef.no>
//          Jostein R. Natvig     <Jostein.R.Natvig@sintef.no>
//          Halvor M. Nilsen      <HalvorMoll.Nilsen@sintef.no>
//          Atgeirr F. Rasmussen  <atgeirr@sintef.no>
//          Xavier Raynaud        <Xavier.Raynaud@sintef.no>
//          Bård Skaflestad       <Bard.Skaflestad@sintef.no>
//
//==========================================================================*/


/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.
  Copyright 2012 Statoil ASA.

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
#include <opm/core/wells.h>
#include <opm/core/well_controls.h>

#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include <stdio.h>


struct WellMgmt {
    int well_cpty;
    int perf_cpty;
};


static void
destroy_well_mgmt(struct WellMgmt *m)
{
    free(m);
}


static struct WellMgmt *
create_well_mgmt(void)
{
    struct WellMgmt *m;

    m = malloc(1 * sizeof *m);

    if (m != NULL) {
        m->well_cpty = 0;
        m->perf_cpty = 0;
    }

    return m;
}




/* ---------------------------------------------------------------------- */
static int
wells_allocate(int nwells, struct Wells *W)
/* ---------------------------------------------------------------------- */
{
    int   ok, np;
    void *type, *depth_ref, *comp_frac;
    void *well_connpos;
    void *ctrls, *name;

    np = W->number_of_phases;

    type      = realloc(W->type     ,  1 * nwells * sizeof *W->type);
    depth_ref = realloc(W->depth_ref,  1 * nwells * sizeof *W->depth_ref);
    comp_frac = realloc(W->comp_frac, np * nwells * sizeof *W->comp_frac);
    ctrls     = realloc(W->ctrls    ,  1 * nwells * sizeof *W->ctrls);
    name      = realloc(W->name     ,  1 * nwells * sizeof *W->name);

    well_connpos = realloc(W->well_connpos,
                           (nwells + 1) * sizeof *W->well_connpos);

    ok = 0;
    if (type         != NULL) { W->type         = type        ; ok++; }
    if (depth_ref    != NULL) { W->depth_ref    = depth_ref   ; ok++; }
    if (comp_frac    != NULL) { W->comp_frac    = comp_frac   ; ok++; }
    if (well_connpos != NULL) { W->well_connpos = well_connpos; ok++; }
    if (ctrls        != NULL) { W->ctrls        = ctrls       ; ok++; }
    if (name         != NULL) { W->name         = name        ; ok++; }

    return ok == 6;
}


/* ---------------------------------------------------------------------- */
static int
perfs_allocate(int nperf, struct Wells *W)
/* ---------------------------------------------------------------------- */
{
    int   ok;
    void *well_cells, *WI;

    well_cells = realloc(W->well_cells, nperf * sizeof *W->well_cells);
    WI         = realloc(W->WI        , nperf * sizeof *W->WI        );

    ok = 0;
    if (well_cells != NULL) { W->well_cells = well_cells; ok++; }
    if (WI         != NULL) { W->WI         = WI        ; ok++; }

    return ok == 2;
}


/* ---------------------------------------------------------------------- */
static int
initialise_new_wells(int nwells, struct Wells *W)
/* ---------------------------------------------------------------------- */
{
    int ok, w, p;

    struct WellMgmt *m;

    m = W->data;

    for (w = m->well_cpty; w < nwells; w++) {
        W->type     [w]        = PRODUCER;
        W->depth_ref[w]        = -1.0;
        W->name     [w]        = NULL;

        for (p = 0; p < W->number_of_phases; ++p) {
            W->comp_frac[W->number_of_phases*w + p] = 0.0;
        }

        W->well_connpos[w + 1] = W->well_connpos[w];
    }

    for (w = m->well_cpty, ok = 1; ok && (w < nwells); w++) {
        W->ctrls[w] = well_controls_create( );

        ok = W->ctrls[w] != NULL;
    }

    if (! ok) {
        for (; w < nwells; w++) {
            W->ctrls[w] = NULL;
        }
    }

    return ok;
}


/* ---------------------------------------------------------------------- */
static void
initialise_new_perfs(int nperf, struct Wells *W)
/* ---------------------------------------------------------------------- */
{
    int k;

    struct WellMgmt *m;

    m = W->data;

    for (k = m->perf_cpty; k < nperf; k++) {
        W->well_cells[k] = -1 ;
        W->WI        [k] = 0.0;
    }
}


/* ---------------------------------------------------------------------- */
static int
wells_reserve(int nwells, int nperf, struct Wells *W)
/* ---------------------------------------------------------------------- */
{
    int ok;

    struct WellMgmt *m;

    m = W->data;

    assert (nwells >= m->well_cpty);
    assert (nperf  >= m->perf_cpty);

    ok = 1;

    if (nwells > m->well_cpty) {
        ok = wells_allocate(nwells, W);

        if (ok) {
            ok = initialise_new_wells(nwells, W);
        }

        if (ok) {
            m->well_cpty = nwells;
        }
    }

    if (ok && (nperf > m->perf_cpty)) {
        ok = perfs_allocate(nperf, W);

        if (ok) {
            initialise_new_perfs(nperf, W);
            m->perf_cpty = nperf;
        }
    }

    return ok;
}


/* ---------------------------------------------------------------------- */
static char *
dup_string(const char *s)
/* ---------------------------------------------------------------------- */
{
    char *t;

    assert (s != NULL);

    t = malloc((strlen(s) + 1) * sizeof *t);

    if (t != NULL) {
        strcpy(t, s);
    }

    return t;
}


/* ======================================================================
 * Public entry points below separator.
 * ====================================================================== */


/* ---------------------------------------------------------------------- */
struct Wells *
create_wells(int nphases, int nwells, int nperf)
/* ---------------------------------------------------------------------- */
{
    int           ok;
    struct Wells *W;

    W = malloc(1 * sizeof *W);

    if (W != NULL) {
        W->number_of_wells = 0;
        W->number_of_phases = nphases;

        W->type            = NULL;
        W->depth_ref       = NULL;
        W->comp_frac       = NULL;

        W->well_connpos    = malloc(1 * sizeof *W->well_connpos);
        W->well_cells      = NULL;
        W->WI              = NULL;

        W->ctrls           = NULL;
        W->name            = NULL;

        W->data            = create_well_mgmt();

        ok = (W->well_connpos != NULL) && (W->data != NULL);
        if (ok) {
            W->well_connpos[0] = 0;

            if ((nwells > 0) || (nperf > 0)) {
                ok = wells_reserve(nwells, nperf, W);
            }
        }

        if (! ok) {
            destroy_wells(W);
            W = NULL;
        }
    }

    return W;
}


/* ---------------------------------------------------------------------- */
void
destroy_wells(struct Wells *W)
/* ---------------------------------------------------------------------- */
{
    int w;

    struct WellMgmt *m;

    if (W != NULL) {
        m = W->data;

        for (w = 0; w < m->well_cpty; w++) {
            well_controls_destroy(W->ctrls[w]);
        }

        for (w = 0; w < m->well_cpty; w++) {
            free(W->name[w]);
        }

        destroy_well_mgmt(m);

        free(W->name);
        free(W->ctrls);
        free(W->WI);
        free(W->well_cells);
        free(W->well_connpos);
        free(W->comp_frac);
        free(W->depth_ref);
        free(W->type);
    }

    free(W);
}


/* ---------------------------------------------------------------------- */
static int
alloc_size(int n, int a, int cpty)
/* ---------------------------------------------------------------------- */
{
    if (cpty < n + a) {
        cpty *= 2;              /* log_2(n) allocations */

        if (cpty < n + a) {     /* Typically for the first few allocs */
            cpty = n + a;
        }
    }

    return cpty;
}


/* ---------------------------------------------------------------------- */
int
add_well(enum WellType  type     ,
         double         depth_ref,
         int            nperf    ,
         const double  *comp_frac, /* Injection fraction or NULL */
         const int     *cells    ,
         const double  *WI       , /* Well index per perf (or NULL) */
         const char    *name     , /* Well name (or NULL) */
         struct Wells  *W        )
/* ---------------------------------------------------------------------- */
{
    int ok, nw, np, nperf_tot, off;
    int nwalloc, nperfalloc;

    struct WellMgmt *m;

    assert (W != NULL);

    nw        = W->number_of_wells;
    nperf_tot = W->well_connpos[nw];

    m = W->data;

    ok = (nw < m->well_cpty) && (nperf_tot + nperf <= m->perf_cpty);

    if (! ok) {
        nwalloc    = alloc_size(nw       , 1    , m->well_cpty);
        nperfalloc = alloc_size(nperf_tot, nperf, m->perf_cpty);

        ok = wells_reserve(nwalloc, nperfalloc, W);
    }

    off = W->well_connpos[nw];

    if (ok && (nperf > 0)) {
        assert (cells != NULL);

        memcpy(W->well_cells + off,
               cells, nperf * sizeof *W->well_cells);

        if (WI != NULL) {
            memcpy(W->WI + off, WI, nperf * sizeof *W->WI);
        }
    }

    if (ok) {
        W->type     [nw] = type     ;
        W->depth_ref[nw] = depth_ref;

        if (name != NULL) {
             /* May return NULL, but that's fine for the current
              * purpose. */
            W->name [nw] = dup_string(name);
        }

        np = W->number_of_phases;
        if (comp_frac != NULL) {
            memcpy(W->comp_frac + np*nw, comp_frac, np * sizeof *W->comp_frac);
        }

        W->well_connpos[nw + 1]  = off + nperf;
        W->number_of_wells      += 1;
    }

    return ok;
}


/* ---------------------------------------------------------------------- */
int
append_well_controls(enum WellControlType type,
                     double               target,
                     const double        *distr,
                     int                  well_index,
                     struct Wells        *W)
/* ---------------------------------------------------------------------- */
{
    struct WellControls    *ctrl;

    assert (W != NULL);
    assert ((0 <= well_index) && (well_index < W->number_of_wells));

    ctrl = W->ctrls[well_index];
    assert (ctrl != NULL);
    
    well_controls_assert_number_of_phases( ctrl , W->number_of_phases);
    return well_controls_add_new(type , target , distr , ctrl);
}


/* ---------------------------------------------------------------------- */
void
set_current_control(int well_index, int current_control, struct Wells *W)
/* ---------------------------------------------------------------------- */
{
    assert (W != NULL);
    assert ((0 <= well_index) && (well_index < W->number_of_wells));
    assert (W->ctrls[well_index] != NULL);
    assert (current_control < well_controls_get_num(W->ctrls[well_index]));
    
    well_controls_set_current(W->ctrls[well_index] , current_control);
}


/* ---------------------------------------------------------------------- */
void
clear_well_controls(int well_index, struct Wells *W)
/* ---------------------------------------------------------------------- */
{
    assert (W != NULL);
    assert ((0 <= well_index) && (well_index < W->number_of_wells));

    if (W->ctrls[well_index] != NULL) 
        well_controls_clear( W->ctrls[well_index] );
}




/* ---------------------------------------------------------------------- */
struct Wells *
clone_wells(const struct Wells *W)
/* ---------------------------------------------------------------------- */
{
    int                  np, nperf, ok, pos, w;
    const int           *cells;
    const double        *WI, *comp_frac;

    struct WellControls *ctrl;
    struct Wells        *newWells;

    if (W == NULL) {
        newWells = NULL;
    }
    else {
        np  = W->number_of_phases;
        newWells = create_wells(W->number_of_phases, W->number_of_wells,
                                W->well_connpos[ W->number_of_wells ]);

        if (newWells != NULL) {
            pos = W->well_connpos[ 0 ];
            ok  = 1;

            for (w = 0; ok && (w < W->number_of_wells); w++) {
                nperf = W->well_connpos[w + 1] - pos;
                cells = W->well_cells + pos;

                WI        = W->WI        != NULL ? W->WI        + pos  : NULL;
                comp_frac = W->comp_frac != NULL ? W->comp_frac + w*np : NULL;

                ok = add_well(W->type[ w ], W->depth_ref[ w ], nperf,
                              comp_frac, cells, WI, W->name[ w ], newWells);

                if (ok) {
                    ok = (ctrl = well_controls_clone(W->ctrls[w])) != NULL;
                }

                if (ok) {
                    /* Destroy control set implied by add_well() */
                    well_controls_destroy(newWells->ctrls[w]);

                    /* Assign complete clone of w's control set */
                    newWells->ctrls[w] = ctrl;
                }

                pos = W->well_connpos[w + 1];
            }

            if (! ok) {
                destroy_wells(newWells);
                newWells = NULL;
            }
        }
    }

    assert (wells_equal(newWells, W, false));

    return newWells;
}

/* ---------------------------------------------------------------------- */
bool
wells_equal(const struct Wells *W1, const struct Wells *W2 , bool verbose)
/* ---------------------------------------------------------------------- */
{
    bool are_equal = true;
    are_equal = (W1->number_of_wells == W2->number_of_wells);
    are_equal = are_equal && (W1->number_of_phases == W2->number_of_phases);
    if (!are_equal) {
        return are_equal;
    }

    for (int i=0; i<W1->number_of_wells; i++) {
        if (are_equal) {
            /*
              The name attribute can be NULL. The comparison is as
              follows:
              
                1. If both names are different from NULL a normal
                   strcmp() is performed.
                2. If both names are NULL they compare as equal.
                3. If one name is NULL and the other is not NULL
                   they are regarded as different.
            */
            if (W1->name[i] && W2->name[i])
                are_equal = are_equal && (strcmp(W1->name[i], W2->name[i]) == 0);
            else 
                are_equal = are_equal && (W1->name[i] == W2->name[i]);  

            if (verbose && !are_equal) 
                printf("Well name[%d] %s and %s are different \n",  i , W1->name[i] ,  W2->name[i]);
        }
        if (W1->type[i] != W2->type[i]) {
            are_equal = false;
            if (verbose)
                printf("Well->type[%d] different %d %d \n",i , W1->type[i] , W2->type[i] );
        }
        if (W1->depth_ref[i] != W2->depth_ref[i]) {
            are_equal = false;
            if (verbose)
                printf("Well->depth_ref[%d] different %g %g \n",i , W1->depth_ref[i] , W2->depth_ref[i] );
        }
        if (!well_controls_equal(W1->ctrls[i], W2->ctrls[i],verbose)) {
            are_equal = false;
            if (verbose)
                printf("Well controls are different for well[%d]:%s \n",i,W1->name[i]);
        }
    }


    {
        struct WellMgmt* mgmt1 = W1->data;
        struct WellMgmt* mgmt2 = W2->data;
        are_equal = are_equal && (mgmt1->perf_cpty == mgmt2->perf_cpty);
        are_equal = are_equal && (mgmt1->well_cpty == mgmt2->well_cpty);
    }

    if (memcmp(W1->comp_frac, W2->comp_frac, W1->number_of_wells * W1->number_of_phases * sizeof *W1->comp_frac ) != 0) {
        are_equal = false;
        if (verbose)
            printf("Component fractions different \n");
    }
    
    if (memcmp(W1->well_connpos, W2->well_connpos, (1 + W1->number_of_wells) * sizeof *W1->well_connpos ) != 0) {
        are_equal = false;
        if (verbose)
            printf("perforation position map difference \n");
    }

    {
        int number_of_perforations = W1->well_connpos[W1->number_of_wells];

        are_equal = are_equal && (memcmp(W1->well_cells, W2->well_cells, number_of_perforations * sizeof *W1->well_cells ) == 0);
        are_equal = are_equal && (memcmp(W1->WI, W2->WI, number_of_perforations * sizeof *W1->WI ) == 0);
    }

    return are_equal;
}

/* ---------------------------------------------------------------------- */
