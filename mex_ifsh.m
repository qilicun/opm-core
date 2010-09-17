function varargout = mex_ifsh(varargin)
%Discretise and solve flow equation using compiled C code.
%
% SYNOPSIS:
%   state = mex_ifsh(state, G, rock, bc, src)
%
% DESCRIPTION:
%   Equivalent to the standard MRST call sequence
%
%   fluid = initSingleFluid('mu', 1, 'rho', 1)
%   S     = computeMimeticIP(G, rock)
%   state = solveIncompFlow(state, G, S, fluid, 'src', src, 'bc', bc)
%
%   Note in particular that the inner products are computed at each call.
%
% PARAMETERS:
%   state   - Reservori state.
%
%   G, rock - Grid and rock data structures, respectivelly
%
%   bc, src - Boundary condition and source data structures as defined by
%             'addBC' and 'addSource', respectively.  Either may be empty.
%
% RETURNS:
%   state   - Updated reservoir state.  Contains new values for
%             'state.pressure' and 'state.flux'.
%
% SEE ALSO:
%   mex_ip_simple, mex_schur_comp_symm, test_mex_schur_comp_symm.

%{
#COPYRIGHT#
%}

% $Date$
% $Revision$

   buildmex CFLAGS="\$CFLAGS -Wall -Wextra -ansi -pedantic           ...
        -Wformat-nonliteral -Wcast-align -Wpointer-arith             ...
        -Wbad-function-cast -Wmissing-prototypes -Wstrict-prototypes ...
        -Wmissing-declarations -Winline -Wundef -Wnested-externs     ...
        -Wcast-qual -Wshadow -Wconversion -Wwrite-strings            ...
        -Wno-conversion -Wchar-subscripts -Wredundant-decls"         ...
   ...
        -O -largeArrayDims         ...
        -I/usr/include/suitesparse ...
   ...
        mex_ifsh.c ifsh.c hybsys.c hybsys_global.c sparse_sys.c ...
        call_umfpack.c flow_bc.c well.c hash_set.c mimetic.c    ...
        mrst_api.c ...
   ...
       -lmwumfpack -lmwamd -lmwlapack -lmwblas

   % Call MEX'ed edition.
   [varargout{1:nargout}] = mex_ifsh(varargin{:});
end
