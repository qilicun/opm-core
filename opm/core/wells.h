/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

  This file is part of the Open Porous Media project (OPM).

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

#ifndef OPM_WELLS_H_INCLUDED
#define OPM_WELLS_H_INCLUDED

#include <stdbool.h>
#include <opm/core/well_controls.h>

/**
 * \file
 *
 * Main OPM-Core well data structure along with functions
 * to create, populate and destroy it.
 */


#ifdef __cplusplus
extern "C" {
#endif

/**
 * Well type indicates desired/expected well behaviour.
 */
enum WellType {
    INJECTOR,  /**< Well is an injector */
    PRODUCER   /**< Well is a producer */
};


/**
 *  Data structure aggregating static information about all wells in a scenario.
 */
struct Wells
{
    int number_of_wells;  /**< Number of wells. */
    int number_of_phases; /**< Number of phases. */

    /**
     * Array of well types.
     */
    enum WellType *type;

    /**
     * Array of well reference depths.
     */
    double *depth_ref;

    /**
     * Component fractions for each well.  Array of size
     * <CODE>number_of_wells * number_of_phases</CODE>.
     * For injection wells, this gives the injected component mix.
     * For production wells the component fractions of the wellbore
     * will vary and cannot be specified a priori, the component mix
     * given here should be considered a default or preferred mix.
     */
    double *comp_frac;

    /**
     * Array of indices into well_cells (and WI).  For a well @c w,
     * <CODE>well_connpos[w]</CODE> and <CODE>well_connpos[w+1]</CODE> are start
     * and one-beyond-end indices into the @c well_cells array for accessing
     * @c w's perforation cell indices.
     */
    int *well_connpos;

    /**
     * Array of perforation cell indices.
     * Size is number of perforations (== well_connpos[number_of_wells]).
     */
    int *well_cells;

    /**
     *  Well productivity index, same size and structure as well_cells.
     */
    double *WI;


    /**
     * Well controls, one set of controls for each well.
     */
    struct WellControls **ctrls;

    /**
     * Well names. One string for each well.
     */
    char **name;

    /**
     * Internal management structure.
     */
    void *data;
};


/**
 * Data structure aggregating dynamic information about all wells in a scenario.
 * All arrays in this structure contain data for each perforation, ordered the
 * same as Wells::well_cells and Wells:WI.  The array sizes are, respectively,
 *
 *     wdp        NP
 *     A          n²*NP (matrix in column-major (i.e., Fortran) order).
 *     phasemob   n*NP
 *
 * in which "n" denotes the number of active fluid phases (and constituent
 * components) and "NP" is the total number of perforations,
 * <CODE>well_connpos[ number_of_wells ]</CODE>.
 */
struct CompletionData
{
    /**
     * Gravity potentials.
     */
    double *wdp;

    /**
     * Volumes to surface-components matrix, A = RB^{-1}.
     */
    double *A;

    /**
     * Phase mobilities for all perforations, stored consecutively with the
     * phase index cycling the most rapidly.
     */
    double *phasemob;
};

/**
 * Construct a Wells object initially capable of managing a given
 * number of wells and total number of well connections
 * (perforations).
 *
 * Function add_well() is used to populate the Wells object.  No
 * reallocation occurs in function add_well() as long as the
 * initially indicated capacities are sufficient.  Call function
 * destroy_wells() to dispose of the Wells object and its allocated
 * memory resources.
 *
 * \param[in] nphases Number of active phases in simulation scenario.
 *
 * \param[in] nwells  Expected number of wells in simulation scenario.
 *                    Pass zero if the total number of wells is unknown.
 *
 * \param[in] nperf   Expected total number of well connections
 *                    (perforations) for all wells in simulation
 *                    scenario.  Pass zero if the total number of well
 *                    connections is unknown.
 *
 * \return A valid Wells object with no wells if successful, and NULL
 * otherwise.
 */
struct Wells *
create_wells(int nphases, int nwells, int nperf);


/**
 * Append a new well to an existing Wells object.
 *
 * Increments W->number_of_wells by one if successful.  The new well
 * does not include operational constraints.  Such information is
 * specified using function append_well_controls().  The current
 * control index is set to -1 (invalid).
 *
 * \param[in] type       Type of well.
 * \param[in] depth_ref  Reference depth for well's BHP.
 * \param[in] nperf      Number of perforations.
 * \param[in] comp_frac  Injection fraction array (size equal to W->number_of_phases) or NULL.
 * \param[in] cells      Grid cells in which well is perforated.  Should
 *                       ideally be track ordered.
 * \param[in] WI         Well production index per perforation, or NULL.
 * \param[in] name       Name of new well. NULL if no name.
 * \param[in,out] W      Existing set of wells to which new well will
 *                       be added.
 *
 * \return Non-zero (true) if successful and zero otherwise.
 */
int
add_well(enum WellType  type     ,
         double         depth_ref,
         int            nperf    ,
         const double  *comp_frac,
         const int     *cells    ,
         const double  *WI       ,
         const char    *name     ,
         struct Wells  *W        );


/**
 * Append operational constraint to an existing well.
 *
 * Increments ctrl->num by one if successful.  Introducing a new
 * operational constraint does not affect the well's notion of the
 * currently active constraint represented by ctrl->current.
 * Note that *_RATE controls now require a phase distribution array
 * to be associated with the control, see WellControls.
 *
 * \param[in] type       Control type.
 * \param[in] target     Target value for the control.
 * \param[in] distr      Array of size W->number_of_phases or NULL.
 * \param[in] well_index Index of well to receive additional control.
 * \param[in,out] W  Existing set of well controls.
 * \return Non-zero (true) if successful and zero (false) otherwise.
 */


int
append_well_controls(enum WellControlType type  ,
                     double               target,
                     const double        *distr,
                     int                  well_index,
                     struct Wells        *W);


/**
 * Set the current/active control for a single well.
 *
 * The new control ID must refer to a previously defined control mode.
 * Total number of defined control modes available through function
 * well_controls_get_num().
 *
 * \param[in]     well_index
 *                   Identity of particular well.  Must be in
 *                   \code [0 .. number_of_wells - 1] \endcode.
 *
 * \param[in]     current_control
 *                   Index of new control mode.
 *
 * \param[in,out] W  Existing set of wells.
 */
void
set_current_control(int well_index, int current_control, struct Wells *W);


/**
 * Clear all controls from a single well.
 *
 * Does not affect the control set capacity.
 *
 * \param[in]     well_index
 *                   Identity of particular well.  Must be in
 *                   \code [0 .. number_of_wells - 1] \endcode.
 *
 * \param[in,out] W  Existing set of wells.
 */
void
clear_well_controls(int well_index, struct Wells *W);


/**
 * Wells object destructor.
 *
 * Disposes of all resources managed by the Wells object.
 *
 * The Wells object must be built using function create_wells() and
 * subsequently populated using function add_well().
 */
void
destroy_wells(struct Wells *W);


/**
 * Create a deep-copy (i.e., clone) of an existing Wells object, including its
 * controls.
 *
 * @param[in] W Existing Wells object.
 * @return Complete clone of the input object.  Dispose of resources using
 * function destroy_wells() when no longer needed.  Returns @c NULL in case of
 * allocation failure.
 */
struct Wells *
clone_wells(const struct Wells *W);


/**
 * Compare well structures for equality.
 *
 * Two sets of wells are equal if all of the following conditions hold
 *  - They have the same number of wells
 *
 *  - They have the same number of completions/connections
 *
 *  - They specify the same number of phases
 *
 *  - Individual wells with corresponding well IDs have the same names
 *    (including both being \c NULL).
 *
 *  - Individual wells with corresponding well IDs have the same
 *    completions
 *
 *  - Individual wells with corresponding well IDs have the same well
 *    types
 *
 *  - Individual wells with corresponding well IDs specify the same
 *    reference depths
 *
 *  - Individual wells with corresponding well IDs have the same set
 *    of defined and active operational constraints as determined by
 *    function well_controls_equal()
 *
 * \param[in] W1      Existing set of wells.
 * \param[in] W2      Existing set of wells.
 *
 * \param[in] verbose Flag for whether or not to report which
 *                    conditions do not hold.  Use \code verbose =
 *                    true \endcode to print transcript to \c stdout.
 *
 * \return Whether or not well structures \c W1 and \c W2 represent
 * the same set of wells.
 */
bool
wells_equal(const struct Wells *W1, const struct Wells *W2 , bool verbose);


#ifdef __cplusplus
}
#endif

#endif /* OPM_WELLS_H_INCLUDED */
