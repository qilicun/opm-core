/*
  Copyright (c) 2013 Uni Research AS

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

#ifndef OPM_OUTPUT_WRITER_HPP
#define OPM_OUTPUT_WRITER_HPP

#include <memory>  // unique_ptr, shared_ptr

struct UnstructuredGrid;

namespace Opm {

// forward declaration
class EclipseState;
namespace parameter { class ParameterGroup; }
class SimulatorState;
class SimulatorTimer;
class WellState;
struct PhaseUsage;

/*!
 * Interface for writing non-compositional (blackoil, two-phase) simulation
 * state to files.
 *
 * Use the create() function to setup a chain of writer based on the
 * configuration values, e.g.
 *
 * \example
 * \code{.cpp}
 *  ParameterGroup params (argc, argv, false);
 *  auto parser = std::make_shared <const Deck> (
 *                      params.get <string> ("deck_filename"));
 *
 *  std::unique_ptr <OutputWriter> writer =
 *          OutputWriter::create (params, parser);
 *
 *  // before the first timestep
 *  writer->writeInit (timer);
 *
 *  // after each timestep
 *  writer->writeTimeStep (timer, state, wellState);
 *
 * \endcode
 */
class OutputWriter {
public:
    /// Allow derived classes to be used in the unique_ptr that is returned
    /// from the create() method. (Every class that should be delete'd should
    /// have a proper constructor, and if the base class isn't virtual then
    /// the compiler won't call the right one when the unique_ptr goes out of
    /// scope).
    virtual ~OutputWriter () { }

    /**
     * Write the static data (grid, PVT curves, etc) to disk.
     *
     * This routine should be called before the first timestep (i.e. when
     * timer.currentStepNum () == 0)
     */
    virtual void writeInit(const SimulatorTimer &timer) = 0;

    /*!
     * \brief Write a blackoil reservoir state to disk for later inspection with
     *        visualization tools like ResInsight
     *
     * \param[in] reservoirState The thermodynamic state of the reservoir
     * \param[in] wellState The production/injection data for all wells
     *
     * This routine should be called after the timestep has been advanced,
     * i.e. timer.currentStepNum () > 0.
     */
    virtual void writeTimeStep(const SimulatorTimer& timer,
                               const SimulatorState& reservoirState,
                               const WellState& wellState) = 0;

    /*!
     * Create a suitable set of output formats based on configuration.
     *
     * @param params Configuration properties. This function will setup a
     *               multiplexer of applicable output formats based on the
     *               desired configuration values.
     *
     * @param deck Input deck used to set up the simulation.
     *
     * @param eclipseState The internalized input deck.
     *
     * @return       Pointer to a multiplexer to all applicable output formats.
     *
     * @see Opm::share_obj
     */
    static std::unique_ptr <OutputWriter>
    create (const parameter::ParameterGroup& params,
            std::shared_ptr <const EclipseState> eclipseState,
            const Opm::PhaseUsage &phaseUsage,
            std::shared_ptr <const UnstructuredGrid> grid);
};

} // namespace Opm

#endif /* OPM_OUTPUT_WRITER_HPP */
