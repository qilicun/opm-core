/*
  Copyright 2015 Statoil ASA.

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

#include <opm/core/grid/CpGrid.hpp>
#include <opm/core/grid.h>
#include <opm/core/grid/GridManager.hpp>

#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>


int main(int argc, char** argv)
{
    using namespace Opm;
    parameter::ParameterGroup param(argc, argv, false);
    Opm::DeckConstPtr deck;
    std::unique_ptr<GridManager> grid;
    ParserPtr parser(new Opm::Parser());
    std::string deck_filename = param.get<std::string>("deck_filename");
    deck = parser->parseFile(deck_filename);
    grid.reset(new GridManager(deck));
    std::vector<Opm::TransHelper::Simplex> cpcell;
    createCpGridCells(*grid->c_grid(), cpcell);

    
    for (size_t c = 0; c < cpcell.size(); ++c) {
        std::cout << "cell " << c << std::endl;
        TransHelper::printVertex(cpcell[c].cell_centroid, "cell_centroid");
        //TransHelper::printVertex(cpcell[c].cell_centroid, "cell_centroid");
        
        for (size_t n = 0; n < cpcell[c].corner.size(); ++n) {
            TransHelper::printVertex(cpcell[c].corner[n], "cell_corner");
        }
    }

    return 0;
}
