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

#ifndef OPM_CPGRID_HEADER_INCLUDED
#define OPM_CPGRID_HEADER_INCLUDED


#include <opm/core/grid/GridHelpers.hpp>
#include <array>
#include <climits>
#include <vector>
#include <array>
#include <algorithm>

namespace Opm {
namespace TransHelper {

    typedef std::array<double, 3> vertex; // per vertex has 3 coordinates.

    struct Simplex {
        Simplex()  {};
        std::vector<vertex> corner;
        vertex cell_centroid;
        std::vector<double> face_areas;
        std::vector<vertex> face_normals;
        std::vector<vertex> face_centroids;
        double cell_volume;
    };

    void createCpGridCells(const UnstructuredGrid& grid,
                           //const std::vector<int>& actnum,
                           std::vector<Simplex>& cpcell)
    {
        Simplex cell;
        const int nc = grid.number_of_cells;
        const int dim = grid.dimensions;
        
        for (int c = 0; c < nc; ++c) {
            std::vector<vertex> nodes;
            const double* cc = Opm::UgGridHelpers::cellCentroid(grid, c);
            std::copy(cc, cc+dim, cell.cell_centroid.begin());
            cell.cell_volume = Opm::UgGridHelpers::cellVolume(grid, c);
            for (int i = grid.cell_facepos[c]; i < grid.cell_facepos[c+1]; ++i) {
                int f = grid.cell_faces[i];
                double area = Opm::UgGridHelpers::faceArea(grid, f);
                cell.face_areas.push_back(area);
                vertex tmp;
                const double* fn = Opm::UgGridHelpers::faceNormal(grid, f);
                std::copy(fn, fn+dim, tmp.begin());
                cell.face_normals.push_back(tmp);

                vertex tmp2;
                const double* fc = Opm::UgGridHelpers::faceCentroid(grid, f);
                std::copy(fc, fc+dim, tmp2.begin());
                cell.face_centroids.push_back(tmp2);
                
                //get the corner coordinates.
                for (int nd = grid.face_nodepos[f]; nd < grid.face_nodepos[f+1]; ++nd) {
                    int node_index = grid.face_nodes[nd];
                    vertex tmp3;
                    const double* vc = Opm::UgGridHelpers::vertexCoordinates(grid, node_index);
                    std::copy(vc, vc+dim, tmp3.begin());
                    nodes.push_back(tmp3);
                }
            }
            // Uniform the normals.
            for (int f = 0; f < 6; ++f) {
                for (int n = 0; n < 3; ++n) {
                    cell.face_normals[f][n] /= cell.face_areas[f];
                }
            }

            // unique the corners.
            for (size_t n = 0; n < nodes.size(); ++n) {
                if (std::count(cell.corner.begin(), cell.corner.end(), nodes[n]) < 1) {
                    cell.corner.push_back(nodes[n]);
                }
            }
            cpcell.push_back(cell);
        }
    }

} // namespace TransHelper.
} // namespace Opm.

#endif // OPM_CPGRID_HEADER_INCLUDED
