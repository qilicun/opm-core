/*
  Copyright 2015, Li Xiang, Reservoir Simulation Lab, COE, PKU.
  Copyright 2015, Statoil ASA.

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

#ifndef OPM_TRANSHELPERS_HEADER_INCLUDED
#define OPM_TRANSHELPERS_HEADER_INCLUDED

#include <opm/core/grid/GridHelpers.hpp>
#include <array>
#include <climits>
#include <vector>
#include <array>
#include <algorithm>

namespace Opm {
namespace TransHelpers {
    
    typedef std::array<double, 3> vertex; // per points has 3 coordinates.

    inline vertex operate=(const vertex& v1, const vertex& v2) {
        return (v1[0] == v1[0] && v1[1] == v2[1] && v1[2] == v2[2]);
    }
        
    struct Simplex {
        std::array<vertex, 8> corner;
        int isActive;
        vertex cell_centroid;
        std::array<double, 6> face_areas;
        std::array<vertex, 6> face_normals;
        std::array<vertex, 6> face_centroids;
        double cell_volume;
    };

    void createCpGridCells(const UnstructuredGrid& grid,
                           const std::vector<int>& actnum,
                           std::vector<Simplex>& cpcell)
    {
        Simplex cell;
        const int nc = grid.number_of_cells;
        const int dim = grid.dimensions;
        
        for (int c = 0; c < nc; ++c) {
            std::vector<vertex> nodes;
            cell.cell_centroid.data() = Opm::UgGridHelpers::cellCentroid(g, c);
            cell.isActive = actnum[c];
            cell.cell_volume = Opm::UgGridHelpers::cellVolume(g, c);
            for (int i = grid.cell_facepos[c]; i < grid.cell_facepos[c+1]; ++i) {
                f = grid.cell_faces[i];
                double area = Opm::UgGridHelpers::faceArea(g, f);
                cell.face_areas.push_back(area);
                vertex tmp;
                tmp.data() = Opm::UgGridHelpers::faceNormal(g, f);
                cell.face_normals.push_back(tmp);
                tmp.clear();
                tmp.data() = Opm::UgGridHelpers::faceCentroid(g, f);
                cell.face_centroids.push_back(tmp);
                //get the corner coordinates.
                for (int nd = grid.face_nodepos[f]; nd < grid.face_nodepos[f+1]; ++nd) {
                    node_index = grid.face_nodes[nd];
                    vertex tmp2;
                    tmp2.data() = Opm::UgGridHelpers::vertexCoordinates(g, node_index);
                    nodes.push_back(tmp2);
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
    
    struct NNCs {
        int i;         // cell index.
        int j;         // cell index.
        double wi;     // trans coffecient.
        double wj;     // trans coffecient.
        double area;   // contract face area.
        int dir;       // directions.
        
        int left() const { return i; };
        int right() const {return j; };
        int dest() const { return j; };

        double alpha_i () const { return wi; };
        double alpha_j () const { return wj; };
        double Area    () const { return area; };
        double Dir     () const { dir; };

        int& left() const { return i; };
        int& right() const {return j; };
        int& dest() const { return j; };
        
        double& alpha_i () const { return wi; };
        double& alpha_j () const { return wj; };
        double& Area    () const { return area; };
        double& Dir     () const { dir; };
    };
    
    static void
    cross(const vertex& v1, const vertex& v2, const vertex& v)
    {
        v[0] = v1[2] * v2[1] - v1[1] * v2[2];
        v[1] = v1[0] * v2[2] - v1[2] * v2[0];
        v[2] = v1[0] * v2[1] - v1[1] * v2[0];
    }
    
    int crossType(const double z1, const double z2, const double z3, const double z4, const double za, const double zb)
    {
        if (za < z1) {
            if      (zb < z3) { return 1; }
            else if (zb < z4) { return 3; }
            else              { return 9; }
        } else if (za < z2) {
            if      (zb < z3) { return 2; }
            else if (zb < z4) { return 4; }
            else              { return 5; }
        } else {
            if      (zb < z3) { return 8; }
            else if (zb < z4) { return 6; }
            else              { return 7; }
        }
    }

    void contactType(int& conType, int& crossType1, int& crossType2, const corner& corner1, const corner& corner2)
    {
        double z1, z5, z3, z7, z0, z2, z4, z6;
        // For corner of cell 1.
        z1 = corner1[1][2];
        z3 = corner1[3][2];
        z5 = corner1[5][2];
        z7 = corner1[7][2];
        //For corner of cell 2.
        z0 = corner2[0][2];
        z2 = corner2[2][2];
        z4 = corner2[4][2];
        z6 = corner2[6][2];

        if (z1==z0 && z5==z4 && z7==z6 && z3==z2) {
            crossType1 = 1;
            crossType2 = 7;
            conType    = 2;
        } else {
            crossType1 = CrossType(z1, z5, z3, z7, z0, z2);
            crossType2 = CrossType(z1, z5, z3, z7, z4, z6);

            if ((crossType1==1 && crossType2==1) || (crossType1==7 && crossType2 ==7)) {
                conType = 0;
            } else {
                conType = 1;
            }
        }
    }

    int zContactType(const corner& corner1, corner& corner2)
    {
        if (corner2[0][2] - corner1[4][2] > std::numeric_limits<double>::min() ||
            corner2[1][2] - corner1[5][2] > std::numeric_limits<double>::min() ||
            corner2[2][2] - corner1[6][2] > std::numeric_limits<double>::min() ||
            corner2[3][2] - corner1[7][2] > std::numeric_limits<double>::min()) {
            
            return 0;
        } else {
            
            return 2;
        }
    }

    double zContactArea(const UnstructuredGrid& g, int index)
    {
        return Opm::UgGridHelpers::faceArea(g, index);
    }

    // Cross point of two lines, assuming they are on the same plane.
    vertex
    intersect(const vertex& endpt0,
              const vertex& endpt1,
              const vertex& endpt2,
              const vertex& endpt3)
    {
        double t, dx1, dx2, dy1, dy2, dz1, dz2, projA;
        vertex crosspoint;

        dx1 = (endpt1[0] - endpt0[0]);
        dx2 = (endpt3[0] - endpt2[0]);
        dy1 = (endpt1[1] - endpt0[1]);
        dy2 = (endpt3[1] - endpt2[1]);
        dz1 = (endpt1[2] - endpt0[2]);
        dz2 = (endpt3[2] - endpt2[2]);

        if (std::fabs(projA = (dy1*dz2 - dz1*dy2)) > std::numeric_limits<double>::min()) {
            t = ((endpt2[1] - endpt0[1])*dz2 - (endpt2[2] - endpt0[2])*dy2) / projA;
        } else if (std::fabs(projA = (dx1*dz2 - dz1*dx2)) > std::numeric_limits<double>::min()) {
            t = ((endpt2[0] - endpt0[0])*dz2 - (endpt2[2] - endpt0[2])*dx2) / projA;
        } else {
            t = ((endpt2[1] - endpt0[1])*dx2 - (endpt2[0] - endpt0[0])*dy2) / (dy1*dx2 - dx1*dy2);
        }

        crosspoint[0] = t*dx1 + endpt0[0];
        crosspoint[1] = t*dy1 + endpt0[1];
        crosspoint[2] = t*dz1 + endpt0[2];

        return crosspoint;
    }
    
    inline vertex
    triProj(const vertex& v1, const vertex& v2, const vertex& v3)
    {
        vertex A(0);

        A[0] = 0.5 * ((v3[2]-v2[2])*(v1[1]-v2[1]) - (v3[1]-v2[1])*(v1[2]-v2[2]));
        A[1] = 0.5 * ((v3[0]-v2[0])*(v1[2]-v2[2]) - (v3[2]-v2[2])*(v1[0]-v2[0]));
        A[2] = 0.5 * ((v3[1]-v2[1])*(v1[0]-v2[0]) - (v3[0]-v2[0])*(v1[1]-v2[1]));

        return A;
    }

    template<const int N>
    vertex polyProj(const std::vector<vertex>& v)
    {
        return polyProj<N-1>(v) + triProj(v[0], v[N-2], v[N-1]);
    }
    
    template<>
    inline vertex
    polyProj<3>(const std::vector<vertex>& v)
    {
        return triProj(v[0], v[1], v[2]);
    }


    vertex
    contactArea(const corner& corner1,
                const corner& corner2,
                const int crossType1,
                const int crossType2)
    {
        vertex A;
        std::vector<double> polyg;

        if (crossType1 == 1) {
            if (crossType2 == 2) { // #1
                polyg.push_back(corner1[1]);
                polyg.push_back(corner2[4]);
                polyg.push_back(intersect(corner1[1], corner1[3], corner2[4], corner2[6]));
                A = polygProj<3>(polyg);
                
            } else if (crossType2 == 3) { // #2
                polyg.push_back(corner1[3]);
                polyg.push_back(intersect(corner1[1], corner1[3], corner2[4], corner2[6]));
                polyg.push_back(corner2[6]);
                A = polyProj<3>(polyg);
                
            } else if (crosstype2 == 4) { // #3
                polyg.push_back(corner1[1]);
                polyg.push_back(corner2[4]);
                polyg.push_back(corner2[6]);
                polyg.push_back(corner1[3]);

                A = polyProj<4>(polyg);
                
            } else if ( crossType2 == 5) { // #4
                polyg.push_back(corner1[1]);
                polyg.push_back(corner2[4]);
                polyg.push_back(intersect(corner1[5], corner1[7], corner2[4], corner2[6]));
                polyg.push_back(corner1[7]);
                polyg.push_back(corner1[3]);

                A = polyProj<5>(polyg);
                
            } else if (crossType2 == 6) { // #5
                polyg.push_back(corner1[1]);
                polyg.push_back(corner1[5]);
                polyg.push_back(intersect(corner1[5], corner1[7], corner2[4], corner2[6]));
                polyg.push_back(corner2[6]);
                polyg.push_back(corner1[3]);

                A = polyProj<5>(polyg);
                
            } else if (crossType2 == 7) { // #6
                polyg.push_back(corner1[1]);
                polyg.push_back(corner1[5]);
                polyg.push_back(corner1[7]);
                polyg.push_back(corner1[3]);
                
                A = polyProj<4>(polyg);
                
            } else if (crossType2 == 8) { // #7
                polyg.push_back(corner1[1]);
                polyg.push_back(corner1[5]);
                polyg.push_back(intersect(corner1[5], corner1[7], corner2[4], corner2[6]));
                polyg.push_back(intersect(corner1[1], corner1[3], corner2[4], corner2[6]));
                
                A = polyProj<4>(polyg);
                
            } else if (crossType == 9) { // #8
                polyg.push_back(intersect(corner1[1], corner1[3], corner2[4], corner2[6]));
                polyg.push_back(intersect(corner1[5], corner1[7], corner2[4], corner2[6]));
                polyg.push_back(corner1[7]);
                polyg.push_back(corner1[3]);

                A = polyProj<4>(polyg);
                
            } else {
                A[0] = A[1] = A[2] = 0.0;
            }
        }

        else if (crossType1 == 2) {

            if (crossType2 == 2) { // #9
                polyg.push_back(corner2[0]);
                polyg.push_back(corner2[4]);
                polyg.push_back(intersect(corner1[1], corner1[3], corner2[4], corner2[6]));
                polyg.push_back(intersect(corner1[1], corner1[3], corner2[0], corner2[2]));

                A = polyProj<4>(polyg);
                
            } else if (crossType2 == 4) { // #10
                polyg.push_back(corner2[0]);
                polyg.push_back(corner2[4]);
                polyg.push_back(corner2[6]);
                polyg.push_back(corner1[3]);
                polyg.push_back(intersect(cornerl1[1], corner1[3], corner2[0], corner2[2]));

                A = polyProj<5>(polyg);
                
            } else if (crossType2 == 5) { // #11
                polyg.push_back(corner2[0]);
                polyg.push_back(corner2[4]);
                polyg.push_back(intersect(corner1[5], corner1[7], corner2[4], corner2[6]));
                polyg.push_back(corner1[7]);
                polyg.push_back(corner1[3]);
                polyg.push_back(intersect(corner1[1], corner1[3], corner2[0], corner2[2]));

                A = polyProj<6>(polyg);
                
            } else if (crossType2 == 6) { // #12
                polyg.push_back(corner2[0]);
                polyg.push_back(corner1[5]);
                polyg.push_back(intersect(corner1[1], corner1[3], corner2[0], corner2[2]));
                polyg.push_back(corner2[6]);
                polyg.push_back(corner1[3]);
                polyg.push_back(intersect(corner1[1], corner1[3], corner2[0], corner2[2]));

                A = polyProj<6>(polyg);
                
            } else if (crossType2 == 7) { //#13
                polyg.push_back(corner2[0]);
                polyg.push_back(corner1[5]);
                polyg.push_back(corner1[7]);
                polyg.push_back(corner1[3]);
                polyg.push_back(intersect(corner1[1], corner1[3], corner2[0], corner2[2]));

                A = polyProj<5>(polyg);
                
            } else if (crossType2 == 8) { // #14
                polyg.push_back(corner2[0]);
                polyg.push_back(corner1[5]);
                polyg.push_back(intersect(corner1[5], corner1[7], corner2[4], corner2[6]));
                polyg.push_back(intersect(corner1[1], corner1[3], corner2[4], corner2[6]));
                polyg.push_back(intersect(corner1[1], corner1[3], corner2[0], corner2[2]));

                A = polyProj<5>(polyg);
                
            } else {
                A[0] = A[1] = A[2] = 0.0;
            }
        }
        
        else if (crossType1 == 3) {
            if (crossType2 == 3) { // #15
                polyg.push_back(intersect(corner1[1], corner1[3], corner2[4], corner2[6]));
                polyg.push_back(corner2[6]);
                polyg.push_back(corner2[2]);
                polyg.push_back(intersect(corner1[1], corner1[3], corner2[0], corner2[2]));

                A = polyProj<4>(polyg);
                
            } else if (crossType2 == 4) { // #16
                polyg.push_back(corner1[1]);
                polyg.push_back(corner2[4]);
                polyg.push_back(corner2[6]);
                polyg.push_back(corner2[2]);
                polyg.push_back(intersect(corner1[1], corner1[3], corner2[0], corner2[2]));

                A = polyProj<5>(polyg);
                
            } else if (corssType2 == 5) { // #17
                polyg.push_back(corner1[1]);
                polyg.push_back(corner2[4]);
                polyg.push_back(intersect(corner1[5], corner1[7], corner2[4], corner2[6]));
                polyg.push_back(corner1[7]);
                polyg.push_back(corner2[2]);
                polyg.push_back(intersect(corner1[1], corner1[3], corner2[0], corner2[2]));

                A = polyProj<6>(polyg);
                
            } else if (crossType2 == 6) { // #18
                polyg.push_back(corner1[1]);
                polyg.push_back(corner1[5]);
                polyg.push_back(intersect(corner1[5], corner1[7], corner2[4], corner2[6]));
                polyg.push_back(corner2[6]);
                polyg.push_back(corner2[2]);
                polyg.push_back(intersect(corner1[1], corner1[3], corner2[0], corner2[2]));

                A = polyProj<6>(polyg);
                
            } else if (crossType2 == 7) { // #19
                polyg.push_back(corner1[1]);
                polyg.push_back(corner1[5]);
                polyg.push_back(corner1[7]);
                polyg.push_back(corner2[2]);
                polyg.push_back(intersect(corner1[1], corner1[3], corner2[0], corner2[2]));

                A = polyProj<5>(polyg);
                
            } else if (crossType2 == 9) { // #20
                polyg.push_back(intersect(corner1[1], corner1[3], corner2[4], corner2[6]));
                polyg.push_back(intersect(corner1[5], corner1[7], corner2[4], corner2[6]));
                polyg.push_back(corner1[7]);
                polyg.push_back(corner2[2]);
                polyg.push_back(intersect(corner1[1], corner1[3], corner2[0], corner2[2]));

                A = polyProj<5>(polyg);
                
            } else {
                A[0] = A[1] = A[2] = 0.0;
            }
        }

        else if (crossType1 == 4) {
            if (crossType2 == 4) { // #21
                polyg.push_back(corner2[0]);
                polyg.push_back(corner2[4]);
                polyg.push_back(corner2[6]);
                polyg.push_back(corner2[2]);

                A = polyProj<4>(polyg);
                
            } else if (crossType2 ==5 ) { // #22
                polyg.push_back(corner2[0]);
                polyg.push_back(corner2[4]);
                polyg.push_back(intersect(corner1[5], corner1[7], corner2[4], corner2[6]));
                polyg.push_back(corner1[7]);
                polyg.push_back(corner2[2]);

                A = polyProj<5>(polyg);
                
            } else if (crossType == 6) { // #23
                polyg.push_back(corner2[0]);
                polyg.push_back(corner1[5]);
                polyg.push_back(intersect(corner1[5], corner1[7], corner2[4], corner2[6]));
                polyg.push_back(corner2[6]);
                polyg.push_back(corner2[2]);

                A = polyProj<5>(polyg);
                
            } else if (crossType2 == 7) { //#24
                polyg.push_back(corner2[0]);
                polyg.push_back(corner1[5]);
                polyg.push_back(corner1[7]);
                polyg.push_back(corner2[2]);

                A = polyProj<4>(polyg);
                
            } else {
                A[0] = A[1] = A[2] = 0.0;
            }
        }

        else if (crossType1 == 5) {
            if (crossType2 == 5) { //#25
                polyg.push_back(corner2[0]);
                polyg.push_back(corner2[4]);
                polyg.push_back(intersect(corner1[5], corner1[7], corner2[4], corner2[6]));
                polyg.push_back(intersect(corner1[5], corner1[7], corner2[0], corner2[2]));

                A = polyProj<4>(polyg);
                
            } else if (crossType == 7) { // #26
                polyg.push_back(corner2[0]);
                polyg.push_back(corner1[5]);
                polyg.push_back(intersect(corner1[5], corner1[7], corner2[4], corner2[6]));

                A = polyProj<3>(polyg);
                
            } else {
                A[0] = A[1] = A[2] = 0.0;
            }
        }

        else if (crossType1 == 6) {
            if (crossType2 == 6) { // #27
                polyg.push_back(intersect(corner1[5], corner1[7], corner2[0], corner2[2]));
                polyg.push_back(intersect(corner1[5], corner1[7], corner2[4], corner2[6]));
                polyg.push_back(corner2[6]);
                polyg.push_back(corner2[2]);

                A = polyProj<4>(polyg);
                
            } else if (crosstype2 == 7) { // #28
                polyg.push_back(intersect(corner1[5], corner1[7], corner2[0], corner2[2]));
                polyg.push_back(corner1[7]);
                polyg.push_back(corner2[2]);

                A = polyProj<3>(polyg);
                
            } else {
                A[0] = A[1] = A[2] = 0.0;
            }
        }

        else if (crossType1 == 8) {
            if (crossType2 == 6) { // #29
                polyg.push_back(intersect(corner1[5], corner1[7], corner2[0], corner2[2]));
                polyg.push_back(intersect(corner1[5], corner1[7], corner2[4], corner2[6]));
                polyg.push_back(corner2[6]);
                polyg.push_back(corner1[3]);
                polyg.push_back(intersect(corner1[1], corner1[3], corner2[0], corner2[2]));

                A = polyProj<5>(polyg);
                
            } else if (crossType2 == 7) { // #30
                polyg.push_back(intersect(corner1[5], corner1[7], corner2[0], corner2[2]));
                polyg.push_back(corner1[7]);
                polyg.push_back(corner1[3]);
                polyg.push_back(intersect(corner1[1], corner1[3], corner2[0], corner2[2]));

                A = polyProj<4>(polyg);
                
            } else if (crossType2 == 8) { // #31
                polyg.push_back(intersect(corner1[5], corner1[7], corner2[0], corner2[2]));
                polyg.push_back(intersect(corner1[5], corner1[7], corner2[4], corner2[6]));
                polyg.push_back(intersect(corner1[1], corner1[3], corner2[4], corner2[6]));
                polyg.push_back(intersect(corner1[1], corner1[3], corner2[0], corner2[2]));
                A = polyProj<4>(polyg);
                
            } else {
                A[0] = A[1] = A[2] = 0.0;
            }
        }

        else if (crosstype1 == 9) {
            if (crossType2 == 5) { // #32
                polyg.push_back(corner1[1]);
                polyg.push_back(corner2[4]);
                polyg.push_back(intersect(corner1[5], corner1[7], corner2[4], corner2[6]));
                polyg.push_back(intersect(corner1[5], corner1[7], corner2[0], corner2[2]));
                polyg.push_back(intersect(corner1[1], corner1[3], corner2[0], corner2[2]));
                A = polyProj<5>(polyg);
                
            } else if (crossType2 == 7) { // #33
                polyg.push_back(corner1[1]);
                polyg.push_back(corner1[5]);
                polyg.push_back(intersect(corner1[5], corner1[7], corner2[0], corner2[2]));
                polyg.push_back(intersect(corner1[1], corner1[3], corner2[0], corner2[2]));
                A = polyProj<4>(polyg);
                
            } else if (crossType2 == 9) { // #34
                polyg.push_back(intersect(corner1[1], corner1[3], corner2[4], corner2[6]));
                polyg.push_back(intersect(corner1[5], corner1[7], corner2[4], corner2[6]));
                polyg.push_back(intersect(corner1[5], corner1[7], corner2[0], corner2[2]));
                polyg.push_back(intersect(corner1[1], corner1[3], corner2[0], corner2[2]));

                A = polyProj<4>(polyg);
                
            } else {
                A[0] = A[1] = A[2] = 0.0;
            }
        }

        else {
            A[0] = A[1] = A[2] = 0.0;
        }

        return A;
    }



    int NewTran(const CpCell* cell, const int NX, const int NY, const int NZ, const std::vector<int>& actnum, std::vector<NNCs>& tpcnc)
    {
        int nPinchOut = 0;
        const int nc = NX * NY * NZ;

        for (int c = 0; c < nc; ++c) {
            if (!actnum[c]) {
                continue;
            }
            int i, j, k;
            getLocalIJK(c, NX, NY, i, j, k);

            if (i+1 < NX) { // x-direction.
                for (int k2 = 0; k2 < NZ; k2++) {
                    int c2 = 1 + (NY*k2 + j) * NX + 1; // global index.
                    int conType, crossType1, crossType2;
                    contactType(conType, crossType1, crossType2, cell[c], cell[c2]);
                    if (conType != 0 && !actnum[c2]) {
                        vertex VA = contactArea(cell[c], cell[c2], crossType1, crossType2);
                        double A = VecNorm(VA);
                        NNCs tmp;
                        if (A > std::numeric_limits<double>::min()) {
                            tmp.left() = c;
                            tmp.right() = c2;
                            tmp.Area = A;
                            tmp.Dir = XDir;
                            tmp.alpha_i() = computeTransCoeff(VA, cell[c].face_centroids);
                            tmp.alpha_j() = computeTransCoeff(VA, cell[c2].face_centroids);
                            if (tmp.alpha_i() > std::numeric_limits<double>::min() && tmp.alpha_j() > std::numeric_limits<double>::min()) {
                                tpcnc.push_back(tmp);
                            }
                        } else {
                            ++nPinchOut;
                        }
                    } 
                    if (conType == 2) break;
                }
            }

            if (j + 1 < NY) { // y-direction.
                for (int k2 = 0; k2 < NZ; k2++) {
                    int c2 = i  + (j + 1 + k2 * NY) * NX;
                    int conType, crossType1, crossType2;
                    contactType(conType, crossType1, crossType2, cell[c], cell[c2]);
                    if (conType != 0 && !actnum[c2]) {
                        vertex VA = contactAres(cell[c], cell[c2], crossType1, crossType2);
                        double A = VecNorm(VA);
                        NNCs tmp;
                        if (A > std::numeric_limits<double>::min()) {
                            tmp.left() = c;
                            tmp.dest() = c2;
                            tmp.Area() = A;
                            tmp.Dir = YDir;
                            tmp.alpha_i() = computeTransCoeff(VA, cell[c].face_centroids);
                            tmp.alpha_j() = computeTransCoeff(VA, cell[c2].face_centroids);
                            if (tmp.alpha_i() > std::numeric_limits<double>::min() && tmp.alpha_j() > std::numeric_limits<double>::min()) {
                                tpcnc.push_back(tmp);
                            }
                        } else {
                            ++nPinchOut;
                        }
                    }
                    if (conType == 2) break;
                }
            }

            if (k + 1 < NZ) { // z-direction.
                int c2 = i + (j + (k + 1) * NY) * NX;
                int conType = zContactType(cell[c], cell[c2]);
                if (conType != 0 && !actnum[c2]) {
                    vertex VA = contactAreaZ(cell[c]);
                    NNCs tmp;
                    tmp.left() = c;
                    tmp.dest() = c2;
                    tmp.Area() = VecNorm(VA);
                    tmp.Dir() = ZDir;

                    tmp.alpha_i() = computeTransCoeff(VA, cell[c].face_centroids);
                    tmp.alpha_j() = computeTransCoeff(VA, cell[c2].face_centroids);
                    if (tmp.alpha_i() > std::numeric_limits<double>::min() && tmp.alpha_j() > std::numeric_limits<double>::min()) {
                        tpcnc.push_back(tmp);
                    }
                }
            }
        }

        return nPinchOut;
    }
} // namespace TransHelpers
} //namespace Opm
#endif // OPM_TRANSHELPERS_HEADER_INCLUDED
