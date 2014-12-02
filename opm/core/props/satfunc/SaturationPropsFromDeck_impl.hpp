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

#ifndef OPM_SATURATIONPROPSFROMDECK_IMPL_HEADER_INCLUDED
#define OPM_SATURATIONPROPSFROMDECK_IMPL_HEADER_INCLUDED


#include <opm/core/utility/UniformTableLinear.hpp>
#include <opm/core/utility/NonuniformTableLinear.hpp>
#include <opm/core/props/phaseUsageFromDeck.hpp>
#include <opm/core/grid.h>
#include <opm/core/grid/GridHelpers.hpp>

#include <opm/parser/eclipse/Utility/EndscaleWrapper.hpp>
#include <opm/parser/eclipse/Utility/ScalecrsWrapper.hpp>

#include <iostream>
#include <map>

namespace Opm
{


    // ----------- Methods of SaturationPropsFromDeck ---------


    /// Default constructor.
    template <class SatFuncSet>
    SaturationPropsFromDeck<SatFuncSet>::SaturationPropsFromDeck()
    {
    }

    /// Initialize from deck.
    template <class SatFuncSet>
    void SaturationPropsFromDeck<SatFuncSet>::init(Opm::DeckConstPtr deck,
                                                   Opm::EclipseStateConstPtr eclipseState,
                                                   const UnstructuredGrid& grid,
                                                   const int samples)
    {
        this->init(deck, eclipseState, grid.number_of_cells,
                   grid.global_cell, grid.cell_centroids,
                   grid.dimensions, samples);
    }

    /// Initialize from deck.
    template <class SatFuncSet>
    template<class T>
    void SaturationPropsFromDeck<SatFuncSet>::init(Opm::DeckConstPtr deck,
                                                   Opm::EclipseStateConstPtr eclipseState,
                                                   int number_of_cells,
                                                   const int* global_cell,
                                                   const T& begin_cell_centroids,
                                                   int dimensions,
                                                   const int samples)
    {
        phase_usage_ = phaseUsageFromDeck(deck);

        // Extract input data.
        // Oil phase should be active.
        if (!phase_usage_.phase_used[Liquid]) {
            OPM_THROW(std::runtime_error, "SaturationPropsFromDeck::init()   --  oil phase must be active.");
        }
        
        // Check SATOPTS status
        bool hysteresis_switch = false;
        if (deck->hasKeyword("SATOPTS")) {
            const std::vector<std::string>& satopts = deck->getKeyword("SATOPTS")->getStringData();
            for (size_t i = 0; i < satopts.size(); ++i) {
                if (satopts[i] == std::string("HYSTER")) {
                    hysteresis_switch = true;
                } else {
                    OPM_THROW(std::runtime_error, "Keyword SATOPTS:  Switch " << satopts[i] << " not supported. ");
                }
            }
        }

        // Obtain SATNUM, if it exists, and create cell_to_func_.
        // Otherwise, let the cell_to_func_ mapping be just empty.
        int satfuncs_expected = 1;
        if (deck->hasKeyword("SATNUM")) {
            const std::vector<int>& satnum = deck->getKeyword("SATNUM")->getIntData();
            satfuncs_expected = *std::max_element(satnum.begin(), satnum.end());
            const int num_cells = number_of_cells;
            cell_to_func_.resize(num_cells);
            const int* gc = global_cell;
            for (int cell = 0; cell < num_cells; ++cell) {
                const int deck_pos = (gc == NULL) ? cell : gc[cell];
                cell_to_func_[cell] = satnum[deck_pos] - 1;
            }
        }

        // Find number of tables, check for consistency.
        enum { Uninitialized = -1 };
        int num_tables = Uninitialized;
        if (phase_usage_.phase_used[Aqua]) {
            num_tables = deck->getKeyword("SWOF")->size();
            if (num_tables < satfuncs_expected) {
                OPM_THROW(std::runtime_error, "Found " << num_tables << " SWOF tables, SATNUM specifies at least " << satfuncs_expected);
            }
        }
        if (phase_usage_.phase_used[Vapour]) {
            int num_sgof_tables = deck->getKeyword("SGOF")->size();
            if (num_sgof_tables < satfuncs_expected) {
                OPM_THROW(std::runtime_error, "Found " << num_tables << " SGOF tables, SATNUM specifies at least " << satfuncs_expected);
            }
            if (num_tables == Uninitialized) {
                num_tables = num_sgof_tables;
            } else if (num_tables != num_sgof_tables) {
                OPM_THROW(std::runtime_error, "Inconsistent number of tables in SWOF and SGOF.");
            }
        }

        // Initialize tables.
        satfuncset_.resize(num_tables);
        for (int table = 0; table < num_tables; ++table) {
            satfuncset_[table].init(eclipseState, table, phase_usage_, samples);
        }
        
        // Check EHYSTR status
        do_hyst_ = false;
        if (hysteresis_switch && deck->hasKeyword("EHYSTR")) {
           const int& relative_perm_hyst = deck->getKeyword("EHYSTR")->getRecord(0)->getItem(1)->getInt(0);
           const std::string& limiting_hyst_flag = deck->getKeyword("EHYSTR")->getRecord(0)->getItem(4)->getString(0);
           if (relative_perm_hyst != int(0)) {
               OPM_THROW(std::runtime_error, "Keyword EHYSTR, item 2: Flag '" << relative_perm_hyst << "' found, only '0' is supported. ");
           }          
           if (limiting_hyst_flag != std::string("KR")) {
               OPM_THROW(std::runtime_error, "Keyword EHYSTR, item 5: Flag '" << limiting_hyst_flag << "' found, only 'KR' is supported. ");
           }                   
           if ( ! deck->hasKeyword("ENDSCALE")) {
               // TODO When use of IMBNUM is implemented, this constraint will be lifted.
               OPM_THROW(std::runtime_error, "Currently hysteris effects is only available through endpoint scaling.");
           }
           do_hyst_ = true;
        } else if (hysteresis_switch) {
           OPM_THROW(std::runtime_error, "Switch HYSTER of keyword SATOPTS is active, but keyword EHYSTR not found.");
        } else if (deck->hasKeyword("EHYSTR")) {
           OPM_THROW(std::runtime_error, "Found keyword EHYSTR, but switch HYSTER of keyword SATOPTS is not set.");
        }

        // Saturation table scaling
        do_eps_ = false;
        do_3pt_ = false;
        if (deck->hasKeyword("ENDSCALE")) {
            Opm::EndscaleWrapper endscale(deck->getKeyword("ENDSCALE"));
            if (endscale.directionSwitch() != std::string("NODIR")) {
                OPM_THROW(std::runtime_error,
                          "SaturationPropsFromDeck::init()   --  ENDSCALE: "
                          "Currently only 'NODIR' accepted.");
            }
            if (!endscale.isReversible()) {
                OPM_THROW(std::runtime_error,
                          "SaturationPropsFromDeck::init()   --  ENDSCALE: "
                          "Currently only 'REVERS' accepted.");
            }
            if (deck->hasKeyword("SCALECRS")) {
                Opm::ScalecrsWrapper scalecrs(deck->getKeyword("SCALECRS"));
                if (scalecrs.isEnabled()) {
                    do_3pt_ = true;
                }
            }
            do_eps_ = true;
                 
            // Make a consistency check of ENDNUM: #regions = NTENDP (ENDSCALE::3, TABDIMS::8)...      
            if (deck->hasKeyword("ENDNUM")) {
                const std::vector<int>& endnum = deck->getKeyword("ENDNUM")->getIntData();
                int endnum_regions = *std::max_element(endnum.begin(), endnum.end());
                if (endnum_regions > endscale.numEndscaleTables()) {
                    OPM_THROW(std::runtime_error,
                        "ENDNUM:  Found " << endnum_regions << 
                         " regions.  Maximum allowed is " << endscale.numEndscaleTables() <<
                         " (confer item 3 of keyword ENDSCALE).");  // TODO  See also item 8 of TABDIMS ...
                }
            }
            // TODO: ENPTVD/ENKRVD: Too few tables gives a cryptical message from parser, 
            //       superfluous tables are ignored by the parser without any warning ...

            const std::vector<std::string> eps_kw{"SWL", "SWU", "SWCR", "SGL", "SGU", "SGCR", "SOWCR",
                "SOGCR", "KRW", "KRG", "KRO", "KRWR", "KRGR", "KRORW", "KRORG", "PCW", "PCG"};
            eps_transf_.resize(number_of_cells);
            initEPS(deck, eclipseState, number_of_cells, global_cell, begin_cell_centroids,
                    dimensions, eps_kw, eps_transf_);

            if (do_hyst_) {
                if (deck->hasKeyword("KRW")
                    || deck->hasKeyword("KRG")
                    || deck->hasKeyword("KRO")
                    || deck->hasKeyword("KRWR")
                    || deck->hasKeyword("KRGR")
                    || deck->hasKeyword("KRORW")
                    || deck->hasKeyword("KRORG")
                    || deck->hasKeyword("ENKRVD")
                    || deck->hasKeyword("IKRG")
                    || deck->hasKeyword("IKRO")
                    || deck->hasKeyword("IKRWR")
                    || deck->hasKeyword("IKRGR")
                    || deck->hasKeyword("IKRORW")
                    || deck->hasKeyword("IKRORG") ) {
                    OPM_THROW(std::runtime_error,"Currently hysteresis and relperm value scaling cannot be combined.");
                }
                
                if (deck->hasKeyword("IMBNUM")) {
                    const std::vector<int>& imbnum = deck->getKeyword("IMBNUM")->getIntData();
                    int imbnum_regions = *std::max_element(imbnum.begin(), imbnum.end());
                    if (imbnum_regions > num_tables) {
                        OPM_THROW(std::runtime_error,
                            "IMBNUM:  Found " << imbnum_regions << 
                            " regions.  Maximum allowed is " << num_tables <<
                            " (number of tables provided by SWOF/SGOF).");
                    }
                    const int num_cells = number_of_cells;
                    cell_to_func_imb_.resize(num_cells);
                    const int* gc = global_cell;
                    for (int cell = 0; cell < num_cells; ++cell) {
                        const int deck_pos = (gc == NULL) ? cell : gc[cell];
                        cell_to_func_imb_[cell] = imbnum[deck_pos] - 1;
                    }
                    // TODO: Make actual use of IMBNUM.  For now we just consider the imbibition curve
                    //       to be a scaled version of the drainage curve (confer Norne model).
                }

                const std::vector<std::string> eps_i_kw{"ISWL", "ISWU", "ISWCR", "ISGL", "ISGU", "ISGCR", "ISOWCR",
                    "ISOGCR", "IKRW", "IKRG", "IKRO", "IKRWR", "IKRGR", "IKRORW", "IKRORG", "IPCW", "IPCG"};
                eps_transf_hyst_.resize(number_of_cells);
                sat_hyst_.resize(number_of_cells);                
                initEPS(deck, eclipseState, number_of_cells, global_cell, begin_cell_centroids,
                        dimensions, eps_i_kw, eps_transf_hyst_);
            }
        }
    }




    /// \return   P, the number of phases.
    template <class SatFuncSet>
    int SaturationPropsFromDeck<SatFuncSet>::numPhases() const
    {
        return phase_usage_.num_phases;
    }




    /// Relative permeability.
    /// \param[in]  n      Number of data points.
    /// \param[in]  s      Array of nP saturation values.
    /// \param[in]  cells  Array of n cell indices to be associated with the s values.
    /// \param[out] kr     Array of nP relperm values, array must be valid before calling.
    /// \param[out] dkrds  If non-null: array of nP^2 relperm derivative values,
    ///                    array must be valid before calling.
    ///                    The P^2 derivative matrix is
    ///                           m_{ij} = \frac{dkr_i}{ds^j},
    ///                    and is output in Fortran order (m_00 m_10 m_20 m01 ...)
    template <class SatFuncSet>
    void SaturationPropsFromDeck<SatFuncSet>::relperm(const int n,
                                          const double* s,
                                          const int* cells,
                                          double* kr,
                                          double* dkrds) const
    {
        assert(cells != 0);

        const int np = phase_usage_.num_phases;
        if (dkrds) {
// #pragma omp parallel for
            for (int i = 0; i < n; ++i) {
                if (do_hyst_) {
                   funcForCell(cells[i]).evalKrDeriv(s + np*i, kr + np*i, dkrds + np*np*i, &(eps_transf_[cells[i]]), &(eps_transf_hyst_[cells[i]]), &(sat_hyst_[cells[i]]));
                } else if (do_eps_) {
                   funcForCell(cells[i]).evalKrDeriv(s + np*i, kr + np*i, dkrds + np*np*i, &(eps_transf_[cells[i]]));
                } else {
                   funcForCell(cells[i]).evalKrDeriv(s + np*i, kr + np*i, dkrds + np*np*i);
                }
            }
        } else {
// #pragma omp parallel for
            for (int i = 0; i < n; ++i) {
                if (do_hyst_) {
                   funcForCell(cells[i]).evalKr(s + np*i, kr + np*i, &(eps_transf_[cells[i]]), &(eps_transf_hyst_[cells[i]]), &(sat_hyst_[cells[i]]));
                } else if (do_eps_) {
                   funcForCell(cells[i]).evalKr(s + np*i, kr + np*i, &(eps_transf_[cells[i]]));
                } else {
                   funcForCell(cells[i]).evalKr(s + np*i, kr + np*i);
                }
            }
        }
    }




    /// Capillary pressure.
    /// \param[in]  n      Number of data points.
    /// \param[in]  s      Array of nP saturation values.
    /// \param[in]  cells  Array of n cell indices to be associated with the s values.
    /// \param[out] pc     Array of nP capillary pressure values, array must be valid before calling.
    /// \param[out] dpcds  If non-null: array of nP^2 derivative values,
    ///                    array must be valid before calling.
    ///                    The P^2 derivative matrix is
    ///                           m_{ij} = \frac{dpc_i}{ds^j},
    ///                    and is output in Fortran order (m_00 m_10 m_20 m01 ...)
    template <class SatFuncSet>
    void SaturationPropsFromDeck<SatFuncSet>::capPress(const int n,
                                           const double* s,
                                           const int* cells,
                                           double* pc,
                                           double* dpcds) const
    {
        assert(cells != 0);

        const int np = phase_usage_.num_phases;
        if (dpcds) {
// #pragma omp parallel for
            for (int i = 0; i < n; ++i) {
                if (do_eps_) {
                   funcForCell(cells[i]).evalPcDeriv(s + np*i, pc + np*i, dpcds + np*np*i, &(eps_transf_[cells[i]]));
                } else {
                   funcForCell(cells[i]).evalPcDeriv(s + np*i, pc + np*i, dpcds + np*np*i);
                }
            }
        } else {
// #pragma omp parallel for
            for (int i = 0; i < n; ++i) {         
                if (do_eps_) {
                   funcForCell(cells[i]).evalPc(s + np*i, pc + np*i, &(eps_transf_[cells[i]]));
                } else {
                   funcForCell(cells[i]).evalPc(s + np*i, pc + np*i);
                }
            }
        }
    }




    /// Obtain the range of allowable saturation values.
    /// \param[in]  n      Number of data points.
    /// \param[in]  cells  Array of n cell indices.
    /// \param[out] smin   Array of nP minimum s values, array must be valid before calling.
    /// \param[out] smax   Array of nP maximum s values, array must be valid before calling.
    template <class SatFuncSet>
    void SaturationPropsFromDeck<SatFuncSet>::satRange(const int n,
                                           const int* cells,
                                           double* smin,
                                           double* smax) const
    {
        assert(cells != 0);
        const int np = phase_usage_.num_phases;
       
        if (do_eps_) {
            const int wpos = phase_usage_.phase_pos[BlackoilPhases::Aqua];
            const int opos = phase_usage_.phase_pos[BlackoilPhases::Liquid];
            const int gpos = phase_usage_.phase_pos[BlackoilPhases::Vapour];
            for (int i = 0; i < n; ++i) {
                smin[np*i + opos] = 1.0;
                smax[np*i + opos] = 1.0;
                if (phase_usage_.phase_used[Aqua]) {
                    smin[np*i + wpos] = eps_transf_[cells[i]].wat.doNotScale ? funcForCell(cells[i]).smin_[wpos]
                                                                             : eps_transf_[cells[i]].wat.smin;
                    smax[np*i + wpos] = eps_transf_[cells[i]].wat.doNotScale ? funcForCell(cells[i]).smax_[wpos]
                                                                             : eps_transf_[cells[i]].wat.smax;
                    smin[np*i + opos] -= smax[np*i + wpos];
                    smax[np*i + opos] -= smin[np*i + wpos];
                }  
                if (phase_usage_.phase_used[Vapour]) {
                    smin[np*i + gpos] = eps_transf_[cells[i]].gas.doNotScale ? funcForCell(cells[i]).smin_[gpos]
                                                                             : eps_transf_[cells[i]].gas.smin;
                    smax[np*i + gpos] = eps_transf_[cells[i]].gas.doNotScale ? funcForCell(cells[i]).smax_[gpos]
                                                                             : eps_transf_[cells[i]].gas.smax;
                    smin[np*i + opos] -= smax[np*i + gpos];
                    smax[np*i + opos] -= smin[np*i + gpos];
                }
                if (phase_usage_.phase_used[Vapour] && phase_usage_.phase_used[Aqua]) {
                    smin[np*i + opos] = std::max(0.0,smin[np*i + opos]);
                }
            }
        } else {
            for (int i = 0; i < n; ++i) {
                for (int p = 0; p < np; ++p) {
                    smin[np*i + p] = funcForCell(cells[i]).smin_[p];
                    smax[np*i + p] = funcForCell(cells[i]).smax_[p];
                }
            }
        }
    }

        
    /// Update saturation state for the hysteresis tracking 
    /// \param[in]  n      Number of data points. 
    /// \param[in]  s      Array of nP saturation values.
    template <class SatFuncSet>
    void SaturationPropsFromDeck<SatFuncSet>::updateSatHyst(const int n,
                                                            const int* cells,
                                                            const double* s)
    {        
        assert(cells != 0);

        const int np = phase_usage_.num_phases;
        if (do_hyst_) {
// #pragma omp parallel for
            for (int i = 0; i < n; ++i) {
                funcForCell(cells[i]).updateSatHyst(s + np*i, &(eps_transf_[cells[i]]), &(eps_transf_hyst_[cells[i]]), &(sat_hyst_[cells[i]]));
            }
        } 
    }


    /// Update capillary pressure scaling according to pressure diff. and initial water saturation.
    /// \param[in]     cell  Cell index.
    /// \param[in]     pcow  P_oil - P_water.
    /// \param[in/out] swat  Water saturation. / Possibly modified Water saturation.
    template <class SatFuncSet>
    void SaturationPropsFromDeck<SatFuncSet>::swatInitScaling(const int cell,
                                                              const double pcow,
                                                              double& swat)
    {
        if (phase_usage_.phase_used[BlackoilPhases::Aqua]) {
            const double pc_low_threshold = 1.0e-8;
            // TODO: Mixed wettability systems - see ecl kw OPTIONS switch 74
            if (swat <= eps_transf_[cell].wat.smin) {
                swat = eps_transf_[cell].wat.smin;
            } else if (pcow < pc_low_threshold) {
                swat = eps_transf_[cell].wat.smax;
            } else {
                const int wpos = phase_usage_.phase_pos[BlackoilPhases::Aqua];
                const int max_np = BlackoilPhases::MaxNumPhases;
                double s[max_np] = { 0.0 };
                s[wpos] = swat;
                double pc[max_np] = { 0.0 };
                funcForCell(cell).evalPc(s, pc, &(eps_transf_[cell]));
                if (pc[wpos] > pc_low_threshold) {
                    eps_transf_[cell].wat.pcFactor *= pcow/pc[wpos];
                }
            }
        } else {
            OPM_THROW(std::runtime_error, "swatInitScaling: no water phase! ");
        }
    }


    // Map the cell number to the correct function set.
    template <class SatFuncSet>
    const typename SaturationPropsFromDeck<SatFuncSet>::Funcs&
    SaturationPropsFromDeck<SatFuncSet>::funcForCell(const int cell) const
    {
        return cell_to_func_.empty() ? satfuncset_[0] : satfuncset_[cell_to_func_[cell]];
    }

    // Initialize saturation scaling parameters
    template <class SatFuncSet>
    template<class T>
    void SaturationPropsFromDeck<SatFuncSet>::initEPS(Opm::DeckConstPtr deck,
                                                      Opm::EclipseStateConstPtr eclipseState,
                                                      int number_of_cells,
                                                      const int* global_cell,
                                                      const T& begin_cell_centroid,
                                                      int dimensions,
                                                      const std::vector<std::string>& eps_kw,
                                                      std::vector<EPSTransforms>& eps_transf)
    {
        std::vector<std::vector<double> > eps_vec(eps_kw.size());
        const std::vector<double> dummy;
        
        for (size_t i = 0; i < eps_kw.size(); ++i) {
            initEPSKey(deck, eclipseState, number_of_cells, global_cell, begin_cell_centroid, dimensions,
                       eps_kw[i], eps_vec[i]);
        }

        const int wpos = phase_usage_.phase_pos[BlackoilPhases::Aqua];
        const int gpos = phase_usage_.phase_pos[BlackoilPhases::Vapour];
        const bool oilWater = phase_usage_.phase_used[Aqua] && phase_usage_.phase_used[Liquid] && !phase_usage_.phase_used[Vapour];
        const bool oilGas = !phase_usage_.phase_used[Aqua] && phase_usage_.phase_used[Liquid] && phase_usage_.phase_used[Vapour];
        const bool threephase = phase_usage_.phase_used[Aqua] && phase_usage_.phase_used[Liquid] && phase_usage_.phase_used[Vapour];

        for (int cell = 0; cell < number_of_cells; ++cell) {
            if (threephase || oilWater) {
                // ### krw
                initEPSParam(cell, eps_transf[cell].wat, false,
                             funcForCell(cell).smin_[wpos],
                             funcForCell(cell).swcr_,
                             funcForCell(cell).smax_[wpos],
                             funcForCell(cell).sowcr_,
                             oilWater ? -1.0 : funcForCell(cell).smin_[gpos],
                             funcForCell(cell).krwr_,
                             funcForCell(cell).krwmax_,
                             funcForCell(cell).pcwmax_,
                             eps_vec[0], eps_vec[2], eps_vec[1], eps_vec[6], eps_vec[3], eps_vec[11], eps_vec[8], eps_vec[15]);
                // ### krow
                initEPSParam(cell, eps_transf[cell].watoil, true,
                             0.0,
                             funcForCell(cell).sowcr_,
                             funcForCell(cell).smin_[wpos],
                             funcForCell(cell).swcr_,
                             oilWater ? -1.0 : funcForCell(cell).smin_[gpos],
                             funcForCell(cell).krorw_,
                             funcForCell(cell).kromax_,
                             0.0,
                             eps_vec[0], eps_vec[6], eps_vec[0], eps_vec[2], eps_vec[3], eps_vec[13], eps_vec[10], dummy);
            }
            if (threephase || oilGas) {
                // ### krg
                initEPSParam(cell, eps_transf[cell].gas, false,
                             funcForCell(cell).smin_[gpos],
                             funcForCell(cell).sgcr_,
                             funcForCell(cell).smax_[gpos],
                             funcForCell(cell).sogcr_,
                             oilGas ? -1.0 : funcForCell(cell).smin_[wpos],
                             funcForCell(cell).krgr_,
                             funcForCell(cell).krgmax_,
                             funcForCell(cell).pcgmax_,
                             eps_vec[3], eps_vec[5], eps_vec[4], eps_vec[7], eps_vec[0], eps_vec[12], eps_vec[9], eps_vec[16]);
                // ### krog
                initEPSParam(cell, eps_transf[cell].gasoil, true,
                             0.0,
                             funcForCell(cell).sogcr_,
                             funcForCell(cell).smin_[gpos],
                             funcForCell(cell).sgcr_,
                             oilGas ? -1.0 : funcForCell(cell).smin_[wpos],
                             funcForCell(cell).krorg_,
                             funcForCell(cell).kromax_,
                             0.0,
                             eps_vec[3], eps_vec[7], eps_vec[3], eps_vec[5], eps_vec[0], eps_vec[14], eps_vec[10], dummy);
            }
        }
    }

    // Initialize saturation scaling parameter
    template <class SatFuncSet>
    template<class T>
    void SaturationPropsFromDeck<SatFuncSet>::initEPSKey(Opm::DeckConstPtr deck,
                                                         Opm::EclipseStateConstPtr eclipseState,
                                                         int number_of_cells,
                                                         const int* global_cell,
                                                         const T& begin_cell_centroid,
                                                         int dimensions,
                                                         const std::string& keyword,
                                                         std::vector<double>& scaleparam)
    { 
        const bool useAqua = phase_usage_.phase_used[Aqua];
        const bool useLiquid = phase_usage_.phase_used[Liquid];
        const bool useVapour = phase_usage_.phase_used[Vapour];
        bool useKeyword = deck->hasKeyword(keyword);
        bool useStateKeyword = eclipseState->hasDoubleGridProperty(keyword);
        const std::map<std::string, int> kw2tab = {
            {"SWL", 1}, {"SWCR", 2}, {"SWU", 3}, {"SGL", 4},
            {"SGCR", 5}, {"SGU", 6}, {"SOWCR", 7}, {"SOGCR", 8},
            {"ISWL", 1}, {"ISWCR", 2}, {"ISWU", 3}, {"ISGL", 4},
            {"ISGCR", 5}, {"ISGU", 6}, {"ISOWCR", 7}, {"ISOGCR", 8}};
        bool hasENPTVD = deck->hasKeyword("ENPTVD");
        bool hasENKRVD = deck->hasKeyword("ENKRVD");
        int itab = 0;
        std::vector<std::vector<double> > param_col;
        std::vector<std::vector<double> > depth_col;
        std::vector<std::string> col_names;

        // Active keyword assigned default values for each cell (in case of possible box-wise assignment)
        if ((keyword[0] == 'S' && (useStateKeyword || hasENPTVD)) || (keyword[1] == 'S' && useStateKeyword) ) {
            if (useAqua && (useStateKeyword || columnIsMasked_(deck, "ENPTVD", kw2tab.find(keyword)->second-1))) {
                itab = kw2tab.find(keyword)->second;
                scaleparam.resize(number_of_cells);
            }
            if (!useKeyword && itab > 0) {
                const auto& enptvdTables = eclipseState->getEnptvdTables();
                int num_tables = enptvdTables.size();
                param_col.resize(num_tables);
                depth_col.resize(num_tables);
                col_names.resize(9);
                for (int table_num=0; table_num<num_tables; ++table_num) {
                    const auto& enptvdTable = enptvdTables[table_num];
                    depth_col[table_num] = enptvdTable.getDepthColumn();
                    param_col[table_num] = enptvdTable.getColumn(itab); // itab=[1-8]: swl swcr swu sgl sgcr sgu sowcr sogcr
                }
            }
        } else if ((keyword[0] == 'K' && (useKeyword || hasENKRVD)) || (keyword[1] == 'K' && useKeyword) ) {
            if (keyword == std::string("KRW") || keyword == std::string("IKRW") ) {
                if (useAqua && (useKeyword || columnIsMasked_(deck, "ENKRVD", 0))) {
                    itab = 1;
                    scaleparam.resize(number_of_cells);
                    for (int i=0; i<number_of_cells; ++i)
                        scaleparam[i] = funcForCell(i).krwmax_;
                }
            } else if (keyword == std::string("KRG") || keyword == std::string("IKRG") ) {
                if (useVapour && (useKeyword || columnIsMasked_(deck, "ENKRVD", 1))) {
                    itab = 2;
                    scaleparam.resize(number_of_cells);
                    for (int i=0; i<number_of_cells; ++i)
                        scaleparam[i] = funcForCell(i).krgmax_;
                }
            } else if (keyword == std::string("KRO") || keyword == std::string("IKRO") ) {
                if (useLiquid && (useKeyword || columnIsMasked_(deck, "ENKRVD", 2))) {
                    itab = 3;
                    scaleparam.resize(number_of_cells);
                    for (int i=0; i<number_of_cells; ++i)
                        scaleparam[i] = funcForCell(i).kromax_;
                }
            } else if (keyword == std::string("KRWR") || keyword == std::string("IKRWR") ) {
                if (useAqua && (useKeyword || columnIsMasked_(deck, "ENKRVD", 3))) {
                    itab = 4;
                    scaleparam.resize(number_of_cells);
                    for (int i=0; i<number_of_cells; ++i)
                        scaleparam[i] = funcForCell(i).krwr_;
                }
            } else if (keyword == std::string("KRGR") || keyword == std::string("IKRGR") ) {
                if (useVapour && (useKeyword || columnIsMasked_(deck, "ENKRVD", 4))) {
                    itab = 5;
                    scaleparam.resize(number_of_cells);
                    for (int i=0; i<number_of_cells; ++i)
                        scaleparam[i] = funcForCell(i).krgr_;
                }
            } else if (keyword == std::string("KRORW") || keyword == std::string("IKRORW") ) {
                if (useAqua && (useKeyword || columnIsMasked_(deck, "ENKRVD", 5))) {
                    itab = 6;
                    scaleparam.resize(number_of_cells);
                    for (int i=0; i<number_of_cells; ++i)
                        scaleparam[i] = funcForCell(i).krorw_;
                }
            } else if (keyword == std::string("KRORG") || keyword == std::string("IKRORG") ) {
                if (useVapour && (useKeyword || columnIsMasked_(deck, "ENKRVD", 6))) {
                    itab = 7;
                    scaleparam.resize(number_of_cells);
                    for (int i=0; i<number_of_cells; ++i)
                        scaleparam[i] = funcForCell(i).krorg_;
                }
            } else {
                OPM_THROW(std::runtime_error, " -- unknown keyword: '" << keyword << "'");
            }
            if (!useKeyword && itab > 0) {
                const auto& enkrvdTables = eclipseState->getEnkrvdTables();
                int num_tables = enkrvdTables.size();
                param_col.resize(num_tables);
                depth_col.resize(num_tables);
                col_names.resize(8);
                for (int table_num=0; table_num<num_tables; ++table_num) {
                    const auto &enkrvdTable = enkrvdTables[table_num];
                    depth_col[table_num] = enkrvdTable.getDepthColumn();
                    param_col[table_num] = enkrvdTable.getColumn(itab); // itab=[1-7]: krw krg kro krwr krgr krorw krorg
                }
            }
        } else if (useKeyword && (keyword[0] == 'P' || keyword[1] == 'P') ) {
             if (useAqua && (keyword == std::string("PCW") || keyword == std::string("IPCW")) ) {
                 scaleparam.resize(number_of_cells);
                 for (int i=0; i<number_of_cells; ++i)
                     scaleparam[i] = funcForCell(i).pcwmax_;
             } else if (useVapour && (keyword == std::string("PCG") || keyword == std::string("IPCG")) ) {
                 scaleparam.resize(number_of_cells);
                 for (int i=0; i<number_of_cells; ++i)
                     scaleparam[i] = funcForCell(i).pcgmax_;
             }
        }

        if (scaleparam.empty()) {
            return;
        }

        if (useKeyword || useStateKeyword) {
            // Keyword values from deck
            std::cout << "--- Scaling parameter '" << keyword << "' assigned." << std::endl;
            const int* gc = global_cell;
            std::vector<double> val;
            if (keyword[0] == 'S' || keyword[1] == 'S') { // Saturation from EclipseState
                val = eclipseState->getDoubleGridProperty(keyword)->getData();
            } else {
                val = deck->getKeyword(keyword)->getSIDoubleData(); //KR and PC directly from deck.
            }
            for (int c = 0; c < int(scaleparam.size()); ++c) {
                const int deck_pos = (gc == NULL) ? c : gc[c];
                scaleparam[c] = val[deck_pos];
            }
        }

        if (itab > 0) {
            const int dim = dimensions;
            std::vector<int> endnum;
            if ( deck->hasKeyword("ENDNUM")) {
                const std::vector<int>& e = 
                    deck->getKeyword("ENDNUM")->getIntData();              
                endnum.resize(number_of_cells);                                   
                const int* gc = global_cell;
                for (int cell = 0; cell < number_of_cells; ++cell) {
                    const int deck_pos = (gc == NULL) ? cell : gc[cell];
                    endnum[cell] = e[deck_pos] - 1; // Deck value zero prevents scaling via ENPTVD/ENKRVD
                }
            }
            else {
                // Default deck value is one
                endnum.assign(number_of_cells, 0);
            }
            if (keyword[0] == 'S' || keyword[1] == 'S') { // From EclipseState
                for (int cell = 0; cell < number_of_cells; ++cell) {
                    if (!std::isfinite(scaleparam[cell]) && endnum[cell] >= 0 && param_col[endnum[cell]][0] >= 0.0) {
                        double zc = UgGridHelpers
                            ::getCoordinate(UgGridHelpers::increment(begin_cell_centroid, cell, dim),
                                           dim-1);
                        if (zc >= depth_col[endnum[cell]].front() && zc <= depth_col[endnum[cell]].back()) { //don't want extrap outside depth interval
                            scaleparam[cell] = linearInterpolation(depth_col[endnum[cell]], param_col[endnum[cell]], zc);
                        }
                    } else if (!std::isfinite(scaleparam[cell]) && endnum[cell] >= 0) {
                        // As of 1/9-2014:  Reflects remaining work on opm/parser/eclipse/EclipseState/Grid/GridPropertyInitializers.hpp ...
                        OPM_THROW(std::runtime_error, " -- Inconsistent EclipseState: '" << keyword << "' (ENPTVD)");
                    }
                }
            } else { //KR and PC from deck.
                for (int cell = 0; cell < number_of_cells; ++cell) {
                    if (endnum[cell] >= 0 && param_col[endnum[cell]][0] >= 0.0) {
                        double zc = UgGridHelpers
                            ::getCoordinate(UgGridHelpers::increment(begin_cell_centroid, cell, dim),
                                           dim-1);
                        if (zc >= depth_col[endnum[cell]].front() && zc <= depth_col[endnum[cell]].back()) { //don't want extrap outside depth interval
                            scaleparam[cell] = linearInterpolation(depth_col[endnum[cell]], param_col[endnum[cell]], zc);
                        }
                    }
                }
            }
        }

//        std::cout << keyword << ":" << std::endl;
//        for (int c = 0; c < int(scaleparam.size()); ++c) {
//                std::cout << c << "    " << scaleparam[c] << std::endl;
//        }

    }

    // Saturation scaling
    template <class SatFuncSet>
    void SaturationPropsFromDeck<SatFuncSet>::initEPSParam(const int cell,
                                                           EPSTransforms::Transform& data,
                                                           const bool oil,          // flag indicating krow/krog calculations
                                                           const double sl_tab,     // minimum saturation (for krow/krog calculations this is normally zero)
                                                           const double scr_tab,    // critical saturation
                                                           const double su_tab,     // maximum saturation (for krow/krog calculations this is minimum water/gas saturation)
                                                           const double sxcr_tab,   // second critical saturation (not used for 2pt scaling)
                                                           const double s0_tab,     // threephase complementary minimum saturation (-1.0 indicates 2-phase)
                                                           const double krsr_tab,   // relperm at displacing critical saturation
                                                           const double krmax_tab,  // relperm at maximum saturation
                                                           const double pcmax_tab,  // cap-pres at maximum saturation (zero => no scaling)
                                                           const std::vector<double>& sl,  // For krow/krog calculations this is not used
                                                           const std::vector<double>& scr,
                                                           const std::vector<double>& su,  // For krow/krog calculations this is SWL/SGL
                                                           const std::vector<double>& sxcr,
                                                           const std::vector<double>& s0,
                                                           const std::vector<double>& krsr,
                                                           const std::vector<double>& krmax,
                                                           const std::vector<double>& pcmax) // For krow/krog calculations this is not used
    {
        if (scr.empty() && su.empty() && (sxcr.empty() || !do_3pt_) && s0.empty()) {
            data.doNotScale = true;
            data.smin = sl_tab;
            if (oil) {
                data.smax = (s0_tab < 0.0) ? 1.0 - su_tab : 1.0 - su_tab - s0_tab;
            } else {
                data.smax = su_tab;
            }
            data.scr = scr_tab;
        } else {
            data.doNotScale = false;
            data.do_3pt = do_3pt_;
            double s_r;
            if (s0_tab < 0.0) { // 2phase
                s_r = 1.0-sxcr_tab;
                if (do_3pt_) data.sr = sxcr.empty() ? s_r : 1.0-sxcr[cell];
            } else { // 3phase
                s_r = 1.0-sxcr_tab-s0_tab;
                if (do_3pt_)data.sr = 1.0 - (sxcr.empty() ? sxcr_tab : sxcr[cell])
                                          - (s0.empty() ? s0_tab : s0[cell]);
            }
            data.scr = scr.empty() ? scr_tab : scr[cell];
            double s_max = su_tab;
            if (oil) {
                data.smin = sl_tab;
                if (s0_tab < 0.0) { // 2phase
                    s_max = 1.0 - su_tab;
                    data.smax = 1.0 - (su.empty() ? su_tab : su[cell]);
                } else { // 3phase
                    s_max = 1.0 - su_tab - s0_tab;
                    data.smax = 1.0 - (su.empty() ? su_tab : su[cell])
                                    - (s0.empty() ? s0_tab : s0[cell]);
                }
            } else {
                data.smin = sl.empty() ? sl_tab : sl[cell];
                data.smax = su.empty() ? su_tab : su[cell];
            }
            if (do_3pt_) {
                data.slope1 = (s_r-scr_tab)/(data.sr-data.scr);
                data.slope2 = (s_max-s_r)/(data.smax-data.sr);
            } else {
                data.slope2 = data.slope1 = (s_max-scr_tab)/(data.smax-data.scr);
                // Inv transform of tabulated critical displacing saturation to prepare for possible value scaling (krwr etc)
                data.sr = data.scr + (s_r-scr_tab)*(data.smax-data.scr)/(s_max-scr_tab);
            }
        }
        
        data.doKrMax = !krmax.empty();
        data.doKrCrit = !krsr.empty();
        data.doSatInterp = false;
        data.krsr = krsr.empty() ? krsr_tab : krsr[cell];
        data.krmax = krmax.empty() ? krmax_tab : krmax[cell];
        data.krSlopeCrit = data.krsr/krsr_tab;
        data.krSlopeMax = data.krmax/krmax_tab;
        if (data.doKrCrit) {
            if (data.sr > data.smax-1.0e-6) {
                //Ignore krsr and do two-point (one might consider combining krsr and krmax linearly between scr and smax ... )
                data.doKrCrit = false;
            } else if (std::fabs(krmax_tab- krsr_tab) > 1.0e-6) { // interpolate kr
                data.krSlopeMax = (data.krmax-data.krsr)/(krmax_tab-krsr_tab);
            } else { // interpolate sat
                data.doSatInterp = true;
                data.krSlopeMax = (data.krmax-data.krsr)/(data.smax-data.sr);
            }
        }

        if (std::fabs(pcmax_tab) < 1.0e-8 || pcmax.empty() || pcmax_tab*pcmax[cell] < 0.0) {
            data.pcFactor = 1.0;
        } else {
            data.pcFactor = pcmax[cell]/pcmax_tab;
        }

    }

} // namespace Opm

#endif // OPM_SATURATIONPROPSFROMDECK_IMPL_HEADER_INCLUDED
