//========================================================================
// FILE: cviewmol.h
//
// This class is a layer to separate core structure manipulations from
// graphic output and GUI controls.
//
// Copyright 2011-2015 by Edmanuel Torres
// email:   eetorres@gmail.com
//
// Lastest update: 10/07/2015
//========================================================================
//  This file is part of xmolview                                       //
//                                                                      //
//  xmolview is free software: you can redistribute it and/or modify    //
//  it under the terms of the GNU General Public License as published by//
//  the Free Software Foundation, either version 3 of the License, or   //
//  (at your option) any later version.                                 //
//                                                                      //
//  xmolview is distributed in the hope that it will be useful,         //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of      //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       //
//  GNU General Public License for more details.                        //
//                                                                      //
//  You should have received a copy of the GNU General Public License   //
//  along with Foobar.  If not, see <http://www.gnu.org/licenses/>.     //
//======================================================================//

#ifndef _CVIEWCAR_H_
#define _CVIEWCAR_H_

#include<cmapcell.h>

class CViewmol : public CMapcell {

public:

  CViewmol(){};
  ~CViewmol(){};

  bool read_files_view(void);
  void initialize_view(void);
  void update_view(void);
  //
  void set_view_active_fragment(const unsigned int);
  void set_view_bounding_box(void);
  //
  void set_view_translation_map_1d(void);
  void set_view_translation_map_2d(void);
  void set_view_twist_map(void);
  void set_view_tilt_map(void);
  void set_view_precession_map(void);
  //
  void set_view_translation_steps_1d(const unsigned int);
  void set_view_translation_steps_2d(const unsigned int);
  void set_view_twist_steps(const unsigned int);
  void set_view_tilt_steps(const unsigned int);
  void set_view_precession_steps(const unsigned int);
  //
  void set_view_axis_twist(const real);
  void set_view_axis_precession(const real);
  void set_view_axis_tilt(const real);
  //
  void set_view_position_u(const real);
  void set_view_position_v(const real);
  void set_view_position_w(const real);
  //
  void set_view_twist1(const real);
  void set_view_twist2(const real);
  void set_view_tilt1(const real);
  void set_view_tilt2(const real);
  void set_view_precession1(const real);
  void set_view_precession2(const real);
  //
  void set_view_position_uvw(const TVector<real>&);
  void set_view_position_uvw1(const TVector<real>&);
  void set_view_position_uvw2(const TVector<real>&);
  void set_view_position_uvw3(const TVector<real>&);
  void set_view_position_rpt(const TVector<real>&);
  void set_view_position_rpt1(const TVector<real>&);
  void set_view_position_rpt2(const TVector<real>&);
  //
  uint get_view_total_atoms(void);
  uint get_view_total_fragments(void);
  uint get_view_active_fragment(void);
  uint get_view_tilt_axis_begin(unsigned int);
  uint get_view_tilt_axis_end(unsigned int);
  uint get_view_file_type(std::string);
  //
  real get_view_axis_precession(void);
  real get_view_axis_tilt(void);
  real get_view_backbone_precession(void);
  real get_view_backbone_tilt(void);
  //
  TVector<std::string> get_view_atomic_labels(void);
  TVector<std::string> get_view_atomic_symbols(void);
  TVector<std::string> get_view_atomic_symbol_table(void);
  //
  TVector<uint> get_view_atom_table(void);
  TVector<uint> get_view_atomic_composition_table(void);
  TVector<uint> get_view_atomic_number_table(void);
  TVector<uint> get_view_fragment_table(void);
  TVector<uint> get_view_atomic_numbers(void);
  //
  TVector<real> get_view_axis_angles(void);
  TVector<real> get_view_basis_direct(void);
  TVector<real> get_view_position_direct(void);
  TVector<real> get_view_position_uvw(void);
  TVector<real> get_view_position_cartesian(void);
  TVector<real> get_view_position_spheric(void);
  TVector<real> get_view_centered_position_cartesian(void);
  TVector<real> get_view_axis_cartesian_origin(void);
  //
  TMatrix<real> get_view_bounding_box(void);
  TMatrix<real> get_view_cartesian(void);
  TMatrix<real> get_view_direct(void);
  // Conditionals
  bool is_view_periodic(void){ return supercell.is_periodic();};
  //
  // pure virtual functions
  //virtual void view_redraw(void)=0;
  //virtual void view_update_coordinates(void)=0;

private:

  uint view_number_of_atoms;
  uint view_file_format;
  TMatrix<real> m_view_bbox;

};

#endif


