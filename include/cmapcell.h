//========================================================================
// FILE: CMapcell.h
//
// This an utility program to manipulate and genarate structure files
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

#ifndef _CMAPMOL_H_
#define _CMAPMOL_H_

#include<csupercell.h>

class CMapcell {

public:

  CMapcell();
  ~CMapcell(){};
  //
  CSupercell supercell;
  //
  void initialize_map(void);
  //
  void eval_map_cartesian(void);
  void eval_map_direct(void);
  void eval_map_properties(void);
  //
  //
  //void initialize_map_translation_1d(void);
  //void initialize_map_translation_2d(void);
  //void initialize_twist_map(void);
  //void initialize_tilt_map(void);
  //void initialize_precession_map(void);
  //
  //void compute_translation_map(void);
  //void compute_twist_map(void);
  //void compute_tilt_map(void);
  //void compute_precession_map(void);
  //void compute_map(void);
  void generate_map(void);
  //
  void map_update_active_fragment(void);
  //
  void set_translation_last_point(const bool);
  void set_map_use_cartesian(const bool);
  void set_preview(const bool);
  void set_translation_map(const bool);
  void set_twist_map(const bool);
  void set_tilt_map(const bool);
  void set_precession_map(const bool);
  //
  void set_map_active_fragment(const unsigned int);
  void set_map_translation_steps_1d(const unsigned int);
  void set_map_translation_steps_2d(const unsigned int);
  void set_map_twist_steps(const unsigned int);
  void set_map_tilt_steps(const unsigned int);
  void set_map_precession_steps(const unsigned int);
  //
  void set_map_fragment_twist(const real);
  void set_map_fragment_precession(const real);
  void set_map_fragment_tilt(const real);
  //
  void set_map_fragment_position_u(const real);
  void set_map_fragment_position_v(const real);
  void set_map_fragment_position_w(const real);
  //
  void set_map_twist1(const real);
  void set_map_twist2(const real);
  void set_map_tilt1(const real);
  void set_map_tilt2(const real);
  void set_map_precession1(const real);
  void set_map_precession2(const real);
  //
  void set_map_fragment_position_uvw(const TVector<real>&);
  void set_map_position_uvw1(const TVector<real>&);
  void set_map_position_uvw2(const TVector<real>&);
  void set_map_position_uvw3(const TVector<real>&);
  void set_map_position_rpt1(const TVector<real>&);
  void set_map_position_rpt2(const TVector<real>&);
  //
  void set_scan_title(const std::string);
  //
  //void set_write_xyz_file(const std::string);
  //void set_write_uvw_file(const std::string);
  //
  unsigned int get_map_total_atoms(void);
  //
  TVector<real> get_map_position_uvw1(uint);
  TVector<real> get_map_position_uvw2(uint);
  //
  real get_translation_distance(void);
  real get_translation_distance_2d(void);
  real get_translation_step(void);
  real get_translation_step_2d(void);
  real get_twist_step(void);
  real get_tilt_step(void);
  real get_precession_step(void);
  //
  real get_map_axis_precession(void);
  real get_map_axis_tilt(void);
  real get_map_backbone_precession(void);
  real get_map_backbone_tilt(void);
  //
  TVector<std::string> get_map_atomic_labels(void);
  TVector<std::string> get_map_atomic_symbols(void);
  TVector<unsigned int> get_map_atomic_numbers(void);
  TVector<unsigned int> get_map_atomic_table(void);
  TVector<real> get_map_axis_angles(void);
  TVector<real> get_map_basis_direct(void);
  TVector<real> get_map_position_direct(void);
  TVector<real> get_map_position_uvw(void);
  //
  TVector<real> get_map_position_cartesian(void);
  TVector<real> get_map_centered_position_cartesian(void);
  TVector<real> get_map_axis_direct_origin(void);
  TVector<real> get_map_axis_cartesian_origin(void);
  //
  //TMatrix<real> get_map_bounding_box(void);
  TMatrix<real> get_map_cartesian(void);
  TMatrix<real> get_map_direct(void);
  //
  // pure virtual functions
  virtual void view_redraw(void)=0;
  virtual void update_data(void)=0;
  virtual void set_update_coordinates(bool)=0;
  virtual void save_wysiwyg_as(std::string)=0;
  virtual void save_wysiwyg_extension(std::string)=0;
  virtual void save_wysiwyg_as(std::string,std::string)=0;

private:

  bool __if_preview;
  bool __if_translation_map;
  bool __if_translation_map_1d;
  bool __if_translation_map_2d;
  bool __if_twist_map;
  bool __if_tilt_map;
  bool __if_precession_map;
  bool __if_translation_last_point;
  bool __if_cartesian_spheric;
  //
  uint __active_fragment;
  uint __translation_map_steps_1d;
  uint __translation_map_steps_2d;
  uint __twist_map_steps;
  uint __tilt_map_steps;
  uint __precession_map_steps;
  //
  real __twist_map_step, __angle_map_twist1, __angle_map_twist2;
  real __tilt_map_step, __angle_map_tilt1, __angle_map_tilt2;
  real __precession_map_step, __theta_step;
  real __angle_map_precession1, __angle_map_precession2;
  real __angle_theta1, __angle_theta2;
  real __translation_step_1d, __translation_step_2d;
  real __translation_distance_1d, __translation_distance_2d;
  real __inv_stp_1d, __inv_stp_2d;
  //
  std::string s_scan_directory;
  std::string __scan_path;
  //
  TVector<std::string>  v_atomic_labels;
  TVector<std::string>  v_atomic_symbols;
  TVector<unsigned int> v_atomic_numbers;
  TVector<unsigned int> v_atomic_table;

  TVector<real> v_position_map_uvw1;
  TVector<real> v_position_map_uvw2;
  TVector<real> v_position_map_uvw3;
  TVector<real> v_position_map_uvw4;
  TVector<real> v_position_map_step_1d;
  TVector<real> v_position_map_step_2d;

  TMatrix<real> m_cartesian;
  TMatrix<real> m_direct;
  //
  void compute_axis_angles(void);
};

#endif

