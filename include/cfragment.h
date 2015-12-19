//========================================================================
// FILE: cfragment.h
//
// This an utility program to manipulate and generate structure files
//
// Copyright 2011-2015 by Edmanuel Torres
// email:   eetorres@gmail.com
//
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

#ifndef _CFRAGMENT_H_
#define _CFRAGMENT_H_

#include<config_debug.h>
#include<string.h>
#include<msmvtl/const.h>
#include<msmvtl/tmatrix.h>
#include<msmvtl/tvmath.h>
#include<msmvtl/linalg.h>
#include<atom_color.h>


class CFragment {

public:
  CFragment();
  ~CFragment();
  //
  void copy(CFragment&, TVector<uint>&);
  void move(CFragment&, TVector<uint>&);
  void size(const uint);
  uint size(void);
  void clear(void);
  //
  void is_initialized(bool);
  void is_pbc(bool);
  //
  void eval_integrity(void);
  // Deprecated
  //void eval_integrity(uint,real);
  void eval_scaled_bond_integrity(uint,real);
  void eval_radial_integrity(uint,real);
  void eval_initial_position(void);
  void eval_initial_orientation(void);
  void eval_position(void);
  void eval_orientation(void);
  //
  void compute_origin_cartesian(void);
  void compute_origin_direct(void);
  void compute_origin_direct_from_direct(void);
  //
  TVector<uint> compute_vdw_fragment(uint,real);
  TVector<uint> compute_radial_fragment(uint,real);
  //
  void set_atomic_label(uint,std::string);
  void set_atomic_symbol(uint,std::string);
  void set_atomic_number(uint,uint);
  //
  void set_backbone_precession(real);
  void set_backbone_tilt(real);
  //
  void set_axis_twist(const real);
  void set_axis_precession(const real);
  void set_axis_tilt(const real);
  //
  void set_origin_u(const real);
  void set_origin_v(const real);
  void set_origin_w(const real);
  //
  void set_axis_index(const TVector<uint>&);
  void set_axis_index(const uint, const uint);
  void set_atomic_label(const TVector<std::string>&);
  void set_atomic_symbol(const TVector<std::string>&);
  void set_atomic_number(const TVector<uint>&);

  void set_direct(uint,const TVector<real>&);
  void set_cartesian(uint,const TVector<real>&);
  void set_position_direct(uint,const TVector<real>&);
  void set_position_cartesian(uint,const TVector<real>&);
  //
  void set_uvw(const TMatrix<real>&);
  void set_direct(const TMatrix<real>&);
  void set_cartesian(const TMatrix<real>&);
  void set_position_direct(const TMatrix<real>&);
  void set_position_cartesian(const TMatrix<real>&);
  void set_centered_direct(const TMatrix<real>&);
  void set_centered_cartesian(const TMatrix<real>&);
  void set_rotation_cartesian(const TMatrix<real>&);
  void set_uvw_to_xyz(const TMatrix<real>&);
  void set_unit_uvw_to_xyz(const TMatrix<real>&);
  //
  bool is_initialized(void);
  std::string get_atomic_label(uint);
  std::string get_atomic_symbol(uint);
  uint get_axis_origin_index(void);
  uint get_atomic_number(uint);
  //
  real get_axis_tilt(void);
  real get_axis_precession(void);
  real get_axis_twist(void);
  real get_backbone_tilt(void);
  real get_backbone_precession(void);
  //
  void show_information(void);
  //
  TVector<real> get_axis_angles(void);
  TVector<real> get_basis_direct(void);
  TVector<real> get_origin_direct(void);
  TVector<real> get_origin_uvw(void);
  TVector<real> get_origin_cartesian(void);
  TVector<real> get_centered_origin_cartesian(void);

  TVector<real> get_direct(uint);
  TVector<real> get_cartesian(uint);
  TVector<real> get_centered_cartesian(uint);

  TMatrix<real> get_direct(void);
  TMatrix<real> get_cartesian(void);
  TMatrix<real> get_centered_cartesian(void);
  //
  // Utilities
  TVector<real> pbc(const TVector<real>&);
  real pbc(const real);
  TVector<std::string> v_atomic_label;
  TVector<std::string> v_atomic_symbol;

private:

  bool __is_show;
  bool __is_initialized;
  bool __is_centered;
  static bool __is_pbc;
  uint __size;
  //
  TVector<uint> v_atomic_number;
  TVector<uint> v_axis_index;
  TVector<real> v_axis_angle;
  TVector<real> v_initial_axis_angle;
  TVector<real> v_backbone_angle;
  TVector<real> v_basis_direct;
  TVector<real> v_pbc_box;
  TVector<real> v_cartesian_axes;
  //
  TVector<real> v_origin_uvw;
  TVector<real> v_origin_direct;
  TVector<real> v_origin_cartesian;
  TVector<real> v_initial_origin_cartesian;
  //
  TMatrix<real> m_uvwTxyz;
  TMatrix<real> m_xyzTuvw;
  TMatrix<real> m_unit_uvwTxyz;
  TMatrix<real> m_unit_xyzTuvw;
  //
  TMatrix<real> m_uvw;
  TMatrix<real> m_direct;
  TMatrix<real> m_cartesian;
  TMatrix<real> m_position_direct;
  TMatrix<real> m_position_cartesian;
  TMatrix<real> m_centered_direct;
  TMatrix<real> m_centered_cartesian;
  TMatrix<real> m_rotation_cartesian;
  // generalized coordinates
  //
};

#endif

