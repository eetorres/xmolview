//========================================================================
// FILE: csupercell.h
//
// This an utility program to manipulate and genarate structure files
//
// Copyright 2006-2015 by Edmanuel Torres
// email:   eetorres@gmail.com
//
// Lastest update: 13/07/2015
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

#ifndef _CFRAGMOL_H_
#define _CFRAGMOL_H_

#include<config_debug.h>
#include<cfragment.h>
#include<cfile.h>

class CSupercell {

public:

  CSupercell();
  ~CSupercell();
  //
  CFile gsf;
  //
  void clear(void);
  bool read_input_file(void);                  // read the input file
  void initialize_fragments(void);             //
  //
  /////////////////////////////////////////////////////////////////////////////////
  // Save functions
  /////////////////////////////////////////////////////////////////////////////////
  void save_input_file(void);                  // save the file
  void save_as_file(std::string,bool,bool,TMatrix<real> m,int,int,int,int,bool);
  void save_as_file(std::string,std::string,bool,bool,TMatrix<real> m,int,int,int,int,bool);
  void save_topmol_file(std::string,std::string);

  /////////////////////////////////////////////////////////////////////////////////
  // Compute functions
  /////////////////////////////////////////////////////////////////////////////////
  void compute_fragmol_cartesian(uint);
  void compute_fragmol_direct(uint);
  void compute_fragmol_all_cartesian(void);
  void compute_fragmol_all_direct(void);
  void compute_fragmol_direct(void);
  void compute_fragmol_position_cartesian(void);
  //
  /////////////////////////////////////////////////////////////////////////////////
  // Evaluation functions
  /////////////////////////////////////////////////////////////////////////////////
  void eval_fragmol_initial_position(void);
  void eval_fragmol_initial_orientation(void);
  void eval_cell_table(void);
  void eval_connections(const TMatrix<uint>&,uint);
  bool eval_new_fragment(const TVector<uint>&);               // Create all posible scaled framgent
  bool eval_scaled_fragment(uint,bool,real);      // Create a framgent using scaled distance
  bool eval_radial_fragment(uint,bool,real);      // Create a framgent using scaled distance
  bool eval_scaled_fragment(uint,real);           // Create single scaled distance fragment
  void eval_scaled_fragments(real);               // Create all posible scaled framgent
  // Begin deprecated
  void eval_vdw_fragments(void);               // Create all vdW framgents
  void eval_atom_fragments(void);               // Create all atom framgents
  //bool eval_atom_fragment(uint);               // Create one atom framgent
  //bool eval_vdw_fragment(uint);                // Create one vdW framgent
  // End deprecated
  //bool eval_vdw_fragment(uint,bool,real);      // Create one framgent using scaling
  bool eval_merge_fragment(uint,bool);         // Merge two fragments
  //
  bool delete_atom(uint u);                    // Delete an atom
  //
  void update_fragmol_cartesian(void);
  void update_fragmol_direct(void);
  //
  void is_fragmol_initialized(bool);
  void apply_pbc(bool);

  /////////////////////////////////////////////////////////////////////////////////
  // Set functions
  /////////////////////////////////////////////////////////////////////////////////
  void set_cartesian(void);
  //
  void set_gsf_modified(const bool);
  //
  void set_fragmol_active_fragment(const uint);
  void set_input_file_type(const uint);
  void set_input_file_units(const uint);
  void set_output_file_type(const uint);
  void set_output_file_format(const uint);
  void set_export_format(const uint);
  //
  void set_fragmol_fragment_twist(const real);
  void set_fragmol_fragment_precession(const real);
  void set_fragmol_fragment_tilt(const real);
  //
  void set_fragmol_fragment_position_u(const real);
  void set_fragmol_fragment_position_v(const real);
  void set_fragmol_fragment_position_w(const real);
  //
  void set_input_file(std::string);
  void set_dir(std::string);
  void set_topmol_directory(std::string);
  //
  void set_fragment_axis(const TVector<uint>&);
  void set_fragment_axis(const uint, const uint);
  void set_fragment_direct(const uint,const TMatrix<real>&);
  void set_fragment_cartesian(const uint,const TMatrix<real>&);
  //
  /////////////////////////////////////////////////////////////////////////////////
  // Is  functions
  /////////////////////////////////////////////////////////////////////////////////
  bool is_direct(void);
  bool is_periodic(void);
  void is_periodic(bool);
  bool is_potmol(void);
  bool is_fragmol_initialized();

  /////////////////////////////////////////////////////////////////////////////////
  //  Get  functions
  /////////////////////////////////////////////////////////////////////////////////
  std::string get_dir(void);
  std::string get_fragment_atomic_label(uint,uint);
  std::string get_fragment_atomic_symbol(uint,uint);
  //
  uint get_fragmol_number_of_fragments(void);
  uint get_fragmol_active_fragment(void);
  uint get_fragment_size(uint);
  uint get_fragment_atomic_number(uint,uint);
  //
  uint get_fragmol_total_atoms(void);
  uint get_species(void);
  uint get_atomic_composition(uint);
  uint get_total_atoms(void);
  uint get_input_file_format(void);
  uint get_output_file_type(void);
  uint get_output_file_format(void);
  //
  real get_fragmol_axis_tilt(void);
  real get_fragmol_axis_precession(void);
  real get_fragmol_backbone_tilt(void);
  real get_fragmol_backbone_precession(void);
  //
  TVector<std::string> get_fragmol_atomic_symbol_table(void);
  //
  TVector<uint> get_fragmol_atomic_composition_table(void);
  TVector<uint> get_fragmol_atomic_number_table(void);
  TVector<uint> get_fragmol_atom_table(void);
  TVector<uint> get_fragmol_fragment_table(void);
  //
  TVector<real> get_fragmol_axis_angles(void);
  TVector<real> get_fragmol_basis_direct(void);
  TVector<real> get_fragmol_position_direct(void);
  TVector<real> get_fragmol_position_uvw(void);
  TVector<real> get_fragmol_position_cartesian(void);
  TVector<real> get_fragmol_centered_position_cartesian(void);
  //
  TVector<real> get_cartesian(uint);
  TVector<real> get_direct(uint);
  TVector<real> get_vector_diff(uint,uint);
  //
  real get_distance(uint,uint);
  real get_angle(uint,uint,uint);
  real get_dihedral(uint,uint,uint,uint);
  //
  TVector<real> get_fragment_direct(uint,uint);
  TVector<real> get_fragment_cartesian(uint,uint);
  TVector<real> get_fragment_centered_cartesian(uint,uint);
  TMatrix<real> get_cartesian(void);
  TMatrix<real> get_direct(void);
  TMatrix<real> get_fragment_direct(uint);
  TMatrix<real> get_fragment_cartesian(uint);
  TMatrix<real> get_fragment_centered_cartesian(uint);
  TMatrix<real> get_uvw_to_xyz(void);
  TMatrix<real> get_unit_uvw_to_xyz(void);
  TMatrix<real> get_bounding_box(void);


private:

  bool __is_direct;
  bool __is_periodic;
  bool __is_potmol;
  bool __is_input_fragment;
  //
  uint __input_format;
  uint __output_format;
  uint __total_atoms;
  uint __atomic_species;
  uint __number_of_fragments;
  uint __active_fragment;
  //
  std::ifstream iposmol;
  std::ifstream itopmol;
  std::ofstream oposmol;
  std::ofstream otopmol;
  //
  TVector<CFragment> v_fragments;
  TVector<uint> v_atomic_composition_table;
  TVector<uint> v_atomic_number_table;
  TVector<uint> v_atomic_numbers;
  TVector<uint> v_fragment_table;
  TVector<uint> v_atom_type_table;
  TVector<uint> v_atom_cell_table;
  TVector<std::string> v_atomic_labels;
  TVector<std::string> v_atomic_symbols;
  TVector<std::string> v_atomic_symbol_table;
  //
  TMatrix<real> m_xyz;
  TMatrix<real> m_uvw;
  TMatrix<real> m_uvw_to_xyz;
  TMatrix<real> m_uvw_to_xyz_u;
  //
  std::string inputfile, output_filename, potmolfile, __sdir;
  //
  void clear_fragments(void);
  void create_initial_fragments(void);     //
  void eval_initial_fragments(void);               //
  //
};

#endif

