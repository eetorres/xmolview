//========================================================================
// FILE: cread.h
//
// This an utility program to manipulate and genarate structure files
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

#ifndef _CFILE_H_
#define _CFILE_H_

#include <config_debug.h>
#include <global.h>
#include <ctopmol.h>
#include <cpotcar.h>
#include <cposcar.h>
#include <cxyz.h>
#include <cgau.h>
#include <czmat.h>
#include <cpdb.h>
#include <cdlp.h>
#include <czmat.h>

class CFile : public CTopmol{

public:

  CFile();
  ~CFile(){};
  //
  void clear(void);
  void clear_xyz(void);
  void clear_gau(void);
  void clear_pdb(void);
  void clear_dlp(void);
  void clear_poscar(void);
  void clear_potcar(void);
  //
  bool read_input_file(void);
  bool read_poscar(void);
  bool read_xyz(void);
  bool read_gau(void);
  bool read_pdb(void);
  bool read_dlp(void);
  bool read_zmt(void);
  bool read_potcar(void);
  //
  void save_input_file(void);
  void save_xyz(void);
  void save_as_file(std::string);
  //void save_as_file(std::string,std::string);
  //void save_poscar_as(std::string);
  //void save_xyz_as(std::string);
  //void save_gau_as(std::string);
  //void save_pdb_as(std::string);
  //void save_dlp_as(std::string);
  void save_as_file(std::string,std::string);
  void save_poscar_as(std::string,std::string,uint);
  void save_xyz_as(std::string,std::string,uint);
  void save_gau_as(std::string,std::string,uint);
  void save_pdb_as(std::string,std::string,uint);
  void save_dlp_as(std::string,std::string,uint);
  void save_zmt_as(std::string,std::string,uint);
  //
  void eval_connections(const TMatrix<uint>&,uint);
  //
  void set_input_file(std::string s);
  void set_output_file_name(std::string s);
  void set_potcar_file(std::string s);
  void set_input_dir(std::string s);
  //
  void set_input_type(uint u);
  void set_input_units(uint u);
  void set_input_format(uint u);
  //
  void set_output_file_type(uint u);
  void set_output_file_format(uint u);
  void set_export_format(uint u);
  //
  void set_fragment_table(TVector<uint> v,uint u);
  //
  void set_bounding_box(bool b);
  void set_labels(bool b);
  void set_numbers(bool b);
  void set_fragments(bool b);
  void set_modified(bool b);
  //
  void set_xyz(TMatrix<real> m);
  void set_xyz_cells(int x, int y, int z, int t);
  //void set_input_file(std::string inputfile, potcarfile, s_actual_dir;
  uint get_total_atoms(void){ return u_total_atoms;}
  uint get_atomic_species(void){ return u_atomic_species;}
  uint get_input_type(void){ return u_input_file_type;}
  uint get_input_units(void){ return u_input_file_units;}
  uint get_input_format(void){ return u_input_format;}
  uint get_output_file_format(void){ return u_output_file_format;}
  uint get_output_file_type(void){ return u_output_file_type;}
  //
  TVector<uint> get_atomic_composition_table(void){ return v_atomic_composition_table;}
  TVector<uint> get_atomic_numbers(void){ return v_atomic_numbers;}
  TVector<uint> get_fragment_table(void){ return v_fragment_table;}
  TVector<uint> get_atom_table(void){ return v_atom_table;}
  TVector<uint> get_atomic_number_table(void){ return v_atomic_number_table;}
  //
  TVector<std::string> get_atomic_labels(void){ return v_atomic_labels;}
  TVector<std::string> get_atomic_symbols(void){ return v_atomic_symbols;}
  TVector<std::string> get_atomic_symbol_table(void){ return v_atomic_symbol_table;}
  //
  bool get_is_direct(void){ return is_direct;}
  bool get_is_periodic(void){ return is_periodic;}
  //
  TMatrix<real> get_xyz(void){ return m_xyz;}
  TMatrix<real> get_uvw(void){ return m_uvw;}
  TMatrix<real> get_uvw_to_xyz_u(void){ return m_uvw_to_xyz_u;}
  TMatrix<real> get_uvw_to_xyz(void){ return m_uvw_to_xyz;}

private:

  bool is_direct;
  bool is_periodic;
  bool is_potcar;
  //bool is_input_fragment;
  bool is_bounding_box;
  bool is_labels;
  bool is_numbers;
  bool is_fragments;
  bool is_charges;
  bool is_modified;
  //
  int x_cells;
  int y_cells;
  int z_cells;
  int total_cells;
  //
  uint u_input_format;
  uint u_input_file_type;
  uint u_input_file_units;
  uint u_output_file_type;
  uint u_output_file_format;
  uint u_export_format;
  uint u_total_atoms;
  uint u_total_fragments;
  uint u_atomic_species;

  std::string inputfile;
  std::string output_filename;
  std::string potcarfile;
  std::string s_actual_dir;

  CPotcar file_potcar;
  CPoscar file_poscar;
  CXyz    file_xyz;
  CGau    file_gau;
  CZmat   file_zmt;
  CPdb    file_pdb;
  CDlp    file_dlp;

  //bool is_potcar(void);
  TVector<std::string> v_atomic_symbol_table;
  TVector<uint> v_atomic_composition_table;
  TVector<uint> v_atomic_number_table;
  TVector<uint> v_fragment_table;
  TVector<uint> v_atom_table;
  TVector<uint> v_atomic_numbers;
  TVector<real> v_atomic_charges;
  TVector<int> v_connections;
  //
  TVector<std::string> v_atomic_symbols;
  TVector<std::string> v_atomic_labels;
  //
  TMatrix<real> m_file_input;
  TMatrix<real> m_xyz;
  TMatrix<real> m_uvw;
  TMatrix<real> m_uvw_to_xyz;
  TMatrix<real> m_uvw_to_xyz_u;
};

#endif
