//========================================================================
// FILE: cposmol.h 
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

#ifndef _POSCAR_H_
#define _POSCAR_H_

//#include <config.h>
#include <config_debug.h>
#include <time.h>
#include <string.h>
#include <msmvtl/tmatrix.h>
#include <msmvtl/tvmath.h>
#include <global.h>

class CPoscar{

public:

  CPoscar();
  ~CPoscar(){};

  // read/write files
  void clear(void);
  // POSCAR VASP file
  bool is_periodic(void);
  bool get_poscar_format(void);
  bool read_file(std::string,std::string);
  ////bool read_potcar(string,string);
  //
  void save_file(void);
  void save_file_as(std::string);
  void save_file_as(std::string,std::string);
  // Coordinate functions
  //void apply_pbc_to_cartesian(void);

  // set functions
  void set_velocities_zero(void);
  void set_accelerations_zero(void);
  void set_format(bool b){ __is_direct=b;};
  void set_sorted(bool b){ __is_sorted=b;};
  void set_centered(bool b){ __is_centered=b;};
  void set_total_atoms(uint u){ __total_atoms=u;};
  void set_total_cells(const uint, const uint, const uint, const uint);
  void set_lattice_constant(real f){__lattice_constant=f;};
  void set_atomic_numbers(TVector<uint>& v){ v_atomic_numbers=v;};
  void set_atomic_number_table(TVector<uint>& v){ v_atomic_number_table=v;};
  void set_atomic_composition_table(TVector<uint>& v){ v_atomic_composition_table=v;};
  void set_atomic_symbol_table(TVector<std::string>& v){ v_atomic_symbol_table=v;};
  void set_uvw_to_xyz(TMatrix<real>& v){ uvwTxyz=v;};
  void set_direct(const TMatrix<real>&);
  void set_cartesian(const TMatrix<real>&);
  // get functions
  //
  bool get_format(void);
  bool is_labels(void){ return __is_labels;};
  bool is_potcar(void){ return __is_potcar;}
  uint get_total_species(void);
  uint get_total_atoms(void);
  //TMatrix<real> get_bounding_box(void);

  ////TVector<string> get_symbols(void);
  TVector<uint> get_composition(void);
  ////string get_symbols(int);
  uint get_composition(uint);

  real get_cartesian(uint,uint);
  real get_direct(uint,uint);

  TVector<std::string> get_atomic_symbols(void){ return v_atomic_symbols;};

  TVector<real> get_direct(uint);
  TVector<real> get_cartesian(uint);
  TVector<real> get_direct_basis(void);

  TMatrix<real> get_xyz_input(void);
  TMatrix<real> get_direct(void);
  TMatrix<real> get_cartesian(void);
  TMatrix<real> get_centered_cartesian(void);

  //TMatrix<real> get_red_uvw_to_xyz(void){ return scl_uvwTxyz;};
  TMatrix<real> get_unit_uvw_to_xyz(void){ return unit_uvwTxyz;};
  TMatrix<real> get_uvw_to_xyz(void){ return uvwTxyz;};
  // vectorial funtions
  void get_rot_angle(void);
  // Extra functions
  void write_copyright(std::ofstream&);
  void write_potcar(std::string);
  void eval_sorting(void);
  void eval_centering(void);
  std::string get_date(void){
    time_t rawtime;
    char * _st, tmp_buffer[256];
    std::string stime;
    time (&rawtime);
    sprintf(tmp_buffer,"%s",ctime(&rawtime));
    _st = strtok(tmp_buffer,"\n");
    stime = _st;
    return stime;
  };

protected:
  //
  char poscar_head_buffer[10][256];
  //char tmp_buffer[256], *pch, dyn;
  std::string __dir, __filename, __system_title;
  std::string __poscar_file;
  std::string __config_filename;
  //
  bool __is_atom_comments;
  bool __is_selected_dynamics;
  bool __is_mit_format;
  bool __is_direct;
  bool __is_periodic;
  bool __is_velocity_reset;
  bool __is_acceleration_reset;
  bool __is_sorted;
  bool __is_centered;
  bool __fix_all_atoms;
  bool __is_labels;
  bool __is_read;
  bool __is_potcar;
  //
  real __lattice_constant;
  uint  __head_lines;
  uint  __atomic_species;
  //uint  __potcar_species;
  uint  __total_atoms, a1;
  uint __total_cells;
  // POTCAR data in MIT format
  TVector<std::string> v_atomic_symbols;
  //
  TVector<real> _vx, _vy, _vz, v_xyz;
  TVector<real> v_cells;
  TVector<uint> v_sorted_position;
  TVector<uint> v_atomic_composition;
  TVector<uint> v_atomic_numbers;
  TVector<uint> v_atomic_number_table;
  TVector<uint> v_atomic_composition_table;
  TVector<std::string> v_atomic_symbol_table;
  //
  TMatrix<char> m_dynamic;
  TMatrix<real> m_xyz_input;
  TMatrix<real> m_coordinates_uvw;
  TMatrix<real> m_coordinates_xyz;
  // transformation matrices
  TMatrix<real> uvwTxyz, uvwTxyz_t, unit_uvwTxyz;
  TMatrix<real> m_xyz, scl_uvwTxyz;
  //
};

#endif
