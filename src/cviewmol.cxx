//========================================================================
// FILE: cviewmol.cxx
//
// Class layer to separate core structure manipulations from
// graphic output and GUI controls.
//
// Copyright 2006-2015 by Edmanuel Torres
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

#include<cviewmol.h>

bool CViewmol::read_files_view(void){
  supercell.gsf.read_topmol();
  return supercell.read_input_file();
}

void CViewmol::initialize_view(void){
  // instialization before the TOPCAR is read
  initialize_map();
  set_view_bounding_box();
}

// this the function in charge of coordinates actualization
void CViewmol::update_view(void){
  //+++supercell.compute_fragmol_position_cartesian();
  //---supercell.update_fragmol_cartesian();
  //+++eval_map_cartesian();
  //---supercell.update_fragmol_direct();
  //+++eval_map_direct();
}

void CViewmol::set_view_active_fragment(const unsigned int ui){
  set_map_active_fragment(ui);
}

void CViewmol::set_view_bounding_box(void){
  m_view_bbox = supercell.get_uvw_to_xyz();
}

void CViewmol::set_view_translation_map_1d(void){
  //initialize_map_translation_1d();
}

void CViewmol::set_view_translation_map_2d(void){
  //initialize_map_translation_2d();
}

void CViewmol::set_view_twist_map(void){
  //initialize_twist_map();
}

void CViewmol::set_view_tilt_map(void){
  //initialize_tilt_map();
}

void CViewmol::set_view_precession_map(void){
  //initialize_precession_map();
}

void CViewmol::set_view_translation_steps_1d(const unsigned int i){
  set_map_translation_steps_1d(i);
}

void CViewmol::set_view_translation_steps_2d(const unsigned int i){
  set_map_translation_steps_2d(i);
}

void CViewmol::set_view_twist_steps(const unsigned int i){
  set_map_twist_steps(i);
}

void CViewmol::set_view_tilt_steps(const unsigned int i){
  set_map_tilt_steps(i);
}

void CViewmol::set_view_precession_steps(const unsigned int i){
  set_map_precession_steps(i);
}

void CViewmol::set_view_axis_twist(const real r){
  set_map_fragment_twist(r);
}

void CViewmol::set_view_axis_precession(const real r){
  set_map_fragment_precession(r);
}

void CViewmol::set_view_axis_tilt(const real r){
  set_map_fragment_tilt(r);
}

void CViewmol::set_view_position_u(const real r){
  set_map_fragment_position_u(r);
}

void CViewmol::set_view_position_v(const real r){
  set_map_fragment_position_v(r);
}

void CViewmol::set_view_position_w(const real r){
  set_map_fragment_position_w(r);
}

void CViewmol::set_view_twist1(const real r){
  set_map_twist1(r);
}

void CViewmol::set_view_twist2(const real r){
  set_map_twist2(r);
}

void CViewmol::set_view_tilt1(const real r){
  set_map_tilt1(r);
}

void CViewmol::set_view_tilt2(const real r){
  set_map_tilt2(r);
}

void CViewmol::set_view_precession1(const real r){
  set_map_precession1(r);
}

void CViewmol::set_view_precession2(const real r){
  set_map_precession2(r);
}

void CViewmol::set_view_position_uvw(const TVector<real>& _v){
  set_map_fragment_position_uvw(_v);
}

void CViewmol::set_view_position_uvw1(const TVector<real>& _v){
  set_map_position_uvw1(_v);
}

void CViewmol::set_view_position_uvw2(const TVector<real>& _v){
  set_map_position_uvw2(_v);
}

void CViewmol::set_view_position_uvw3(const TVector<real>& _v){
  set_map_position_uvw3(_v);
}

void CViewmol::set_view_position_rpt(const TVector<real>& _v){
  TVector<real> _xyz = v_rpt_to_xyz(_v);
  set_map_fragment_position_uvw(_xyz);
}

void CViewmol::set_view_position_rpt1(const TVector<real>& _v){
  set_map_position_rpt1(_v);
}

void CViewmol::set_view_position_rpt2(const TVector<real>& _v){
  set_map_position_rpt2(_v);
}

uint CViewmol::get_view_total_atoms(void){
  return get_map_total_atoms();
}

uint CViewmol::get_view_total_fragments(void){
  return supercell.get_fragmol_number_of_fragments();
}

uint CViewmol::get_view_active_fragment(void){
  return supercell.get_fragmol_active_fragment();
}

uint CViewmol::get_view_tilt_axis_begin(unsigned int i){
  return supercell.gsf.get_topology_axis_index_origin(i);
}

uint CViewmol::get_view_tilt_axis_end(unsigned int i){
  return supercell.gsf.get_topology_axis_index_direction(i);
}

uint CViewmol::get_view_file_type(std::string s){
#ifdef _VIEWMOL_MESSAGES_
  std::cout<<"######################### VIEWMOL new file #########################"<<std::endl;
#endif
  if(strstr(s.c_str(),"POSCAR")!=0 || strstr(s.c_str(),"CONTCAR")!=0 || strstr(s.c_str(),"vasp")!=0 || strstr(s.c_str(),"vsp")!=0){
#ifdef _VIEWMOL_MESSAGES_
    std::cout<<" VIEWMOL: POSCAR found"<<std::endl;
#endif
    return INPUT_FILE_TYPE_VSP;
  }else if(strstr(s.c_str(),"xyz")!=0 || strstr(s.c_str(),"XYZ")!=0){
#ifdef _VIEWMOL_MESSAGES_
    std::cout<<" VIEWMOL: XYZ found"<<std::endl;
#endif
    return INPUT_FILE_TYPE_XYZ;
  }else if(strstr(s.c_str(),"gau")!=0){
#ifdef _VIEWMOL_MESSAGES_
    std::cout<<" VIEWMOL: GAUSSIAN found"<<std::endl;
#endif
    return INPUT_FILE_TYPE_GAU;
  }else if(strstr(s.c_str(),"pdb")!=0 || strstr(s.c_str(),"PDB")!=0){
#ifdef _VIEWMOL_MESSAGES_
    std::cout<<" VIEWMOL: PDB found"<<std::endl;
#endif
    return INPUT_FILE_TYPE_PDB;
  }else if(strstr(s.c_str(),"dlp")!=0 || strstr(s.c_str(),"CONFIG")!=0 || strstr(s.c_str(),"CFGMIN")!=0 || strstr(s.c_str(),"REVCON")!=0){
#ifdef _VIEWMOL_MESSAGES_
    std::cout<<" VIEWMOL: DLP found"<<std::endl;
#endif
    return INPUT_FILE_TYPE_DLP;
  }else if(strstr(s.c_str(),"zmat")!=0 || strstr(s.c_str(),"gmat")!=0){
#ifdef _VIEWMOL_MESSAGES_
    std::cout<<" VIEWMOL: ZMAT found"<<std::endl;
#endif
    return INPUT_FILE_TYPE_ZMT;
  }
  std::cout<<" unknown file type"<<std::endl;
  return INPUT_FILE_TYPE_UNKNOWN;
}

real CViewmol::get_view_axis_precession(void){
  return get_map_axis_precession();
}

real CViewmol::get_view_axis_tilt(void){
  return get_map_axis_tilt();
}

real CViewmol::get_view_backbone_precession(void){
  return get_map_backbone_precession();
}

real CViewmol::get_view_backbone_tilt(void){
  return get_map_backbone_tilt();
}

TVector<std::string> CViewmol::get_view_atomic_labels(void){
   return get_map_atomic_labels();
}

TVector<std::string> CViewmol::get_view_atomic_symbols(void){
   return get_map_atomic_symbols();
}

TVector<std::string> CViewmol::get_view_atomic_symbol_table(void){
   return supercell.get_fragmol_atomic_symbol_table();
}

TVector<uint> CViewmol::get_view_atom_table(void){
  return supercell.get_fragmol_atom_table();
}

TVector<uint> CViewmol::get_view_atomic_composition_table(void){
   return supercell.get_fragmol_atomic_composition_table();
}

TVector<uint> CViewmol::get_view_atomic_number_table(void){
   return supercell.get_fragmol_atomic_number_table();
}

TVector<uint> CViewmol::get_view_atomic_numbers(void){
  return get_map_atomic_numbers();
}

TVector<uint> CViewmol::get_view_fragment_table(void){
  return supercell.get_fragmol_fragment_table();
}

TVector<real> CViewmol::get_view_axis_angles(void){
  return get_map_axis_angles();
}

TVector<real> CViewmol::get_view_basis_direct(void){
  //return get_map_basis_direct();
  return supercell.get_fragmol_basis_direct();
}

TVector<real> CViewmol::get_view_position_direct(void){
  return get_map_position_direct();
}

TVector<real> CViewmol::get_view_position_uvw(void){
  return get_map_position_uvw();
}

TVector<real> CViewmol::get_view_position_spheric(void){
  TVector<real> _rpt;
  TVector<real> _xyz = get_map_position_uvw();
  std:: cout<<" cartesian = "<<_xyz<<std::endl;
  _rpt = v_xyz_to_rpt(_xyz);
  std:: cout<<" spherical = "<<_rpt<<std::endl;
  return _rpt;
}

TVector<real> CViewmol::get_view_position_cartesian(void){
  return get_map_position_cartesian();
}

TVector<real> CViewmol::get_view_centered_position_cartesian(void){
  return get_map_centered_position_cartesian();
}

TMatrix<real> CViewmol::get_view_bounding_box(void){
  // becausee of the symetry only half is needed for the visualization
  return 0.5*m_view_bbox;
}

// This function return the coordinates showed by the rendering
TMatrix<real> CViewmol::get_view_cartesian(void){
 return get_map_cartesian();
}

TMatrix<real> CViewmol::get_view_direct(void){
  return get_map_direct();
}

