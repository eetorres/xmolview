//========================================================================
// FILE: cmapcell.cxx
//
// Operations and modification to the Suprecell.
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

#include<cmapcell.h>

CMapcell::CMapcell(){
  __if_preview=true;
  __if_translation_map=false;
  __if_translation_map_1d=false;
  __if_translation_map_2d=false;
  __if_twist_map=false;
  __if_tilt_map=false;
  __if_precession_map=false;
  __if_translation_last_point=true;
  __if_cartesian_spheric=true;
  __active_fragment=0;
  s_scan_directory="scan";
}

void CMapcell::initialize_map(void){
#ifdef _debugging_messages_
  std::cout<<"MAPMOL: initialize_map [begin]"<<std::endl;
#endif
  supercell.initialize_fragments();
  eval_map_cartesian();
  eval_map_direct();
  eval_map_properties();
#ifdef _debugging_messages_
  std::cout<<"MAPMOL: initialize_map [end]"<<std::endl;
#endif
}

TVector<real> CMapcell::get_map_position_uvw1(uint u){
  real ang;
  TVector<real> _v, _uvw1, _uvw2;
  if(!__if_cartesian_spheric){
    ang = (u*__theta_step);
    _uvw1 = vrot_z(v_position_map_uvw1,ang);

    _uvw2 = vrot_z(v_position_map_uvw2,ang);
    _v=(_uvw2-_uvw1);
    v_position_map_step_2d=(__inv_stp_2d*_v);
  }else{
    _uvw1=v_position_map_uvw1+(real)u*v_position_map_step_1d;
  }
  return _uvw1;
}

TVector<real> CMapcell::get_map_position_uvw2(uint u){
  return (real)u*v_position_map_step_2d;
}

void CMapcell::eval_map_cartesian(void){
  supercell.update_fragmol_cartesian();
}

void CMapcell::eval_map_direct(void){
  supercell.update_fragmol_direct();
}

void CMapcell::eval_map_properties(void){
  TVector<real> _v;
  TVector<unsigned int> v_l;
  unsigned int _n, _s, _a;
  _a = supercell.get_fragmol_total_atoms();
  v_atomic_labels.resize(_a);
  v_atomic_symbols.resize(_a);
  v_atomic_numbers.resize(_a);
  _n = supercell.get_fragmol_number_of_fragments();
#ifdef _debugging_messages
  std::cout<<"MAPMOL: total atoms: "<<_a<<std::endl;
#endif
  for(unsigned int i=0; i<_n; i++){
    _s = supercell.get_fragment_size(i);
    v_l = supercell.gsf.get_topology_atoms(i);
    for(unsigned int j=0;j<_s;j++){
      v_atomic_labels[v_l[j]] = supercell.get_fragment_atomic_label(i,j);
      v_atomic_symbols[v_l[j]] = supercell.get_fragment_atomic_symbol(i,j);
      v_atomic_numbers[v_l[j]] = supercell.get_fragment_atomic_number(i,j);
    }
  }
#ifdef _debugging_messages
  std::cout<<"MAPMOL: Properties ready"<<std::endl;
  std::cout<<"MAPMOL: v_atomic_symbols: "<<v_atomic_symbols;
  std::cout<<"MAPMOL: v_atomic_numbers: "<<v_atomic_numbers;
#endif
}

// This function have to be called everytime the axis are changed
void CMapcell::set_map_active_fragment(const unsigned int i){
  __active_fragment=i;
#ifdef _debugging_messages_
  std::cout<<"MAPMOL: active fragment ="<<__active_fragment<<std::endl;
#endif
  supercell.set_fragmol_active_fragment(i);
  if(!supercell.is_fragmol_initialized()){
#ifdef _debugging_messages_
    std::cout<<"MAPMOL: initializing the fragment"<<std::endl;
#endif
    supercell.eval_fragmol_initial_position();
    // this is now done inside cfragment
    // compute_fragmol_axis_angles();
    supercell.eval_fragmol_initial_orientation();
    supercell.is_fragmol_initialized(true);
#ifdef _debugging_messages_
    std::cout<<"MAPMOL: fragment initialized"<<std::endl;
#endif
  }
  supercell.compute_fragmol_position_cartesian();
  eval_map_cartesian();
  eval_map_direct();
}

void CMapcell::map_update_active_fragment(void){
  supercell.compute_fragmol_position_cartesian();
  eval_map_cartesian();
  eval_map_direct();
}

void CMapcell::set_translation_last_point(const bool b){
  __if_translation_last_point=b;
}

void CMapcell::set_map_use_cartesian(const bool b){
  __if_cartesian_spheric=b;
#ifdef _debugging_messages_
  if(__if_cartesian_spheric)
    std::cout<<" MAPMOL: using cartesian"<<std::endl;
  else
    std::cout<<" MAPMOL: using spherical"<<std::endl;
#endif
}

void CMapcell::set_preview(const bool b){
  __if_preview=b;
}

void CMapcell::set_translation_map(const bool b){
  __if_translation_map=b;
}

void CMapcell::set_twist_map(const bool b){
  __if_twist_map=b;
}

void CMapcell::set_tilt_map(const bool b){
  __if_tilt_map=b;
}

void CMapcell::set_precession_map(const bool b){
  __if_precession_map=b;
}

void CMapcell::set_map_translation_steps_1d(const unsigned int i){
  __translation_map_steps_1d=i;
}

void CMapcell::set_map_translation_steps_2d(const unsigned int i){
  __translation_map_steps_2d=i;
}

void CMapcell::set_map_twist_steps(const unsigned int i){
  __twist_map_steps=i;
}

void CMapcell::set_map_tilt_steps(const unsigned int i){
  __tilt_map_steps=i;
}

void CMapcell::set_map_precession_steps(const unsigned int i){
  __precession_map_steps=i;
}

void CMapcell::set_map_fragment_twist(const real r){
  supercell.set_gsf_modified(true);
  supercell.set_fragmol_fragment_twist(r);
}

void CMapcell::set_map_fragment_precession(const real r){
  supercell.set_gsf_modified(true);
  supercell.set_fragmol_fragment_precession(r);
}

void CMapcell::set_map_fragment_tilt(const real r){
  supercell.set_gsf_modified(true);
  supercell.set_fragmol_fragment_tilt(r);
}

void CMapcell::set_map_fragment_position_u(const real r){
  supercell.set_gsf_modified(true);
  supercell.set_fragmol_fragment_position_u(r);
}

void CMapcell::set_map_fragment_position_v(const real r){
  supercell.set_gsf_modified(true);
  supercell.set_fragmol_fragment_position_v(r);
}

void CMapcell::set_map_fragment_position_w(const real r){
  supercell.set_gsf_modified(true);
  supercell.set_fragmol_fragment_position_w(r);
}

void CMapcell::set_map_twist1(const real r){
  __angle_map_twist1=r;
}

void CMapcell::set_map_twist2(const real r){
  __angle_map_twist2=r;
}

void CMapcell::set_map_tilt1(const real r){
  __angle_map_tilt1=r;
}

void CMapcell::set_map_tilt2(const real r){
  __angle_map_tilt2=r;
}

void CMapcell::set_map_precession1(const real r){
  __angle_map_precession1=r;
}

void CMapcell::set_map_precession2(const real r){
  __angle_map_precession2=r;
}

void CMapcell::set_map_fragment_position_uvw(const TVector<real>& _v){
  supercell.set_gsf_modified(true);
  supercell.set_fragmol_fragment_position_u(_v[0]);
  supercell.set_fragmol_fragment_position_v(_v[1]);
  supercell.set_fragmol_fragment_position_w(_v[2]);
}

void CMapcell::set_map_position_uvw1(const TVector<real>& _v){
  v_position_map_uvw1=_v;
}

void CMapcell::set_map_position_uvw2(const TVector<real>& _v){
  v_position_map_uvw2=_v;
}

void CMapcell::set_map_position_uvw3(const TVector<real>& _v){
  v_position_map_uvw3=_v;
}

void CMapcell::set_map_position_rpt1(const TVector<real>& _v){
  __angle_theta1=_v[2];
  TVector<real> _xyz = v_rpt_to_xyz(_v);
  set_map_position_uvw1(_xyz);
}

void CMapcell::set_map_position_rpt2(const TVector<real>& _v){
  __angle_theta2=_v[2];
  TVector<real> _xyz = v_rpt_to_xyz(_v);
  set_map_position_uvw2(_xyz);
}

void CMapcell::set_scan_title(const std::string s){
  s_scan_directory=s;
}

unsigned int CMapcell::get_map_total_atoms(void){
  return supercell.get_fragmol_total_atoms() ;
}

real CMapcell::get_translation_distance(void){
  return __translation_distance_1d;
}

real CMapcell::get_translation_distance_2d(void){
  return __translation_distance_2d;
}

real CMapcell::get_translation_step(void){
  return __translation_step_1d;
}

real CMapcell::get_translation_step_2d(void){
  return __translation_step_2d;
}

real CMapcell::get_twist_step(void){
  return __twist_map_step;
}

real CMapcell::get_tilt_step(void){
  return __tilt_map_step;
}

real CMapcell::get_precession_step(void){
  return __precession_map_step;
}

real CMapcell::get_map_axis_precession(void){
  return supercell.get_fragmol_axis_precession();
}

real CMapcell::get_map_axis_tilt(void){
  return supercell.get_fragmol_axis_tilt();
}

real CMapcell::get_map_backbone_precession(void){
  return supercell.get_fragmol_backbone_precession();
}

real CMapcell::get_map_backbone_tilt(void){
  return supercell.get_fragmol_backbone_tilt();
}

TVector<std::string> CMapcell::get_map_atomic_labels(void){
  return v_atomic_labels;
}

TVector<std::string> CMapcell::get_map_atomic_symbols(void){
  return v_atomic_symbols;
}

TVector<unsigned int> CMapcell::get_map_atomic_numbers(void){
  return v_atomic_numbers;
}

TVector<unsigned int> CMapcell::get_map_atomic_table(void){
  return v_atomic_table;
}

TVector<real> CMapcell::get_map_axis_angles(void){
  return supercell.get_fragmol_axis_angles();
}

TVector<real> CMapcell::get_map_basis_direct(void){
  return supercell.get_fragmol_basis_direct();
}

TVector<real> CMapcell::get_map_position_direct(void){
  return supercell.get_fragmol_position_direct();
}

TVector<real> CMapcell::get_map_position_uvw(void){
  return supercell.get_fragmol_position_uvw();
}

TVector<real> CMapcell::get_map_position_cartesian(void){
  return supercell.get_fragmol_position_cartesian();
}

TVector<real> CMapcell::get_map_centered_position_cartesian(void){
  return supercell.get_fragmol_centered_position_cartesian();
}

TMatrix<real> CMapcell::get_map_cartesian(void){
  //return m_cartesian;
  return supercell.get_cartesian();
}

TMatrix<real> CMapcell::get_map_direct(void){
  //return m_direct;
  return supercell.get_direct();
}

