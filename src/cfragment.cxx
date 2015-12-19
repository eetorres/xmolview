//========================================================================
// FILE: cfragment.h
//
// This an utility program to manipulate and genarate structure files
//
// Copyright 2011-2015 by Edmanuel Torres
// email:   eetorres@gmail.com
//
// Latest update: 2013
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

#include<cfragment.h>

bool CFragment::__is_pbc=true;

CFragment::CFragment(){
  __is_show = true;
  //__is_pbc = true;
  __is_initialized = false;
  __is_centered = false;
  v_axis_angle.resize(3);
  v_backbone_angle.resize(2);
  v_origin_direct.resize(3);
  v_origin_cartesian.resize(3);
  v_basis_direct.resize(3);
}

CFragment::~CFragment(){};

void CFragment::clear(void){
  __size=0;
  __is_initialized = false;
  v_atomic_label.clear();
  v_atomic_symbol.clear();
  v_atomic_number.clear();
  m_direct.clear();
  m_cartesian.clear();
  m_position_direct.clear();
  m_position_cartesian.clear();
  m_centered_cartesian.clear();
  m_rotation_cartesian.clear();
}

void CFragment::copy(CFragment& f, TVector<uint>& v){
  uint s = v.size();
  TMatrix<real> m_tmp(s,3);
  TVector<std::string> v_str(s);
  TVector<uint> v_uint(s);
  for(uint i=0; i<s; i++){
    m_tmp[i]=m_uvw[v[i]];
  }
  f.set_uvw(m_tmp);
  //
  for(uint i=0; i<s; i++){
    m_tmp[i]=m_direct[v[i]];
  }
  f.set_direct(m_tmp);
  //
  for(uint i=0; i<s; i++){
    m_tmp[i]=m_cartesian[v[i]];
  }
  f.set_cartesian(m_tmp);
  //
  for(uint i=0; i<s; i++){
    m_tmp[i]=m_position_direct[v[i]];
  }
  f.set_position_direct(m_tmp);
  //
  for(uint i=0; i<s; i++){
    m_tmp[i]=m_position_cartesian[v[i]];
  }
  f.set_position_cartesian(m_tmp);
  //
  for(uint i=0; i<s; i++){
    m_tmp[i]=m_centered_direct[v[i]];
  }
  f.set_centered_direct(m_tmp);
  //
  for(uint i=0; i<s; i++){
    m_tmp[i]=m_centered_cartesian[v[i]];
  }
  f.set_centered_cartesian(m_tmp);
  //
  for(uint i=0; i<s; i++){
    m_tmp[i]=m_rotation_cartesian[v[i]];
  }
  f.set_rotation_cartesian(m_tmp);
  //
  for(uint i=0; i<s; i++){
    v_str[i]=v_atomic_label[v[i]];
  }
  f.set_atomic_label(v_str);
  //
  for(uint i=0; i<s; i++){
    v_str[i]=v_atomic_symbol[v[i]];
  }
  f.set_atomic_symbol(v_str);
  //
  for(uint i=0; i<s; i++){
    v_uint[i]=v_atomic_number[v[i]];
  }
  f.set_atomic_number(v_uint);
  //
  //
  f.set_unit_uvw_to_xyz(m_unit_uvwTxyz);
  f.set_uvw_to_xyz(m_uvwTxyz);
  //
}

// Move the list of atoms v to a new fragemnet f
void CFragment::move(CFragment& f, TVector<uint>& v){
#ifdef _FRAGMENT_DEBUG_MOVE_
  std::cout<<" FRAGMENT: start move"<<std::endl;
#endif
#ifdef _FRAGMENT_DEBUG_MOVE_
  std::cout<<" FRAGMENT: v_remove = "<<v;
#endif
  uint s = v.size();
  __size -= s;
  //TVector<uint> vs = v;
  //vs.sort_max();
#ifdef _FRAGMENT_DEBUG_MESSAGES_
  std::cout<<" FRAGMENT: size = "<<__size<<std::endl;
  std::cout<<" FRAGMENT: s = "<<s<<std::endl;
  std::cout<<" FRAGMENT: v = "<<v;
  //std::cout<<" FRAGMENT: vs = "<<vs;
#endif
  TMatrix<real> m_tmp(s,3);
  TVector<std::string> v_str(s);
  TVector<uint> v_uint(s);
#ifdef _FRAGMENT_DEBUG_MOVE_messages_
  std::cout<<" FRAGMENT: start move 0"<<std::endl;
#endif
  for(uint i=0; i<s; i++){
    //std::cout<<" FRAGMENT: i="<<i<<" v[i]="<<v[i]<<std::endl;
    //std::cout<<" FRAGMENT: start move 01"<<std::endl;
    m_tmp[i]=m_uvw[v[i]];
    //std::cout<<" FRAGMENT: start move 02"<<std::endl;
    //std::cout<<" FRAGMENT: m_uvw="<<m_uvw<<std::endl;
    //std::cout<<" FRAGMENT: m_tmp="<<m_tmp<<std::endl;
    //m_uvw.remove_row(v[i]);
  }
  for(uint i=0; i<s; i++){
    m_uvw.remove_row(v[i]);
    //std::cout<<" FRAGMENT: m_uvw="<<m_uvw<<std::endl;
  }
#ifdef _FRAGMENT_DEBUG_MOVE_messages_
  std::cout<<" FRAGMENT: start move 03"<<std::endl;
#endif
  f.set_uvw(m_tmp);
  //
#ifdef _FRAGMENT_DEBUG_MOVE_messages_
  std::cout<<" FRAGMENT: move 1"<<std::endl;
#endif
  for(uint i=0; i<s; i++){
    m_tmp[i]=m_direct[v[i]];
    //m_direct.remove_row(v[i]);
  }
  for(uint i=0; i<s; i++){
    m_direct.remove_row(v[i]);
    //std::cout<<" FRAGMENT: m_uvw="<<m_uvw<<std::endl;
  }
  f.set_direct(m_tmp);
  //
#ifdef _FRAGMENT_DEBUG_MOVE_MESSAGES_
  std::cout<<" FRAGMENT: move 2"<<std::endl;
#endif
  for(uint i=0; i<s; i++){
    m_tmp[i]=m_cartesian[v[i]];
    //m_cartesian.remove_row(v[i]);
  }
  for(uint i=0; i<s; i++){
    m_cartesian.remove_row(v[i]);
    //std::cout<<" FRAGMENT: m_uvw="<<m_uvw<<std::endl;
  }
  f.set_cartesian(m_tmp);
  //
#ifdef _FRAGMENT_DEBUG_MOVE_MESSAGES_
  std::cout<<" FRAGMENT: move 3"<<std::endl;
#endif
  for(uint i=0; i<s; i++){
    m_tmp[i]=m_position_direct[v[i]];
    //m_position_direct.remove_row(v[i]);
  }
  for(uint i=0; i<s; i++){
    m_position_direct.remove_row(v[i]);
    //std::cout<<" FRAGMENT: m_uvw="<<m_uvw<<std::endl;
  }
  f.set_position_direct(m_tmp);
  //
#ifdef _FRAGMENT_DEBUG_MOVE_MESSAGES_
  std::cout<<" FRAGMENT: move 4"<<std::endl;
#endif
  for(uint i=0; i<s; i++){
    m_tmp[i]=m_position_cartesian[v[i]];
    //m_position_cartesian.remove_row(v[i]);
  }
  for(uint i=0; i<s; i++){
    m_position_cartesian.remove_row(v[i]);
    //std::cout<<" FRAGMENT: m_uvw="<<m_uvw<<std::endl;
  }
  f.set_position_cartesian(m_tmp);
  //
#ifdef _FRAGMENT_DEBUG_MOVE_MESSAGES_
  std::cout<<" FRAGMENT: move 5"<<std::endl;
#endif
  for(uint i=0; i<s; i++){
    m_tmp[i]=m_centered_direct[v[i]];
    //m_centered_direct.remove_row(v[i]);
  }
  for(uint i=0; i<s; i++){
    m_centered_direct.remove_row(v[i]);
    //std::cout<<" FRAGMENT: m_uvw="<<m_uvw<<std::endl;
  }
  f.set_centered_direct(m_tmp);
  //
#ifdef _FRAGMENT_DEBUG_MOVE_MESSAGES_
  std::cout<<" FRAGMENT: move 6"<<std::endl;
#endif
  for(uint i=0; i<s; i++){
    m_tmp[i]=m_centered_cartesian[v[i]];
    //m_centered_cartesian.remove_row(v[i]);
  }
  for(uint i=0; i<s; i++){
    m_centered_cartesian.remove_row(v[i]);
    //std::cout<<" FRAGMENT: m_uvw="<<m_uvw<<std::endl;
  }
  f.set_centered_cartesian(m_tmp);
  //
#ifdef _FRAGMENT_DEBUG_MOVE_MESSAGES_
  std::cout<<" FRAGMENT: move 7"<<std::endl;
#endif
  for(uint i=0; i<s; i++){
    m_tmp[i]=m_rotation_cartesian[v[i]];
    //m_rotation_cartesian.remove_row(v[i]);
  }
  for(uint i=0; i<s; i++){
    m_rotation_cartesian.remove_row(v[i]);
    //std::cout<<" FRAGMENT: m_uvw="<<m_uvw<<std::endl;
  }
  f.set_rotation_cartesian(m_tmp);
  //
#ifdef _FRAGMENT_DEBUG_MOVE_MESSAGES_
  std::cout<<" FRAGMENT: move 8"<<std::endl;
#endif
  for(uint i=0; i<s; i++){
    //std::cout<<" FRAGMENT: start move 01"<<std::endl;
    v_str[i]=v_atomic_label[v[i]];
    //v_atomic_label.remove(v[i]);
  }
  //std::cout<<" FRAGMENT: start move 02"<<std::endl;
  for(uint i=0; i<s; i++){
    //std::cout<<" FRAGMENT: start move 03"<<std::endl;
    //std::cout<<" FRAGMENT: i="<<i<<" v[i]="<<v[i]<<std::endl;
    //std::cout<<" FRAGMENT: v_atomic_label="<<v_atomic_label<<std::endl;
    v_atomic_label.remove(v[i]);
    //std::cout<<" FRAGMENT: v_atomic_label="<<v_atomic_label<<std::endl;
    //std::cout<<" FRAGMENT: m_uvw="<<m_uvw<<std::endl;
  }
  //std::cout<<" FRAGMENT: start move 04"<<std::endl;
  f.set_atomic_label(v_str);
  //
#ifdef _FRAGMENT_DEBUG_MOVE_MESSAGES_
  std::cout<<" FRAGMENT: move 9"<<std::endl;
#endif
  for(uint i=0; i<s; i++){
    v_str[i]=v_atomic_symbol[v[i]];
    //v_atomic_symbol.remove(v[i]);
  }
  for(uint i=0; i<s; i++){
    v_atomic_symbol.remove(v[i]);
    //std::cout<<" FRAGMENT: m_uvw="<<m_uvw<<std::endl;
  }
  f.set_atomic_symbol(v_str);
  //
#ifdef _FRAGMENT_DEBUG_MOVE_MESSAGES_
  std::cout<<" FRAGMENT: move 10"<<std::endl;
#endif
#ifdef _FRAGMENT_DEBUG_MOVE_
  std::cout<<" FRAGMENT: v_atomic_number = "<<v_atomic_number;
#endif
  for(uint i=0; i<s; i++){
#ifdef _FRAGMENT_DEBUG_MOVE_
    std::cout<<" FRAGMENT: remove i="<<i<<", v[i]="<<v[i]<<std::endl;
#endif
    v_uint[i]=v_atomic_number[v[i]];
    //v_atomic_number.remove(v[i]);
  }
  for(uint i=0; i<s; i++){
    v_atomic_number.remove(v[i]);
    //std::cout<<" FRAGMENT: m_uvw="<<m_uvw<<std::endl;
  }
  f.set_atomic_number(v_uint);
#ifdef _FRAGMENT_DEBUG_MOVE_
  std::cout<<" FRAGMENT: v_atomic_number = "<<v_atomic_number;
  std::cout<<" FRAGMENT: v_uint = "<<v_uint;
#endif
  //
  //
  f.set_unit_uvw_to_xyz(m_unit_uvwTxyz);
  f.set_uvw_to_xyz(m_uvwTxyz);
  //
}

void CFragment::size(const uint i){
  __size=i;
  __is_initialized = false;
  v_atomic_label.resize(i);
  v_atomic_symbol.resize(i);
  v_atomic_number.resize(i);
  m_direct.resize(i,3);
  m_cartesian.resize(i,3);
  m_position_direct.resize(i,3);
  m_position_cartesian.resize(i,3);
  m_centered_cartesian.resize(i,3);
  m_rotation_cartesian.resize(i,3);
}

uint CFragment::size(){
  return __size;
}

void CFragment::is_initialized(bool b){
  __is_initialized=b;
}

void CFragment::is_pbc(bool b){
  __is_pbc=b;
}

// The function check whether a molecule is splited due periodic boundary conditions
void CFragment::eval_integrity(void){
#ifdef _FRAGMENT_DEBUG_INTEGRITY_
  std::cout<<" FRAGMENT: analysing the integrity"<<std::endl;
#endif
  eval_scaled_bond_integrity(get_axis_origin_index(),1.1);
  //__origin=get_axis_origin_index();
}

// Find the list atoms conected through bonds to the _o atom
// The function check whether a bond is splited due periodic boundary conditions
void CFragment::eval_scaled_bond_integrity(uint _o, real _s){
#ifdef _FRAGMENT_DEBUG_INTEGRITY_
  std::cout<<" FRAGMENT: evaluating integrity"<<std::endl;
#endif
  uint __length=__size-1, __origin;
  uint _p, _c, _f;
  real _r1, _r2, _rr, _n;
  bool _is_pbc=false;
  TVector<real> _v3, _u1, _u2, _u3;
  TVector<uint> v_i(__size);
  TVector<uint> v_l(__size);
  __origin=_o;
  v_i[0]=__origin;
#ifdef _FRAGMENT_DEBUG_INTEGRITY_
  //std::cout<<" FRAGMENT: initial list="<<v_i;
  std::cout<<" FRAGMENT: m_uvw="<<m_uvw;
  std::cout<<" FRAGMENT: basis direct="<<v_basis_direct;
  std::cout<<" FRAGMENT: v_pbc_box="<<v_pbc_box;
#endif
  for(uint k=0;k<__size;k++)
    v_l[k]=k;
  // in the case the origin is not the last one.
  if(__origin!=(__size-1))
    v_l[__origin]=v_l.last();
  //cout<<" atom list="<<v_l;
  _c=0;
  for(uint k=0;k<__size;k++){
    _p=v_i[k];
    _u1 = m_uvw[_p];
    _r1 = atom_rrgb[v_atomic_number[_p]][0];
#ifdef _FRAGMENT_DEBUG_INTEGRITY_
    std::cout<<" FRAGMENT: atomic number="<<v_atomic_number;
#endif
#ifdef _FRAGMENT_DEBUG_INTEGRITY_DATA_
    std::cout<<" FRAGMENT: for - "<<_p<<" - _u1="<<_u1;
#endif
    for(uint l=0;l<__length;l++){
      _f=v_l[l];
#ifdef _FRAGMENT_DEBUG_INTEGRITY_DATA_
      std::cout<<std::endl<<" FRAGMENT: check "<<_f<<" - ";
#endif
      if(_f!=_p){
        _is_pbc=false;
        _u2 = m_uvw[_f];
        _r2 = atom_rrgb[v_atomic_number[_f]][0];
        _u3 = _u2-_u1;
#ifdef _FRAGMENT_DEBUG_INTEGRITY_DATA_
        std::cout<<" FRAGMENT: for - "<<_p<<" - _u2="<<_u2;
        std::cout<<" FRAGMENT: no-pbc _u3="<<_u3;
#endif
        for(uint j=0; j<3; j++){
          if(_u3[j]>v_pbc_box[j]){
#ifdef _FRAGMENT_DEBUG_INTEGRITY_DATA_
            std::cout<<" FRAGMENT: apply pcb (-) ";
#endif
            _u2[j]-=v_basis_direct[j];
            _is_pbc=true;
          }else if(_u3[j]<-v_pbc_box[j]){
#ifdef _FRAGMENT_DEBUG_INTEGRITY_DATA_
            std::cout<<" FRAGMENT: apply pcb (+) ";
#endif
            _u2[j]+=v_basis_direct[j];
            _is_pbc=true;
          }
        }
        if(_is_pbc) _u3 = _u2-_u1;
#ifdef _FRAGMENT_DEBUG_INTEGRITY_DATA_
        std::cout<<" FRAGMENT: _u3="<<_u3;
#endif
        _v3 = m_unit_uvwTxyz*_u3;
        _n = _v3.norm();
#ifdef _FRAGMENT_DEBUG_INTEGRITY_DATA_
        std::cout<<" FRAGMENT: _v3="<<_v3;
        std::cout<<" FRAGMENT: norm="<<_n;
#endif
        _rr = _r1+_r2;
        _rr=_s*(_rr*_rr);
#ifdef _FRAGMENT_DEBUG_INTEGRITY_DATA_
        std::cout<<" FRAGMENT: bond="<<_rr;
#endif
        if(_n<_rr){
          _c++;
          if(_is_pbc)  // Fix PBC bonds
            m_uvw[_f]=_u2;
#ifdef _FRAGMENT_DEBUG_INTEGRITY_
          std::cout<<" FRAGMENT: #found "<<_f<<" #";
#endif
          if(l!=(__length-1)){
            v_l[l]=v_l[__length-1];
            l--;
          }
          v_i[_c]=_f;
          __length--;
        }
      }
    }
    //v_i.sort();
#ifdef _FRAGMENT_DEBUG_INTEGRITY_
    std::cout<<std::endl<<" FRAGMENT: updated atom list="<<v_l;
    std::cout<<std::endl<<" FRAGMENT: partial list="<<v_i<<std::endl;
#endif
  }
#ifdef _FRAGMENT_DEBUG_INTEGRITY_
  std::cout<<" FRAGMENT: integrity list="<<v_i;
  std::cout<<" FRAGMENT: integrity evaluated"<<std::endl;
  std::cout<<" -------------------------------------------------"<<std::endl;
#endif
}

// Find the list atoms conected through bonds to the _o atom
// The function check whether a bond is splited due periodic boundary conditions
void CFragment::eval_radial_integrity(uint _o, real _r){
#ifdef _FRAGMENT_DEBUG_INTEGRITY_
  std::cout<<" FRAGMENT: evaluating integrity"<<std::endl;
#endif
  uint __length=__size-1, __origin;
  uint _p, _c, _f;
  real _rr, _n;
  bool _is_pbc=false;
  TVector<real> _v3, _u1, _u2, _u3;
  TVector<uint> v_i(__size);
  TVector<uint> v_l(__size);
  __origin=_o;
  v_i[0]=__origin;
#ifdef _FRAGMENT_DEBUG_INTEGRITY_
  //std::cout<<" FRAGMENT: initial list="<<v_i;
  std::cout<<" FRAGMENT: m_uvw="<<m_uvw;
  std::cout<<" FRAGMENT: basis direct="<<v_basis_direct;
  std::cout<<" FRAGMENT: v_pbc_box="<<v_pbc_box;
#endif
  //_rr=_r+_r;
  _rr=_r*_r;
  for(uint k=0;k<__size;k++)
    v_l[k]=k;
  // in the case the origin is not the last one.
  if(__origin!=(__size-1))
    v_l[__origin]=v_l.last();
  //cout<<" atom list="<<v_l;
  _c=0;
  for(uint k=0;k<__size;k++){
    _p=v_i[k];
    _u1 = m_uvw[_p];
    //_r1 = atom_rrgb[v_atomic_number[_p]][0];
#ifdef _FRAGMENT_DEBUG_INTEGRITY_
    std::cout<<" FRAGMENT: atomic number="<<v_atomic_number;
#endif
#ifdef _FRAGMENT_DEBUG_INTEGRITY_DATA_
    std::cout<<" FRAGMENT: for - "<<_p<<" - _u1="<<_u1;
#endif
    for(uint l=0;l<__length;l++){
      _f=v_l[l];
#ifdef _FRAGMENT_DEBUG_INTEGRITY_DATA_
      std::cout<<std::endl<<" FRAGMENT: check "<<_f<<" - ";
#endif
      if(_f!=_p){
        _is_pbc=false;
        _u2 = m_uvw[_f];
        //_r2 = atom_rrgb[v_atomic_number[_f]][0];
        _u3 = _u2-_u1;
#ifdef _FRAGMENT_DEBUG_INTEGRITY_DATA_
        std::cout<<" FRAGMENT: for - "<<_p<<" - _u2="<<_u2;
        std::cout<<" FRAGMENT: no-pbc _u3="<<_u3;
#endif
        for(uint j=0; j<3; j++){
          if(_u3[j]>v_pbc_box[j]){
#ifdef _FRAGMENT_DEBUG_INTEGRITY_DATA_
            std::cout<<" FRAGMENT: apply pcb (-) ";
#endif
            _u2[j]-=v_basis_direct[j];
            _is_pbc=true;
          }else if(_u3[j]<-v_pbc_box[j]){
#ifdef _FRAGMENT_DEBUG_INTEGRITY_DATA_
            std::cout<<" FRAGMENT: apply pcb (+) ";
#endif
            _u2[j]+=v_basis_direct[j];
            _is_pbc=true;
          }
        }
        if(_is_pbc) _u3 = _u2-_u1;
#ifdef _FRAGMENT_DEBUG_INTEGRITY_DATA_
        std::cout<<" FRAGMENT: _u3="<<_u3;
#endif
        _v3 = m_unit_uvwTxyz*_u3;
        _n = _v3.norm();
#ifdef _FRAGMENT_DEBUG_INTEGRITY_DATA_
        std::cout<<" FRAGMENT: _v3="<<_v3;
        std::cout<<" FRAGMENT: norm="<<_n;
#endif

#ifdef _FRAGMENT_DEBUG_INTEGRITY_DATA_
        std::cout<<" FRAGMENT: bond="<<_rr;
#endif
        if(_n<_rr){
          _c++;
          if(_is_pbc)  // Fix PBC bonds
            m_uvw[_f]=_u2;
#ifdef _FRAGMENT_DEBUG_INTEGRITY_
          std::cout<<" FRAGMENT: #found "<<_f<<" #";
#endif
          if(l!=(__length-1)){
            v_l[l]=v_l[__length-1];
            l--;
          }
          v_i[_c]=_f;
          __length--;
        }
      }
    }
    //v_i.sort();
#ifdef _FRAGMENT_DEBUG_INTEGRITY_
    std::cout<<std::endl<<" FRAGMENT: updated atom list="<<v_l;
    std::cout<<std::endl<<" FRAGMENT: partial list="<<v_i<<std::endl;
#endif
  }
#ifdef _FRAGMENT_DEBUG_INTEGRITY_
  std::cout<<" FRAGMENT: integrity list="<<v_i;
  std::cout<<" FRAGMENT: integrity evaluated"<<std::endl;
  std::cout<<" -------------------------------------------------"<<std::endl;
#endif
}


// split fragment using the atom u as the seed for the new vdW fragment
TVector<uint> CFragment::compute_vdw_fragment(uint u, real _s){
  uint __length=__size, atom_seed;
  uint _p, _c, _f;
  real _r1, _r2, _rr, _norm;
  //bool _is_pbc=false;
  TVector<real> _v3, _u1, _u2, _u3;
  //TVector<uint> v_i(__size);
  // New fragment atom list
  TVector<uint> v_i;
  TVector<uint> v_l(__size);
#ifdef _FRAGMENT_DEBUG_NEW_FRAGMENT_
  std::cout<<" +++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  std::cout<<" FRAGMENT: create_new_fragment for "<<u<<std::endl;
#endif
  atom_seed=u;
  //eval_integrity(atom_seed);
  v_i.push_back(atom_seed);
  // atom_seed=get_axis_origin_index();
  // v_i[0]=atom_seed;
  //v_i[0]=u;
  //m_uvw = m_cartesian*m_unit_xyzTuvw;
#ifdef _FRAGMENT_DEBUG_DATA_
  //std::cout<<" FRAGMENT: initial list="<<v_i;
  std::cout<<" FRAGMENT: m_uvw="<<m_uvw;
  std::cout<<" FRAGMENT: basis="<<v_basis_direct;
#endif
  for(uint k=0;k<__size;k++)
    v_l[k]=k;
  // in the case the origin is not the last one.
  //if(atom_seed!=(__size-1))
    //v_l[atom_seed]=v_l.last();
#ifdef _FRAGMENT_DEBUG_NEW_FRAGMENT_
  std::cout<<" atom list="<<v_l;
#endif
  _c=0;
  //for(uint k=0;k<__size;k++){
  for(uint k=0;k<v_i.size();k++){
    _p=v_i[k];
#ifdef _FRAGMENT_DEBUG_NEW_FRAGMENT_data_
    std::cout<<" FRAGMENT: checking for: ["<<_p<<"]"<<std::endl;
#endif
    _u1 = m_uvw[_p];
    _r1 = atom_rrgb[v_atomic_number[_p]][0];
#ifdef _FRAGMENT_DEBUG_NEW_FRAGMENT_
    std::cout<<" FRAGMENT: atomic number="<<v_atomic_number;
#endif
#ifdef _FRAGMENT_DEBUG_NEW_FRAGMENT_data_
    std::cout<<" FRAGMENT: checking for: "<<_p<<" - _u1="<<_u1;
#endif
    for(uint l=0;l<__length;l++){
      _f=v_l[l];
#ifdef _FRAGMENT_DEBUG_NEW_FRAGMENT_
      std::cout<<std::endl<<" FRAGMENT: check "<<_f<<" - ";
#endif
      if( (_f!=_p) && (_f!=atom_seed) ){
        //_is_pbc=false;
        _u2 = m_uvw[_f];
        _r2 = atom_rrgb[v_atomic_number[_f]][0];
        _u3 = _u2-_u1;
#ifdef _FRAGMENT_DEBUG_NEW_FRAGMENT_data_
        std::cout<<" FRAGMENT: _u2="<<_u2;
        std::cout<<" FRAGMENT: _u3="<<_u3;
#endif
        _v3 = m_unit_uvwTxyz*_u3;
        _norm = _v3.norm();
#ifdef _FRAGMENT_DEBUG_NEW_FRAGMENT_data_
        std::cout<<" FRAGMENT: _v3="<<_v3;
        std::cout<<" FRAGMENT: norm="<<_n;
#endif
        _rr = _r1+_r2;
        _rr=_s*(_rr*_rr);
#ifdef _FRAGMENT_DEBUG_NEW_FRAGMENT_data_
        std::cout<<" FRAGMENT: bond="<<_rr;
#endif
        if(_norm<_rr){
          _c++;
#ifdef _FRAGMENT_DEBUG_NEW_FRAGMENT_
          std::cout<<" FRAGMENT: # connected "<<_f<<" #";
#endif
          if(l!=(__length-1)){
            v_l[l]=v_l[__length-1];
            l--;
          }
          //v_i[_c]=_f;
          v_i.push_back(_f);
          __length--;
        }
      }
    }
#ifdef _FRAGMENT_DEBUG_NEW_FRAGMENT_
    std::cout<<std::endl<<" FRAGMENT: updated atom list="<<v_l;
    std::cout<<std::endl<<" FRAGMENT: partial list="<<v_i<<std::endl;
#endif
  }
  v_i.sort();
#ifdef _FRAGMENT_DEBUG_NEW_FRAGMENT_
  std::cout<<" FRAGMENT: new fragment list="<<v_i;
  std::cout<<" +++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
#endif
  return v_i;
}

// split fragment using the atom u as the seed for the new vdW fragment
TVector<uint> CFragment::compute_radial_fragment(uint _u, real _r){
  uint __length=__size, atom_seed;
  uint _p, _c, _f;
  real _rr, _norm;
  //bool _is_pbc=false;
  TVector<real> _v3, _u1, _u2, _u3;
  // New fragment atom list
  TVector<uint> v_i;
  TVector<uint> v_l(__size);
  atom_seed=_u;
  v_i.push_back(atom_seed);
  for(uint k=0;k<__size;k++)
    v_l[k]=k;
  _c=0;
  //for(uint k=0;k<__size;k++){
  //for(uint k=0;k<v_i.size();k++){
    _rr=(_r*_r);
    _p=v_i[0];
    _u1 = m_uvw[_p];
    //_r1 = atom_rrgb[v_atomic_number[_p]][0];
    for(uint l=0;l<__size;l++){
      _f=v_l[l];
      if( (_f!=_p) && (_f!=atom_seed) ){
        //_is_pbc=false;
        _u2 = m_uvw[_f];
        //_r2 = atom_rrgb[v_atomic_number[_f]][0];
        _u3 = _u2-_u1;
        _v3 = m_unit_uvwTxyz*_u3;
        _norm = _v3.norm();
        if(_norm<_rr){
          _c++;
          if(l!=(__length-1)){
            v_l[l]=v_l[__length-1];
            l--;
          }
          v_i.push_back(_f);
          //__length--;
        }
      }
    }
  //}
  v_i.sort();
  return v_i;
}

void CFragment::compute_origin_direct(void){
  m_position_direct = (m_position_cartesian*m_xyzTuvw);
#ifdef _FRAGMENT_DEBUG_DATA_
  std::cout<<" FRAGMENT: m_position_direct "<<m_position_direct;
#endif
}

void CFragment::compute_origin_direct_from_direct(void){
  for(uint f=0;f<__size;f++){
    m_position_direct[f] = m_direct[f] + v_origin_direct;
  }
}

void CFragment::compute_origin_cartesian(void){
  if(is_initialized()){
    //std::cout<<" Edited"<<std::endl;
    eval_orientation();
    eval_position();
    //m_position_direct = (m_position_cartesian*m_xyzTuvw);
    compute_origin_direct();
    //m_position_direct=pcb(m_position_direct);
    if(__is_pbc){
      for(uint i=0; i<__size;i++){
#ifdef _FRAGMENT_DEBUG_PBC_MESSAGES_
        std::cout<<" apply pbc m_position_direct[i] = "<<m_position_direct[i]<<std::endl;
#endif
        m_position_direct[i]=pbc(m_position_direct[i]);
#ifdef _FRAGMENT_DEBUG_PBC_MESSAGES_
        std::cout<<" done pbc m_position_direct[i] = "<<m_position_direct[i]<<std::endl;
#endif
      }
    }
  }
  m_centered_direct=(m_position_direct-0.5);
  m_centered_cartesian=(m_centered_direct*m_uvwTxyz);
  //v_origin_cartesian=m_centered_cartesian[atom_seed];
#ifdef _FRAGMENT_DEBUG_DATA_
  std::cout<<" FRAGMENT: m_position_direct "<<m_position_direct;
  std::cout<<" FRAGMENT: m_centered_direct "<<m_centered_direct;
  std::cout<<" FRAGMENT: m_centered_cartesian "<<m_centered_cartesian;
  std::cout<<" FRAGMENT: m_uvwTxyz "<<m_uvwTxyz;
#endif
}

void CFragment::eval_initial_position(void){
#ifdef _FRAGMENT_DEBUG_MESSAGES_
  std::cout<<" FRAGMENT: Moving the cluster to (0,0,0) reference"<<std::endl;
#endif
  v_origin_direct = m_position_direct[get_axis_origin_index()];  //
  m_position_cartesian = m_position_direct*m_uvwTxyz;
  v_origin_cartesian = m_position_cartesian[get_axis_origin_index()];
  v_origin_uvw = v_origin_cartesian*m_unit_xyzTuvw;
  //v_origin_cartesian = m_unit_xyzTuvw*v_origin_uvw;
  v_initial_origin_cartesian=v_origin_cartesian;
  //if(!__is_centered){
    for(uint f=0;f<__size;f++){
      m_cartesian[f] = m_position_cartesian[f] - v_origin_cartesian;
      m_direct[f] = m_position_direct[f] - v_origin_direct;
    }
    //__is_centered=true;
  //}
  // Note (10/07/2015): It seems it is only need at the begenin
  if(!is_initialized()){
    compute_origin_direct_from_direct();
  }
  compute_origin_cartesian();
  m_uvw = m_cartesian*m_unit_xyzTuvw; // <--------------------------- test
  eval_integrity(); // <--------------------------- ???
#ifdef _FRAGMENT_DEBUG_DATA_
  std::cout<<" FRAGMENT: direct: "<<m_direct;
  std::cout<<" FRAGMENT: cartesian: "<<m_cartesian;
  std::cout<<" FRAGMENT: position direct = "<<v_origin_direct;
  std::cout<<" FRAGMENT: position cartesian = "<<v_origin_cartesian;
  std::cout<<" FRAGMENT: atomic symbols = "<<v_atomic_symbol;
  std::cout<<" FRAGMENT: atomic number = "<<v_atomic_number;
#endif
}

void CFragment::eval_initial_orientation(void){
  real p, t1, t2;
  // Moving the molecule to (0,0,0) angles
#ifdef _FRAGMENT_DEBUG_MESSAGES_
  std::cout<<" FRAGMENT: Rotating the molecule to (0,0,0) angles"<<std::endl;
#endif
  if(__size>1){
    TVector<real> _v1, _v2, _v3, _u1;
    real _prec, _tilt, _twist;
    real lz, lxy, _l;
    uint i,j;
    i=v_axis_index[0];
    j=v_axis_index[1];
    _v1 = get_centered_cartesian(i);
    _v2 = get_centered_cartesian(j);
    _v3 = _v2-_v1;
    _u1 = _v3/_v3.magnitude();
    lz = _v3[2];
    _v3[2] = 0;
    lxy = _v3.magnitude();
    // precession angle
    _prec = atan2(_v3[1],_v3[0]);
    // tilt angle
    _tilt = atan2(lxy,lz);
    //v_angle[1] = atan2(lxy,lz);
    set_axis_precession(_prec);
    set_axis_tilt(_tilt);
    // Backbone plane direction
    i=v_axis_index[2];
    j=v_axis_index[3];
    _v1 = get_centered_cartesian(i);
    _v2 = get_centered_cartesian(j);
    _v3 = _v2-_v1;
    _l = _v3*_u1;
    _v2 = _l*_u1;
    _v1 = _v3-_v2;
    lz = _v1[2];
    _v1[2] = 0;
    lxy = _v1.magnitude();
    _v3=vrot_z(_v3,-_prec);
    _v3=vrot_y(_v3,-_tilt);
#ifdef _FRAGMENT_DEBUG_DATA_
    //cout<<" FRAGMENT: origin= "<<_v1<<std::endl;
    std::cout<<" FRAGMENT: twist direction= "<<_v3<<std::endl;
#endif
    // axis twist angle
    _twist = atan2(_v3[1],_v3[0]);
    set_axis_twist(_twist);
    // backnode precession angle
    _prec = atan2(_v1[1],_v1[0]);
    // backbone tilt angle
    _tilt = atan2(lxy,lz);
    //v_angle[1] = atan2(lxy,lz);
    set_backbone_precession(_prec);
    set_backbone_tilt(_tilt);
#ifdef _FRAGMENT_DEBUG_MESSAGES_
    std::cout<<" FRAGMENT: precession, tilt, twist= ";
    std::cout<<RAD_DEG*_prec<<", "<<RAD_DEG*_tilt<<", "<<RAD_DEG*_twist<<std::endl;
#endif
    v_initial_axis_angle=v_axis_angle;
    p=get_axis_precession();
    t1=get_axis_tilt();
    t2=get_axis_twist();
    for(uint f=0;f<__size;f++){
      _v1=m_cartesian[f];
      // removing precession
      _v1=vrot_z(_v1,-p);
      // removing tilt
      _v1=vrot_y(_v1,-t1);
      // removing twist
      _v1=vrot_z(_v1,-t2);
      m_cartesian[f]=_v1;
    }
  }
  m_direct = (m_cartesian*m_xyzTuvw);
  compute_origin_direct_from_direct();
  compute_origin_cartesian();
}

void CFragment::eval_position(void){
  for(uint f=0;f<__size;f++){
    m_position_cartesian[f] = m_rotation_cartesian[f] + v_origin_cartesian;
  }
  v_origin_direct = m_position_direct[get_axis_origin_index()];
}

void CFragment::eval_orientation(void){
  real p, t1, t2;
  // Moving the molecule to (0,0,0) angles
#ifdef _FRAGMENT_DEBUG_MESSAGES_
  std::cout<<" FRAGMENT: Rotating the molecule to (0,0,0) angles"<<std::endl;
#endif
  TVector<real> _v;
  TVector<real> _u(3);

  p=get_axis_precession();
  t1=get_axis_tilt();
  t2=get_axis_twist();
  for(uint f=0;f<__size;f++){
    _v=m_cartesian[f];
    // add twist
    _v=vrot_z(_v,t2);
    // add tilt
    _v=vrot_y(_v,t1);
    // add precession
    _v=vrot_z(_v,p);
    m_rotation_cartesian[f]=_v;
  }
  _u.constant(1);
  _u[1]=0;
  _u[2]=0;
  _u=_u/_u.magnitude();
  _u=vrot_z(_u,t2);
  _u=vrot_y(_u,t1);
  _u=vrot_z(_u,p);
  real lz = _u[2];
  _u[2] = 0;
  real lxy = _u.magnitude();
  // precession angle
  v_backbone_angle[1] = atan2(_u[1],_u[0]);
  // tilt angle
  v_backbone_angle[0] = atan2(lxy,lz);
}

void CFragment::set_atomic_label(uint i,std::string _s){
  v_atomic_label[i]=_s;
}

void CFragment::set_atomic_symbol(uint i,std::string _s){
  v_atomic_symbol[i]=_s;
}

void CFragment::set_atomic_number(uint i, uint _a){
  v_atomic_number[i]=_a;
}

void CFragment::set_backbone_tilt(real _r){
  v_backbone_angle[0]=_r;
}

void CFragment::set_backbone_precession(const real _r){
  v_backbone_angle[1]=_r;
}

void CFragment::set_axis_tilt(const real _r){
  v_axis_angle[0]=_r;
}

void CFragment::set_axis_precession(const real _r){
  v_axis_angle[1]=_r;
}

void CFragment::set_axis_twist( const real _r){
  v_axis_angle[2]=_r;
}

void CFragment::set_origin_u(const real r){
  v_origin_uvw[0]=r;
  v_origin_cartesian = v_origin_uvw*m_unit_uvwTxyz;
}

void CFragment::set_origin_v(const real r){
  v_origin_uvw[1]=r;
  v_origin_cartesian = v_origin_uvw*m_unit_uvwTxyz;
}

void CFragment::set_origin_w(const real r){
  v_origin_uvw[2]=r;
  v_origin_cartesian = v_origin_uvw*m_unit_uvwTxyz;
}

void CFragment::set_axis_index(const TVector<uint>& _v){
  v_axis_index=_v;
}

void CFragment::set_axis_index(const uint idx, const uint val){
  v_axis_index[idx]=val;
}

void CFragment::set_atomic_label(const TVector<std::string>& _v){
  v_atomic_label=_v;
}

void CFragment::set_atomic_symbol(const TVector<std::string>& _v){
  v_atomic_symbol=_v;
}

void CFragment::set_atomic_number(const TVector<uint>& _v){
  v_atomic_number=_v;
}

void CFragment::set_direct(uint i,const TVector<real>& _v){
  m_direct[i]=_v;
}

void CFragment::set_cartesian(uint i,const TVector<real>& _v){
  m_cartesian[i]=_v;
}

void CFragment::set_position_direct(uint i,const TVector<real>& _v){
  m_position_direct[i]=_v;
}

void CFragment::set_position_cartesian(uint i,const TVector<real>& _v){
  m_position_cartesian[i]=_v;
}

void CFragment::set_uvw(const TMatrix<real>& _m){
   m_uvw=_m;
}

void CFragment::set_direct(const TMatrix<real>& _m){
  m_direct=_m;
}

void CFragment::set_cartesian(const TMatrix<real>& _m){
  m_cartesian=_m;
}

void CFragment::set_position_direct(const TMatrix<real>& _m){
  m_position_cartesian=_m;
}

void CFragment::set_position_cartesian(const TMatrix<real>& _m){
  m_position_cartesian=_m;
}

void CFragment::set_centered_direct(const TMatrix<real>& _m){
  m_centered_direct=_m;
}

void CFragment::set_centered_cartesian(const TMatrix<real>& _m){
  m_centered_cartesian=_m;
}

void CFragment::set_rotation_cartesian(const TMatrix<real>& _m){
  m_rotation_cartesian=_m;
}

void CFragment::set_unit_uvw_to_xyz(const TMatrix<real>& _m){
  m_unit_uvwTxyz=_m;
#ifdef _FRAGMENT_DEBUG_MATRICES_
  std::cout<<" FRAGMENT: normalized uvwTxyz="<<m_unit_uvwTxyz;
#endif
  // HERE THERE SEEMS TO BE A PROBLEM WITH THE INVERSE MATRIX COMPUTATION
  m_unit_xyzTuvw=m_unit_uvwTxyz.inverse();
#ifdef _FRAGMENT_DEBUG_MATRICES_
  std::cout<<" FRAGMENT: normalized xyzTuvw="<<m_unit_xyzTuvw;
#endif
}

void CFragment::set_uvw_to_xyz(const TMatrix<real>& _m){
  m_uvwTxyz=_m;
#ifdef _FRAGMENT_DEBUG_MATRICES_
  std::cout<<" FRAGMENT: uvwTxyz="<<m_uvwTxyz;
#endif
  m_xyzTuvw=m_uvwTxyz.inverse();
#ifdef _FRAGMENT_DEBUG_MATRICES_
  std::cout<<" FRAGMENT: xyzTuvw="<<m_xyzTuvw;
#endif
  v_basis_direct[0]=m_uvwTxyz[0].magnitude();
  v_basis_direct[1]=m_uvwTxyz[1].magnitude();
  v_basis_direct[2]=m_uvwTxyz[2].magnitude();
#ifdef _FRAGMENT_DEBUG_MATRICES_
  std::cout<<" FRAGMENT: direct axes="<<v_basis_direct;
#endif
  v_pbc_box=v_basis_direct/2.0;
#ifdef _FRAGMENT_DEBUG_MATRICES_
  std::cout<<" FRAGMENT: pbc box"<<v_pbc_box;
#endif
}

bool CFragment::is_initialized(void){
  return __is_initialized;
}

std::string CFragment::get_atomic_label(uint i){
   return v_atomic_label[i];
}

std::string CFragment::get_atomic_symbol(uint i){
   return v_atomic_symbol[i];
}

uint CFragment::get_axis_origin_index(void){
  return v_axis_index[0];
}

uint CFragment::get_atomic_number(uint i){
  return v_atomic_number[i];
}

real CFragment::get_axis_tilt(void){
  return v_axis_angle[0];
}

real CFragment::get_axis_precession(void){
  return v_axis_angle[1];
}

real CFragment::get_axis_twist(void){
  return v_axis_angle[2];
}

real CFragment::get_backbone_tilt(void){
  return v_backbone_angle[0];
}

real CFragment::get_backbone_precession(void){
  return v_backbone_angle[1];
}

TVector<real> CFragment::get_axis_angles(void){
  return v_axis_angle;
}

TVector<real> CFragment::get_basis_direct(void){
  return v_basis_direct;
}

TVector<real> CFragment::get_origin_direct(void){
  return v_origin_direct;
}

TVector<real> CFragment::get_origin_uvw(void){
  return v_origin_uvw;
}

TVector<real> CFragment::get_origin_cartesian(void){
  return v_origin_cartesian;
}

TVector<real> CFragment::get_centered_origin_cartesian(void){
  //return v_origin_cartesian;
  return m_centered_cartesian[get_axis_origin_index()];
}

TVector<real> CFragment::get_direct(uint i){
  return m_position_direct[i];
}

TVector<real> CFragment::get_cartesian(uint i){
  return m_position_cartesian[i];
}

TVector<real> CFragment::get_centered_cartesian(uint i){
  return m_centered_cartesian[i];
}

TMatrix<real> CFragment::get_direct(void){
  return m_position_direct;
}

TMatrix<real> CFragment::get_cartesian(void){
  return m_position_cartesian;
}

TMatrix<real> CFragment::get_centered_cartesian(void){
  return m_centered_cartesian;
}

// PBC funtions

TVector<real> CFragment::pbc(const TVector<real>& _v){
  TVector<real> _vpbc(3);
#ifdef _FRAGMENT_DEBUG_PBC_MESSAGES_
  std::cout<<" debug pbc"<<std::endl;
#endif
  _vpbc[0]= pbc(_v[0]);
  _vpbc[1]= pbc(_v[1]);
  _vpbc[2]= pbc(_v[2]);
#ifdef _FRAGMENT_DEBUG_PBC_MESSAGES_
  std::cout<<"pbc results = "<<_vpbc;
#endif
  return _vpbc;
};

real CFragment::pbc(const real r){
  real _r = r;
  if(_r<0){ //_r += pbc;
    while(_r<0) _r += 1.0;
  }else if(_r>=1.0){ //_r -= pbc;
    while(_r>=1.0) _r -= 1.0;
  }
  return _r;
}

void CFragment::show_information(void){
  std::cout<<" FRAGMENT: ================================================================ "<<std::endl;
  std::cout<<" FRAGMENT: v_atomic_number="<<v_atomic_number;
  std::cout<<" FRAGMENT: v_atomic_label="<<v_atomic_label;
  std::cout<<" FRAGMENT: v_atomic_symbol="<<v_atomic_symbol;
  std::cout<<" FRAGMENT: v_axis_index="<<v_axis_index;
  std::cout<<" FRAGMENT: v_axis_angle="<<v_axis_angle;
  std::cout<<" FRAGMENT: v_initial_axis_angle="<<v_initial_axis_angle;
  std::cout<<" FRAGMENT: v_backbone_angle="<<v_backbone_angle;
  std::cout<<" FRAGMENT: v_basis_direct="<<v_basis_direct;
  std::cout<<" FRAGMENT: v_pbc_box="<<v_pbc_box;
  std::cout<<" FRAGMENT: v_cartesian_axes="<<v_cartesian_axes;
  std::cout<<" FRAGMENT: v_origin_uvw="<<v_origin_uvw;
  std::cout<<" FRAGMENT: v_origin_direct="<<v_origin_direct;
  std::cout<<" FRAGMENT: v_origin_cartesian="<<v_origin_cartesian;
  std::cout<<" FRAGMENT: v_initial_origin_cartesian="<<v_initial_origin_cartesian;
  //
  std::cout<<" FRAGMENT: m_uvwTxyz="<<m_uvwTxyz;
  std::cout<<" FRAGMENT: m_xyzTuvw="<<m_xyzTuvw;
  std::cout<<" FRAGMENT: m_unit_uvwTxyz="<<m_unit_uvwTxyz;
  std::cout<<" FRAGMENT: m_unit_xyzTuvw="<<m_unit_xyzTuvw;
  //
  std::cout<<" FRAGMENT: m_uvw: "<<m_uvw;
  std::cout<<" FRAGMENT: m_direct: "<<m_direct;
  std::cout<<" FRAGMENT: m_cartesian: "<<m_cartesian;
  std::cout<<" FRAGMENT: m_position_direct: "<<m_position_direct;
  std::cout<<" FRAGMENT: m_position_cartesian: "<<m_position_cartesian;
  std::cout<<" FRAGMENT: m_centered_direct: "<<m_centered_direct;
  std::cout<<" FRAGMENT: m_centered_cartesian: "<<m_centered_cartesian;
  //
  std::cout<<" FRAGMENT: v_origin_direct = "<<v_origin_direct;
  std::cout<<" FRAGMENT: v_origin_cartesian = "<<v_origin_cartesian;
  std::cout<<" FRAGMENT: v_atomic_symbols = "<<v_atomic_symbol;
  std::cout<<" FRAGMENT: v_atomic_number = "<<v_atomic_number;
  std::cout<<" FRAGMENT: ================================================================ "<<std::endl;
}

/// Deprecated functions

/*
void CFragment::eval_integrity(void){
#ifdef _FRAGMENT_DEBUG_INTEGRITY_
  std::cout<<" FRAGMENT: analysing the integrity"<<std::endl;
#endif
  uint __length=__size-1, __origin;
  uint _p, _c, _f;
  real _r1, _r2, _rr, _norm;
  bool _is_pbc=false;
  TVector<real> _v3, _u1, _u2, _u3;
  TVector<uint> v_i(__size);
  TVector<uint> v_l(__size);
  __origin=get_axis_origin_index();
  v_i[0]=__origin;
#ifdef _FRAGMENT_DEBUG_INTEGRITY_
  std::cout<<" FRAGMENT: initial list="<<v_i;
  std::cout<<" FRAGMENT: m_uvw="<<m_uvw;
  std::cout<<" FRAGMENT: basis direct="<<v_basis_direct;
#endif
  for(uint k=0;k<__size;k++)
    v_l[k]=k;
  // in the case the origin is not the last one.
  if(__origin!=(__size-1))
    v_l[__origin]=v_l.last();
  //cout<<" atom list="<<v_l;
  _c=0;
  for(uint k=0;k<__size;k++){
    _p=v_i[k];
    _u1 = m_uvw[_p];
    _r1 = atom_rrgb[v_atomic_number[_p]][0];
#ifdef _FRAGMENT_DEBUG_INTEGRITY_
    std::cout<<" FRAGMENT: atomic number="<<v_atomic_number;
#endif
#ifdef _FRAGMENT_DEBUG_INTEGRITY_data_
    std::cout<<" FRAGMENT: for - "<<_p<<" - _u1="<<_u1;
#endif
    for(uint l=0;l<__length;l++){
      _f=v_l[l];
#ifdef _FRAGMENT_DEBUG_INTEGRITY_data_
      std::cout<<std::endl<<" FRAGMENT: check "<<_f<<" - ";
#endif
      if(_f!=_p){
        _is_pbc=false;
        _u2 = m_uvw[_f];
        _r2 = atom_rrgb[v_atomic_number[_f]][0];
        _u3 = _u2-_u1;
#ifdef _FRAGMENT_DEBUG_INTEGRITY_data_
        std::cout<<" FRAGMENT: no-pbc _u3="<<_u3;
#endif
        for(uint j=0; j<3; j++){
          if(_u3[j]>v_pbc_box[j]){
#ifdef _FRAGMENT_DEBUG_INTEGRITY_data_
            std::cout<<" FRAGMENT: apply pcb (-) ";
#endif
            _u2[j]-=v_basis_direct[j];
            _is_pbc=true;
          }else if(_u3[j]<-v_pbc_box[j]){
#ifdef _FRAGMENT_DEBUG_INTEGRITY_data_
            std::cout<<" FRAGMENT: apply pcb (+) ";
#endif
            _u2[j]+=v_basis_direct[j];
            _is_pbc=true;
          }
        }
        if(_is_pbc){
          //m_uvw[_f]=_u2;
          _u3 = _u2-_u1;
        }
#ifdef _FRAGMENT_DEBUG_INTEGRITY_data_
        std::cout<<" FRAGMENT: _u3="<<_u3;
#endif
        _v3 = m_unit_uvwTxyz*_u3;
        _norm = _v3.norm();
#ifdef _FRAGMENT_DEBUG_INTEGRITY_data_data_
        std::cout<<" FRAGMENT: _v3="<<_v3;
        std::cout<<" FRAGMENT: norm="<<_n;
#endif
        _rr = _r1+_r2;
        _rr=1.01*(_rr*_rr);
#ifdef _FRAGMENT_DEBUG_INTEGRITY_data_
        std::cout<<" FRAGMENT: bond="<<_rr;
#endif
        if(_norm<_rr){
          _c++;
          if(_is_pbc) // Fix PBC bonds
            m_uvw[_f]=_u2;
#ifdef _FRAGMENT_DEBUG_INTEGRITY_
          std::cout<<" FRAGMENT: #found "<<_f<<" #";
#endif
          if(l!=(__length-1)){
            v_l[l]=v_l[__length-1];
            l--;
          }
          v_i[_c]=_f;
          __length--;
        }
      }
    }
    //v_i.sort();
#ifdef _FRAGMENT_DEBUG_INTEGRITY_
    std::cout<<std::endl<<" FRAGMENT: updated atom list="<<v_l;
    std::cout<<std::endl<<" FRAGMENT: partial list="<<v_i<<std::endl;
#endif
  }
#ifdef _FRAGMENT_DEBUG_INTEGRITY_
  std::cout<<" FRAGMENT: integrity list="<<v_i;
  std::cout<<" FRAGMENT: integrity analysed"<<std::endl;
  std::cout<<" -------------------------------------------------"<<std::endl;
#endif
}

// Find the list atoms conected through bonds to the _o atom
// The function check whether a bond is splited due periodic boundary conditions
void CFragment::eval_integrity(uint _o, real _s){
#ifdef _FRAGMENT_DEBUG_INTEGRITY_
  std::cout<<" FRAGMENT: evaluating integrity"<<std::endl;
#endif
  uint __length=__size-1, __origin;
  uint _p, _c, _f;
  real _r1, _r2, _rr, _n;
  bool _is_pbc=false;
  TVector<real> _v3, _u1, _u2, _u3;
  TVector<uint> v_i(__size);
  TVector<uint> v_l(__size);
  __origin=_o;
  v_i[0]=__origin;
#ifdef _FRAGMENT_DEBUG_INTEGRITY_
  //std::cout<<" FRAGMENT: initial list="<<v_i;
  std::cout<<" FRAGMENT: m_uvw="<<m_uvw;
  std::cout<<" FRAGMENT: basis direct="<<v_basis_direct;
  std::cout<<" FRAGMENT: v_pbc_box="<<v_pbc_box;
#endif
  for(uint k=0;k<__size;k++)
    v_l[k]=k;
  // in the case the origin is not the last one.
  if(__origin!=(__size-1))
    v_l[__origin]=v_l.last();
  //cout<<" atom list="<<v_l;
  _c=0;
  for(uint k=0;k<__size;k++){
    _p=v_i[k];
    _u1 = m_uvw[_p];
    _r1 = atom_rrgb[v_atomic_number[_p]][0];
#ifdef _FRAGMENT_DEBUG_INTEGRITY_
    std::cout<<" FRAGMENT: atomic number="<<v_atomic_number;
#endif
#ifdef _FRAGMENT_DEBUG_INTEGRITY_DATA_
    std::cout<<" FRAGMENT: for - "<<_p<<" - _u1="<<_u1;
#endif
    for(uint l=0;l<__length;l++){
      _f=v_l[l];
#ifdef _FRAGMENT_DEBUG_INTEGRITY_DATA_
      std::cout<<std::endl<<" FRAGMENT: check "<<_f<<" - ";
#endif
      if(_f!=_p){
        _is_pbc=false;
        _u2 = m_uvw[_f];
        _r2 = atom_rrgb[v_atomic_number[_f]][0];
        _u3 = _u2-_u1;
#ifdef _FRAGMENT_DEBUG_INTEGRITY_DATA_
        std::cout<<" FRAGMENT: for - "<<_p<<" - _u2="<<_u2;
        std::cout<<" FRAGMENT: no-pbc _u3="<<_u3;
#endif
        for(uint j=0; j<3; j++){
          if(_u3[j]>v_pbc_box[j]){
#ifdef _FRAGMENT_DEBUG_INTEGRITY_DATA_
            std::cout<<" FRAGMENT: apply pcb (-) ";
#endif
            _u2[j]-=v_basis_direct[j];
            _is_pbc=true;
          }else if(_u3[j]<-v_pbc_box[j]){
#ifdef _FRAGMENT_DEBUG_INTEGRITY_DATA_
            std::cout<<" FRAGMENT: apply pcb (+) ";
#endif
            _u2[j]+=v_basis_direct[j];
            _is_pbc=true;
          }
        }
        if(_is_pbc) _u3 = _u2-_u1;
#ifdef _FRAGMENT_DEBUG_INTEGRITY_DATA_
        std::cout<<" FRAGMENT: _u3="<<_u3;
#endif
        _v3 = m_unit_uvwTxyz*_u3;
        _n = _v3.norm();
#ifdef _FRAGMENT_DEBUG_INTEGRITY_DATA_
        std::cout<<" FRAGMENT: _v3="<<_v3;
        std::cout<<" FRAGMENT: norm="<<_n;
#endif
        _rr = _r1+_r2;
        _rr=_s*(_rr*_rr);
#ifdef _FRAGMENT_DEBUG_INTEGRITY_DATA_
        std::cout<<" FRAGMENT: bond="<<_rr;
#endif
        if(_n<_rr){
          _c++;
          if(_is_pbc)  // Fix PBC bonds
            m_uvw[_f]=_u2;
#ifdef _FRAGMENT_DEBUG_INTEGRITY_
          std::cout<<" FRAGMENT: #found "<<_f<<" #";
#endif
          if(l!=(__length-1)){
            v_l[l]=v_l[__length-1];
            l--;
          }
          v_i[_c]=_f;
          __length--;
        }
      }
    }
    //v_i.sort();
#ifdef _FRAGMENT_DEBUG_INTEGRITY_
    std::cout<<std::endl<<" FRAGMENT: updated atom list="<<v_l;
    std::cout<<std::endl<<" FRAGMENT: partial list="<<v_i<<std::endl;
#endif
  }
#ifdef _FRAGMENT_DEBUG_INTEGRITY_
  std::cout<<" FRAGMENT: integrity list="<<v_i;
  std::cout<<" FRAGMENT: integrity evaluated"<<std::endl;
  std::cout<<" -------------------------------------------------"<<std::endl;
#endif
}
*/

// END

