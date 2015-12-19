//========================================================================
// TOPMOL UTILITY
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

// #include<const.h>
#include <config_debug.h>
#include <ctopmol.h>
//#define _debugging_messages_

CTopmol::CTopmol(void){
  __filename = "TOPCAR"; // default name
#ifdef _debugging_messages_
  std::cout<<"CTopmol " <<__filename<<std::endl;
#endif
  __if_topmol = false;
  __number_of_topologies=0;
}

CTopmol::~CTopmol(void){
}

bool CTopmol::read_topmol(void){
  bool res;
#ifdef _debugging_messages_
  std::cout<<"read_topmol " <<__filename<<std::endl;
#endif
  res=read_topmol(__dir, __filename);
  //if(!res){
    //res=read_topmol(__dir,"TOPCAR");
    //if(res)
      //__filename = "TOPCAR";
  //}
  return res;
}

void CTopmol::clear_topmol(void){
  for(uint i=0; i<v_topology_list.size(); i++){
    v_topology_list[i].clear();
  }
  v_topology_list.clear();
  v_topology_names.clear();
}

bool CTopmol::read_topmol(std::string d, std::string f){
  uint _c=0, _s, _t;
  char _ch;
  clear_topmol();
  __number_of_topologies=0;
  //
  std::string topmol_file = d+"/"+f;
  std::ifstream topmol;
  topmol.open(topmol_file.c_str());
  if(!topmol.is_open()){
    __if_topmol = false;
#ifdef _debugging_messages_
    std::cout<<"WARNING: The topology file: "<<topmol_file<<std::endl;
    std::cout<<"         !!!...was not found...!!!  :(\n"<<std::endl;
#endif
    //__number_of_topologies=1;
    return false;
  }else{
    __if_topmol = true;
#ifdef _debugging_messages_
    printf("The topology was successful opened :)\n");
#endif
    _ch = topmol.peek();
    while(_ch == '#'){
      topmol.ignore(1024,NLN);
      _ch = topmol.peek();
    }
    topmol>>__number_of_topologies;
    uint _data;
    for(uint i=0; i<__number_of_topologies;i++){
      __new_topology.clear();
      //
      //topmol.ignore(1024,NLN);
      /*
      _ch = topmol.peek();
      while(_ch == '#' || _ch == NLN || _ch == SPC){
        topmol.ignore(1024,NLN);
        _ch = topmol.peek();
      }*/
      get_next(topmol);
      //
      topmol>>_t;
#ifdef _debugging_messages_
      std::cout<<"Topology type: "<<_t<<std::endl;
#endif
      //
      //topmol.ignore(1024,NLN);
      /*
      _ch = topmol.peek();
      while(_ch == '#' || _ch == NLN || _ch == SPC){
        topmol.ignore(1024,NLN);
        _ch = topmol.peek();
      }*/
      get_next(topmol);
      //
      topmol>>_s;
      for(uint j=0; j<_s;j++){
        topmol>>_data;
        __new_topology.v_atoms.push_back(_data);
        _c++;
      }
      get_next(topmol);
      if(_t>1){
        if(_t==2) _s=2;
        else _s=4;
        for(uint j=0; j<_s;j++){
          topmol>>_data;
          __new_topology.v_axis_index[j]=_data;
        }
      }
      __new_topology.type=_t;
      v_topology_list.push_back(__new_topology);
    }
#ifdef _debugging_messages_
    display_topology_summary();
#endif
    topmol.close();
  }
  __total_atoms=_c;
  return true;
}

bool CTopmol::get_next(std::ifstream& _s){
  char _ch;
  _ch = _s.peek();
  while(_ch == '#' || _ch == NLN || _ch == SPC){
    _s.ignore(1024,NLN);
    _ch = _s.peek();
  }
  return true;
}

bool CTopmol::save_topmol(void){
#ifdef _debugging_messages_
  std::cout<<" save_topmol() "<<__filename<<std::endl;
#endif
  return save_topmol(__dir, __filename);
}

bool CTopmol::save_topmol(std::string d, std::string f){
  //uint _c=0, _s, _t;
  //char _ch;
  uint _s;
  set_topmol_filename(d,f);
  //
#ifdef _debugging_messages_
  std::cout<<" save_topmol(d,f) "<<f<<std::endl;
#endif
  std::string topmol_file = __dir+"/"+__filename;
  std::ofstream topmol;
  topmol.open(topmol_file.c_str());
  if(!topmol.is_open()){
    __if_topmol = false;
#ifdef _debugging_messages_
    std::cout<<"WARNING: The TOPCAR: "<<topmol_file<<std::endl;
    std::cout<<"         !!!...was not created...!!!  :(\n"<<std::endl;
#endif
    return false;
  }else{
#ifdef _debugging_messages_
    std::cout<<"The TOPCAR was successful created :)"<<std::endl;
    std::cout<<__dir<<std::endl;
    std::cout<<__filename<<std::endl;
#endif
    topmol<<"# Edit this files if you know what you are doing"<<std::endl;
    topmol<<"# Number of Topologies:"<<std::endl;
    topmol<<__number_of_topologies<<std::endl;
#ifdef _debugging_messages_
    std::cout<<"Number of Topologies: "<<__number_of_topologies<<std::endl;
#endif
    uint _data;
    for(uint i=0; i<__number_of_topologies;i++){
      //
      topmol<<"# Topology type:"<<std::endl;
      topmol<<v_topology_list[i].type<<std::endl;
#ifdef _debugging_messages_
      std::cout<<"Topology type: "<<v_topology_list[i].type<<std::endl;
#endif
      _data = v_topology_list[i].size();
      topmol<<"# Number of atoms:"<<std::endl;
      topmol<<_data<<std::endl;
      //
      //
      for(uint j=0; j<_data; j++){
        topmol<<v_topology_list[i].v_atoms[j];
        if(j<(_data-1)) topmol<<" ";
        else topmol<<std::endl;
      }
      _data = v_topology_list[i].type;
      topmol<<"# Axis definition:"<<std::endl;
      if(_data>1){
        if(_data==2) _s=2;
        else _s=4;
        for(uint j=0; j<_s;j++){
          topmol<<v_topology_list[i].v_axis_index[j];
          if(j == 1 || j == 3) topmol<<std::endl;
          else topmol<<" ";
        }
      }
    }
    topmol<<"# end of topologies"<<std::endl;
#ifdef _debugging_messages_
    display_topology_summary();
#endif
    topmol.close();
  }
  //__total_atoms=_c;
  return true;
}

void CTopmol::set_topmol_multi_topology(TVector<uint>& v){
  //uint _s, _t;
  clear_topmol();
  //v_topology_list.clear();
  __number_of_topologies=0;
  uint u = v.size();
  for(uint j=0; j<u; j++){
    if(is_name(v[j])){
      v_topology_names.push_back(v[j]);
    }
  }
  __number_of_topologies=v_topology_names.size();
  //
  for(uint i=0; i<__number_of_topologies; i++){
    __new_topology.clear();
    for(uint j=0; j<u; j++){
      if(v[j] == v_topology_names[i]){
        __new_topology.v_atoms.push_back(j);
      }
    }
    for(uint j=0; j<2;j++){
    __new_topology.v_axis_index[j]=j;
    }
    __new_topology.type=2;
    v_topology_list.push_back(__new_topology);
  }
#ifdef _debugging_messages_
  std::cout<<" TOPCAR: topology names found:  "<<v_topology_names;
#endif
}

bool CTopmol::is_name(uint u){
  for(uint i=0; i<v_topology_names.size(); i++){
    if(v_topology_names[i]==u)
     return false;
  }
  return true;
}

void CTopmol::set_topmol_axis(uint frag, uint idx, uint val){
  v_topology_list[frag].set_index(idx,val);
}

void CTopmol::set_topmol_single_topology(uint u){
  clear_topmol();
  __number_of_topologies=0;
  set_default_topology(u,true);
}

void CTopmol::set_complete_topology(uint u){
  set_default_topology(u,false);
}

void CTopmol::set_default_topology(uint u, bool b){
  __new_topology.clear();
  for(uint j=0; j<u; j++){
    if(b){
      __new_topology.v_atoms.push_back(j);
    }else if(!if_topology(j)){
      __new_topology.v_atoms.push_back(j);
    }
  }
  for(uint j=0; j<2;j++){
    __new_topology.v_axis_index[j]=j;
  }
  __new_topology.type=2;
  v_topology_list.push_back(__new_topology);
  __number_of_topologies++;
#ifdef _debugging_messages_
  display_topology_summary();
#endif
}

void CTopmol::set_new_topology(TVector<uint>& v, uint u){
}

void CTopmol::eval_topmol_delete_atoms(TVector<uint>& v, uint u){
  v_topology_list[u].delete_atoms(v);
}

TVector<uint> CTopmol::get_topmol_atoms(TVector<uint>& v, uint u){
  return v_topology_list[u].get_atoms(v);
}

void CTopmol::add_topmol_atoms(TVector<uint>& v, uint u){
  v_topology_list[u].add_atoms(v);
}

void CTopmol::remove_topmol_topology(uint u){
  v_topology_list.remove(u);
  __number_of_topologies--;
}

void CTopmol::add_topmol_topology(CTopology& t){
  v_topology_list.push_back(t);
  __number_of_topologies++;
}

void CTopmol::set_topology(uint comp_begin, uint comp_end, uint axis_begin, uint axis_end, std::string name){
  //
#ifdef _debugging_messages_
  std::cout<<"fragment name: "<<name<<std::endl;
  std::cout<<"first atom: "<<comp_begin<<std::endl;
  std::cout<<"last atom: "<<comp_end<<std::endl;
  std::cout<<"axis begin: "<<axis_begin<<std::endl;
  std::cout<<"axis end: "<<axis_end<<std::endl;
#endif
  //
  __new_topology.name = name;
  //__new_fragment.ordered = true;
  __new_topology.v_atoms.push_back(comp_begin);
  __new_topology.v_atoms.push_back(comp_end);
  __new_topology.v_axis_index.push_back(axis_begin);
  __new_topology.v_axis_index.push_back(axis_end);
  //
  v_topology_list.push_back(__new_topology);
  //std::cout<<"atomic species: "<<ato_spe<<std::endl;
  //std::cout<<"atomic composition: "<<ato_comp<<std::endl;
}

void CTopmol::topmol_filename(std::string s){
#ifdef _debugging_messages_
  std::cout<<"topmol_filename "<<s<<std::endl;
#endif
  if(strstr(s.c_str(),"POSCAR")!=0 || strstr(s.c_str(),"CONTCAR")!=0 || strstr(s.c_str(),"vasp")!=0){
    __filename = "TOPCAR";
  }else{
    size_t pos = strcspn(s.c_str(),".");
    __filename = s.substr(0,pos);
    __filename += ".top";
  }
#ifdef _debugging_messages_
  std::cout<<"topmol_filename "<<__filename<<std::endl;
#endif
}

void CTopmol::topmol_dir(std::string s){
  __dir = s;
}

void CTopmol::set_topmol_filename(std::string _p, std::string _f){
  topmol_dir(_p);
  topmol_filename(_f);
#ifdef _debugging_messages_
  std::cout<<"set_topmol_filename " <<__filename<<std::endl;
#endif
}

bool CTopmol::if_topmol(void){
  return __if_topmol;
}

bool CTopmol::if_topology(const uint t){
  bool _if_exist=false;
  for(uint i=0;i<__number_of_topologies;i++){
    for(uint j=0;j<v_topology_list[i].size();j++){
      _if_exist=(t==v_topology_list[i].v_atoms[j]);
      if(_if_exist) break;
    }
  }
  return _if_exist;
}

/*
void CTopmol::is_complet_topologies(bool b){
  return __is_complete_topology=b;
}*/

uint CTopmol::get_number_of_topologies(void){
  return __number_of_topologies;
}

uint CTopmol::get_total_topology_atoms(void){
  return __total_atoms;
}

uint CTopmol::get_topology_size(uint i){
  return v_topology_list[i].size();
}

uint CTopmol::get_topology_axis_index_origin(uint i){
  return v_topology_list[i].axis_index_origin();
}

uint CTopmol::get_topology_axis_index_direction(uint i){
  return v_topology_list[i].axis_index_direction();
}

TVector<uint> CTopmol::get_topology_axis(const uint i){
  return v_topology_list[i].get_axis_index();
}

TVector<uint> CTopmol::get_topology_atoms(const uint i){
  return v_topology_list[i].v_atoms;
}

///int CTopmol::get_topology_twist_axis_begin(int i){
//  return v_topology_list[i].twist_axis_begin();
//}

//int CTopmol::get_topology_twist_axis_end(int i){
//  return v_topology_list[i].twist_axis_end();
//}

void CTopmol::display_topology_summary(void){
  std::cout<<"######################### TOPCAR summary #########################"<<std::endl;
  std::cout<<"Number of topologies: "<<__number_of_topologies<<std::endl;
  for(uint i=0; i<__number_of_topologies;i++){
    std::cout<<"Atoms in topology ["<<i<<"]: "<<v_topology_list[i].size()<<std::endl;
    std::cout<<"Topology type: "<<v_topology_list[i].type<<std::endl;
    std::cout<<"Atom list: "<<v_topology_list[i].v_atoms<<std::endl;
    if(v_topology_list[i].type>1)
      std::cout<<"Axis list: "<<v_topology_list[i].v_axis_index<<std::endl;
  }
  std::cout<<"######################### TOPCAR summary #########################"<<std::endl;
}
//

