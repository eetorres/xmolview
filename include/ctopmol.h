//========================================================================
// TOPMOL
//
// Define the atomic structure
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

#ifndef _CTOPMOL_H_
#define _CTOPMOL_H_

#include <config_debug.h>
#include <string.h>
#include <msmvtl/tmatrix.h>
#include <msmvtl/tvmath.h>

//#define _topmol_debugging_messages_


//#include<string>

class CTopology{

public:
  CTopology(){
    v_axis_index.resize(4);
    clear();
  };

  ~CTopology(){
     clear();
     v_axis_index.clear();
  };

  std::string name;

  uint type;
  TVector<uint> v_atoms;
  TVector<uint> v_axis_index;

  void clear(void){
    name="topology name";
    type=-1;
    v_atoms.clear();
    //v_axis_index[0]=0;
    set_index(0,0);
    //v_axis_index[1]=0;
    set_index(1,0);
    //v_axis_index[2]=0;
    set_index(2,0);
    //v_axis_index[3]=0;
    set_index(3,0);
    //v_axis_index.clear();
  }

  void set_index(uint idx, uint val){
    v_axis_index[idx]=val;
  }

  void delete_atoms(TVector<uint>& v_l){
#ifdef _debugging_messages_
    std::cout<<" TOPMOL_H: old topology = "<<v_atoms<<std::endl;
#endif
    for(uint i=0; i<v_l.size(); i++){
      for(uint j=0; j<v_atoms.size(); j++){
        if(v_l[i]==v_atoms[j]){
          if(j != v_atoms.size()-1){
            v_atoms[j]=v_atoms[v_atoms.size()-1];
            v_atoms.resize(v_atoms.size()-1);
          }else{
            v_atoms.resize(v_atoms.size()-1);
          }
        }
      }
    }
#ifdef _debugging_messages_
    std::cout<<" TOPMOL_H: new topology = "<<v_atoms<<std::endl;
#endif
  }

  void add_atoms(TVector<uint>& v_l){
    for(uint i=0; i<v_l.size(); i++){
      v_atoms.push_back(v_l[i]);
    }
  }

  TVector<uint> get_atoms(TVector<uint>& v_l){
    TVector<uint> _v(v_l.size());
    for(uint i=0; i<v_l.size(); i++){
      _v[i]=v_atoms[v_l[i]];
    }
    return _v;
  }

  void size(uint i){ v_atoms.resize(i);};
  uint size(void){return v_atoms.size();};
  uint axis_index_origin(void){ return v_axis_index[0];};
  uint axis_index_direction(void){ return v_axis_index[1];};
  TVector<uint> get_axis_index(void){ return v_axis_index;};
  //int twist_axis_end(void){ return v_axis[2];};

};

class CTopmol {

public:

  CTopmol();
  ~CTopmol();
  //
  void topmol_filename(std::string);
  void topmol_dir(std::string);
  void set_topmol_filename(std::string,std::string);
  //
  bool if_topmol(void);
  bool read_topmol(void);                    // read a TOPCAR file
  void clear_topmol(void);                   // read a TOPCAR file
  bool save_topmol(void);                    // svae a TOPCAR file
  bool is_name(uint);
  bool read_topmol(std::string,std::string); // read a TOPCAR file
  bool save_topmol(std::string,std::string); // save a TOPCAR file
  //
  void set_topology(uint,uint,uint,uint,std::string);
  //void set_topology_list(void);
  void set_topmol_multi_topology(TVector<uint>&);
  void set_topmol_axis(uint, uint, uint);
  void set_topmol_single_topology(uint);
  void set_complete_topology(uint);
  void set_default_topology(uint,bool);
  void set_new_topology(TVector<uint>&,uint);
  //
  void eval_topmol_delete_atoms(TVector<uint>&,uint);
  TVector<uint> get_topmol_atoms(TVector<uint>&,uint);
  void add_topmol_atoms(TVector<uint>&,uint);
  void remove_topmol_topology(uint u);
  void add_topmol_topology(CTopology&);
  //
  //void is_complete_topology(bool);
  //
  uint get_total_topology_atoms(void);
  uint get_number_of_topologies(void);
  uint get_topology_size(void);
  uint get_topology_size(uint);
  //
  uint get_topology_axis_index_origin(uint);
  uint get_topology_axis_index_direction(uint);
  //int get_topology_twist_axis_begin(int);
  //int get_topology_twist_axis_end(int);
  TVector<uint> get_topology_axis(const uint);
  TVector<uint> get_topology_atoms(const uint);

private:

  bool __if_topmol, __if_complete_topology;
  uint __number_of_topologies, __total_atoms;
  //
  CTopology __new_topology;
  //
  TVector<CTopology> v_topology_list;
  TVector<uint> v_topology_names;
  //
  std::string __filename, __dir;
  //
  void display_topology_summary(void);
  bool if_topology(const uint);
  bool get_next(std::ifstream& _s);

};

#endif

