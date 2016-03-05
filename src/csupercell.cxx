//========================================================================
// FILE: csupercell.cxx -> csupercell
//
// Abstraction layer of the supercell and fragments.
// The manipulation to individual fragaments is hidden.
//
// Copyright 2011-2015 by Edmanuel Torres
// email:   eetorres@gmail.com
//
// Comment: change this class/file names to CSupercell/csupercell
//
// Lastest update: Tue Jul 14 16:22:34 EDT 2015
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

#include <config_debug.h>
#include <csupercell.h>

CSupercell::CSupercell(){
  __is_potmol=false;
}

CSupercell::~CSupercell(void){
  clear_fragments();
}

void CSupercell::clear(void){
  __number_of_fragments=0;
  __total_atoms=0;
  __is_potmol=false;
  __active_fragment=0;
}

void CSupercell::clear_fragments(void){
  for(unsigned int i=0; i<v_fragments.size(); i++)
    v_fragments[i].clear();
  v_fragments.clear();
}

bool CSupercell::read_input_file(void){
  bool res=false;
  clear();
  // check if a POTCAR is available
  gsf.read_potcar();
  // Read the input structure file
  res = gsf.read_input_file();
  if(res){
    __total_atoms=gsf.get_total_atoms();
    __atomic_species=gsf.get_atomic_species();
    __is_direct=gsf.get_is_direct();
    __is_periodic=gsf.get_is_periodic();
    __input_format=gsf.get_input_format();
#ifdef _FRAGMOL_DEBUG_MESSAGES_
    std::cout<<" FRAGMOL: total atoms = "<<__total_atoms<<std::endl;
    std::cout<<" FRAGMOL: total species = "<<__atomic_species<<std::endl;
    std::cout<<" FRAGMOL: is direct? "<<__is_direct<<std::endl;
    std::cout<<" FRAGMOL: is periodic? "<<__is_periodic<<std::endl;
    std::cout<<" FRAGMOL: input format = "<<__input_format<<std::endl;
#endif
    //
    v_atomic_composition_table=gsf.get_atomic_composition_table();
    v_atomic_number_table=gsf.get_atomic_number_table();
    v_atomic_symbol_table=gsf.get_atomic_symbol_table();
    v_atom_type_table=gsf.get_atom_table();
    v_atomic_labels=gsf.get_atomic_labels();
    v_atomic_symbols=gsf.get_atomic_symbols();
    v_atomic_numbers=gsf.get_atomic_numbers();
    if(__input_format==2 || __input_format==4){
      v_fragment_table=gsf.get_fragment_table();
    }
    v_atom_cell_table.resize(__total_atoms);
#ifdef _FRAGMOL_DEBUG_MESSAGES_
    std::cout<<" FRAGMOL: Composition ="<<v_atomic_composition_table;
    std::cout<<" FRAGMOL: Z-Number Table ="<<v_atomic_number_table;
    std::cout<<" FRAGMOL: Symbol Table ="<<v_atomic_symbol_table;
    std::cout<<" FRAGMOL: Type ="<<v_atom_type_table;
    std::cout<<" FRAGMOL: Labels ="<<v_atomic_labels;
    std::cout<<" FRAGMOL: Symbols ="<<v_atomic_symbols;
    std::cout<<" FRAGMOL: Z-Numbers ="<<v_atomic_numbers;
#endif
    //
    //m_xyz_centered=gsf.get_xyz();
    m_xyz=gsf.get_xyz();
    m_uvw=gsf.get_uvw();
    m_uvw_to_xyz_u=gsf.get_uvw_to_xyz_u();
    m_uvw_to_xyz=gsf.get_uvw_to_xyz();
#ifdef _FRAGMOL_DATA_MESSAGES_
    std::cout<<" FRAGMOL: XYZ = "<<m_xyz;
    std::cout<<" FRAGMOL: UVW = "<<m_uvw;
    std::cout<<" FRAGMOL: uvwTxyz U = "<<m_uvw_to_xyz_u;
    std::cout<<" FRAGMOL: uvwTxyz  = "<<m_uvw_to_xyz;
#endif
    return res;
  }
  return res;
}

bool CSupercell::delete_atom(uint u){
  bool res=false;
  //clear();
  //
  /*
  __total_atoms=gsf.get_total_atoms();
  __atomic_species=gsf.get_atomic_species();
  //
  v_atomic_composition_table=gsf.get_atomic_composition_table();
  v_atomic_number_table=gsf.get_atomic_number_table();
  v_atomic_symbol_table=gsf.get_atomic_symbol_table();
  v_atom_type_table=gsf.get_atom_table();
  v_atomic_labels=gsf.get_atomic_labels();
  v_atomic_symbols=gsf.get_atomic_symbols();
  v_atomic_numbers=gsf.get_atomic_numbers();
  v_fragment_table=gsf.get_fragment_table();
  v_atom_cell_table.resize(__total_atoms);
#ifdef _FRAGMOL_DEBUG_MESSAGES_
  std::cout<<" FRAGMOL: atomic symbols="<<v_atomic_symbols;
  std::cout<<" FRAGMOL: atomic labels="<<v_atomic_labels;
#endif
  //
  m_xyz=gsf.get_xyz();
  m_uvw=gsf.get_uvw();
  m_uvw_to_xyz_u=gsf.get_uvw_to_xyz_u();
  m_uvw_to_xyz=gsf.get_uvw_to_xyz();
  */
  return res;
}

// the code below should be migrated to a new cread.h file API
void CSupercell::save_input_file(void){
  gsf.save_input_file();
}

void CSupercell::save_as_file(std::string _f, bool _l, bool _n,TMatrix<real> m, int xc, int yc, int zc, int tc, bool _bb){
  if(__output_format == OUTPUT_FORMAT_ATM_FRG || __output_format == OUTPUT_FORMAT_NAT_FRG )
    gsf.set_fragment_table(v_fragment_table,__number_of_fragments);
  gsf.set_bounding_box(_bb);
  gsf.set_labels(_l);
  gsf.set_numbers(_n);
  gsf.set_xyz_cells(xc,yc,zc,tc);
  gsf.set_xyz(m);
  gsf.save_as_file(_f);
}

void CSupercell::save_as_file(std::string _p, std::string _f, bool _l, bool _n, TMatrix<real> m, int xc, int yc, int zc, int tc, bool _bb){
  if(__output_format == OUTPUT_FORMAT_ATM_FRG || __output_format == OUTPUT_FORMAT_NAT_FRG )
    gsf.set_fragment_table(v_fragment_table,__number_of_fragments);
  gsf.set_bounding_box(_bb);
  gsf.set_labels(_l);
  gsf.set_numbers(_n);
  gsf.set_xyz_cells(xc,yc,zc,tc);
  gsf.set_xyz(m);
  gsf.save_as_file(_p,_f);
}

// create fragments from the input topology
void CSupercell::create_initial_fragments(void){
  // making space for the fragments
  __number_of_fragments = gsf.get_number_of_topologies();
  v_fragments.resize(get_fragmol_number_of_fragments());
  for(uint i=0;i<__number_of_fragments;i++){
    v_fragments[i].size(gsf.get_topology_size(i));
  }
}

void CSupercell::eval_initial_fragments(void){
  uint _s=0;//, _c=0;
  TVector<uint> v_l, v_n, v_i;
  TVector<std::string> v_s;
  TVector<real> _v;
  for(uint i=0;i<get_fragmol_number_of_fragments();i++){
    _s=v_fragments[i].size();
    v_l = gsf.get_topology_atoms(i);
    v_i = gsf.get_topology_axis(i);
    v_fragments[i].set_axis_index(v_i);
#ifdef _FRAGMOL_DEBUG_MESSAGES_
    if(is_direct())
      std::cout<<" FRAGMOL: set direct"<<std::endl;
    else
      std::cout<<" FRAGMOL: set cartesian"<<std::endl;
#endif
    // the posmol indicates the coordinate input format
    // (1) direct, (0) cartesian
    for(uint j=0;j<_s;j++){
      // store the position of the current atom inside the fragment as a table
      v_atom_cell_table[v_l[j]]=j;
      if(is_direct()){
        _v = get_direct(v_l[j]);
        v_fragments[i].set_position_direct(j,_v);
      }else{
        _v = get_cartesian(v_l[j]);
        v_fragments[i].set_position_cartesian(j,_v);
      }
      v_fragment_table[v_l[j]]=i+1;
      v_fragments[i].set_atomic_label(j,v_atomic_labels[(v_l[j])]);
      v_fragments[i].set_atomic_symbol(j,v_atomic_symbols[(v_l[j])]);
      v_fragments[i].set_atomic_number(j,v_atomic_numbers[(v_l[j])]);
      v_fragments[i].is_pbc(__is_periodic);
    }
  }
  // (1) direct, (0) cartesian
  if(is_direct()){
    compute_fragmol_all_cartesian();
#ifdef _FRAGMOL_DEBUG_MESSAGES_
    std::cout<<" FRAGMOL: the cartesian were computed"<<std::endl;
#endif
  }else{
    // It is needed when the input is in cartesian coordinates
    compute_fragmol_all_direct();
    //compute_direct();
#ifdef _FRAGMOL_DEBUG_MESSAGES_
    std::cout<<" FRAGMOL: the direct were computed"<<std::endl;
#endif
  }
}

void CSupercell::initialize_fragments(void){
#ifdef _FRAGMOL_DEBUG_MESSAGES_
  std::cout<<" FRAGCAR: initialize fragments"<<std::endl;
#endif
  clear_fragments();
  __total_atoms=get_total_atoms();
  if(!gsf.if_topmol()){
#ifdef _FRAGMOL_DEBUG_MESSAGES_
    std::cout<<" FRAGCAR: !!! No TOPCAR fragment file to load !!!"<<std::endl;
    std::cout<<" FRAGCAR: Build the best topology for the system"<<std::endl;
#endif
    if(__input_format==2 || __input_format==4){
#ifdef _FRAGMOL_DEBUG_MESSAGES_
      std::cout<<" FRAGCAR: Topology based on the fragment numbers"<<std::endl;
#endif
      gsf.set_topmol_multi_topology(v_fragment_table);
    }else{
#ifdef _FRAGMOL_DEBUG_MESSAGES_
      std::cout<<" FRAGCAR: !!! No topology definition to use !!!"<<std::endl;
      std::cout<<" FRAGCAR: Build a single fragment system"<<std::endl;
#endif
      gsf.set_topmol_single_topology(__total_atoms);
    }
  }else if(__total_atoms > gsf.get_total_topology_atoms()){
#ifdef _FRAGMOL_DEBUG_MESSAGES_
    std::cout<<" FRAGCAR: !!! Fragments not for all atoms !!!"<<std::endl;
    std::cout<<" FRAGCAR: Build a fragment with remaining atoms"<<std::endl;
#endif
    gsf.set_complete_topology(__total_atoms);
  }
#ifdef _FRAGMOL_DEBUG_MESSAGES_
  else{
    std::cout<<" FRAGCAR: !!! All atoms in fragments !!!"<<std::endl;
  }
#endif
  v_fragment_table.resize(__total_atoms);
  create_initial_fragments();
  eval_initial_fragments();
#ifdef _FRAGMOL_DEBUG_MESSAGES_
  std::cout<<" FRAGCAR:  ["<<get_fragmol_number_of_fragments()<<"] fragments loaded"<<std::endl;
#endif
}

void CSupercell::eval_cell_table(void){
  TVector<uint> v_l;
  uint _s;
  for(uint i=0;i<get_fragmol_number_of_fragments();i++){
    _s=v_fragments[i].size();
    v_l = gsf.get_topology_atoms(i);
    for(uint j=0;j<_s;j++){
      // store the position of the current atom inside the fragment as a table
      v_atom_cell_table[v_l[j]]=j;
    }
  }
#ifdef _FRAGMOL_DEBUG_MESSAGES_
  std::cout<<" FRAGMOL: updated atom cell table: "<<v_atom_cell_table;
#endif
}

void CSupercell::eval_connections(const TMatrix<uint>& _m, uint u){
  gsf.eval_connections(_m,u);
}

// Construct a new fragment give a list of atoms
bool CSupercell::eval_new_fragment(const TVector<uint>& _iv){
  bool res=false;
  TVector<uint> new_topology_atoms, v_l, v_i;;
  CFragment new_frag;
  CTopology new_top;
  TVector<uint> new_fragment_atoms = _iv;
  // sort the atom list
  new_fragment_atoms.sort_max();
  // if the new atom list is smaller than the active fragment, create a new fragment
  if( (new_fragment_atoms.size() > 0) &&  (new_fragment_atoms.size() < v_fragments[__active_fragment].size()) ){
    new_topology_atoms=gsf.get_topmol_atoms(new_fragment_atoms,__active_fragment);
    // it works now
    gsf.eval_topmol_delete_atoms(new_topology_atoms,__active_fragment);
    new_top.v_atoms = new_topology_atoms;
    new_top.type=2;
    __number_of_fragments++;
    gsf.add_topmol_topology(new_top);
    //
    new_frag.size(new_topology_atoms.size());
    // BUG BELOW ???
    v_fragments[__active_fragment].move(new_frag,new_fragment_atoms);
    //new_frag.show_information();
    v_fragments.push_back(new_frag);
    //v_fragments[__active_fragment].show_information();
    // update the fragment table
    v_l = gsf.get_topology_atoms(__number_of_fragments-1);
    for(uint j=0;j<v_fragments[__number_of_fragments-1].size();j++){
      v_fragment_table[v_l[j]]=__number_of_fragments;
    }
    v_i = gsf.get_topology_axis(__number_of_fragments-1);
    v_fragments[__number_of_fragments-1].set_axis_index(v_i);
    // the new fragmented part should be updated
    // there is a bug when an initially structure ???
    // is moved and then fragmented ???
    // see water case:/home/etorres/src/utils/xmol/test/xyz/WaterWater
    res = true;
  }
#ifdef _FRAGMOL_DEBUG_MESSAGES_
  else{
    std::cout<<" FRAGMOL: !!! No more atoms can be fragmented !!!"<<std::endl;;
  }
#endif
  //std::cout<<" ******************************************"<<std::endl;
  return res;
}

// param u: initial selected atom
// param sw: atom list switch
// param _s: distance scaling factor
bool CSupercell::eval_scaled_fragment(uint _u, bool _sw, real _scale){
  bool res=false;
  TVector<uint> new_fragment_atoms;
  uint atom_seed;
  if(_sw) // the atom number given by the fragment table
    atom_seed=v_atom_cell_table[_u];
  else   // the atom number in the fragment
    atom_seed=_u;
  // check if the structure is splited due PBC
  v_fragments[__active_fragment].eval_scaled_bond_integrity(atom_seed,_scale);
  // get the linked list of atoms
  new_fragment_atoms=v_fragments[__active_fragment].compute_vdw_fragment(atom_seed,_scale);
  // Create the new fragment
  res = eval_new_fragment(new_fragment_atoms); //<--------------------------
  // (1) direct, (0) cartesian
  if(is_direct()){
    //compute_fragmol_all_cartesian();
    compute_fragmol_cartesian(__active_fragment);
    compute_fragmol_cartesian(__number_of_fragments-1);
#ifdef _FRAGMOL_DEBUG_MESSAGES_
    std::cout<<" FRAGMOL: the cartesian were computed"<<std::endl;
#endif
  }else{
    // It is needed when the input structure is given in cartesian coordinates
    compute_fragmol_direct(__active_fragment);
    compute_fragmol_direct(__number_of_fragments-1);
#ifdef _FRAGMOL_DEBUG_MESSAGES_
    std::cout<<" FRAGMOL: the direct were computed"<<std::endl;
#endif
  }
#ifdef _FRAGMOL_DEBUG_MESSAGES_
  new_fragment_atoms=gsf.get_topology_atoms(__active_fragment);
  std::cout<<" FRAGMOL: new active topology atoms = "<<new_fragment_atoms;
#endif
  return res;
}

// param u: initial selected atom
// param sw: atom list switch
// param _s: distance scaling factor
bool CSupercell::eval_radial_fragment(uint _u, bool _sw, real _scale){
  bool res=false;
  TVector<uint> new_fragment_atoms;
  uint atom_seed;
  if(_sw) // the atom number given by the fragment table
    atom_seed=v_atom_cell_table[_u];
  else   // the atom number in the fragment
    atom_seed=_u;
  // check if the structure is splited due PBC
  v_fragments[__active_fragment].eval_radial_integrity(atom_seed,_scale);
  // get the linked list of atoms
  new_fragment_atoms=v_fragments[__active_fragment].compute_radial_fragment(atom_seed,_scale);
  // Create the new fragment
  res = eval_new_fragment(new_fragment_atoms); //<--------------------------
  // (1) direct, (0) cartesian
  if(is_direct()){
    //compute_fragmol_all_cartesian();
    compute_fragmol_cartesian(__active_fragment);
    compute_fragmol_cartesian(__number_of_fragments-1);
#ifdef _FRAGMOL_DEBUG_MESSAGES_
    std::cout<<" FRAGMOL: the cartesian were computed"<<std::endl;
#endif
  }else{
    // It is needed when the input structure is given in cartesian coordinates
    compute_fragmol_direct(__active_fragment);
    compute_fragmol_direct(__number_of_fragments-1);
#ifdef _FRAGMOL_DEBUG_MESSAGES_
    std::cout<<" FRAGMOL: the direct were computed"<<std::endl;
#endif
  }
#ifdef _FRAGMOL_DEBUG_MESSAGES_
  new_fragment_atoms=gsf.get_topology_atoms(__active_fragment);
  std::cout<<" FRAGMOL: new active topology atoms = "<<new_fragment_atoms;
#endif
  return res;
}

bool CSupercell::eval_scaled_fragment(uint u, real _s){
  //std::cout<<" FRAGMOL: m_cartesian: "<<m_xyz;
#ifdef _FRAGMOL_DATA_MESSAGES_
  v_fragments[__active_fragment].show_information();
#endif
  bool is_new_frag=true;
  // Use the atom number
  is_new_frag=eval_scaled_fragment(u,true,_s);
  /////////////////////////set_map_active_fragment(__active_fragment);
  eval_cell_table();
  update_fragmol_cartesian();
  update_fragmol_direct();
  /////////////////////initialize_map();
  //std::cout<<" FRAGMOL: m_cartesian: "<<m_xyz;
  return is_new_frag;
}

// auto search for fragments separated by van der Waals radius distances
void CSupercell::eval_scaled_fragments(real _s){
  //std::cout<<" FRAGMOL: m_cartesian: "<<m_xyz;
#ifdef _FRAGMOL_DATA_MESSAGES_
  v_fragments[__active_fragment].show_information();
  std::cout<<" FRAGMOL: number of fragmesnt = "<<__number_of_fragments<<std::endl;
#endif
  bool is_new_frag=true;
  while(is_new_frag){
    // the new fragmented part should be updated
    // there is a bug if the initial structure
    // is moved and then fragmented
    // see water case:/home/etorres/src/utils/xmol/test/xyz/WaterWater
    // Use the atom number inside the fragment
    is_new_frag=eval_scaled_fragment(0,false,_s);
    ////////////////set_map_active_fragment(__active_fragment);
    v_fragments[__active_fragment].is_initialized(false);
    v_fragments[__number_of_fragments-1].eval_initial_position();
    v_fragments[__number_of_fragments-1].eval_initial_orientation();
    v_fragments[__number_of_fragments-1].compute_origin_cartesian();
  }
  eval_cell_table();
  update_fragmol_cartesian();
  update_fragmol_direct();
  //////////////initialize_map();
  //std::cout<<" FRAGMOL: m_cartesian: "<<m_xyz;
}

// auto search for fragments separated by van der Waals radius distances
void CSupercell::eval_vdw_fragments(void){
  //std::cout<<" FRAGMOL: m_cartesian: "<<m_xyz;
#ifdef _FRAGMOL_DATA_MESSAGES_
  v_fragments[__active_fragment].show_information();
  std::cout<<" FRAGMOL: number of fragmesnt = "<<__number_of_fragments<<std::endl;
#endif
  bool new_frag=true;
  while(new_frag){
    // the new fragmented part should be updated
    // there is a bug if the initial structure
    // is moved and then fragmented
    // see water case:/home/etorres/src/utils/xmol/test/xyz/WaterWater
    // Use the atom number inside the fragment
    new_frag=eval_scaled_fragment(0,false,1.1);
    ////////////////set_map_active_fragment(__active_fragment);
    v_fragments[__active_fragment].is_initialized(false);
    v_fragments[__number_of_fragments-1].eval_initial_position();
    v_fragments[__number_of_fragments-1].eval_initial_orientation();
    v_fragments[__number_of_fragments-1].compute_origin_cartesian();
  }
  eval_cell_table();
  update_fragmol_cartesian();
  update_fragmol_direct();
  //////////////initialize_map();
  //std::cout<<" FRAGMOL: m_cartesian: "<<m_xyz;
}

void CSupercell::eval_atom_fragments(void){
  //std::cout<<" FRAGMOL: m_cartesian: "<<m_xyz;
#ifdef _FRAGMOL_DATA_MESSAGES_
  v_fragments[__active_fragment].show_information();
  std::cout<<" FRAGMOL: number of fragmesnt = "<<__number_of_fragments<<std::endl;
#endif
  bool new_frag=true;
  while(new_frag){
    // the new fragmented part should be updated
    // there is a bug if the initial structure
    // is moved and then fragmented
    // see water case:/home/etorres/src/utils/xmol/test/xyz/WaterWater
    // Use the atom number inside the fragment
    new_frag=eval_scaled_fragment(0,false,1.1);
    ////////////////set_map_active_fragment(__active_fragment);
    v_fragments[__active_fragment].is_initialized(false);
    v_fragments[__number_of_fragments-1].eval_initial_position();
    v_fragments[__number_of_fragments-1].eval_initial_orientation();
    v_fragments[__number_of_fragments-1].compute_origin_cartesian();
  }
  eval_cell_table();
  update_fragmol_cartesian();
  update_fragmol_direct();
  //////////////initialize_map();
  //std::cout<<" FRAGMOL: m_cartesian: "<<m_xyz;
}

// merge the fragment with the atom (u) with the active fragment
bool CSupercell::eval_merge_fragment(uint u, bool sw){
  TVector<uint> fragment_atoms;
  uint atom_seed, atom_fragment;
  if(sw)
    atom_seed=v_atom_cell_table[u]; // atom number inside the fracment
  else
    atom_seed=u; // atom number from the atom table
  // get the fragment number for a given atom from the table
  atom_fragment=v_fragment_table[atom_seed]-1;
  #ifdef _FRAGMOL_DEBUG_MESSAGES_
  std::cout<<" selected atom = "<<u<<std::endl;
  std::cout<<" atom seed = "<<atom_seed<<std::endl;
  std::cout<<" atom fragment = "<<atom_fragment<<std::endl;
#endif
  if( atom_fragment != __active_fragment ){
    fragment_atoms=gsf.get_topology_atoms(atom_fragment);
#ifdef _FRAGMOL_DEBUG_MESSAGES_
    std::cout<<" FRAGMOL: actual fragment table = "<<v_fragment_table;
    std::cout<<" FRAGMOL: fragment atoms = "<<fragment_atoms;
#endif
    // Update the cartesian coordinates
    set_cartesian();
    // Add the list of merged atom to the active topology
    gsf.add_topmol_atoms(fragment_atoms,__active_fragment);
    v_fragments[__active_fragment].size(gsf.get_topology_size(__active_fragment));
    // Decreae the fragment list by one.
    __number_of_fragments--;
    // Remove the merged frament and topology
    v_fragments.remove(atom_fragment);
    gsf.remove_topmol_topology(atom_fragment);
    // re-evaluate the fragments using the new topologies
    eval_initial_fragments();
#ifdef _FRAGMOL_DEBUG_MESSAGES_
    std::cout<<" FRAGMOL: new number of fragments = "<<get_fragmol_number_of_fragments()<<std::endl;;
#endif
    // Set the active fragment to the fragment containing the merged atom
    if( __active_fragment >= __number_of_fragments )
      __active_fragment=atom_fragment;
#ifdef _FRAGMOL_DEBUG_MESSAGES_
    fragment_atoms=gsf.get_topology_atoms(__active_fragment);
    std::cout<<" FRAGMOL: new active topology atoms = "<<fragment_atoms;
#endif
    return true;
  }else{
    return false;
  }
}

void CSupercell::compute_fragmol_cartesian(uint u){
  TMatrix<real> _U = get_unit_uvw_to_xyz();
  TMatrix<real> _T = get_uvw_to_xyz();
  //for(uint i=0;i<__number_of_fragments;i++){
    // this funtion should be declared static
#ifdef _FRAGMOL_DEBUG_MESSAGES_
    std::cout<<" -----------------------------------------"<<std::endl;
    std::cout<<" FRAGMOL: eval cartesian for fragment "<<u<<" [begin]"<<std::endl;
#endif
    v_fragments[u].set_unit_uvw_to_xyz(_U);
    v_fragments[u].set_uvw_to_xyz(_T);
    v_fragments[u].compute_origin_cartesian();
#ifdef _FRAGMOL_DEBUG_MESSAGES_
    std::cout<<" FRAGMOL: eval cartesian for fragment "<<u<<" [end]"<<std::endl;
    std::cout<<" -----------------------------------------"<<std::endl;
#endif
  //}
}

void CSupercell::compute_fragmol_direct(uint u){
  TMatrix<real> _U = get_unit_uvw_to_xyz();
  TMatrix<real> _T = get_uvw_to_xyz();
  //for(uint i=0;i<__number_of_fragments;i++){
    // this function should be declared static
#ifdef _FRAGMOL_DEBUG_MESSAGES_
    std::cout<<" -----------------------------------------"<<std::endl;
    std::cout<<" FRAGMOL: eval direct for fragment "<<u<<" [begin]"<<std::endl;
#endif
    v_fragments[u].set_unit_uvw_to_xyz(_U);
    v_fragments[u].set_uvw_to_xyz(_T);
    v_fragments[u].compute_origin_direct();
    v_fragments[u].compute_origin_cartesian();
#ifdef _FRAGMOL_DEBUG_MESSAGES_
    std::cout<<" FRAGMOL: eval direct for fragment "<<u<<" [end]"<<std::endl;
    std::cout<<" -----------------------------------------"<<std::endl;
#endif
  //}
}

void CSupercell::compute_fragmol_all_cartesian(void){
  TMatrix<real> _U = get_unit_uvw_to_xyz();
  TMatrix<real> _T = get_uvw_to_xyz();
  for(uint i=0;i<__number_of_fragments;i++){
    // this funtion should be declared static
    compute_fragmol_cartesian(i);
#ifdef _FRAGMOL_DEBUG_MESSAGES_
    std::cout<<" -----------------------------------------"<<std::endl;
    std::cout<<" FRAGMOL: eval cartesian for fragment "<<i<<" [begin]"<<std::endl;
#endif
//    v_fragments[i].set_unit_uvw_to_xyz(_U);
//    v_fragments[i].set_uvw_to_xyz(_T);
//    v_fragments[i].compute_origin_cartesian();
#ifdef _FRAGMOL_DEBUG_MESSAGES_
    std::cout<<" FRAGMOL: cartesian for fragment "<<i<<" [end]"<<std::endl;
    std::cout<<" -----------------------------------------"<<std::endl;
#endif
  }
}

void CSupercell::compute_fragmol_all_direct(void){
  TMatrix<real> _U = get_unit_uvw_to_xyz();
  TMatrix<real> _T = get_uvw_to_xyz();
  for(uint i=0;i<__number_of_fragments;i++){
    // this function should be declared static
    compute_fragmol_direct(i);
#ifdef _FRAGMOL_DEBUG_MESSAGES_
    std::cout<<" -----------------------------------------"<<std::endl;
    std::cout<<" FRAGMOL: eval direct for fragment "<<i<<" [begin]"<<std::endl;
#endif
//    v_fragments[i].set_unit_uvw_to_xyz(_U);
//    v_fragments[i].set_uvw_to_xyz(_T);
//    v_fragments[i].compute_origin_direct();
//    v_fragments[i].compute_origin_cartesian();
#ifdef _FRAGMOL_DEBUG_MESSAGES_
    std::cout<<" FRAGMOL: direct for fragment "<<i<<" [end]"<<std::endl;
    std::cout<<" -----------------------------------------"<<std::endl;
#endif
  }
}

void CSupercell::eval_fragmol_initial_position(void){
#ifdef _FRAGMOL_DEBUG_MESSAGES_
  std::cout<<" FRAGMOL: eval_fragmol_initial_position"<<std::endl;
#endif
  v_fragments[__active_fragment].eval_initial_position();
}

void CSupercell::eval_fragmol_initial_orientation(void){
#ifdef _FRAGMOL_DEBUG_MESSAGES_
  std::cout<<" FRAGMOL: eval_fragmol_initial_orientation"<<std::endl;
#endif
  v_fragments[__active_fragment].eval_initial_orientation();
}

void CSupercell::set_cartesian(void){
  TVector<real> _v;
  TVector<unsigned int> v_l;
  unsigned int _n, _s;
  _n = get_fragmol_number_of_fragments();
#ifdef _FRAGMOL_DEBUG_MESSAGES_
  std::cout<<"FRAGMOL: TOTAL FRAGMENTS: "<<_n<<std::endl;
#endif
  for(unsigned int i=0; i<_n; i++){
    _s = get_fragment_size(i);
    v_l = gsf.get_topology_atoms(i);
    for(unsigned int j=0;j<_s;j++){
      _v = get_fragment_cartesian(i,j);
      // put each atomic coordinate in the right place inside the matrix
      m_xyz[v_l[j]]=_v;
    }
  }
#ifdef _FRAGMOL_FORMAT_MESSAGES_
  std::cout<<" FRAGMOL: Cartesian coordinates ready"<<std::endl;
#endif
#ifdef _FRAGMOL_DATA_MESSAGES_
  std::cout<<" FRAGMOL: m_cartesian: "<<m_xyz;
#endif
}

void CSupercell::update_fragmol_cartesian(void){
  TVector<real> _v;
  TVector<unsigned int> v_l;
  unsigned int _n, _s;
  _n = get_fragmol_number_of_fragments();
#ifdef _FRAGMOL_DEBUG_MESSAGES_
  std::cout<<"FRAGMOL: TOTAL FRAGMENTS: "<<_n<<std::endl;
#endif
  for(unsigned int i=0; i<_n; i++){
    _s = get_fragment_size(i);
    v_l = gsf.get_topology_atoms(i);
    for(unsigned int j=0;j<_s;j++){
      _v = get_fragment_centered_cartesian(i,j);
      // put each atomic coordinate in the right place inside the matrix
      m_xyz[v_l[j]]=_v;
    }
  }
#ifdef _FRAGMOL_FORMAT_MESSAGES_
  std::cout<<" FRAGMOL: Cartesian coordinates ready"<<std::endl;
#endif
#ifdef _FRAGMOL_DATA_MESSAGES_
  std::cout<<" FRAGMOL: centered cartesian: "<<m_xyz;
#endif
}

void CSupercell::update_fragmol_direct(void){
  TVector<real> _v;
  TVector<unsigned int> v_l;
  unsigned int _n, _s;
  _n = get_fragmol_number_of_fragments();
#ifdef _FRAGMOL_DEBUG_MESSAGES_
  std::cout<<"FRAGMOL: TOTAL FRAGMENTS: "<<_n<<std::endl;
#endif
  for(unsigned int i=0; i<_n; i++){
    _s = get_fragment_size(i);
    v_l = gsf.get_topology_atoms(i);
    for(unsigned int j=0;j<_s;j++){
      _v = get_fragment_direct(i,j);
      // put each atomic coordinate in the right place inside the matrix
      m_uvw[v_l[j]]=_v;
    }
  }
#ifdef _FRAGMOL_FORMAT_MESSAGES_
  std::cout<<" FRAGMOL: Direct coordinates ready"<<std::endl;
#endif
#ifdef _FRAGMOL_DATA_MESSAGES_
  std::cout<<" FRAGMOL: m_direct: "<<m_uvw;
#endif
}

void CSupercell::compute_fragmol_position_cartesian(void){
  v_fragments[__active_fragment].compute_origin_cartesian();
}

void CSupercell::is_fragmol_initialized(bool b){
  v_fragments[__active_fragment].is_initialized(b);
}

void CSupercell::apply_pbc(bool b){
  v_fragments[__active_fragment].is_pbc(b);
}

void CSupercell::set_gsf_modified(const bool b){
#ifdef _FRAGMOL_DEBUG_MESSAGES_
  std::cout<<" FRAGMOL: the gsf was modified"<<std::endl;
#endif
  gsf.set_modified(b);
}

void CSupercell::set_fragmol_active_fragment(const uint i){
  __active_fragment=i;
#ifdef _FRAGMOL_DEBUG_MESSAGES_
  std::cout<<" FRAGMOL: active fragment= "<<__active_fragment<<std::endl;
#endif
}

void CSupercell::set_input_file_type(const uint u){
  //__input_type=u;
  gsf.set_input_type(u);
}

void CSupercell::set_input_file_units(const uint u){
  //__input_type=u;
  gsf.set_input_units(u);
}

void CSupercell::set_output_file_type(const uint u){
  //__output_file_type=u;
  gsf.set_output_file_type(u);
}

void CSupercell::set_output_file_format(const uint u){
  __output_format=u;
  gsf.set_output_file_format(u);
}

void CSupercell::set_export_format(const uint u){
  //__export_format=u;
//#ifdef _FRAGMOL_DEBUG_MESSAGES_
//  std::cout<<" FRAGMOL: export format = "<<__export_format<<std::endl;
//#endif
  gsf.set_export_format(u);
}

void CSupercell::set_fragment_axis(const TVector<uint>& _v){
  is_fragmol_initialized(false);
  v_fragments[__active_fragment].set_axis_index(_v);
}

void CSupercell::set_fragment_axis(const uint idx, const uint val){
  is_fragmol_initialized(false);
  v_fragments[__active_fragment].set_axis_index(idx,v_atom_cell_table[val]);
}

void CSupercell::set_fragmol_fragment_twist(const real r){
  v_fragments[__active_fragment].set_axis_twist(r);
}

void CSupercell::set_fragmol_fragment_precession(const real r){
  v_fragments[__active_fragment].set_axis_precession(r);
}

void CSupercell::set_fragmol_fragment_tilt(const real r){
  v_fragments[__active_fragment].set_axis_tilt(r);
}

void CSupercell::set_fragmol_fragment_position_u(const real r){
  v_fragments[__active_fragment].set_origin_u(r);
}

void CSupercell::set_fragmol_fragment_position_v(const real r){
  v_fragments[__active_fragment].set_origin_v(r);
}

void CSupercell::set_fragmol_fragment_position_w(const real r){
  v_fragments[__active_fragment].set_origin_w(r);
}

void CSupercell::set_input_file(std::string s){
  //inputfile = s;
  gsf.set_input_file(s);
  gsf.topmol_filename(s);
}

//void CSupercell::set_output_file_name(std::string s){
  //output_filename = s;
//  gsf.set_output_file_name(s);
//}

void CSupercell::set_dir(std::string s){
  // use__sdir in case you would like to change the default directory
  __sdir = s;
  gsf.set_input_dir(s);
}


void CSupercell::set_topmol_directory(std::string s){
  // use__sdir in case you would like to change the default directory
  //__sdir = s;
  gsf.topmol_dir(s);
}

void CSupercell::save_topmol_file(std::string _d,std::string _f){
    gsf.save_topmol(_d,_f);
}

void CSupercell::set_fragment_direct(const uint i, const TMatrix<real>& _m){
  v_fragments[i].set_position_direct(_m);
}
void CSupercell::set_fragment_cartesian(const uint i, const TMatrix<real>& _m){
  v_fragments[i].set_position_cartesian(_m);
}

real CSupercell::get_fragmol_backbone_tilt(void){
  return v_fragments[__active_fragment].get_backbone_tilt();
}

real CSupercell::get_fragmol_backbone_precession(void){
  return v_fragments[__active_fragment].get_backbone_precession();
}

real CSupercell::get_fragmol_axis_tilt(void){
  return v_fragments[__active_fragment].get_axis_tilt();
}

real CSupercell::get_fragmol_axis_precession(void){
  return v_fragments[__active_fragment].get_axis_precession();
}

bool CSupercell::is_direct(void){
  return __is_direct;
}

bool CSupercell::is_periodic(void){
  return __is_periodic;
}

void CSupercell::is_periodic(bool b){
  __is_periodic = b;
}

bool CSupercell::is_potmol(void){
  return __is_potmol;
}

bool CSupercell::is_fragmol_initialized(void){
  return v_fragments[__active_fragment].is_initialized();
}

std::string CSupercell::get_dir(void){
  return __sdir;
}

std::string CSupercell::get_fragment_atomic_label(uint i, uint j){
  return v_fragments[i].get_atomic_label(j);
}

std::string CSupercell::get_fragment_atomic_symbol(uint i, uint j){
  return v_fragments[i].get_atomic_symbol(j);
}

uint CSupercell::get_fragmol_number_of_fragments(void){
  return __number_of_fragments;
}

uint CSupercell::get_fragmol_active_fragment(void){
  return __active_fragment;
}

uint CSupercell::get_fragment_size(uint i){
  return v_fragments[i].size();
}

uint CSupercell::get_fragment_atomic_number(uint i,uint j){
  return v_fragments[i].get_atomic_number(j);
}

uint CSupercell::get_fragmol_total_atoms(void){
  return __total_atoms;
}

TVector<real> CSupercell::get_fragmol_axis_angles(void){
  return v_fragments[__active_fragment].get_axis_angles();
}

TVector<real> CSupercell::get_fragmol_basis_direct(void){
  return v_fragments[__active_fragment].get_basis_direct();
}

TVector<real> CSupercell::get_fragmol_position_direct(void){
  return v_fragments[__active_fragment].get_origin_direct();
}

TVector<real> CSupercell::get_fragmol_position_uvw(void){
  return v_fragments[__active_fragment].get_origin_uvw();
}

TVector<real> CSupercell::get_fragmol_position_cartesian(void){
  return v_fragments[__active_fragment].get_origin_cartesian();
}

TVector<real> CSupercell::get_fragmol_centered_position_cartesian(void){
  return v_fragments[__active_fragment].get_centered_origin_cartesian();
}

TVector<real> CSupercell::get_fragment_direct(uint i,uint j){
  return v_fragments[i].get_direct(j);
}

TVector<real> CSupercell::get_fragment_cartesian(uint i,uint j){
  return v_fragments[i].get_cartesian(j);
  //return v_fragments[i].get_centered_cartesian(j);
}

TVector<real> CSupercell::get_fragment_centered_cartesian(uint i,uint j){
  return v_fragments[i].get_centered_cartesian(j);
}

TMatrix<real> CSupercell::get_fragment_direct(uint i){
  return v_fragments[i].get_direct();
}

TMatrix<real> CSupercell::get_fragment_cartesian(uint i){
  return v_fragments[i].get_cartesian();
}

TMatrix<real> CSupercell::get_fragment_centered_cartesian(uint i){
  return v_fragments[i].get_centered_cartesian();
}

// STRUCTURE

uint CSupercell::get_species(void){
  return __atomic_species;
}

uint CSupercell::get_total_atoms(void){
  return __total_atoms;
}

uint CSupercell::get_input_file_format(void){
  return __input_format;
}

uint CSupercell::get_output_file_type(void){
  return gsf.get_output_file_type();
}

uint CSupercell::get_output_file_format(void){
  return gsf.get_output_file_format();
}

uint CSupercell::get_atomic_composition(uint i){
  return v_atomic_composition_table[i];
}

TVector<real> CSupercell::get_cartesian(uint i){
  return m_xyz[i];
}

TVector<real> CSupercell::get_direct(uint i){
  return m_uvw[i];
}

real CSupercell::get_distance(uint i, uint j){
  //TVector<real> v1, v2, v3;
  TVector<real> v3;
  //v1 = get_cartesian(i);
  //v2 = get_cartesian(j);
  // v3 = v2-v1;
  v3 = get_vector_diff(i,j);
  return v3.magnitude();
}

TVector<real> CSupercell::get_vector_diff(uint i, uint j){
  TVector<real> v1, v2, v3;
  v1 = get_cartesian(i);
  v2 = get_cartesian(j);
  v3 = v2-v1;
  return v3;
}

real CSupercell::get_angle(uint i, uint j, uint k){
  TVector<real> v1, v2, v3;
  v1 = get_cartesian(i);
  v2 = get_cartesian(j);
  v3 = get_cartesian(k);
  v1 = (v2-v1);
  v3 = (v2-v3);
  real res = (v1*v3)/(v1.magnitude()*v3.magnitude());
  res = acos(res)*RAD_DEG;
  return res;
}

real CSupercell::get_dihedral(uint i, uint j, uint k, uint l){
  TVector<real> v1, v2, v3, v4, vd;
  real d12, d13, d14, d23, d24, d34;
  v1 = get_cartesian(i);
  v2 = get_cartesian(j);
  v3 = get_cartesian(k);
  v4 = get_cartesian(l);
  vd = v2-v1;
  d12 = vd.magnitude();
  vd = v3-v1;
  d13 = vd.magnitude();
  vd = v4-v1;
  d14 = vd.magnitude();
  vd = v3-v2;
  d23 = vd.magnitude();
  vd = v4-v2;
  d24 = vd.magnitude();
  vd = v4-v3;
  d34 = vd.magnitude();
  real P = SQ(d12) * ( SQ(d23)+SQ(d34)-SQ(d24)) +
            SQ(d23) * (-SQ(d23)+SQ(d34)+SQ(d24)) +
            SQ(d13) * ( SQ(d23)-SQ(d34)+SQ(d24)) -
            2 * SQ(d23) * SQ(d14);
  real Q = (d12 + d23 + d13) * ( d12 + d23 - d13) *
            (d12 - d23 + d13) * (-d12 + d23 + d13 ) *
            (d23 + d34 + d24) * ( d23 + d34 - d24 ) *
            (d23 - d34 + d24) * (-d23 + d34 + d24 );
  //v1 = (v1-v2);
  //v3 = (v4-v3);
  real res = P/sqrt(Q);
  //real res = (v1*v3)/(v1.magnitude()*v3.magnitude());
  res = acos(res)*RAD_DEG;
  return res;
}

//////
TVector<std::string> CSupercell::get_fragmol_atomic_symbol_table(void){
  return v_atomic_symbol_table;
}

TVector<uint> CSupercell::get_fragmol_atomic_composition_table(void){
  return v_atomic_composition_table;
}

TVector<uint> CSupercell::get_fragmol_atomic_number_table(void){
  return v_atomic_number_table;
}

TVector<uint> CSupercell::get_fragmol_atom_table(void){
  return v_atom_type_table;
}
//////

TVector<uint> CSupercell::get_fragmol_fragment_table(void){
  return v_fragment_table;
}

TMatrix<real> CSupercell::get_unit_uvw_to_xyz(void){
  return m_uvw_to_xyz_u;
}

TMatrix<real> CSupercell::get_uvw_to_xyz(void){
  return m_uvw_to_xyz;
}

TMatrix<real> CSupercell::get_bounding_box(void){
  return m_uvw_to_xyz;
}

TMatrix<real> CSupercell::get_cartesian(void){
  return m_xyz;
}

TMatrix<real> CSupercell::get_direct(void){
  return m_uvw;
}

/* Deprecated
bool CSupercell::eval_vdw_fragment(uint u){
  //std::cout<<" FRAGMOL: m_cartesian: "<<m_xyz;
#ifdef _FRAGMOL_DATA_MESSAGES_
  v_fragments[__active_fragment].show_information();
#endif
  bool is_new_frag=true;
  // Use the atom number
  is_new_frag=eval_scaled_fragment(u,true,1.1);
  /////////////////////////set_map_active_fragment(__active_fragment);
  eval_cell_table();
  update_fragmol_cartesian();
  update_fragmol_direct();
  /////////////////////initialize_map();
  //std::cout<<" FRAGMOL: m_cartesian: "<<m_xyz;
  return is_new_frag;
}
  Deprecated
bool CSupercell::eval_atom_fragment(uint u){
  //std::cout<<" FRAGMOL: m_cartesian: "<<m_xyz;
#ifdef _FRAGMOL_DATA_MESSAGES_
  v_fragments[__active_fragment].show_information();
#endif
  bool new_frag=true;
  new_frag=eval_scaled_fragment(u,true,0.1);
  /////////////////////////set_map_active_fragment(__active_fragment);
  eval_cell_table();
  update_fragmol_cartesian();
  update_fragmol_direct();
  /////////////////////initialize_map();
  //std::cout<<" FRAGMOL: m_cartesian: "<<m_xyz;
  return new_frag;
}
*/

// FRAGMOL END
