//========================================================================
// FILE: cread.cxx
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

#include <config_debug.h>
#include <cfile.h>

//#define _DEBUGING_MESSAGES_

CFile::CFile(){
  potcarfile = "POTCAR"; // default potcar file
  s_actual_dir  = "./";     // default directory
  clear();
  // 0: standart format
  u_input_file_type=0;
  u_input_file_units=0;
  u_input_format=0;
  u_output_file_format=OUTPUT_FORMAT_STD;
  u_export_format=0;
}

void CFile::clear(void){
  //clear_xyz();
  //clear_poscar();
  //clear_potcar();
  //__number_of_fragments=0;
  //__total_atoms=0;
  is_periodic=false;
  is_potcar=false;
  is_bounding_box=false;
  is_labels=false;
  is_numbers=false;
  is_charges=false;
  is_modified=false;
  v_connections.clear();
}

void CFile::clear_xyz(void){
  file_xyz.clear();
}

void CFile::clear_gau(void){
  file_gau.clear();
}

void CFile::clear_pdb(void){
  file_pdb.clear();
}

void CFile::clear_dlp(void){
  file_dlp.clear();
}

void CFile::clear_poscar(void){
  file_poscar.clear();
}

void CFile::clear_potcar(void){
  file_potcar.clear();
}

bool CFile::read_potcar(void){
  clear_potcar();
  //bool _res;
  is_potcar=file_potcar.read_potcar(s_actual_dir, potcarfile);
  if(is_potcar)
    v_atomic_symbol_table=file_potcar.get_symbols();
  //v_atomic_numbers=file_potcar.get_atomic_numbers();
  return is_potcar;
}

bool CFile::read_input_file(void){
  bool res=false;
  clear();
#ifdef _DEBUGING_MESSAGES_
  std::cout<<" FILE: Read input file type "<<u_input_file_type<<"  "<<std::endl;
#endif
  try{
    switch(u_input_file_type){
        case INPUT_FILE_TYPE_VSP: // POSCAR
        res=read_poscar();
        break;
      case INPUT_FILE_TYPE_XYZ: // XYZ
        res=read_xyz();
        break;
      case INPUT_FILE_TYPE_GAU: // Gaussian
        res=read_gau();
        break;
      case INPUT_FILE_TYPE_PDB: // PDB
         res=read_pdb();
        break;
      case INPUT_FILE_TYPE_DLP: // DL_POLY
        res=read_dlp();
        break;
      case INPUT_FILE_TYPE_ZMT: // Gaussian ZMat
        res=read_zmt();
        break;
      }
  }catch(...){
    std::cout<<"Hi";
  }
#ifdef _FILE_DEBUGING_MESSAGES_
  if(!res)
    std::cout<<" FILE: Wrong file format"<<std::endl;
  else
    std::cout<<" FILE: Correct file format"<<std::endl;
#endif
  return res;
}

void CFile::save_input_file(void){
  switch(u_input_file_type){
    case OUTPUT_FILE_TYPE_VSP:
      //save_poscar();
      save_xyz();
    break;
    case OUTPUT_FILE_TYPE_XYZ:
      save_xyz();
    break;
    case OUTPUT_FILE_TYPE_GAU:
      save_xyz();
    break;
    case OUTPUT_FILE_TYPE_PDB:
      save_xyz();
    break;
    case OUTPUT_FILE_TYPE_DLP:
      save_xyz();
    break;
  }
}


void CFile::save_as_file(std::string _f){
  save_as_file(s_actual_dir,_f);
  set_topmol_filename(s_actual_dir,_f);
#ifdef _debugging_messages_
  std::cout<<"save_as_file(_f) "<<_f<<std::endl;
#endif
}

void CFile::save_as_file(std::string _p, std::string _f){
  set_topmol_filename(_p,_f);
#ifdef _debugging_messages_
  std::cout<<"save_as_file(_p,_f) "<<_p<<"  "<<_f<<std::endl;
#endif
  switch(u_output_file_type){
    case OUTPUT_FILE_TYPE_VSP:
      save_poscar_as(_p,_f,VASP_FORMAT_CARTES);
    break;
    case OUTPUT_FILE_TYPE_XYZ:
      save_xyz_as(_p,_f,u_output_file_format);
    break;
    case OUTPUT_FILE_TYPE_GAU:
      save_gau_as(_p,_f,u_output_file_format);
    break;
    case OUTPUT_FILE_TYPE_PDB: // reserved for xmatrix
      save_pdb_as(_p,_f,u_output_file_format);
    break;
    case OUTPUT_FILE_TYPE_DLP:
      save_dlp_as(_p,_f,u_output_file_format);
    break;
    case OUTPUT_FILE_TYPE_ZMT:
      save_zmt_as(_p,_f,u_output_file_format);
    break;
  }
}

bool CFile::read_poscar(void){
  TVector<uint> v_n;
  TVector<std::string> v_s;
  uint _c=0; // _s=0;
  bool res;
  clear_poscar();
  res=file_poscar.read_file(s_actual_dir,inputfile);
  if(res){
    u_total_atoms=file_poscar.get_total_atoms();
    u_atomic_species=file_poscar.get_total_species();
    v_atomic_composition_table=file_poscar.get_composition();
    //
    v_atomic_numbers.resize(u_total_atoms);
    v_atomic_symbols.resize(u_total_atoms);
    v_atom_table.resize(u_total_atoms);
    if(file_poscar.is_potcar()){
      file_potcar.set_symbols(file_poscar.get_atomic_symbols());
      file_potcar.eval_atomic_numbers();
    }
    //////////////////////////////////////////////////////////////////////////
    // Get POTCAR information
    v_atomic_symbol_table = file_potcar.get_symbols();
    v_atomic_number_table = file_potcar.get_atomic_numbers();
    //////////////////////////////////////////////////////////////////////////
#ifdef _FILE_DEBUGING_MESSAGES_
    std::cout<<" FILE: potcar atomic symbols: "<<v_atomic_symbol_table;
    std::cout<<" FILE: potcar atomic numbers: "<<v_atomic_number_table;
#endif
    u_input_format=file_poscar.get_format();
    // This is only for POSCAR files
    for(uint i=0;i<u_atomic_species;i++){
      for(uint j=0;j<v_atomic_composition_table[i];j++){
         v_atomic_numbers[_c]=v_atomic_number_table[i];
         v_atomic_symbols[_c]=v_atomic_symbol_table[i];
         v_atom_table[_c]=i;
         _c++;
      }
    }
    v_atomic_labels=v_atomic_symbols;
    is_direct=file_poscar.get_format();
    is_periodic=file_poscar.is_periodic();
    //
    m_file_input=file_poscar.get_xyz_input();
    m_xyz=file_poscar.get_cartesian();
    m_uvw=file_poscar.get_direct();
    m_uvw_to_xyz_u=file_poscar.get_unit_uvw_to_xyz();
    m_uvw_to_xyz=file_poscar.get_uvw_to_xyz();
  }
  return res;
}

bool CFile::read_xyz(void){
  bool res;
  clear_xyz();
  // depreacted code: Wed Jan 18 10:07:27 MST 2012
  // CXyz is more robust and can manage this situation
  // this must be done before read the xyz file
  // in order to know the species in the potcar file
  //if(is_potcar){ // no necesary any longer
    //file_xyz.set_file_format(u_input_format);
    //file_xyz.set_atomic_symbol_table(get_potcar_symbols());
  //}
  // read the xyz file
  // depreacted code: Wed Jan 18 10:07:27 MST 2012
  // not necessary, CXyz autodetect the format
  //file_xyz.set_file_format(u_input_file_type);
  file_xyz.set_input_units(u_input_file_units);
  res=file_xyz.read_file(s_actual_dir,inputfile);
  if(res){
    u_total_atoms=file_xyz.get_total_atoms();
    u_atomic_species=file_xyz.get_total_species();
    v_atomic_composition_table=file_xyz.get_atomic_composition_table();
    v_atomic_number_table=file_xyz.get_atomic_number_table();
    v_atomic_symbol_table=file_xyz.get_atomic_symbol_table();
    v_fragment_table=file_xyz.get_fragment_table();
    v_atom_table=file_xyz.get_atom_table();
    v_atomic_symbols=file_xyz.get_atomic_symbols();
    //v_atomic_labels=v_atomic_symbols;
    v_atomic_labels=file_xyz.get_atomic_labels();
    v_atomic_numbers=file_xyz.get_atomic_numbers();
#ifdef _FILE_DEBUGING_MESSAGES_
    std::cout<<" READ: atomic symbols: "<<v_atomic_symbols<<std::endl;
    std::cout<<" READ: atomic labels: "<<v_atomic_labels<<std::endl;
    std::cout<<" READ: atomic symbol table: "<<v_atomic_symbol_table<<std::endl;
    std::cout<<" READ: atomic numbers: "<<v_atomic_numbers<<std::endl;
    std::cout<<" READ: atomic number table: "<<v_atomic_number_table<<std::endl;
#endif
    is_direct=file_xyz.is_direct();
    if(file_xyz.is_charges()){
      is_charges=true;
      v_atomic_charges=file_xyz.get_charges();
    }else{
      is_charges=false;
    }
    //if(file_xyz.is_fragments()){
    //v_fragment_table=file_xyz.get_fragment_table();
    //}
    is_periodic=file_xyz.is_periodic();
    //// adding format autodetection
    u_input_format=file_xyz.get_format();
#ifdef _FILE_DEBUGING_MESSAGES_
    std::cout<<" READ: Reading file of format "<<u_input_format<<"  "<<std::endl;
#endif
    ////
    //  if(u_input_file_units==0)
    //    m_xyz=file_xyz.get_cartesian();
    //  else
    //    m_xyz=0.52917720859*file_xyz.get_cartesian();
    //
    m_file_input=file_xyz.get_xyz_input();
    m_xyz=file_xyz.get_cartesian();
    m_uvw=file_xyz.get_direct();
    m_uvw_to_xyz_u=file_xyz.get_unit_uvw_to_xyz();
    m_uvw_to_xyz=file_xyz.get_uvw_to_xyz();
  }
  return res;
}

bool CFile::read_gau(void){
  bool res;
  clear_gau();
  res=file_gau.read_file(s_actual_dir,inputfile);
  if(res){
    u_total_atoms=file_gau.get_total_atoms();
    u_atomic_species=file_gau.get_total_species();
    v_atomic_composition_table=file_gau.get_atomic_composition_table();
    v_atomic_number_table=file_gau.get_atomic_number_table();
    v_atomic_symbol_table=file_gau.get_atomic_symbol_table();
    v_fragment_table=file_gau.get_fragment_table();
    v_atom_table=file_gau.get_atom_table();
    v_atomic_symbols=file_gau.get_atomic_symbols();
    v_atomic_labels=v_atomic_symbols;
    v_atomic_numbers=file_gau.get_atomic_numbers();
    is_direct=file_gau.is_direct();
    is_periodic=file_gau.is_periodic();
    //// adding format autodetection
    u_input_format=file_gau.get_format();
#ifdef _FILE_DEBUGING_MESSAGES_
    std::cout<<" READ: Reading file of format "<<u_input_format<<"  "<<std::endl;
#endif
    ////
    m_file_input=file_gau.get_xyz_input();
    m_xyz=file_gau.get_cartesian();
    m_uvw=file_gau.get_direct();
    m_uvw_to_xyz_u=file_gau.get_unit_uvw_to_xyz();
    m_uvw_to_xyz=file_gau.get_uvw_to_xyz();
  }
  return res;
}

bool CFile::read_zmt(void){
  bool res;
  //clear_gau();
  res=file_zmt.read_file(s_actual_dir,inputfile);
  if(res){
    u_total_atoms=file_zmt.get_total_atoms();
    u_atomic_species=file_zmt.get_total_species();
    v_atomic_composition_table=file_zmt.get_atomic_composition_table();
    v_atomic_number_table=file_zmt.get_atomic_number_table();
    v_atomic_symbol_table=file_zmt.get_atomic_symbol_table();
    v_fragment_table=file_zmt.get_fragment_table();
    v_atom_table=file_zmt.get_atom_table();
    v_atomic_symbols=file_zmt.get_atomic_symbols();
    v_atomic_labels=v_atomic_symbols;
    v_atomic_numbers=file_zmt.get_atomic_numbers();
    is_direct=file_zmt.is_direct();
    is_periodic=file_zmt.is_periodic();
    //// adding format autodetection
    u_input_format=file_zmt.get_format();
#ifdef _FILE_DEBUGING_MESSAGES_
    std::cout<<" READ: Reading file of format "<<u_input_format<<"  "<<std::endl;
#endif
    ////
    m_file_input=file_zmt.get_xyz_input();
    m_xyz=file_zmt.get_cartesian();
    m_uvw=file_zmt.get_direct();
    m_uvw_to_xyz_u=file_zmt.get_unit_uvw_to_xyz();
    m_uvw_to_xyz=file_zmt.get_uvw_to_xyz();
  }
  return res;
}

bool CFile::read_pdb(void){
  bool res;
  clear_pdb();
  res=file_pdb.read_file(s_actual_dir,inputfile);
  if(res){
    u_total_atoms=file_pdb.get_total_atoms();
    u_atomic_species=file_pdb.get_total_species();
    v_atomic_composition_table=file_pdb.get_atomic_composition_table();
    v_atomic_number_table=file_pdb.get_atomic_number_table();
    v_atomic_symbol_table=file_pdb.get_atomic_symbol_table();
    v_fragment_table=file_pdb.get_fragment_table();
    v_atom_table=file_pdb.get_atom_table();
    v_atomic_symbols=file_pdb.get_atomic_symbols();
    v_atomic_labels=file_pdb.get_atomic_labels();
    v_atomic_numbers=file_pdb.get_atomic_numbers();
    //
#ifdef _FILE_DEBUGING_MESSAGES_
    std::cout<<" READ: atomic symbols: "<<v_atomic_symbols<<std::endl;
    std::cout<<" READ: atomic labels: "<<v_atomic_labels<<std::endl;
    std::cout<<" READ: atomic symbol table: "<<v_atomic_symbol_table<<std::endl;
    std::cout<<" READ: atomic numbers: "<<v_atomic_numbers<<std::endl;
    std::cout<<" READ: atomic number table: "<<v_atomic_number_table<<std::endl;
#endif
    //
    is_direct=file_pdb.is_direct();
    if(file_pdb.is_charges()){
	  is_charges=true;
	  v_atomic_charges=file_pdb.get_charges();
	}else{
	  is_charges=false;
	}
    is_periodic=file_pdb.is_periodic();
    //// adding format autodetection
    u_input_format=file_pdb.get_format();
#ifdef _FILE_DEBUGING_MESSAGES_
    std::cout<<" READ: Reading file of format "<<u_input_format<<"  "<<std::endl;
#endif
    ////
    m_file_input=file_pdb.get_xyz_input();
    m_xyz=file_pdb.get_cartesian();
    m_uvw=file_pdb.get_direct();
    m_uvw_to_xyz_u=file_pdb.get_unit_uvw_to_xyz();
    m_uvw_to_xyz=file_pdb.get_uvw_to_xyz();
  }
  return res;
}

bool CFile::read_dlp(void){
  bool res;
  clear_dlp();
  res=file_dlp.read_file(s_actual_dir,inputfile);
  if(res){
    u_total_atoms=file_dlp.get_total_atoms();
    u_atomic_species=file_dlp.get_total_species();
    v_atomic_composition_table=file_dlp.get_atomic_composition_table();
    v_atomic_number_table=file_dlp.get_atomic_number_table();
    v_atomic_symbol_table=file_dlp.get_atomic_symbol_table();
    v_fragment_table=file_dlp.get_fragment_table();
    v_atom_table=file_dlp.get_atom_table();
    v_atomic_symbols=file_dlp.get_atomic_symbols();
    v_atomic_labels=file_dlp.get_atomic_labels();
    v_atomic_numbers=file_dlp.get_atomic_numbers();
    is_direct=file_dlp.is_direct();
    is_periodic=file_dlp.is_periodic();
    //// adding format autodetection
    u_input_format=file_dlp.get_format();
#ifdef _FILE_DEBUGING_MESSAGES_
    std::cout<<" READ: Reading file of format "<<u_input_format<<"  "<<std::endl;
#endif
    ////
    m_file_input=file_dlp.get_xyz_input();
    m_xyz=file_dlp.get_cartesian();
    m_uvw=file_dlp.get_direct();
    m_uvw_to_xyz_u=file_dlp.get_unit_uvw_to_xyz();
    m_uvw_to_xyz=file_dlp.get_uvw_to_xyz();
  }
  return res;
}

void CFile::save_poscar_as(std::string _p, std::string _f, uint u){
  //clear_xyz();
#ifdef _FILE_DEBUGING_MESSAGES_
  std::cout<<"CFILE: saving poscar"<<std::endl;
#endif
  file_poscar.set_format(u);
  file_poscar.set_total_atoms(u_total_atoms);
  file_poscar.set_sorted(false);
  file_poscar.set_centered(false);
  file_poscar.set_atomic_composition_table(v_atomic_composition_table);
  file_poscar.set_atomic_symbol_table(v_atomic_symbol_table);
  file_poscar.set_atomic_numbers(v_atomic_numbers);
  file_poscar.set_atomic_number_table(v_atomic_number_table);
  file_poscar.set_uvw_to_xyz(m_uvw_to_xyz);
  //if(file_poscar.format())
  file_poscar.set_total_cells(x_cells,y_cells,z_cells,total_cells);
  file_poscar.set_cartesian(m_xyz);
  //else
  file_poscar.set_direct(m_uvw);
  //file_poscar.set_lattice_constant(1.0);
  file_poscar.save_file_as(_p,_f);
}

void CFile::save_xyz(void){
  //clear_xyz();
  file_xyz.set_cartesian(m_xyz);
  file_xyz.set_direct(m_uvw);
  file_xyz.write_file();
}

void CFile::save_xyz_as(std::string _p, std::string _f, uint u){
  //clear_xyz();
  file_xyz.set_export_format(u);
  file_xyz.set_total_atoms(u_total_atoms);
  //std::cout<<"u_export_format="<<u<<std::endl;
  //std::cout<<"v_fragment_table"<<v_fragment_table;
  if(u==OUTPUT_FORMAT_ATM_FRG || u==OUTPUT_FORMAT_NAT_FRG){
    //std::cout<<"v_fragment_table"<<v_fragment_table;
    file_xyz.set_fragments(v_fragment_table);
  }
  file_xyz.set_total_cells(x_cells,y_cells,z_cells,total_cells);
  //file_xyz.set_cartesian(m_xyz);
  if(is_modified)
    file_xyz.set_cartesian(m_xyz);
  else
    file_xyz.set_cartesian(m_file_input);
  //file_xyz.set_atomic_symbols(v_atomic_symbols);
  if(is_labels){
    file_xyz.set_atomic_symbols(v_atomic_symbols);
  }else{
    file_xyz.set_atomic_symbols(v_atomic_labels);
  }
  //file_xyz.set_direct(m_uvw);
  file_xyz.write_file(_p,_f);
}

void CFile::save_gau_as(std::string _p, std::string _f, uint u){
  //clear_xyz();
  file_gau.set_export_format(u);
  file_gau.set_total_atoms(u_total_atoms);
  if(u==OUTPUT_FORMAT_ATM_FRG || u==OUTPUT_FORMAT_NAT_FRG)
    file_gau.set_fragments(v_fragment_table);
  file_gau.set_total_cells(x_cells,y_cells,z_cells,total_cells);
  if(is_modified)
    file_gau.set_cartesian(m_xyz);
  else
    file_gau.set_cartesian(m_file_input);
  file_gau.set_atomic_symbols(v_atomic_symbols);
  //file_xyz.set_direct(m_uvw);
  file_gau.write_file(_p,_f);
}

void CFile::save_pdb_as(std::string _p, std::string _f, uint u){
  //clear_xyz();
  file_pdb.set_export_format(u);
  file_pdb.set_total_atoms(u_total_atoms);
  if(u==OUTPUT_FORMAT_ATM_FRG || u==OUTPUT_FORMAT_NAT_FRG){
    file_pdb.set_fragments(v_fragment_table,u_total_fragments);
  }
  file_pdb.set_total_cells(x_cells,y_cells,z_cells,total_cells);
  //file_pdb.set_cartesian(m_xyz);
  if(is_modified)
    file_pdb.set_cartesian(m_xyz);
  else
    file_pdb.set_cartesian(m_file_input);
  //file_pdb.set_atomic_symbols(v_atomic_symbols);
  if(is_labels){
    file_pdb.set_atomic_symbols(v_atomic_symbols);
  }else{
    file_pdb.set_atomic_symbols(v_atomic_labels);
  }
  if(is_charges){
	file_pdb.set_charges(v_atomic_charges);
  }
  file_pdb.set_atomic_numbers(is_numbers);
  //file_xyz.set_direct(m_uvw);
  file_pdb.write_file(_p,_f);
}

void CFile::save_dlp_as(std::string _p, std::string _f, uint u){
  //clear_xyz();
  file_dlp.set_export_format(u);
  file_dlp.set_labels_format(is_labels);
  file_dlp.set_periodic(is_bounding_box);
  file_dlp.set_total_atoms(u_total_atoms);
  if(u==OUTPUT_FORMAT_ATM_FRG || u==OUTPUT_FORMAT_ATM_NFR){
    file_dlp.set_connections(v_connections);
  }else{
    file_dlp.is_connections(false);
  }
  if(u==OUTPUT_FORMAT_ATM_FRG || u==OUTPUT_FORMAT_NAT_FRG){
    file_dlp.set_fragments(v_fragment_table,u_total_fragments);
    //std::cout<<" using fragments "<<std::endl;
  }
  file_dlp.set_total_cells(x_cells,y_cells,z_cells,total_cells);
  //file_dlp.set_cartesian(m_xyz);
  if(is_modified)
    file_dlp.set_cartesian(m_xyz);
  else
    file_dlp.set_cartesian(m_file_input);
  file_dlp.set_labels_format(is_labels);
  if(is_labels){
    file_dlp.set_atomic_symbols(v_atomic_symbols);
  }else{
    file_dlp.set_atomic_symbols(v_atomic_labels);
  }
  if(is_charges){
	file_dlp.set_charges(v_atomic_charges);
  }
  file_dlp.set_atomic_numbers(is_numbers);
  //file_xyz.set_direct(m_uvw);
  file_dlp.set_uvw_to_xyz(m_uvw_to_xyz);
  file_dlp.write_file(_p,_f);
}

void CFile::save_zmt_as(std::string _p, std::string _f, uint u){
  //clear_xyz();
  file_zmt.set_export_format(u);
  file_zmt.set_total_atoms(u_total_atoms);
  if(u==OUTPUT_FORMAT_ATM_FRG || u==OUTPUT_FORMAT_NAT_FRG)
    file_zmt.set_fragments(v_fragment_table);
  file_zmt.set_total_cells(x_cells,y_cells,z_cells,total_cells);
  //file_zmt.set_cartesian(m_xyz);
  if(is_modified)
    file_zmt.set_cartesian(m_xyz);
  else
    file_zmt.set_cartesian(m_file_input);
  file_zmt.set_atomic_symbols(v_atomic_symbols);
  //file_xyz.set_direct(m_uvw);
  file_zmt.write_file(_p,_f);
}

void CFile::eval_connections(const TMatrix<uint>& _m, uint nb){
  //std::cout<<" eval_connections "<<std::endl;
  //std::cout<<" Total Atoms = "<<u_total_atoms<<std::endl;
  //std::cout<<" Total Bonds = "<<nb<<std::endl;
  v_connections.clear();
  for(uint u=0; u<u_total_atoms; u++){
	//std::cout<<" search for "<<u<<"]"<<std::endl;
	v_connections.push_back(u);
	for(uint b=0; b<nb; b++){
      if(u == _m[b][0]){
        //std::cout<<" Bonds found = "<<_m[b][1]<<std::endl;
        v_connections.push_back(_m[b][1]);
      }else if(u == _m[b][1]){
        //std::cout<<" Bonds found = "<<_m[b][0]<<std::endl;
        v_connections.push_back(_m[b][0]);
      }
    }
    v_connections.push_back(-1);
  }
  //std::cout<<" v_connections = "<<v_connections<<std::endl;
  //std::cout<<" done "<<std::endl;
}

void CFile::set_input_file(std::string s){
  inputfile=s;
}

void CFile::set_output_file_name(std::string s){
  output_filename=s;
}

void CFile::set_potcar_file(std::string s){
  potcarfile=s;
}

void CFile::set_input_dir(std::string s){
  s_actual_dir=s;
}

void CFile::set_input_type(uint u){
  u_input_file_type=u;
}

void CFile::set_input_units(uint u){
  u_input_file_units=u;
}

void CFile::set_input_format(uint u){
  u_input_format=u;
}

void CFile::set_output_file_type(uint u){
  u_output_file_type=u;
}

void CFile::set_output_file_format(uint u){
  u_output_file_format=u;
}

void CFile::set_export_format(uint u){
  u_export_format=u;
}

void CFile::set_fragment_table(TVector<uint> v,uint u){
  v_fragment_table=v;
  u_total_fragments=u;
}

void CFile::set_bounding_box(bool b){
  is_bounding_box=b;
}

void CFile::set_labels(bool b){
  is_labels = b;
}

void CFile::set_numbers(bool b){
  is_numbers = b;
}

void CFile::set_fragments(bool b){
  is_fragments = b;
}

void CFile::set_modified(bool b){
  is_modified = b;
}

void CFile::set_xyz(TMatrix<real> m){
  m_xyz = m;
}

void CFile::set_xyz_cells(int x, int y, int z, int t){
  x_cells=x;
  y_cells=y;
  z_cells=z;
  total_cells=t;
}

// END
