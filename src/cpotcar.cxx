//========================================================================
// FILE: cpotcar.cxx
//
// This an utility program to manipulate and generate structure files
//
// Copyright 2006-2015 by Edmanuel Torres
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

//#include<const.h>
#include<cpotcar.h>

//#define _information_messages_   1
//#define _POTCAR_DEBUG_MESSAGES_     1
//#define _display_data_           1

CPotcar::CPotcar(){
  clear();
}

void CPotcar::clear(void){
  __potcar_species=0;
  v_atomic_symbols.resize(0);
  v_atomic_numbers.resize(0);
}

// load the POTCAR file
bool CPotcar::read_potcar(std::string d, std::string f){
  char tmp_buffer[256], *pch;
  std::string atom_symbol;
  std::ifstream potcar;
  // Reading atomic species
  std::string potcar_file = d+"/"+f;
  potcar.open(potcar_file.c_str());
  if(!potcar.is_open()){
#ifdef _POTCAR_DEBUG_MESSAGES_
    std::cout<<"The POTCAR file was not found"<<std::endl;
#endif
    return false;
  }else{
#ifdef _POTCAR_DEBUG_MESSAGES_
    std::cout<<"The POTCAR was successful opened"<<std::endl;
#endif
  }
  while(!potcar.eof()){
    potcar.getline(tmp_buffer,256);
    //pruintf("%s\n",tmp_buffer);
    pch=strstr(tmp_buffer,"VRHFIN");
    if(pch!=NULL)
    {
      pch=strrchr(tmp_buffer,'=');
      pch = strtok (pch," =:");
#ifdef _POTCAR_DEBUG_MESSAGES_
      std::cout<<"found "<<pch<<std::endl;
#endif
      atom_symbol = pch;
      __potcar_species++;
      v_atomic_symbols.push_back(atom_symbol);
    }
  }

  if(!potcar.is_open()){
   return false;
  }
  potcar.close();
#ifdef _POTCAR_DEBUG_MESSAGES_
  std::cout<<"The POTCAR was successful closed"<<std::endl;
#endif
/*
    if(__atomic_species!=__potcar_species){
     cout<<"WARNING:"<<endl;
     cout<<" > The number of species in the POSCAR and POTCAR mismatch"<<endl;
     cout<<"Number of species POSCAR: "<<__atomic_species<<std::endl;
     cout<<"Number of species POTCAR: "<<__potcar_species<<std::endl;
     cout<<"Please check..."<<std::endl;
    }*/
  eval_atomic_numbers();
#ifdef _POTCAR_DEBUG_MESSAGES_
    std::cout<<"POTCAR: atomic numbers: "<<v_atomic_numbers;
#endif
  return false;
}

void CPotcar::set_symbols(TVector<std::string> _v){
  v_atomic_symbols=_v;
  __potcar_species=_v.size();
#ifdef _POTCAR_DEBUG_MESSAGES_
  std::cout<<"POTCAR: potcar species: "<<__potcar_species<<std::endl;
#endif
}

void CPotcar::eval_atomic_numbers(void){
  v_atomic_numbers.resize(0);
  for(uint i=0; i<__potcar_species; i++){
    // including dummy atoms X
#ifdef _POTCAR_DEBUG_MESSAGES_
    std::cout<<"atomic symbol "<<v_atomic_symbols[i]<<std::endl;
#endif
    for(uint j=0; j<periodic_table_atoms; j++){
#ifdef _POTCAR_DEBUG_MESSAGES_
      std::cout<<" POTCAR: symbol "<<symbol[j]<<std::endl;
#endif
      if(!strcmp(v_atomic_symbols[i].c_str(),symbol[j].c_str())){
#ifdef _POTCAR_DEBUG_MESSAGES_
        std::cout<<" POTCAR: Found atomic number ("<<i<<"): "<<j+1<<std::endl;
#endif
        v_atomic_numbers.push_back(j);
        // here a break may be implemented
      }
    }
  }
}

uint CPotcar::get_atomic_number(uint i){
  return  v_atomic_numbers[i];
}

std::string CPotcar::get_symbol(uint i){
  return  v_atomic_symbols[i];
}

TVector<uint> CPotcar::get_atomic_numbers(void){
  return  v_atomic_numbers;
}

TVector<std::string> CPotcar::get_symbols(void){
  return  v_atomic_symbols;
}

uint CPotcar::get_total_species(void){
  return  __potcar_species;
}


