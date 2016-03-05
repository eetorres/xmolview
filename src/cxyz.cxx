//========================================================================
// FILE: cxyz.cxx
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

#include <config_debug.h>
#include <global.h>
#include <cxyz.h>
#include <ctype.h>
#include <algorithm>

//#define _XYZ_INFO_MESSAGES_
//#define _XYZ_DEBUG_MESSAGES_

// formats
// 1 standar xyz (i.e., avogadro)
// 2 standar including fragments
// 3 standar without number of atoms
// 4 with fragments without number of atoms

CXyz::CXyz(){
  clear();
}

void CXyz::clear(void){
  //__xyz_filename = "xyz";
  __system_title = "System Title";
  __export_format=0;
  __file_format=0;
  //__is_selected_dynamics=false;
  __is_direct=false;
  __is_periodic=false;
  __is_symbol=true;
  __is_charges=false;
  __is_fragments=false;
  __is_sort=true;
  __is_dummy=true;
  __is_labels=false;
  __is_symbol_table=false;
  //__fix_all_atoms = false;
  __lattice_constant=1.0;
  __atomic_species=0;
  __total_atoms=0;
  a1=0;
  //
  v_atomic_composition_table.resize(0);
  v_atomic_symbol_table.resize(0);
  v_atomic_number_table.resize(0);
  v_xyz.resize(3);
  //
  uvwTxyz.resize(3,3);
  unit_uvwTxyz.resize(3,3);
  scl_uvwTxyz.resize(3,3);
}

bool CXyz::get_xyz_format(void){
  int res1=-1, res2=-1, res3=-1, symb=-1, k;
  float a, b, c, f, alat;
  float xx, xy, xz, yx, yy, yz, zx, zy, zz;
  char str[20];
  std::string text_line;
  //char text_line[256];
  std::ifstream xyz;
  __head_lines = 0;
  try{
    xyz.open(__filename.c_str());
    if(!xyz.is_open()){
#ifdef _XYZ_DEBUG_MESSAGES_
      std::cout<<"CXYZ: The XYZ file: "<<__filename<<std::endl;
      std::cout<<"CXYZ: was not found\n"<<std::endl;
#endif
      return false;
    }//else{
#ifdef _XYZ_DEBUG_MESSAGES_
      std::cout<<"CXYZ: Get the format of the XYZ file"<<std::endl;
      std::cout<<"CXYZ: Opening the XYZ file: "<<__filename<<std::endl;
#endif
    clear();
    while(res1==-1 && !xyz.eof()){
      std::getline(xyz,text_line);
#ifdef _XYZ_DEBUG_MESSAGES_
      std::cout<<"CXYZ: text_line: "<<text_line<<std::endl;
#endif
      //xyz.getline(text_line,256);
      // Format of the first line
      res1 = std::sscanf((const char*)text_line.c_str(),"%s %f %f %f %i %*s",str,&a,&b,&c,&k);
      //res1 = std::sscanf((const char*)text_line,"%s %f %f %f %f",str,&a,&b,&c,&f);
      if(res1 < 4) __head_lines++;
    }
#ifdef _XYZ_DEBUG_MESSAGES_
    std::cout<<"CXYZ: res1: "<<res1<<std::endl;
#endif
    if(res1==1){
      res1 = std::sscanf((const char*)text_line.c_str(),"%i %*s",&k);
      __total_atoms=k;
      //__head_lines++;
      if(res1==0)
        return false;
    }else if(res1==2){
      // This code may be moved to "read_file"
      res3 = std::sscanf((const char*)text_line.c_str(),"%*s %f %*s %i",&alat,&k);
      if( res3 == 2){
        __total_atoms=k;
      }else{
        res3 = std::sscanf((const char*)text_line.c_str(),"%*s %f",&alat);
      }
#ifdef _XYZ_DEBUG_MESSAGES_
      std::cout<<"CXYZ: res1: "<<res1<<" res3: "<<res3<<std::endl;
#endif
      std::getline(xyz,text_line);
      __is_periodic = true;
      res3 = std::sscanf((const char*)text_line.c_str(),"%f %f %f %f %f %f %f %f %f",&xx,&xy,&xz,&yx,&yy,&yz,&zx,&zy,&zz);
      scl_uvwTxyz[0][0] = xx;
      scl_uvwTxyz[0][1] = xy;
      scl_uvwTxyz[0][2] = xz;
      scl_uvwTxyz[1][0] = yx;
      scl_uvwTxyz[1][1] = yy;
      scl_uvwTxyz[1][2] = yz;
      scl_uvwTxyz[2][0] = zx;
      scl_uvwTxyz[2][1] = zy;
      scl_uvwTxyz[2][2] = zz;
      uvwTxyz=alat*scl_uvwTxyz;
      for(uint i=0; i<3; i++){
        v_xyz[i]=alat*scl_uvwTxyz[i].magnitude();
        unit_uvwTxyz[i]=scl_uvwTxyz[i]/scl_uvwTxyz[i].magnitude();
      }
#ifdef _XYZ_DEBUG_MESSAGES_
    std::cout<<"CXYZ: res1: "<<res1<<" res3: "<<res3<<std::endl;
#endif
#ifdef _XYZ_DEBUG_MESSAGES_
          std::cout<<"scl_uvwTxyz="<<scl_uvwTxyz;
          std::cout<<"unit_uvwTxyz="<<unit_uvwTxyz;
          std::cout<<"uvwTxyz="<<uvwTxyz;
          std::cout<<"v_xyz="<<v_xyz;
#endif
      __head_lines++;
#ifdef _XYZ_DEBUG_MESSAGES_
      std::cout<<"CXYZ: alat = "<<alat<<std::endl;
#endif
    }
#ifdef _XYZ_DEBUG_MESSAGES_
    std::cout<<"CXYZ: res1: "<<res1<<std::endl;
    std::cout<<"CXYZ: search for the initial atom parameters"<<std::endl;
#endif
    //std::getline(xyz,text_line);
    // Search for the intial coordinates
    while(res2 < 4 && !xyz.eof()){
#ifdef _XYZ_DEBUG_MESSAGES_
      std::cout<<"CXYZ: start"<<std::endl;
#endif
      text_line.clear();
      std::getline(xyz,text_line,'\n');
#ifdef _XYZ_DEBUG_MESSAGES_
      std::cout<<"CXYZ: text_line: "<<text_line<<std::endl;
#endif
      //res1 = std::sscanf((const char*)text_line.c_str(),"%s %f %f %f %i",str,&a,&b,&c,&k);
      res2 = std::sscanf((const char*)text_line.c_str(),"%s %f %f %f %s %*s",str,&a,&b,&c,str);
      //res2 = std::sscanf((const char*)text_line,"%s %f %f %f %f",str,&a,&b,&c,&f);
#ifdef _XYZ_DEBUG_MESSAGES_
      std::cout<<"CXYZ: res1: "<<res1<<std::endl;
      std::cout<<"CXYZ: res2: "<<res2<<std::endl;
      std::cout<<"CXYZ: str: "<<str<<std::endl;
#endif
      symb = std::sscanf((const char*)text_line.c_str(),"%i",&k);
      //symb = std::sscanf((const char*)text_line,"%i",&k);
#ifdef _XYZ_DEBUG_MESSAGES_
      std::cout<<"CXYZ: symb: "<<symb<<std::endl;
#endif
      if(res2 < 4) __head_lines++;
#ifdef _XYZ_DEBUG_MESSAGES_
      std::cout<<"CXYZ: res2: "<<res2<<std::endl;
#endif
    }
    if(res2==5){
      std::string str_m=str;
      if(str_m.find('.')!=std::string::npos){
        __is_charges=true;
#ifdef _XYZ_DEBUG_MESSAGES_
        std::cout<<"CXYZ: read charges"<<std::endl;
#endif
      }else{
        __is_fragments=true;
#ifdef _XYZ_DEBUG_MESSAGES_
        std::cout<<"CXYZ: read fragments"<<std::endl;
#endif
      }
    }
    if(symb==0) __is_symbol=true;
    else __is_symbol=false;
#ifdef _XYZ_DEBUG_MESSAGES_
    std::cout<<"CXYZ: res2: "<<res2<<std::endl;
#endif
    if(__total_atoms > 0){           // standard format
      __file_format = INPUT_FORMAT_ATM_NFR;
      if(res2==5){         // using fragments
        if(__is_fragments)
          __file_format = INPUT_FORMAT_ATM_FRG;
        else
          __file_format = INPUT_FORMAT_ATM_CHG;
      }
    }else if(res1==2)      // standard format but not number of atoms
      __file_format = INPUT_FORMAT_NAT_NFR;
    else if(res1==4)      // standard format but not number of atoms
      __file_format = INPUT_FORMAT_NAT_NFR;
    else if(res1==5){       // fragments but not number of atoms
      if(__is_fragments)
          __file_format = INPUT_FORMAT_NAT_FRG;
        else
          __file_format = INPUT_FORMAT_NAT_CHG;
    }
    if(xyz.is_open()){
      xyz.close();
    }
#ifdef _XYZ_DEBUG_MESSAGES_
    std::cout<<"CXYZ: file format: "<<__file_format<<std::endl;
#endif
    // if the atoms need to be counted
    if(__file_format==INPUT_FORMAT_NAT_NFR || __file_format==INPUT_FORMAT_NAT_FRG || __file_format==INPUT_FORMAT_NAT_CHG ){
#ifdef _XYZ_DEBUG_MESSAGES_
      std::cout<<"CXYZ: atoms will be counted"<<std::endl;
      std::cout<<"CXYZ: skiping head lines: "<<__head_lines<<std::endl;
#endif
      xyz.open(__filename.c_str());
      for(uint i=0; i<__head_lines; i++){
        //xyz.ignore(0,'\n');
        std::getline(xyz,text_line);
      }
      __total_atoms=0;
      std::getline(xyz,text_line);
#ifdef _XYZ_DEBUG_MESSAGES_
      std::cout<<"CXYZ: text_line: "<<text_line<<std::endl;
#endif
      res1 = std::sscanf((const char*)text_line.c_str(),"%s %f %f %f %f %*s",str,&a,&b,&c,&f);
      //res1 = std::sscanf((const char*)text_line,"%s %f %f %f %f",str,&a,&b,&c,&f);
      while((res1 == 4 || res1 == 5) && !xyz.eof()){
        __total_atoms++;
        std::getline(xyz,text_line);
        //xyz.getline(text_line,256);
        res1 = std::sscanf((const char*)text_line.c_str(),"%s %f %f %f %f %*s",str,&a,&b,&c,&f);
        //res1 = std::sscanf((const char*)text_line,"%s %f %f %f %f",str,&a,&b,&c,&f);
      }
      xyz.close();
#ifdef _XYZ_DEBUG_MESSAGES_
      std::cout<<"CXYZ: number of atoms found "<<__total_atoms<<std::endl;
#endif
    }
#ifdef _XYZ_DEBUG_MESSAGES_
    else{
      std::cout<<"CXYZ: atoms will !!!NOT!!! be counted"<<std::endl;
    }
#endif
#ifdef _XYZ_DEBUG_MESSAGES_
    std::cout<<"CXYZ: total number of atoms: "<<__total_atoms<<std::endl;
    std::cout<<"CXYZ: <--!!!The XYZ file was successful closed!!!-->"<<std::endl;
#endif
    if(__file_format==INPUT_FORMAT_NAT_NFR || __file_format==INPUT_FORMAT_NAT_FRG){
      if(__total_atoms>0)
        return true;
      else
        return false;
    }else{
      return true;
    }
  }catch(...){
    return false;
  }
}

// load the XYZ file
bool CXyz::read_file(std::string d, std::string f){
  std::string atomic_symbol;
  real data, frag;
  uint unum;
  TVector<real> v1, v2, v3;
  __filename = d+"/"+f;
  std::ifstream xyz;
  try{
    if(!get_xyz_format()){
      return false;
    }
    //
    xyz.open(__filename.c_str());
    //if(__file_format == INPUT_FORMAT_ATM_NFR || __file_format == INPUT_FORMAT_ATM_FRG || __file_format == INPUT_FORMAT_ATM_CHG){
      //xyz>>data;
      //__total_atoms=data;
    //}
#ifdef _XYZ_DEBUG_MESSAGES_
    std::cout<<"CXYZ: Total atoms: "<<__total_atoms<<std::endl;
#endif
    v_atomic_symbols.resize(__total_atoms);
    v_atom_table.resize(__total_atoms);
    m_xyz.resize(__total_atoms,3);
    m_xyz_input.resize(__total_atoms,3);
    m_uvw.resize(__total_atoms,3);
    //if(__file_format == 4 || __file_format == 2){
    v_fragment_table.resize(__total_atoms);
    //v_atomic_charges.resize(__total_atoms);
    //}else if(__file_format == 5 || __file_format == 6){
    if(__is_charges)
      v_atomic_charges.resize(__total_atoms);
    //}
    //while(xyz.peek()=='\n')
      //xyz.ignore(0,'\n');
    for(uint i=0; i<__head_lines; i++){
      xyz.ignore(1024,'\n');
    }
    //xyz.ignore(0,'\n');
    //xyz.ignore(0,'\n');
#ifdef _XYZ_INFO_MESSAGES_
    std::cout<<"CXYZ: loading atomic coordinates"<<std::endl;
#endif
    for(uint f=0;f<__total_atoms;f++){
      if(__is_symbol){
        xyz>>atomic_symbol;
      }else{
        xyz>>data;
        atomic_symbol=symbol[(uint)(data-1)];
      }
//#ifdef _XYZ_INFO_MESSAGES_
      //std::cout<<"make sure it uses the standard atomic symbols"<<std::endl;
//#endif
      //transform(atomic_symbol.begin(), atomic_symbol.begin()+1, atomic_symbol.begin(), toupper);
      //transform(atomic_symbol.begin()+1, atomic_symbol.end(), atomic_symbol.begin()+1, tolower);
      v_atomic_symbols[f]=atomic_symbol;
#ifdef _XYZ_INFO_MESSAGES_
      std::cout<<"atomic_symbol "<<atomic_symbol<<std::endl;
#endif
      for(uint c=0;c<3;c++){  // read coordinates
        xyz>>data;
        m_xyz_input[f][c]=data;
      }
      if(__is_fragments){   // read fragments
        xyz>>unum;
        v_fragment_table[f]=unum;
      }
      if(__is_charges){
//#ifdef _XYZ_INFO_MESSAGES_
      //std::cout<<"read charges"<<std::endl;
//#endif
        xyz>>frag;
        v_atomic_charges[f]=frag;
      }
      if(xyz.peek()!='\n')
        xyz.ignore(1024,'\n');
    }
    if(u_input_file_units==1){
      m_xyz_input=(0.52917720859*m_xyz_input);
    }
    v_atomic_labels=v_atomic_symbols;
    eval_atomic_composition();
    if(is_periodic()){
      m_xyz=m_xyz_input;
    }else{
      center_coordinates();
    }
    //eval_atomic_number_table();

    //m_uvw=m_uvw-0.5;
    //////compute_general();
#ifdef _XYZ_DEBUG_MESSAGES_
    std::cout<<"--- Summary of the XYZ file ---"<<std::endl;
    std::cout<<"atom table ="<<v_atom_table;
    std::cout<<"atomic symbol table ="<<v_atomic_symbol_table;
    std::cout<<"atomic number table ="<<v_atomic_number_table;
    std::cout<<"atomic composition table ="<<v_atomic_composition_table;
    std::cout<<"atomic symbols ="<<v_atomic_symbols;
    std::cout<<"atomic numbers ="<<v_atomic_numbers;
    std::cout<<"uvwTxyz'="<<uvwTxyz;
    std::cout<<"--- Summary of the XYZ file ---"<<std::endl;
#endif

#ifdef _display_data_
    std::cout<<"Cartesian'="<<m_xyz;
#endif
    if(xyz.is_open()){
      xyz.close();
    }
#ifdef _XYZ_DEBUG_MESSAGES_
    std::cout<<"CXYZ: The XYZ was successful closed!"<<std::endl;
#endif
    return true;
  }catch(...){
    return false;
  }
}

void CXyz::write_file(void){
  write_file(__filename);
}

void CXyz::write_file(std::string d, std::string f){
  std::string _f=d+'/'+f;
  write_file(_f);
}

//  save the XYZ file
void CXyz::write_file(std::string _fn){
  std::string symbol;
  real data;
  TVector<real> v_scale(3);
  std::ofstream oxyz;
  oxyz.open(_fn.c_str());
  //oxyz.width(17);
  if(oxyz.bad())
    std::cout<<" CXYZ: The file was not created"<<std::endl;
#ifdef _XYZ_INFO_MESSAGES_
  std::cout<<" CXYZ: file format: "<<__export_format<<std::endl;
#endif
  if(__export_format==OUTPUT_FORMAT_ATM_NFR || __export_format==OUTPUT_FORMAT_ATM_FRG){
    oxyz<<__total_cells*__total_atoms;
    oxyz<<std::endl;
    oxyz<<std::endl;
  }
#ifdef _XYZ_INFO_MESSAGES_
  std::cout<<"CXYZ: total atoms: "<<__total_cells*__total_atoms<<std::endl;
#endif
  oxyz.precision(16);
#ifdef _XYZ_INFO_MESSAGES_
  std::cout<<"CXYZ: v_fragment_table "<<v_fragment_table;
  std::cout<<"CXYZ: m_xyz "<<m_xyz;
#endif
  for(uint cell=0; cell<__total_cells; cell++){ // repetition cells
    for(uint f=0;f<__total_atoms;f++){
      symbol=v_atomic_symbols[f];
      if(strcmp(symbol.c_str(),"X") || __is_dummy){
        oxyz.width(3);
        oxyz<<std::fixed<<std::left<<symbol;
        for(uint c=0;c<3;c++){
          data=m_xyz[f+cell*__total_atoms][c];
          oxyz.width(22);
          oxyz<<std::fixed<<std::right<<data;
        }
        // include fragments
        if(__export_format==OUTPUT_FORMAT_ATM_FRG || __export_format==OUTPUT_FORMAT_NAT_FRG){
#ifdef _XYZ_INFO_MESSAGES_
        std::cout<<"CXYZ: include fragments"<<std::endl;
#endif
          oxyz.width(5);
          oxyz<<std::fixed<<std::right<<v_fragment_table[f];
        }
        oxyz<<std::endl;
      }
    }
  }
  //if(__export_format==0)
  //write_copyright(oxyz);
  oxyz.close();
}

void CXyz::eval_atomic_composition(void){
  std::string _tstr, tmp_symbol;
  size_t _sw;
  uint _an;
  v_atomic_symbol_table.resize(0);
  v_atomic_number_table.resize(0);
  TVector<real> _tv;
  __vacuum_space=0;
  for(uint f=0;f<__total_atoms;f++){
    tmp_symbol=v_atomic_symbols[f];
    transform(tmp_symbol.begin(), tmp_symbol.begin()+1, tmp_symbol.begin(), toupper);
    transform(tmp_symbol.begin()+1, tmp_symbol.end(), tmp_symbol.begin()+1, tolower);
    for(uint k=0; k<periodic_table_atoms; k++){
      if(k < 95){ _sw = 2;
      }else{ _sw = 1;}
      if(strncmp(tmp_symbol.c_str(),symbol_sorted[k].c_str(),_sw)==0){
#ifdef _XYZ_DEBUG_COMPOSITION_
        std::cout<<" XYZ: k,sw: "<<k<<","<<_sw<<std::endl;
        std::cout<<" XYZ: tmp_symbol: "<<tmp_symbol<<std::endl;
        std::cout<<" XYZ: symbol_sorted[k]: "<<symbol_sorted[k]<<std::endl;
#endif
        v_atomic_symbols[f]=symbol_sorted[k];
        break;
      }
    }
  }
  if(!__is_symbol_table){
    for(uint j=0; j<periodic_table_atoms; j++){
      //if(j < 95) _sw = 2;
      //else _sw = 1;
      for(uint f=0;f<__total_atoms;f++){
        if(!strcmp(v_atomic_symbols[f].c_str(),symbol[j].c_str())){
        //if(strncmp(v_atomic_symbols[f].c_str(),symbol_sorted[j].c_str(),_sw)==0){
          v_atomic_symbol_table.push_back(symbol[j]);
          _an = get_atomic_number(symbol[j]);
          v_atomic_number_table.push_back(_an);
          __vacuum_space=maxi(__vacuum_space,atom_rrgb[_an][0]);
#ifdef _XYZ_INFO_MESSAGES_
          std::cout<<" CXYZ: vacuum space: "<<__vacuum_space<<std::endl;
#endif
          break;
        }
      }
    }
    __atomic_species=v_atomic_symbol_table.size();
    v_atomic_composition_table.resize(__atomic_species);
  }
  // this function check if the atomic species are sorted
  // by atomic symbol, if not then it does it.
  uint _f=0, _c;
  v_atomic_numbers.resize(__total_atoms);
  for(uint i=0; i<__atomic_species;i++){
    _c=0;
    for(uint j=_f; j<__total_atoms;j++){
      v_atomic_numbers[j]=get_atomic_number(v_atomic_symbols[j]);
      //v_atomic_numbers[j]=get_atomic_number(v_atomic_symbols[j],v_atomic_symbol_size[i]);
      if(v_atomic_symbol_table[i]==v_atomic_symbols[j]){
        v_atom_table[j]=i;
        _c++;
      }

    }
    v_atomic_composition_table[i]=_c;
  }
}

// SET FUNCTIONS
void CXyz::set_file_format(const uint u){
  __file_format=u;
}

void CXyz::set_export_format(const uint u){
  __export_format=u;
}

void CXyz::set_total_atoms(const uint u){
  __total_atoms=u;
}

void CXyz::set_total_cells(const uint x, const uint y, const uint z, const uint t){
  __x_cells=x;
  __y_cells=y;
  __z_cells=z;
  __total_cells=t;
}

void CXyz::set_fragments(const TVector<uint>& _v){
  v_fragment_table=_v;
}

void CXyz::set_atomic_symbol_table(const TVector<std::string>& _v){
  v_atomic_symbol_table=_v;
  __atomic_species=v_atomic_symbol_table.size();
  v_atomic_composition_table.resize(__atomic_species);
  __is_symbol_table=true;
}

/*
void CXyz::eval_atomic_number_table(void){
  uint z;
  v_atomic_number_table.resize(0);
  for( uint i=0; i<v_atomic_symbol_table.size(); i++){
    z=get_atomic_number(v_atomic_symbol_table[i]);
    v_atomic_number_table.push_back(z);
  }
}*/

void CXyz::set_cartesian(const TMatrix<real>& _m){
   m_xyz=_m;
}

void CXyz::set_direct(const TMatrix<real>& _m){
  m_uvw=_m;
}

//////////////////////////////
//   Coordinate functions  //
//////////////////////////////

uint CXyz::get_atomic_number(std::string s){
  //transform(symbol.begin(), symbol.end(), symbol.begin(), toupper);
  //transform(s.begin(), s.begin()+1, s.begin(), toupper);
  for(uint j=0; j<periodic_table_atoms; j++){
    //if(!strcmp(s.c_str(),symbol[j].c_str())){
    if(strcmp(s.c_str(),symbol[j].c_str())==0){
//#ifdef _XYZ_DEBUG_MESSAGES_
//      std::cout<<"CXYZ: atomic number "<<j+1<<std::endl;
//#endif
      return j;
    }
  }
  return 118;
}

uint CXyz::get_atomic_composition(uint i){
  return  v_atomic_composition_table[i];
}

bool CXyz::is_direct(){
  return __is_direct;
}

bool CXyz::is_periodic(){
  return __is_periodic;
}

bool CXyz::is_charges(void){
  return __is_charges;
}

bool CXyz::is_fragments(void){
  return __is_fragments;
}

uint CXyz::get_total_species(void){
  return __atomic_species;
}

uint CXyz::get_total_atoms(void){
  return m_xyz.rows();
}

uint CXyz::get_format(void){
  return __file_format;
}

real CXyz::get_cartesian(uint x, uint y){
  return  m_xyz[x][y];
}

real CXyz::get_direct(uint x, uint y){
  return  m_uvw[x][y];
}

TVector<std::string> CXyz::get_atomic_symbols(void){
  return  v_atomic_symbols;
}

TVector<std::string> CXyz::get_atomic_labels(void){
  return  v_atomic_labels;
}

TVector<std::string> CXyz::get_atomic_symbol_table(void){
  return  v_atomic_symbol_table;
}

TVector<uint> CXyz::get_atom_table(void){
  return  v_atom_table;
}

TVector<uint> CXyz::get_atomic_composition_table(void){
  return  v_atomic_composition_table;
}

TVector<uint> CXyz::get_atomic_number_table(void){
  return  v_atomic_number_table;
}

TVector<uint> CXyz::get_atomic_numbers(void){
  return  v_atomic_numbers;
}

TVector<uint> CXyz::get_fragment_table(void){
  return  v_fragment_table;
}

TVector<real> CXyz::get_direct(uint x){
  return  m_uvw[x];
}

TVector<real> CXyz::get_cartesian(uint x){
  return  m_xyz[x];
}

TVector<real> CXyz::get_direct_basis(void){
  return v_xyz;
}

TVector<real> CXyz::get_charges(void){
  return v_atomic_charges;
}

TMatrix<real> CXyz::get_xyz_input(void){
  return  m_xyz_input;
}

TMatrix<real> CXyz::get_cartesian(void){
  return  m_xyz;
}

TMatrix<real> CXyz::get_direct(void){
  return  m_uvw;
}

// Extra functions
void CXyz::write_copyright(std::ofstream& _o){
  time_t rawtime;
  char * _st, tmp_buffer[256];
  std::string stime;
  time (&rawtime);
  sprintf(tmp_buffer,"%s",ctime(&rawtime));
  _st = strtok(tmp_buffer,"\n");
  stime = _st;
  _o<<std::endl<<std::endl<<std::endl;
  _o<<"# ------------------------------------------------------------#"<<std::endl;
  _o<<"# XYZ structure file generated by XMolView (www.xmol.org)     #"<<std::endl;
  _o<<"# ------------------------------------------------------------#"<<std::endl;
  _o<<"# Notes:                                                      #"<<std::endl;
  _o<<"#                                                             #"<<std::endl;
  _o<<"#                                                             #"<<std::endl;
  _o<<"#                                                             #"<<std::endl;
  _o<<"#                                                             #"<<std::endl;
  _o<<"# ------------------------------------------------------------#"<<std::endl;
  _o<<"# Date of generation: "<<stime<<std::endl;
  //_o<<"# ------------------------------------------------------------#"<<endl;
}

void CXyz::center_coordinates(void){
  if(__is_periodic){
    //m_xyz = m_xyz_input;
#ifdef _XYZ_DEBUG_DATA_
    std::cout<<" CXYZ: is_periodic"<<std::endl;
#endif
    for(uint i=0;i<__total_atoms;i++){
      for(uint j=0; j<3; j++){
        m_xyz[i]+=(0.5*uvwTxyz[j]);
      }
    }
  }else{
    TVector<real> min_xyz(3);
    TVector<real> max_xyz(3);
    real _vacuum = +2.2*__vacuum_space;
    min_xyz[0]=m_xyz_input.col_min(0);
    min_xyz[1]=m_xyz_input.col_min(1);
    min_xyz[2]=m_xyz_input.col_min(2);
    #ifdef _XYZ_DEBUG_DATA_
    std::cout<<" CXYZ: initial m_xyz = "<<m_xyz;
#endif
    // set the coordinate reference to (0,0,0)
    for(uint i=0;i<__total_atoms;i++){
      m_xyz[i]= m_xyz_input[i] - min_xyz;
    }
    max_xyz[0]=m_xyz.col_max(0)+_vacuum;
    max_xyz[1]=m_xyz.col_max(1)+_vacuum;
    max_xyz[2]=m_xyz.col_max(2)+_vacuum;
#ifdef _XYZ_DEBUG_DATA_
    std::cout<<" CXYZ: positive m_xyz = "<<m_xyz;
#endif
    if(max_xyz[0]>0.01)
      uvwTxyz[0][0]=max_xyz[0];
    else{
      uvwTxyz[0][0]=_vacuum;
      max_xyz[0]=_vacuum;
    }
    if(max_xyz[1]>0.01)
      uvwTxyz[1][1]=max_xyz[1];
    else{
      uvwTxyz[1][1]=_vacuum;
      max_xyz[1]=_vacuum;
    }
    if(max_xyz[2]>0.01)
      uvwTxyz[2][2]=max_xyz[2];
    else{
      uvwTxyz[2][2]=_vacuum;
      max_xyz[2]=_vacuum;
    }
    for(uint i=0;i<__total_atoms;i++){
      m_xyz[i]=(m_xyz[i]+1.1*__vacuum_space);
    }
#ifdef _XYZ_DEBUG_DATA_
    std::cout<<" CXYZ: boxed m_xyz = "<<m_xyz;
#endif
    unit_uvwTxyz[0]=uvwTxyz[0]/uvwTxyz[0].magnitude();
    unit_uvwTxyz[1]=uvwTxyz[1]/uvwTxyz[1].magnitude();
    unit_uvwTxyz[2]=uvwTxyz[2]/uvwTxyz[2].magnitude();
  }
  //
#ifdef _XYZ_DEBUG_DATA_
  std::cout<<" CXYZ: min_xyz = "<<min_xyz;
  std::cout<<" CXYZ: max_xyz = "<<max_xyz;
#endif
}

// END

