//========================================================================
// FILE: cgau.cxx
//
// This an utility program to manipulate and generate structure files
//
// Copyright 2012-2015 by Edmanuel Torres
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
#include<cgau.h>
#include<ctype.h>
#include<algorithm>

// formats
// 1 gaussian no fragments
// 2 gaussian including fragments

CGau::CGau(){
  clear();
}

void CGau::clear(void){
  //__xyz_filename = "xyz";
  __system_title = "System Title";
  __export_format=0;
  __file_format=0;
  __header_lines=0;
  //__is_selected_dynamics=false;
  __is_direct=false;
  __is_periodic=false;
  __is_symbol=true;
  __is_sort=true;
  __is_dummy=true;
  __is_labels=false;
  __is_symbol_table=false;
  __is_read = false;
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
  //scl_uvwTxyz.resize(3,3);
}

bool CGau::get_gau_format(void){
  int res1, symb;
  int k, l;
  float a, b, c, f;
  char str[20];
  std::string text_line;
  std::ifstream gau;
  try{
    gau.open(__filename.c_str());
    if(!gau.is_open()){
#ifdef _GAU_DEBUG_MESSAGES_
      std::cout<<"CGAU: The GAU file: "<<__filename<<std::endl;
      std::cout<<"CGAU: was not found\n"<<std::endl;
#endif
      return false;
    }//else{
    res1=0;
    __header_lines=0;
    while(res1!=2 && !gau.eof()){
      //gau.ignore(256,'\n');
      std::getline(gau,text_line);
      std::replace(text_line.begin(),text_line.end(), ',', ' ');
      res1 = std::sscanf((const char*)text_line.c_str(),"%i %i %*s",&k,&l);
      //res2 = std::sscanf((const char*)text_line.c_str(),"%i, %i %*s",&k,&l);
      __header_lines++;
#ifdef _GAU_DEBUG_MESSAGES_
      std::cout<<"[TEXT] "<<text_line<<std::endl;
      std::cout<<"CGAU: res1 = "<<res1<<std::endl;
      //std::cout<<"CGAU: res2 = "<<res2<<std::endl;
#endif
    }
    // count the line 6
    std::getline(gau,text_line);
    std::replace(text_line.begin(),text_line.end(), ',', ' ');
    res1 = std::sscanf((const char*)text_line.c_str(),"%s %f %f %f %f",str,&a,&b,&c,&f);
    symb = std::sscanf((const char*)text_line.c_str(),"%i",&k);
#ifdef _GAU_DEBUG_MESSAGES_
    std::cout<<"CGAU: res1 = "<<res1<<std::endl;
    std::cout<<"CGAU: symb = "<<symb<<std::endl;
#endif
    if(symb==0) __is_symbol=true;
    else __is_symbol=false;
#ifdef _GAU_DEBUG_MESSAGES_
    std::cout<<" symbols: "<<__is_symbol<<std::endl;
    std::cout<<" res1: "<<res1<<std::endl;
#endif
    if(res1==4){               // standar format
      __file_format = 1;
    }else if(res1==5){         // using fragments
      __file_format = 2;
    }
#ifdef _GAU_DEBUG_MESSAGES_
    std::cout<<"CGAU: file format: "<<__file_format<<std::endl;
#endif
    // atoms need to be counted
    //gau.open(__filename.c_str());
    __total_atoms=0;
    while((res1 == 4 || res1 == 5) && !gau.eof()){
      __total_atoms++;
      std::getline(gau,text_line);
      std::replace(text_line.begin(),text_line.end(), ',', ' ');
      res1 = std::sscanf((const char*)text_line.c_str(),"%s %f %f %f %f",str,&a,&b,&c,&f);
#ifdef _GAU_DEBUG_MESSAGES_
      std::cout<<"[TEXT] "<<text_line<<std::endl;
      std::cout<<"CGAU: res1 = "<<res1<<std::endl;
#endif
    }
    //
    if(gau.is_open()){
      gau.close();
    }
#ifdef _GAU_DEBUG_MESSAGES_
    std::cout<<"CGAU: number of atoms found "<<__total_atoms<<std::endl;
    std::cout<<"CGAU: The GAU file was successful opened"<<std::endl;
#endif
    if(__total_atoms>0)
      return true;
    else
      return false;
  }catch(...){
    return false;
  }
}

// load the GAU file
bool CGau::read_file(std::string d, std::string f){
  std::string atomic_symbol;
  real data;
  //uint unum;
  //int res1
  int g;
  float a, b, c;
  char str[20];
  TVector<real> v1, v2, v3;
  __filename = d+"/"+f;
  std::string text_line;
  std::ifstream gau;
  //
  if(!get_gau_format()){
    return false;
  }
  //
  try{
    gau.open(__filename.c_str());
#ifdef _GAU_DEBUG_MESSAGES_
    std::cout<<"CGAU: Total atoms: "<<__total_atoms<<std::endl;
#endif
    v_atomic_symbols.resize(__total_atoms);
    v_fragment_table.resize(__total_atoms);
    v_atom_table.resize(__total_atoms);
    m_xyz_input.resize(__total_atoms,3);
    m_xyz.resize(__total_atoms,3);
    m_uvw.resize(__total_atoms,3);
    v_file_header.resize(__header_lines);
    for(uint i=0; i<__header_lines; i++){
	std::getline(gau,text_line);
	v_file_header[i] =  text_line;
      //gau.ignore(256,'\n');
    }
#ifdef _GAU_INFO_MESSAGES_
    std::cout<<"CGAU: loading atomic coordinates"<<std::endl;
#endif
#ifdef _GAU_DEBUG_MESSAGES_
    int res1
#endif
    for(uint f=0;f<__total_atoms;f++){
      std::getline(gau,text_line);
      std::replace(text_line.begin(),text_line.end(), ',', ' ');
      if(__file_format==2){   // read fragments
#ifdef _GAU_DEBUG_MESSAGES_
        res1 = std::sscanf((const char*)text_line.c_str(),"%s %f %f %f %i",str,&a,&b,&c,&g);
#else
        std::sscanf((const char*)text_line.c_str(),"%s %f %f %f %i",str,&a,&b,&c,&g);
#endif
        v_fragment_table[f]=g;
      }else{
#ifdef _GAU_DEBUG_MESSAGES_
        res1 = std::sscanf((const char*)text_line.c_str(),"%s %f %f %f %*s",str,&a,&b,&c);
#else
        std::sscanf((const char*)text_line.c_str(),"%s %f %f %f %*s",str,&a,&b,&c);
#endif
      }
      //gau>>symbol;
      if(__is_symbol){
        atomic_symbol=str;
      }else{
        data=atoi(str);
        atomic_symbol=symbol[(uint)(data-1)];
      }
      // make sure it uses the standar atom symbols
      transform(atomic_symbol.begin(), atomic_symbol.begin()+1, atomic_symbol.begin(), toupper);
      transform(atomic_symbol.begin()+1, atomic_symbol.end(), atomic_symbol.begin()+1, tolower);
      v_atomic_symbols[f]=atomic_symbol;
      m_xyz_input[f][0]=a;
      m_xyz_input[f][1]=b;
      m_xyz_input[f][2]=c;
      //for(uint c=0;c<3;c++){  // read coordinates
        //gau>>data;
        //m_xyz_input[f][c]=data;
      //}
      //if(__file_format==2){   // read fragments
        //gau>>unum;
        //v_fragment_table[f]=unum;
      //}
      if(gau.peek()!='\n')
        gau.ignore(0,'\n');
    }
    v_atomic_labels=v_atomic_symbols;
    eval_atomic_composition();
    center_coordinates();
    eval_atomic_number_table();
    //m_uvw=m_uvw-0.5;
    //////compute_general();
#ifdef _GAU_DEBUG_MESSAGES_
    std::cout<<"--- Summary of the GAU file ---"<<std::endl;
    std::cout<<"atom table ="<<v_atom_table;
    std::cout<<"atomic symbol table ="<<v_atomic_symbol_table;
    std::cout<<"atomic number table ="<<v_atomic_number_table;
    std::cout<<"atomic composition table ="<<v_atomic_composition_table;
    std::cout<<"atomic symbols ="<<v_atomic_symbols;
    std::cout<<"atomic numbers ="<<v_atomic_numbers;
    std::cout<<"uvwTxyz'="<<uvwTxyz;
    std::cout<<"--- Summary of the GAU file ---"<<std::endl;
#endif

#ifdef _GAU_DEBUG_DATA_
    std::cout<<"Cartesian'="<<m_xyz;
#endif
    if(gau.is_open()){
      gau.close();
    }
#ifdef _GAU_DEBUG_MESSAGES_
    std::cout<<"CGAU: The GAU was successful closed!"<<std::endl;
#endif
    __is_read = true;
    return true;
  }catch(...){
    return false;
  }
}

void CGau::eval_atomic_composition(void){
  std::string _tstr;
  uint _an;
  size_t _sw;
  std::string tmp_symbol;
  v_atomic_symbol_table.resize(0);
  TVector<real> _tv;
  TVector<bool> v_not_found(__total_atoms);
  for(uint f=0;f<__total_atoms;f++)
    v_not_found[f] = true;
  __vacuum_space=0;
  if(!__is_symbol_table){
    for(uint j=0; j<periodic_table_atoms; j++){
      if(j < 95) _sw = 2;
      else _sw = 1;
      for(uint f=0;f<__total_atoms;f++){
        if(v_not_found[f]){
          tmp_symbol=v_atomic_symbols[f];
          //std::cout<<" DLP: tmp_symbol ["<<_sw<<"] "<<tmp_symbol<<" <=> "<<symbol_sorted[j]<<" - ";
          transform(tmp_symbol.begin(), tmp_symbol.begin()+1, tmp_symbol.begin(), toupper);
          transform(tmp_symbol.begin()+1, tmp_symbol.end(), tmp_symbol.begin()+1, tolower);
          if(strncmp(tmp_symbol.c_str(),symbol_sorted[j].c_str(),_sw)==0){
            v_atomic_symbol_table.push_back(symbol_sorted[j]);
            //if(strcmp(v_atomic_symbols[f].c_str(),symbol_sorted[j].c_str()))
            v_atomic_symbols[f]=symbol_sorted[j];
            //std::cout<<" YES"<<std::endl;
            v_not_found[f] = false;
            //break;
          }
          //std::cout<<" NO"<<std::endl;
          //std::cout<<" DLP: tmp_symbol: "<<tmp_symbol<<std::endl;
        }
      }
      //std::cout<<" END"<<std::endl;
    }
    __atomic_species=v_atomic_symbol_table.size();
    v_atomic_composition_table.resize(__atomic_species);
  }
  //std::cout<<" v_atomic_symbols = "<< v_atomic_symbols;
  // this function check if the atomic species are sorted
  // by atomic symbol, if not then it does it.
  uint _f=0, _c;
  v_atomic_numbers.resize(__total_atoms);
  for(uint i=0; i<__atomic_species;i++){
    _c=0;
    for(uint j=_f; j<__total_atoms;j++){
      _an = get_atomic_number(v_atomic_symbols[j]);
      v_atomic_numbers[j]=_an;
      __vacuum_space=maxi(__vacuum_space,atom_rrgb[_an][0]);
      if(v_atomic_symbol_table[i]==v_atomic_symbols[j]){
        v_atom_table[j]=i;
        _c++;
      }

    }
    v_atomic_composition_table[i]=_c;
  }
}

void CGau::write_file(void){
  write_file(__filename);
}

void CGau::write_file(std::string d, std::string f){
  std::string _f=d+'/'+f;
  write_file(_f);
}

//  save the GAU file
void CGau::write_file(std::string _fn){
  std::string symbol;
  real data;
  TVector<real> v_scale(3);
  std::ofstream ogau;
  ogau.open(_fn.c_str());
  //ogau.width(17);
  if(ogau.bad())
    std::cout<<" CGAU: The file was not created"<<std::endl;
#ifdef _GAU_DEBUG_MESSAGES_
  std::cout<<" CGAU: file format: "<<__file_format<<std::endl;
#endif
  if(__export_format==OUTPUT_FORMAT_ATM_NFR || __export_format==OUTPUT_FORMAT_ATM_FRG){
    ogau<<__total_atoms;
    ogau<<std::endl;
    ogau<<std::endl;
  }
  std::string stime = get_date();
  if( __is_read ){
	for(uint i=0; i<__header_lines; i++){
	  ogau<<v_file_header[i]<<std::endl;
    }
  }else{
	ogau<<"Gaussian file generated by XMolView (www.xmol.org) ["<<stime<<"]"<<std::endl;
	ogau<<"0 1"<<std::endl;
  }
  //ogau<<"0 1"<<std::endl;
  ogau.precision(16);
  for(uint cell=0; cell<__total_cells; cell++){ // repetition cells
    for(uint f=0;f<__total_atoms;f++){
      symbol=v_atomic_symbols[f];
      if(strcmp(symbol.c_str(),"X") || __is_dummy){
        ogau<<symbol;
        for(uint c=0;c<3;c++){
          data=m_xyz[f+cell*__total_atoms][c];
          ogau.width(22);
          ogau<<std::fixed<<std::right<<data;
        }
        // export format
        // 1: Gaussian
        if(__export_format==OUTPUT_FORMAT_ATM_FRG || __export_format==OUTPUT_FORMAT_NAT_FRG){
          ogau.width(5);
          ogau<<std::fixed<<std::right<<v_fragment_table[f];
        }
        ogau<<std::endl;
      }
    }
  }
  if(__export_format==0)
    write_copyright(ogau);
  ogau.close();
}

// SET FUNCTIONS

void CGau::set_file_format(const uint u){
  __file_format=u;
}

void CGau::set_export_format(const uint u){
  __export_format=u;
}

void CGau::set_total_cells(const uint x, const uint y, const uint z, const uint t){
  __x_cells=x;
  __y_cells=y;
  __z_cells=z;
  __total_cells=t;
}

void CGau::set_fragments(const TVector<uint>& _v){
  v_fragment_table=_v;
#ifdef _GAU_DEBUG_MESSAGES_
  std::cout<<" CGAU: fragments: "<<v_fragment_table;
#endif
}

void CGau::set_atomic_symbol_table(const TVector<std::string>& _v){
  v_atomic_symbol_table=_v;
  __atomic_species=v_atomic_symbol_table.size();
  v_atomic_composition_table.resize(__atomic_species);
  __is_symbol_table=true;
}

void CGau::eval_atomic_number_table(void){
  uint z;
  v_atomic_number_table.resize(0);
  for( uint i=0; i<v_atomic_symbol_table.size(); i++){
    z=get_atomic_number(v_atomic_symbol_table[i]);
    v_atomic_number_table.push_back(z);
  }
}

void CGau::set_cartesian(const TMatrix<real>& _m){
   m_xyz=_m;
#ifdef _GAU_DEBUG_MESSAGES_
  std::cout<<" CGAU: coordinates: "<<m_xyz;
#endif
}

void CGau::set_atomic_symbols(const TVector<std::string>& v){
  v_atomic_symbols=v;
#ifdef _GAU_DEBUG_MESSAGES_
  std::cout<<" CGAU: symbols: "<<v_atomic_symbols;
#endif
}

//void CGau::set_direct(const TMatrix<real>& _m){
//  m_uvw=_m;
//}

//////////////////////////////
//   Coordinate functions  //
//////////////////////////////

uint CGau::get_atomic_number(std::string s){
  //transform(symbol.begin(), symbol.end(), symbol.begin(), toupper);
  //transform(s.begin(), s.begin()+1, s.begin(), toupper);
  for(uint j=0; j<periodic_table_atoms; j++){
    if(!strcmp(s.c_str(),symbol[j].c_str())){
//#ifdef _GAU_DEBUG_MESSAGES_
//      std::cout<<"CGAU: atomic number "<<j+1<<std::endl;
//#endif
      return j;
    }
  }
  return 118;
}

uint CGau::get_atomic_composition(uint i){
  return  v_atomic_composition_table[i];
}

bool CGau::is_direct(){
  return __is_direct;
}

bool CGau::is_periodic(){
  return __is_periodic;
}

uint CGau::get_total_species(void){
  return __atomic_species;
}

uint CGau::get_total_atoms(void){
  return m_xyz.rows();
}

uint CGau::get_format(void){
  return __file_format;
}

real CGau::get_cartesian(uint x, uint y){
  return  m_xyz[x][y];
}

//real CGau::get_direct(uint x, uint y){
//  return  m_uvw[x][y];
//}

TVector<std::string> CGau::get_atomic_symbols(void){
  return  v_atomic_symbols;
}

TVector<std::string> CGau::get_atomic_labels(void){
  return  v_atomic_labels;
}

TVector<std::string> CGau::get_atomic_symbol_table(void){
  return  v_atomic_symbol_table;
}

TVector<uint> CGau::get_atom_table(void){
  return  v_atom_table;
}

TVector<uint> CGau::get_atomic_composition_table(void){
  return  v_atomic_composition_table;
}

TVector<uint> CGau::get_atomic_number_table(void){
  return  v_atomic_number_table;
}

TVector<uint> CGau::get_atomic_numbers(void){
  return  v_atomic_numbers;
}

TVector<uint> CGau::get_fragment_table(void){
  return  v_fragment_table;
}

//TVector<real> CGau::get_direct(uint x){
//  return  m_uvw[x];
//}

TVector<real> CGau::get_cartesian(uint x){
  return  m_xyz[x];
}

TVector<real> CGau::get_direct_basis(void){
  return v_xyz;
}

TMatrix<real> CGau::get_xyz_input(void){
  return  m_xyz_input;
}

TMatrix<real> CGau::get_cartesian(void){
  return  m_xyz;
}

TMatrix<real> CGau::get_direct(void){
  return  m_uvw;
}

// Extra functions
void CGau::write_copyright(std::ofstream& _o){
  time_t rawtime;
  char * _st, tmp_buffer[256];
  std::string stime;
  time (&rawtime);
  sprintf(tmp_buffer,"%s",ctime(&rawtime));
  _st = strtok(tmp_buffer,"\n");
  stime = _st;
  _o<<std::endl<<std::endl<<std::endl;
  _o<<"# ------------------------------------------------------------#"<<std::endl;
  _o<<"# Gaussian input file generated by XMolView (www.xmol.org)    #"<<std::endl;
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

void CGau::center_coordinates(void){
  TVector<real> min_xyz(3);
  TVector<real> max_xyz(3);
  min_xyz[0]=m_xyz_input.col_min(0);
  min_xyz[1]=m_xyz_input.col_min(1);
  min_xyz[2]=m_xyz_input.col_min(2);
  real _vacuum=2.2*__vacuum_space;
#ifdef _GAU_DEBUG_DATA_
  std::cout<<" CXYZ: initial m_xyz = "<<m_xyz;
#endif
  for(uint i=0;i<__total_atoms;i++){
    m_xyz[i] = m_xyz_input[i] - min_xyz;
  }
#ifdef _GAU_DEBUG_DATA_
  std::cout<<" CXYZ: positive m_xyz = "<<m_xyz;
#endif
  max_xyz[0]=m_xyz.col_max(0)+_vacuum;
  max_xyz[1]=m_xyz.col_max(1)+_vacuum;
  max_xyz[2]=m_xyz.col_max(2)+_vacuum;
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
    //m_xyz[i]+=(0.5*max_xyz);
    m_xyz[i]=(m_xyz[i]+1.1*__vacuum_space);
  }
#ifdef _GAU_DEBUG_DATA_
  std::cout<<" CXYZ: boxed m_xyz = "<<m_xyz;
#endif
  unit_uvwTxyz[0]=uvwTxyz[0]/uvwTxyz[0].magnitude();
  unit_uvwTxyz[1]=uvwTxyz[1]/uvwTxyz[1].magnitude();
  unit_uvwTxyz[2]=uvwTxyz[2]/uvwTxyz[2].magnitude();
  //
#ifdef _GAU_DEBUG_DATA_
  std::cout<<" CXYZ: min_xyz = "<<min_xyz;
  std::cout<<" CXYZ: max_xyz = "<<max_xyz;
#endif
  //std::cout<<" xyz = "<<m_xyz; 
  //uvwTxyz[0][0]=1.1*max_xyz[0];
  //uvwTxyz[1][1]=1.1*max_xyz[1];
  //uvwTxyz[2][2]=1.1*max_xyz[2];
}

// END

