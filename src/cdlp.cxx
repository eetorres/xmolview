//========================================================================
// FILE: CDlp.cxx
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

#include<global.h>
#include<cdlp.h>
#include<ctype.h>
#include<algorithm>

// formats
// 1 standard dlp (i.e., avogadro)
// 2 standard including fragments
// 3 standard without number of atoms
// 4 with fragments without number of atoms

CDlp::CDlp(){
  clear();
  __lattice_constant=1.0;
  __head_lines=5;
  v_dlp.resize(3);
  v_xyz.resize(3);
  v_cells.resize(3);
  uvwTxyz.resize(3,3);
  unit_uvwTxyz.resize(3,3);
  scl_uvwTxyz.resize(3,3);
  uvwTxyz.resize(3,3);
}

void CDlp::clear(void){
  //__dlp_filename = "dlp";
  __system_title = "System Title";
  __export_format=0;
  __file_format=0;
  __is_direct=false;
  __is_symbol=true;
  __is_sort=true;
  __is_dummy=true;
  __is_labels=true;
  __is_numbers=false;
  __is_fragments=false;
  __is_charges=false;
  __is_symbol_table=false;
  __is_periodic=true;
  __is_read=false;
  __write_connections=false;
  __atomic_species=0;
  __total_atoms=0;
  // DL_POLY
  __level_of_configuration=0;
  // DL_POLY
  a1=0;
  //
  v_atomic_composition_table.resize(0);
  v_atomic_symbol_table.resize(0);
  v_atomic_number_table.resize(0);
}

bool CDlp::get_dlp_format(void){
  int res1=-1, res2=-1;//, symb=-1;
  int x1, x2, x3;
  float a, b, c, f;
  //float c, f;
  char str[20];
  std::string text_line;
  std::ifstream dlp;
  __head_lines = 0;
  try{
    dlp.open(__filename.c_str());
    if(!dlp.is_open()){
#ifdef _DLP_DEBUG_MESSAGES_
      std::cout<<"DLP: The dlp file: "<<__filename<<std::endl;
      std::cout<<"DLP: was not found\n"<<std::endl;
#endif
      return false;
    }//else{
    dlp.ignore(128,'\n');
    while(res1==-1 && !dlp.eof()){
      std::getline(dlp,text_line);
      res1 = std::sscanf((const char*)text_line.c_str(),"%i %i %i %f %f",&x1,&x2,&x3,&c,&f);
#ifdef _DLP_DEBUG_MESSAGES_
      std::cout<<" flags: "<<x1<<" "<<x2<<" "<<x3<<std::endl;
#endif
    }
    dlp.close();
    if(x2 == 0){
      __head_lines=2;
      __is_periodic=false;
      __vacuum_space=0;
#ifdef _DLP_DEBUG_MESSAGES_
      std::cout<<" non periodic system"<<std::endl;
#endif
    }else{
      __head_lines=5;
      __is_periodic=true;
#ifdef _DLP_DEBUG_MESSAGES_
      std::cout<<" periodic system"<<std::endl;
#endif
    }
    __coordinates_format=0;
    __level_of_configuration=x1;
    // number of atoms
    if(res1 < 3)
      __total_atoms=0;
    else
      __total_atoms=x3;
#ifdef _DLP_DEBUG_MESSAGES_
    std::cout<<" head lines: "<<__head_lines<<std::endl;
    std::cout<<" value1: "<<res1<<std::endl;
    if(res1 < 3)
      std::cout<<" the atoms need to be counted"<<std::endl;
#endif
    if(res1==2)      // standard format but not number of atoms
      __file_format = INPUT_FORMAT_NAT_NFR;
#ifdef _DLP_DEBUG_MESSAGES_
    std::cout<<"DLP: file format: "<<__file_format<<std::endl;
#endif
    // if the atoms need to be counted
    if(__file_format==INPUT_FORMAT_NAT_NFR || __file_format==INPUT_FORMAT_NAT_FRG){
      dlp.open(__filename.c_str());
      for(uint i=0; i<__head_lines; i++){
        dlp.ignore(1024,'\n');
      }
      __total_atoms=0;
      res1=2;
      res2=3;
#ifdef _DLP_DEBUG_MESSAGES_
      std::cout<<" value1: "<<res1<<" value2: "<<res2<<std::endl;
#endif
      while((res1 > 1 || res2 == 3) && !dlp.eof()){
        std::getline(dlp,text_line);
        //std::cout<<" line: "<<text_line<<std::endl;
        res1 = std::sscanf((const char*)text_line.c_str(),"%s %i %f",str,&x1,&a);
        dlp.ignore(0,'\n');
        std::getline(dlp,text_line);
#ifdef _DLP_DEBUG_MESSAGES_
        std::cout<<" line: "<<text_line<<std::endl;
#endif
        res2 = std::sscanf((const char*)text_line.c_str(),"%f %f %f",&a,&b,&c);
        dlp.ignore(0,'\n');
#ifdef _DLP_DEBUG_MESSAGES_
        std::cout<<" value1: "<<res1<<" value2: "<<res2<<std::endl;
#endif
        if(__level_of_configuration>0){
          dlp.ignore(1024,'\n');
        }
        if(__level_of_configuration>1){
          dlp.ignore(1024,'\n');
        }
        if(res1 > 1 || res2 == 3)
          __total_atoms++;
        //std::getline(dlp,text_line);
        //res1 = std::sscanf((const char*)text_line.c_str(),"%s %f %f %f %f",str,&a,&b,&c,&f);
      }
      dlp.close();
#ifdef _DLP_DEBUG_MESSAGES_
      std::cout<<"DLP: number of atoms found "<<__total_atoms<<std::endl;
#endif
    }
#ifdef _DLP_DEBUG_MESSAGES_
    std::cout<<"DLP: number of atoms: "<<__total_atoms<<std::endl;
    std::cout<<"DLP: The dlp file was successful closed"<<std::endl;
#endif
#ifdef _DLP_DEBUG_MESSAGES_
    //if(__total_atoms == 0)
      //exit(0);
#endif
    if(__total_atoms>0)
      return true;
    else
      return false;
  }catch(...){
    return false;
  }
}

// load the dlp file
bool CDlp::read_file(std::string d, std::string f){
  std::string atomic_symbol;
  //char tmp_buffer[256], *pch;
  real data;
  //uint unum;
  TVector<real> v1, v2, v3;
  double x=0, y=0, z=0;
  __filename = d+"/"+f;
  std::ifstream dlp;
  //v_config_format.resize(0);
  if(!get_dlp_format()){
    return false;
  }
  try{
    dlp.open(__filename.c_str());
    if(!dlp.is_open()){
#ifdef _DLP_DEBUG_MESSAGES_
      printf("The DL_POLY was not found\n");
#endif
      return false;
    }else{
#ifdef _DLP_DEBUG_MESSAGES_
      printf("The DL_POLY was successful opened\n");
#endif
      //clear();
    }

#ifdef _DLP_INFO_MESSAGES_
    std::cout<<" [Head file] "<<std::endl;
    std::cout<<" ================================================= "<<std::endl;
#endif
    for(uint _i=0;_i<__head_lines;_i++){
      dlp.getline(dlp_head_buffer[_i],256,'\n');
#ifdef _DLP_INFO_MESSAGES_
      std::cout<<"  "<<dlp_head_buffer[_i]<<std::endl;
#endif
    }
#ifdef _DLP_INFO_MESSAGES_
    std::cout<<" ================================================= "<<std::endl;
#endif
    __system_title = dlp_head_buffer[0];
#ifdef _DLP_DEBUG_MESSAGES_
    std::cout<<" DLP: coordinate format = "<<__coordinates_format<<std::endl;
    std::cout<<" DLP: configuration level = "<<__level_of_configuration<<std::endl;
#endif
    // Read latice vectors
    if(__is_periodic){
      for(uint i=2; i<5; i++){
        sscanf(dlp_head_buffer[i],"%lf %lf %lf",&x,&y,&z);
        //printf("V%i %Lf %Lf %Lf\n",i,x,y,z);
        scl_uvwTxyz[i-2][0] = x;
        scl_uvwTxyz[i-2][1] = y;
        scl_uvwTxyz[i-2][2] = z;
        //uvwTxyz[i-2]=m_uvw[i-2];
        v_xyz[i-2] = scl_uvwTxyz[i-2].magnitude();
        unit_uvwTxyz[i-2]=scl_uvwTxyz[i-2]/scl_uvwTxyz[i-2].magnitude();
      }
#ifdef _DLP_DEBUG_MESSAGES_
      //cout<<"uvwTxyz="<<uvwTxyz;
      std::cout<<"v_xyz="<<v_xyz;
      std::cout<<"Basis vectors = "<<scl_uvwTxyz;
#endif
      uvwTxyz = scl_uvwTxyz;
      //m_xyz = __lattice_constant*scl_uvwTxyz;
#ifdef _DLP_DEBUG_MESSAGES_
      std::cout<<"uvwTxyz = "<<uvwTxyz;
      //cout<<"Basis vectors XYZ = "<<uvwTxyz;
#endif
    }
    //if(__file_format == INPUT_FORMAT_ATM_NFR || __file_format == INPUT_FORMAT_ATM_FRG){
      //dlp>>data;
      //__total_atoms=data;
    //}
#ifdef _DLP_INFO_MESSAGES_
    std::cout<<" [ Configuration ] "<<std::endl;
    std::cout<<" ================================================= "<<std::endl;
    //std::cout<<"  Atomic configuration = "<<v_config_format;
    std::cout<<"  Coordinate format: "<<__coordinates_format<<std::endl;
    std::cout<<"  Level of configuration: "<<__level_of_configuration<<std::endl;
    std::cout<<"  Total atoms: "<<__total_atoms<<std::endl;
    std::cout<<" ================================================= "<<std::endl;
#endif

#ifdef _DLP_DEBUG_MESSAGES_
    std::cout<<"DLP: Total atoms: "<<__total_atoms<<std::endl;
#endif
    v_atomic_symbols.resize(__total_atoms);
    v_fragment_table.resize(__total_atoms);
    v_atom_table.resize(__total_atoms);
    m_xyz_input.resize(__total_atoms,3);
    m_xyz.resize(__total_atoms,3);
    m_uvw.resize(__total_atoms,3);
    if(__level_of_configuration>0){
      m_velocities_xyz.resize(__total_atoms,3);
    }
    if(__level_of_configuration>1){
      m_accelerations_xyz.resize(__total_atoms,3);
    }
    //while(dlp.peek()=='\n')
      //dlp.ignore(0,'\n');
#ifdef _DLP_INFO_MESSAGES_
    std::cout<<"DLP: loading atomic coordinates"<<std::endl;
#endif
    int n;
#ifdef _DLP_INFO_MESSAGES_
    std::cout<<" [loading atomic coordinates] "<<std::endl;
#endif
    for(uint f=0; f<__total_atoms; f++){
      dlp>>atomic_symbol;
#ifdef _DLP_INFO_MESSAGES_
      std::cout<<"DLP: atomic symbol "<<atomic_symbol<<std::endl;
#endif
      // make sure it uses the standar atom symbols
      //transform(atomic_symbol.begin(), atomic_symbol.begin()+1, atomic_symbol.begin(), toupper);
      //transform(atomic_symbol.begin()+1, atomic_symbol.end(), atomic_symbol.begin()+1, tolower);
      v_atomic_symbols[f]=atomic_symbol;
      dlp>>n;
      if(dlp.peek()!='\n')
          dlp.ignore(1024,'\n');
      //printf("%i-",n);
      for(int c=0; c<3; c++){
        dlp>>data;
        //if(data<0) data = (data+1.0);
        //printf("%f-",data);
        m_xyz_input[f][c]=data;
      }
      if(__level_of_configuration>0){
        if(dlp.peek()!='\n')
          dlp.ignore(128,'\n');
        for(int c=0; c<3; c++){
          dlp>>data;
          //m_velocities_xyz[f][c]=data;
        }
      }
      if(__level_of_configuration>1){
        if(dlp.peek()!='\n')
          dlp.ignore(128,'\n');
        for(int c=0; c<3; c++){
          dlp>>data;
          //m_accelerations_xyz[f][c]=data;
        }
      }
      if(dlp.peek()!='\n')
        dlp.ignore(128,'\n');
    }
    v_atomic_labels=v_atomic_symbols;
    //std::cout<<v_atomic_labels;
    eval_atomic_composition();
    center_coordinates();
    //eval_atomic_number_table();
    //m_uvw=m_uvw-0.5;
    //////compute_general();
#ifdef _DLP_DEBUG_MESSAGES_
    std::cout<<"--- Summary of the dlp file ---"<<std::endl;
    std::cout<<"atom table ="<<v_atom_table;
    std::cout<<"atomic symbol table ="<<v_atomic_symbol_table;
    std::cout<<"atomic number table ="<<v_atomic_number_table;
    std::cout<<"atomic composition table ="<<v_atomic_composition_table;
    std::cout<<"atomic labels ="<<v_atomic_labels;
    std::cout<<"atomic symbols ="<<v_atomic_symbols;
    std::cout<<"atomic numbers ="<<v_atomic_numbers;
    std::cout<<"uvwTxyz'="<<uvwTxyz;
    std::cout<<"--- Summary of the dlp file ---"<<std::endl;
#endif

#ifdef _DLP_DEBUG_DATA_
    std::cout<<"Cartesian'="<<m_xyz;
#endif
    if(dlp.is_open()){
      dlp.close();
    }
#ifdef _DLP_DEBUG_MESSAGES_
    std::cout<<"DLP: The dlp was successful closed!"<<std::endl;
#endif
    __is_read=true;
    return true;
  }catch(...){
    return false;
  }
}

void CDlp::eval_atomic_composition(void){
  std::string _tstr;
  uint cont=0;
  size_t _sw;
  bool _found = false;
  uint _an;
  std::string tmp_symbol;
  v_atomic_symbol_table.resize(0);
  v_atomic_symbol_size.resize(0);
  v_atomic_number_table.resize(0);
  TVector<real> _tv;
  TVector<bool> v_found(__total_atoms);
  v_found.constant(false);
  if(!__is_symbol_table){
    for(uint j=0; j<periodic_table_atoms; j++){
      if(j < 95) _sw = 2;
      else _sw = 1;
      _found = false;
      for(uint f=0;f<__total_atoms;f++){
        if(!v_found[f]){
          tmp_symbol=v_atomic_symbols[f];
#ifdef _DLP_COMPOSITON_MESSAGES_
          std::cout<<" DLP: tmp_symbol: "<<tmp_symbol<<std::endl;
#endif
          transform(tmp_symbol.begin(), tmp_symbol.begin()+1, tmp_symbol.begin(), toupper);
          transform(tmp_symbol.begin()+1, tmp_symbol.end(), tmp_symbol.begin()+1, tolower);
          if(strncmp(tmp_symbol.c_str(),symbol_sorted[j].c_str(),_sw)==0){
            v_atomic_symbols[f]=symbol_sorted[j];
            v_atomic_symbol_size.push_back(_sw);
            if(!_found){
              v_atomic_symbol_table.push_back(symbol_sorted[j]);
              _an = get_atomic_number(symbol_sorted[j],_sw);
              v_atomic_number_table.push_back(_an);
              if(!__is_periodic){
                __vacuum_space=maxi(__vacuum_space,atom_rrgb[_an][0]);
              }
              _found=true;
            }
            v_found[f]=true;
            cont++;
            if(cont==__total_atoms) break;
          }
          if(cont==__total_atoms) break;
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
      v_atomic_numbers[j]=get_atomic_number(v_atomic_symbols[j],v_atomic_symbol_size[i]);
      if(strncmp(v_atomic_symbols[j].c_str(),v_atomic_symbol_table[i].c_str(),v_atomic_symbol_size[i])==0){
        v_atom_table[j]=i;
        _c++;
      }
    }
    v_atomic_composition_table[i]=_c;
  }
}

void CDlp::write_file(void){
  write_file(__filename);
}

void CDlp::write_file(std::string d, std::string f){
  std::string _f=d+'/'+f;
  write_file(_f);
  if(__write_connections)
    write_structure_file(d,f);
}

//  save the dlp file
void CDlp::write_file(std::string _fn){
  std::string symbol;
  real data, val;
  TVector<real> v_scale(3);
  std::ofstream odlp;
  odlp.open(_fn.c_str());
  //odlp.width(17);
  if(odlp.bad())
    std::cout<<" DLP: The file was not created"<<std::endl;
#ifdef _DLP_INFO_MESSAGES_
  std::cout<<" DLP: file format: "<<__file_format<<std::endl;
#endif
  //if(__export_format==OUTPUT_FORMAT_ATM_NFR || __export_format==OUTPUT_FORMAT_ATM_FRG){
  std::string stime = get_date();
  if(__is_read){
	odlp<<dlp_head_buffer[0]<<std::endl;
  }else{
    odlp<<"DL_POLY structure file generated by XMolView (www.xmol.org) ["<<stime<<"]"<<std::endl;
  }
  odlp<<"         0         ";
  if(__is_periodic){
    odlp<<"2         "<<__total_cells*__total_atoms<<std::endl;
    odlp.precision(12);
    for(uint f=0; f<3; f++){
    for(uint c=0; c<3; c++){
        //odlp<<" ";
        val=uvwTxyz[f][c]*v_cells[f];
        //if(val>=0) odlp<<" ";
        odlp.width(20);
        odlp<<std::fixed<<std::right<<val;
        //oposcar<<" ";
      }
      //if(f<2) odlp<<"\n ";
      odlp<<"\n";
    }
  }else{
	odlp<<"0         "<<__total_cells*__total_atoms<<std::endl;
  }
  //odlp<<std::endl;
  //}
  char _buff[5];
  TVector<uint> vcont(__total_fragments);
  odlp.precision(12);
  for(uint cell=0; cell<__total_cells; cell++){ // repetition cells
    for(uint f=0;f<__total_atoms;f++){
      //if(__labels_format)
      symbol=v_atomic_symbols[f];
      if(strcmp(symbol.c_str(),"X") || __is_dummy){
        odlp.width(8);
        if(__is_numbers){
          if(__is_fragments){
            //v_fragment_table=_v;
            //__total_fragments=u;
            vcont[v_fragment_table[f]-1]+= 1;
            sprintf(_buff,"%i",vcont[v_fragment_table[f]-1]);
          }else{
            sprintf(_buff,"%i",1+f+cell*__total_atoms);
          }
          //sprintf(_buff,"%i",f);
          symbol = symbol + _buff;
        }
        odlp<<std::fixed<<std::left<<symbol;
        odlp.width(10);
        odlp<<std::fixed<<std::right<<1+f+cell*__total_atoms<<std::endl;
        for(uint c=0;c<3;c++){
          data=m_xyz[f+cell*__total_atoms][c];
          odlp.width(17);
          odlp<<std::fixed<<std::right<<data<<"    ";
        }
      // include fragments
//      if(__export_format==OUTPUT_FORMAT_ATM_FRG || __export_format==OUTPUT_FORMAT_NAT_FRG){
//#ifdef _DLP_INFO_MESSAGES_
//        std::cout<<"DLP: include fragments"<<std::endl;
//#endif
//        odlp.width(5);
//        odlp<<std::fixed<<std::right<<v_fragment_table[f];
//      }
        odlp<<std::endl;
      }
    }
  }
  //if(__export_format==0)
    //write_copyright(odlp);
  odlp.close();
}

void CDlp::write_structure_file(std::string d, std::string f){
  bool cnn=true;
  std::string symbol;
  std::ofstream osf;
  std::string _fn;
  char _buff[5];
  size_t pos = strcspn(f.c_str(),".");
  _fn = f.substr(0,pos);
  _fn += ".sf";
  std::string _sf = d+'/'+_fn;
  osf.open(_sf.c_str());
  std::cout<<"DLP: write connections to "<<_sf<<std::endl;
  osf<<std::endl;
  osf<<"# DL_FIELD Structure"<<std::endl;
  osf<<"MOLECULE Molecule_Name "<<__total_atoms<<" 0.0"<<std::endl;
  for(uint f=0;f<__total_atoms;f++){
    //if(strcmp(symbol.c_str(),"X")){
      symbol = v_atomic_symbols[f];
      if(__is_numbers){
        sprintf(_buff,"%i",1+f);
        symbol = symbol + _buff;
      }
      osf.width(6);
      osf<<std::left<<symbol<<" "<<symbol<<"_type  ";
      if(__is_charges){
        osf.width(8);
        osf.precision(4);
        osf<<std::fixed<<std::right<<v_atomic_charges[f]<<std::endl;
      }else osf<<"0.0"<<std::endl;
    //}
  }
  int s = v_connection.size();
  for(int i=0; i<s; i++){
    if(cnn) osf<<"CONNECT"<<" ";
    if(v_connection[i]<0){
	  osf<<std::endl;
	  cnn=true;
	}else{
	  if(cnn){
	    osf.width(3);
	    symbol = v_atomic_symbols[v_connection[i]];
	    if(__is_numbers){
		  sprintf(_buff,"%i",1+v_connection[i]);
		  symbol = symbol + _buff;
	    }
		osf<<std::fixed<<std::left<<symbol<<" ";
	  }else{
		symbol = v_atomic_symbols[v_connection[i]];
		if(__is_numbers){
		  sprintf(_buff,"%i",1+v_connection[i]);
		  symbol = symbol + _buff;
	    }
		osf<<" "<<symbol;
	  }
	  if(cnn) osf<<"> ";
	  cnn=false;
	}
  }
  osf<<"END MOLECULE Molecule_Name"<<std::endl;
  osf.close();
}

// SET FUNCTIONS

void CDlp::set_file_format(const uint u){
  __file_format=u;
}

void CDlp::set_export_format(const uint u){
  __export_format=u;
}

void CDlp::set_labels_format(const bool b){
  __labels_format=b;
}

void CDlp::set_periodic(const bool b){
  __is_periodic=b;
}

void CDlp::set_total_cells(const uint x, const uint y, const uint z, const uint t){
  v_cells[0]=x;
  v_cells[1]=y;
  v_cells[2]=z;
  __total_cells=t;
}

void CDlp::set_fragments(const TVector<uint>& _v, uint u){
  v_fragment_table=_v;
  __total_fragments=u;
  __is_fragments=true;
}

void CDlp::set_charges(TVector<real>& _v){
  __is_charges=true;
  v_atomic_charges=_v;
}

void CDlp::is_connections(bool b){
  __write_connections=b;
}

void CDlp::set_connections(const TVector<int>& _v){
  v_connection=_v;
  is_connections(true);
}

void CDlp::set_atomic_symbol_table(const TVector<std::string>& _v){
  v_atomic_symbol_table=_v;
  __atomic_species=v_atomic_symbol_table.size();
  v_atomic_composition_table.resize(__atomic_species);
  __is_symbol_table=true;
}

/* Fri May 18 11:39:33 MDT 2012
   deprecated, this is done in eval_composition
void CDlp::eval_atomic_number_table(void){
  uint z;
  v_atomic_number_table.resize(0);
  for( uint i=0; i<v_atomic_symbol_table.size(); i++){
    z=get_atomic_number(v_atomic_symbol_table[i]);
    v_atomic_number_table.push_back(z);
  }
}*/

void CDlp::set_cartesian(const TMatrix<real>& _m){
   m_xyz=_m;
}

void CDlp::set_direct(const TMatrix<real>& _m){
  m_uvw=_m;
}

//////////////////////////////
//   Coordinate functions  //
//////////////////////////////

uint CDlp::get_atomic_number(std::string s, size_t z){
  for(uint j=0; j<periodic_table_atoms; j++){
    if(strncmp(s.c_str(),symbol[j].c_str(),z)==0){
//#ifdef _DLP_DEBUG_MESSAGES_
//      std::cout<<"DLP: atomic number "<<j+1<<std::endl;
//#endif
      return j;
    }
  }
  return 118;
}

uint CDlp::get_atomic_composition(uint i){
  return  v_atomic_composition_table[i];
}

bool CDlp::is_direct(){
  return __is_direct;
}

bool CDlp::is_periodic(){
  return __is_periodic;
}

bool CDlp::is_labels(void){
  return __is_labels;
}

bool CDlp::is_charges(void){
  return __is_charges;
}

uint CDlp::get_total_species(void){
  return __atomic_species;
}

uint CDlp::get_total_atoms(void){
  return m_xyz.rows();
}

uint CDlp::get_format(void){
  return __file_format;
}

real CDlp::get_cartesian(uint x, uint y){
  return  m_xyz[x][y];
}

real CDlp::get_direct(uint x, uint y){
  return  m_uvw[x][y];
}

TVector<std::string> CDlp::get_atomic_symbols(void){
  return  v_atomic_symbols;
}

TVector<std::string> CDlp::get_atomic_labels(void){
  return  v_atomic_labels;
}

TVector<std::string> CDlp::get_atomic_symbol_table(void){
  return  v_atomic_symbol_table;
}

TVector<uint> CDlp::get_atom_table(void){
  return  v_atom_table;
}

TVector<uint> CDlp::get_atomic_composition_table(void){
  return  v_atomic_composition_table;
}

TVector<uint> CDlp::get_atomic_number_table(void){
  return  v_atomic_number_table;
}

TVector<uint> CDlp::get_atomic_numbers(void){
  return  v_atomic_numbers;
}

TVector<uint> CDlp::get_fragment_table(void){
  return  v_fragment_table;
}

TVector<real> CDlp::get_direct(uint x){
  return  m_uvw[x];
}

TVector<real> CDlp::get_cartesian(uint x){
  return  m_xyz[x];
}

TVector<real> CDlp::get_direct_basis(void){
  return v_dlp;
}

TMatrix<real> CDlp::get_xyz_input(void){
  return  m_xyz_input;
}

TMatrix<real> CDlp::get_cartesian(void){
  return  m_xyz;
}

TMatrix<real> CDlp::get_direct(void){
  return  m_uvw;
}

// Extra functions
void CDlp::write_copyright(std::ofstream& _o){
  time_t rawtime;
  char * _st, tmp_buffer[256];
  std::string stime;
  time (&rawtime);
  sprintf(tmp_buffer,"%s",ctime(&rawtime));
  _st = strtok(tmp_buffer,"\n");
  stime = _st;
  _o<<std::endl<<std::endl<<std::endl;
  _o<<"# ------------------------------------------------------------#"<<std::endl;
  _o<<"# DL_POLY input file generated by xmolview (www.xmol.org)     #"<<std::endl;
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

void CDlp::center_coordinates(void){
  TVector<real> min_dlp(3);
  TVector<real> max_dlp(3);
  real _vacuum = 2.2*__vacuum_space;
  if(__is_periodic){
	m_xyz = m_xyz_input;
    for(uint i=0;i<__total_atoms;i++){
      for(uint j=0; j<3; j++){
        m_xyz[i] += 0.5*(uvwTxyz[j]);
      }
    }
  }else{
    min_dlp[0]=m_xyz_input.col_min(0);
    min_dlp[1]=m_xyz_input.col_min(1);
    min_dlp[2]=m_xyz_input.col_min(2);
#ifdef _DLP_DEBUG_DATA_
    std::cout<<" DLP: initial m_xyz = "<<m_xyz;
#endif
    for(uint i=0;i<__total_atoms;i++){
        m_xyz[i] = m_xyz_input[i] - min_dlp;
    }
#ifdef _DLP_DEBUG_DATA_
    std::cout<<" DLP: positive m_xyz = "<<m_xyz;
#endif
    max_dlp[0]=m_xyz.col_max(0)+_vacuum;
    max_dlp[1]=m_xyz.col_max(1)+_vacuum;
    max_dlp[2]=m_xyz.col_max(2)+_vacuum;
    if(max_dlp[0]>0.01)
      uvwTxyz[0][0]=max_dlp[0];
    else{
      uvwTxyz[0][0]=+_vacuum;;
      max_dlp[0]=+_vacuum;;
    }
    if(max_dlp[1]>0.01)
      uvwTxyz[1][1]=max_dlp[1];
    else{
      uvwTxyz[1][1]=+_vacuum;;
      max_dlp[1]=+_vacuum;;
    }
    if(max_dlp[2]>0.01)
      uvwTxyz[2][2]=max_dlp[2];
    else{
      uvwTxyz[2][2]=+_vacuum;;
      max_dlp[2]=+_vacuum;;
    }
    for(uint i=0;i<__total_atoms;i++){
      //m_xyz[i]+=(0.5*max_dlp);
      m_xyz[i]=(m_xyz[i]+1.1*__vacuum_space);
    }
#ifdef _DLP_DEBUG_DATA_
    std::cout<<" DLP: boxed m_xyz = "<<m_xyz;
#endif
    unit_uvwTxyz[0]=uvwTxyz[0]/uvwTxyz[0].magnitude();
    unit_uvwTxyz[1]=uvwTxyz[1]/uvwTxyz[1].magnitude();
    unit_uvwTxyz[2]=uvwTxyz[2]/uvwTxyz[2].magnitude();
#ifdef _DLP_DEBUG_DATA_
    std::cout<<" DLP: min_dlp = "<<min_dlp;
    std::cout<<" DLP: max_dlp = "<<max_dlp;
#endif
  }
}

// END

