//========================================================================
// FILE: cpdb.cxx
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

//#include<const.h>
#include<cpdb.h>
#include<ctype.h>
#include<algorithm>

// formats
// 1 standard PDB no fragments
// 2 including fragments ?

CPdb::CPdb(){
  clear();
}

void CPdb::clear(void){
  //__xyz_filename = "xyz";
  __system_title = "System Title";
  __export_format=0;
  __file_format=0;
  __header_lines=0;
  //__is_selected_dynamics=false;
  __is_direct=false;
  __is_periodic=false;
  __is_sort=true;
  __is_dummy=false;
  __is_labels=false;
  __is_numbers=false;
  __is_fragments=false;
  __is_charges=false;
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

bool CPdb::get_pdb_format(void){
  int res1, res2, res3, res4, res5;
  int k, l;
  //float a, b, c, d;
  float a, b, c, e;
  float alpha, beta, gamma;
  char str[80], stra[80], strb[80], str2[80];
  std::string text_line;
  std::ifstream pdb;
  str2[0]='X';
  try{
    pdb.open(__filename.c_str());
    if(!pdb.is_open()){
#ifdef _PDB_DEBUG_MESSAGES_
      std::cout<<"CPDB: The PDB file: "<<__filename<<std::endl;
      std::cout<<"CPDB: was not found\n"<<std::endl;
#endif
      return false;
    }//else{
    clear();
    res1=0;
    res2=0;
    __header_lines=0;
    while(strcmp("ATOM",str2) && strcmp("HETATM",str2) && !pdb.eof()){
      //pdb.ignore(256,'\n');
      std::getline(pdb,text_line);
      res1 = std::sscanf((const char*)text_line.c_str(),"%s",str2);
      if(!strcmp("CRYST1",str2)){
#ifdef _PDB_DEBUG_MESSAGES_
        std::cout<<" CRYST1 found"<<std::endl;
        std::cout<<" res1: "<<res1<<std::endl;
#endif
        res1 = std::sscanf((const char*)text_line.c_str(),"%*s %f %f %f %f %f %f",&a,&b,&c,&alpha,&beta,&gamma);
#ifdef _PDB_DEBUG_MESSAGES_
        std::cout<<" (a,b,c) =("<<a<<","<<b<<","<<c<<")"<<std::endl;
        std::cout<<" (alpha,beta,gamma) =("<<alpha<<","<<beta<<","<<gamma<<")"<<std::endl;
#endif
        if ( a > 0 && b > 0 && c > 0){
          // Compute lattice vectors coordinates in a cartesian frame.
          // Vector a is along the x axes, vector b is in the (x, y) plane.
          // conversion en radian
          __is_periodic = true;
          float alphar = alpha * DEG_RAD;
          //float betar  = beta  * DEG_RAD;
          float gammar = gamma * DEG_RAD;
          // calcul des vecteurs
           // lattice vector 1
          scl_uvwTxyz[0][0] = a;
          // lattice vector 2
          if(gammar != 90){
            scl_uvwTxyz[1][0] = b*cos(gammar);
            scl_uvwTxyz[1][1] = b*sin(gammar);
          }else{
            scl_uvwTxyz[1][1] = b;
          }
          // lattice vector 3
          if(beta != 90 || gamma != 90){
            //float yy = (cos(alphar) - cos(gammar) * cos(betar)) / (sin(gammar) * sin(betar));
            scl_uvwTxyz[2][0] = c*cos(alphar)*cos(gammar);
            scl_uvwTxyz[2][1] = c*cos(alphar)*sin(gammar);
            scl_uvwTxyz[2][2] = c*sin(alphar);
            //scl_uvwTxyz[2][0] = c*cos(betar);
            //real cy = (cos(alphar)-cos(gammar)*cos(betar))/sin(gammar);
            //scl_uvwTxyz[2][1] = c*cy;
            //real cz = sqrt((sin(betar)*sin(betar))-cy*cy);
            //scl_uvwTxyz[2][2] = c*cz;
          }else{
            scl_uvwTxyz[2][2] = c;
          }
          uvwTxyz=scl_uvwTxyz;
          for(uint i=0; i<3; i++){
            v_xyz[i]=scl_uvwTxyz[i].magnitude();
            unit_uvwTxyz[i]=scl_uvwTxyz[i]/v_xyz[i];
          }
#ifdef _PDB_DEBUG_MESSAGES_
          std::cout<<"scl_uvwTxyz="<<scl_uvwTxyz;
          std::cout<<"unit_uvwTxyz="<<unit_uvwTxyz;
          std::cout<<"uvwTxyz="<<uvwTxyz;
          std::cout<<"v_xyz="<<v_xyz;
#endif
          /*
          self._vecc = [0.,0.,0.]
          self._vecc[0] = self._c * cos( betar )
          cy = ( cos(alphar) - cos(gammar) * cos(betar) ) / sin(gammar)
          self._vecc[1] = self._c * cy
          cz = sqrt( ( sin( betar ) )**2 - cy**2 )
          self._vecc[2] = self._c * cz
          */
        }
      }
      __header_lines++;
#ifdef _PDB_DEBUG_MESSAGES_
      std::cout<<" str2 -> "<<str2<<std::endl;
#endif
    }
    while(pdb.eof()){
      return false;
    }
    //while(res1!=6){
    //pdb.ignore(256,'\n');
    //std::getline(pdb,text_line);
    //std::cout<<" text_line: "<<text_line<<std::endl;
    res1 = std::sscanf((const char*)text_line.c_str(),"%*s %i %s %*s %i %f %f %f %e",&k,str,&l,&a,&b,&c,&e);
#ifdef _PDB_DEBUG_MESSAGES_
    std::cout<<" res1: "<<res1<<std::endl;
#endif
    res2 = std::sscanf((const char*)text_line.c_str(),"%*s %i %s %i %f %f %f",&k,str,&l,&a,&b,&c);
#ifdef _PDB_DEBUG_MESSAGES_
    std::cout<<" res2: "<<res2<<std::endl;
#endif
    res3 = std::sscanf((const char*)text_line.c_str(),"%*s %i %s %f %f %f",&k,str,&a,&b,&c);
#ifdef _PDB_DEBUG_MESSAGES_
    std::cout<<" res3: "<<res3<<std::endl;
#endif
    res4 = std::sscanf((const char*)text_line.c_str(),"%*s %i %s %s %i %f %f %f %e",&k,str,strb,&l,&a,&b,&c,&e);
    //res5 = std::sscanf((const char*)text_line.c_str(),"%*s %i %s %s %i %f %f %f %e",&k,str,strb,&l,&a,&b,&c,&e);
#ifdef _PDB_DEBUG_MESSAGES_
    std::cout<<" res4: "<<res4<<std::endl;
#endif
    res5 = std::sscanf((const char*)text_line.c_str(),"%*s %i %s %s %s %i %f %f %f %e",&k,str,stra,strb,&l,&a,&b,&c,&e);
    //res4 = std::sscanf((const char*)text_line.c_str(),"%*s %i %s %s %i %f %f %f %e",&k,str,stra,strb,&l,&a,&b,&c,&e);
    //__header_lines++;
#ifdef _PDB_DEBUG_MESSAGES_
    std::cout<<" res5: "<<res5<<std::endl;
#endif
    //}
    //
    //#ifdef _PDB_DEBUG_MESSAGES_
    //  std::cout<<" res1: "<<res1<<std::endl;
    //  std::cout<<" res2: "<<res2<<std::endl;
    //#endif
    if(res1==6){                 // standard format
      if(res2 == 2 && res3 == 2) __file_format = 1;
      if(res2 == 6 && res3 == 5) __file_format = 2;
    }else if(res1==7){
      if(res2 == 6){
        __file_format = 2;
      }else if(res4==8){
        __file_format = 4;
        __is_charges=true;
      }
    }else if(res1==2 && res2 == 6 && res3 == 5){           // standard format
      __file_format = 3;
    }else if(res5==9){           // full format
      __file_format = 5;
    }else{
      return false;
    }
#ifdef _PDB_DEBUG_MESSAGES_
    std::cout<<"CPDB: file format: "<<__file_format<<std::endl;
    std::cout<<"CPDB: atoms need to be counted"<<std::endl;
#endif
    // atoms need to be counted
    //pdb.open(__filename.c_str());
    __total_atoms=0;
    if(__file_format >= 1){
      while((((!strcmp("ATOM",str2) || !strcmp("HETATM",str2))) && !pdb.eof() && res1>0) && !pdb.eof()){
        __total_atoms++;
        std::getline(pdb,text_line);
        str2[0]='X';
        res1 = std::sscanf((const char*)text_line.c_str(),"%s",str2);
#ifdef _PDB_DEBUG_MESSAGES_
        std::cout<<" atoms: "<<__total_atoms<<std::endl;
        std::cout<<" res1: "<<res1<<std::endl;
        std::cout<<" str2: "<<str2<<std::endl;
#endif
      }
    }
    //
    if(pdb.is_open())
      pdb.close();
#ifdef _PDB_DEBUG_MESSAGES_
    std::cout<<"CPDB: number of atoms found "<<__total_atoms<<std::endl;
#endif
#ifdef _PDB_DEBUG_MESSAGES_
    std::cout<<"CPDB: The PDB file was successful closed"<<std::endl;
#endif
    if(__total_atoms>0)
      return true;
    else
      return false;
  }catch(...){
    return false;
  }
}

// load the PDB file
bool CPdb::read_file(std::string d, std::string f){
  std::string symbol;
  //int res1;
  int k, l;
  float a, b, c, e;
  char str[16];
  std::string text_line;
  //real data;
  //uint unum;
  TVector<real> v1, v2, v3;
  __filename = d+"/"+f;
  std::ifstream pdb;
  //
  try{
    if(!get_pdb_format()){
      return false;
    }
    //
    pdb.open(__filename.c_str());
    //if(__file_format == 1 || __file_format == 2){
      //pdb>>data;
      //__total_atoms=data;
    //}
#ifdef _PDB_DEBUG_MESSAGES_
    std::cout<<"CPDB: Total atoms: "<<__total_atoms<<std::endl;
#endif
    v_atomic_symbols.resize(__total_atoms);
    v_fragment_table.resize(__total_atoms);
    v_atom_table.resize(__total_atoms);
    m_xyz_input.resize(__total_atoms,3);
    m_xyz.resize(__total_atoms,3);
    m_uvw.resize(__total_atoms,3);
    if(__file_format == 4){
      v_atomic_charges.resize(__total_atoms);
    }
    //while(xyz.peek()=='\n')
    //xyz.ignore(0,'\n');
    for(uint i=0; i<__header_lines-1; i++){
      //pdb.ignore(256,'\n');
      std::getline(pdb,text_line);
#ifdef _PDB_DEBUG_MESSAGES_
      std::cout<<"CPDB: head line: "<<text_line<<std::endl;
#endif
      //res1 = std::sscanf((const char*)text_line.c_str(),"%i %i %*s",str,&k,&l);
      //__header_lines++;
    }
    //pdb.ignore(0,'\n');
    //pdb.ignore(0,'\n');
#ifdef _PDB_INFO_MESSAGES_
    std::cout<<"CPDB: loading atomic coordinates"<<std::endl;
#endif
    //int res1;
    for(uint f=0;f<__total_atoms;f++){
      std::getline(pdb,text_line);
#ifdef _PDB_SHOW_DATA_
      std::cout<<"CPDB: "<<text_line<<std::endl;
#endif
      if(__file_format == 1){
        //res1 = std::sscanf((const char*)text_line.c_str(),"%*s %i %s %*s %i %f %f %f %*s",&k,str,&l,&a,&b,&c);
        std::sscanf((const char*)text_line.c_str(),"%*s %i %s %*s %i %f %f %f %*s",&k,str,&l,&a,&b,&c);
#ifdef _PDB_SHOW_DATA_
        std::cout<<"CPDB: "<<k<<" "<<str<<" "<<l<<" "<<a<<" "<<b<<" "<<c<<std::endl;
#endif
      }else if(__file_format == 2){
        //res1 = std::sscanf((const char*)text_line.c_str(),"%*s %i %s %i %f %f %f %*s",&k,str,&l,&a,&b,&c);
        std::sscanf((const char*)text_line.c_str(),"%*s %i %s %i %f %f %f %*s",&k,str,&l,&a,&b,&c);
#ifdef _PDB_SHOW_DATA_
        std::cout<<"CPDB: "<<k<<" "<<str<<" "<<l<<" "<<a<<" "<<b<<" "<<c<<std::endl;
#endif
      }else if(__file_format == 3){
        std::sscanf((const char*)text_line.c_str(),"%*s %i %s %f %f %f",&k,str,&a,&b,&c);
        //std::sscanf((const char*)text_line.c_str(),"%*s %i %s %f %f %f",&k,str,&a,&b,&c);
      }else if(__file_format == 4){
        //res1 = std::sscanf((const char*)text_line.c_str(),"%*s %i %s %*s %*s %f %f %f %f",&k,str,&a,&b,&c,&e);
        std::sscanf((const char*)text_line.c_str(),"%*s %i %s %*s %*s %f %f %f %f",&k,str,&a,&b,&c,&e);
#ifdef _PDB_SHOW_DATA_
        std::cout<<"CPDB: "<<k<<" str="<<str<<" a="<<a<<" b="<<b<<" c="<<c<<" e="<<e<<std::endl;
#endif
      }else if(__file_format == 5){
        //res1 = std::sscanf((const char*)text_line.c_str(),"%*s %i %s %*s %*s %*s %f %f %f %f",&k,str,&a,&b,&c,&e);
        std::sscanf((const char*)text_line.c_str(),"%*s %i %s %*s %*s %*s %f %f %f %f",&k,str,&a,&b,&c,&e);
#ifdef _PDB_SHOW_DATA_
        std::cout<<"CPDB: "<<k<<" str="<<str<<" a="<<a<<" b="<<b<<" c="<<c<<" e="<<e<<std::endl;
#endif
      }
      //pdb>>symbol;
      symbol = str;
      m_xyz_input[f][0]=a;
      m_xyz_input[f][1]=b;
      m_xyz_input[f][2]=c;
      if(__is_charges){
        v_atomic_charges[f]=e;
      }
      // make sure it uses the standard atom symbols
      //transform(symbol.begin(), symbol.begin()+1, symbol.begin(), toupper);
      //transform(symbol.begin()+1, symbol.end(), symbol.begin()+1, tolower);
      v_atomic_symbols[f]=symbol;
#ifdef _PDB_SHOW_DATA_
        std::cout<<"CPDB: "<<k<<" symbol="<<str<<" a="<<a<<" b="<<b<<" c="<<c<<std::endl;
#endif
      //for(uint c=0;c<3;c++){  // read coordinates
        //pdb>>data;
        //m_xyz[f][c]=data;
      //}
      //if(__file_format==2){   // read fragments
        //pdb>>unum;
        //v_fragment_table[f]=unum;
      //}
      if(pdb.peek()!='\n')
        pdb.ignore(0,'\n');
#ifdef _PDB_SHOW_DATA_
        std::cout<<"CPDB: loop end"<<std::endl;
#endif
    }
    v_atomic_labels=v_atomic_symbols;
    eval_atomic_composition();
    center_coordinates();
    //eval_atomic_number_table();
    //m_uvw=m_uvw-0.5;
    //////compute_general();
#ifdef _PDB_DEBUG_MESSAGES_
    std::cout<<" PDB: atomic labels="<<v_atomic_labels;
    std::cout<<"--- Summary of the PDB file ---"<<std::endl;
    std::cout<<"atom table ="<<v_atom_table;
    std::cout<<"atomic symbol table ="<<v_atomic_symbol_table;
    std::cout<<"atomic number table ="<<v_atomic_number_table;
    std::cout<<"atomic composition table ="<<v_atomic_composition_table;
    std::cout<<"atomic symbols ="<<v_atomic_symbols;
    std::cout<<"atomic numbers ="<<v_atomic_numbers;
    std::cout<<"uvwTxyz="<<uvwTxyz;
    std::cout<<"--- Summary of the PDB file ---"<<std::endl;
#endif

#ifdef _PDB_SHOW_DATA_
    std::cout<<"Cartesian'="<<m_xyz;
#endif
    if(pdb.is_open()){
      pdb.close();
    }
#ifdef _PDB_DEBUG_MESSAGES_
    std::cout<<"CPDB: The PDB was successful closed!"<<std::endl;
#endif
    if((v_atomic_labels.size()>0) && (v_atom_table.size()>0) && (v_atomic_symbols.size()>0) && (v_atomic_numbers.size()>0))
      return true;
    else
      return false;
  }catch(...){
    return false;
  }
}

void CPdb::eval_atomic_composition(void){
  std::string _tstr;
  uint cont = 0;
  size_t _sw;
  bool _found=false;
  uint _an;
  std::string symbol;
  v_atomic_symbol_table.resize(0);
  v_atomic_symbol_size.resize(0);
  v_atomic_number_table.resize(0);
  TVector<real> _tv;
  TVector<bool> v_found(__total_atoms);
  __vacuum_space=0;
  v_found.constant(false);
  // check this problem in Gau and XYZ
  if(!__is_symbol_table){
    for(uint j=0; j<periodic_table_atoms; j++){
      if(j < 95) _sw = 2;
      else _sw = 1;
      _found=false;
      for(uint f=0;f<__total_atoms;f++){
        if(!v_found[f]){
		  // tested with dl_filed example files
          symbol=v_atomic_symbols[f];
          transform(symbol.begin(), symbol.end(), symbol.begin(), tolower);
          std::size_t found = symbol.find_first_of("abcdefghijklmnopqrstuvwxyz");
          std::string str2 = symbol.substr (found,2);
          //std::size_t found = symbol.find(symbol_sorted[j]);
          //if(found==std::string::npos){
          //transform(str2.begin(), str2.begin()+1, str2.begin(), toupper);
          str2[0]=toupper(str2[0]);
          //std::cout<<str2<<std::endl;
          //transform(str2.begin()+1, str2.end(), str2.begin()+1, tolower);
          //found = symbol.find(symbol_sorted[j]);
	      //}
          if(strncmp(str2.c_str(),symbol_sorted[j].c_str(),_sw)==0){

		  //if(found!=std::string::npos){
            v_atomic_symbols[f]=symbol_sorted[j];
            v_atomic_symbol_size.push_back(_sw);
            if(!_found){
              v_atomic_symbol_table.push_back(symbol_sorted[j]);
              _an = get_atomic_number(symbol_sorted[j],_sw);
              v_atomic_number_table.push_back(_an);
              __vacuum_space=maxi(__vacuum_space,atom_rrgb[_an][0]);
              _found=true;
            }
            v_found[f]=true;
            cont++;
            if(cont==__total_atoms) break;
            //break;
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

void CPdb::write_file(void){
  write_file(__filename);
}

void CPdb::write_file(std::string d, std::string f){
  std::string _f=d+'/'+f;
  write_file(_f);
}

//  save the PDB file
void CPdb::write_file(std::string _fn){
  std::string symbol;
  real data;
  TVector<real> v_scale(3);
  std::ofstream opdb;
  opdb.open(_fn.c_str());
  //opdb.width(17);
  if(opdb.bad())
    std::cout<<" CPDB: The file was not created"<<std::endl;
#ifdef _PDB_DEBUG_MESSAGES_
  std::cout<<" CPDB: file format: "<<__file_format<<std::endl;
#endif
  //if(__export_format==0){
    //opdb<<__total_atoms;
    //opdb<<std::endl;
    //opdb<<std::endl;
  //}
  char _buff[5];
  std::string stime = get_date();
  TVector<uint> vcont(__total_fragments);
  opdb<<"HEADER PDB structure file generated by XMolView (www.xmol.org) ["<<stime<<"]"<<std::endl;
  //opdb<<"CRYST1   15.638   15.638   42.547  90.00  90.00  90.00 P 1           1"<<std::endl;
  //opdb.precision(16);
  for(uint cell=0; cell<__total_cells; cell++){ // repetition cells
    for(uint f=0;f<__total_atoms;f++){
      symbol=v_atomic_symbols[f];
      if(strcmp(symbol.c_str(),"X") || __is_dummy){
        if(__is_numbers){
          if(__is_fragments){
            //v_fragment_table=_v;
            //__total_fragments=u;
              vcont[v_fragment_table[f]-1]+= 1;
              sprintf(_buff,"%i",vcont[v_fragment_table[f]-1]);
                }else{
                  sprintf(_buff,"%i",1+f+cell*__total_atoms);
                }
                symbol = symbol + _buff;
        }
        opdb<<"ATOM";
        opdb.width(7);
        opdb<<std::fixed<<std::right<<1+f+cell*__total_atoms<<"  ";
        opdb.width(4);
        opdb<<std::fixed<<std::left<<symbol<<"LIG     1   ";
        opdb.precision(3);
        for(uint c=0;c<3;c++){
          data=m_xyz[f+cell*__total_atoms][c];
          opdb.width(8);
          opdb<<std::fixed<<std::right<<data;
        }
        //opdb<<"  ";
        //1.00";
        if(__is_charges){
          opdb.width(9);
          opdb.precision(4);
          opdb<<std::fixed<<std::right<<v_atomic_charges[f];
        }else opdb<<"  1.00";
        opdb<<"  0.00           "<<v_atomic_symbols[f];
        // export format
        // 1: Gaussian
        //if(__export_format==1 || __file_format==2){
          //opdb.width(5);
          //opdb<<std::fixed<<std::right<<v_fragment_table[f];
        //}
        opdb<<std::endl;
      }
    }
  }
  //if(__export_format==0)
    write_copyright(opdb);
  opdb.close();
}

// SET FUNCTIONS

void CPdb::set_file_format(const uint u){
  __file_format=u;
}

void CPdb::set_export_format(const uint u){
  __export_format=u;
}

void CPdb::set_total_cells(const uint x, const uint y, const uint z, const uint t){
  __x_cells=x;
  __y_cells=y;
  __z_cells=z;
  __total_cells=t;
}

void CPdb::set_fragments(const TVector<uint>& _v, uint u){
  v_fragment_table=_v;
  __total_fragments=u;
  __is_fragments=true;
}

void CPdb::set_charges(TVector<real>& _v){
  __is_charges=true;
  v_atomic_charges=_v;
}

void CPdb::set_atomic_symbol_table(const TVector<std::string>& _v){
  v_atomic_symbol_table=_v;
  __atomic_species=v_atomic_symbol_table.size();
  v_atomic_composition_table.resize(__atomic_species);
  __is_symbol_table=true;
}

/* Fri May 18 11:39:33 MDT 2012
   deprecated, this is done in eval_composition
void CPdb::eval_atomic_number_table(void){
  uint z;
  v_atomic_number_table.resize(0);
  for( uint i=0; i<v_atomic_symbol_table.size(); i++){
    z=get_atomic_number(v_atomic_symbol_table[i]);
    v_atomic_number_table.push_back(z);
  }
}*/

void CPdb::set_cartesian(const TMatrix<real>& _m){
   m_xyz=_m;
}

void CPdb::set_direct(const TMatrix<real>& _m){
  m_uvw=_m;
}

void CPdb::set_atomic_symbols(TVector<std::string>& v){
  v_atomic_symbols =v;
}

void CPdb::set_atomic_numbers(bool b){
  __is_numbers=b;
}

//////////////////////////////
//   Coordinate functions  //
//////////////////////////////

uint CPdb::get_atomic_number(std::string s, size_t z){
  //transform(symbol.begin(), symbol.end(), symbol.begin(), toupper);
  //transform(s.begin(), s.begin()+1, s.begin(), toupper);
  for(uint j=0; j<periodic_table_atoms; j++){
    //if(strncmp(s.c_str(),symbol[j].c_str(),symbol_t[j],z)==0){
    if(strncmp(s.c_str(),symbol[j].c_str(),z)==0){
//#ifdef _PDB_DEBUG_MESSAGES_
//      std::cout<<"CPDB: atomic number "<<j+1<<std::endl;
//#endif
      return j;
    }
  }
  return 118;
}

uint CPdb::get_atomic_composition(uint i){
  return  v_atomic_composition_table[i];
}

bool CPdb::is_direct(){
  return __is_direct;
}

bool CPdb::is_periodic(){
  return __is_periodic;
}

bool CPdb::is_labels(void){
  return __is_labels;
}

bool CPdb::is_charges(void){
  return __is_charges;
}

uint CPdb::get_total_species(void){
  return __atomic_species;
}

uint CPdb::get_total_atoms(void){
  return m_xyz.rows();
}

uint CPdb::get_format(void){
  return __file_format;
}

real CPdb::get_cartesian(uint x, uint y){
  return  m_xyz[x][y];
}

real CPdb::get_direct(uint x, uint y){
  return  m_uvw[x][y];
}

TVector<std::string> CPdb::get_atomic_symbols(void){
  return  v_atomic_symbols;
}

TVector<std::string> CPdb::get_atomic_labels(void){
  return  v_atomic_labels;
}

TVector<std::string> CPdb::get_atomic_symbol_table(void){
  return  v_atomic_symbol_table;
}

TVector<uint> CPdb::get_atom_table(void){
  return  v_atom_table;
}

TVector<uint> CPdb::get_atomic_composition_table(void){
  return  v_atomic_composition_table;
}

TVector<uint> CPdb::get_atomic_number_table(void){
  return  v_atomic_number_table;
}

TVector<uint> CPdb::get_atomic_numbers(void){
  return  v_atomic_numbers;
}

TVector<uint> CPdb::get_fragment_table(void){
  return  v_fragment_table;
}

TVector<real> CPdb::get_direct(uint x){
  return  m_uvw[x];
}

TVector<real> CPdb::get_cartesian(uint x){
  return  m_xyz[x];
}

TVector<real> CPdb::get_direct_basis(void){
  return v_xyz;
}

TVector<real> CPdb::get_charges(void){
  return v_atomic_charges;
}

TMatrix<real> CPdb::get_xyz_input(void){
  return  m_xyz_input;
}

TMatrix<real> CPdb::get_cartesian(void){
  return  m_xyz;
}

TMatrix<real> CPdb::get_direct(void){
  return  m_uvw;
}

// Extra functions
void CPdb::write_copyright(std::ofstream& _o){
  //time_t rawtime;
  //char * _st, tmp_buffer[256];
  //std::string stime;
  //time (&rawtime);
  //sprintf(tmp_buffer,"%s",ctime(&rawtime));
  //_st = strtok(tmp_buffer,"\n");
  //stime = _st;
  std::string stime = get_date();
  _o<<std::endl<<std::endl<<std::endl;
  _o<<"# ------------------------------------------------------------#"<<std::endl;
  _o<<"# PDB structure file generated by XMolView (www.xmol.org)     #"<<std::endl;
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

void CPdb::center_coordinates(void){
  TVector<real> min_xyz(3);
  TVector<real> max_xyz(3);
  min_xyz[0]=m_xyz_input.col_min(0);
  min_xyz[1]=m_xyz_input.col_min(1);
  min_xyz[2]=m_xyz_input.col_min(2);
#ifdef _PDB_DEBUG_DATA_
  std::cout<<" CPDB: initial m_xyz = "<<m_xyz;
  std::cout<<" CPDB: min_xyz = "<<min_xyz;
#endif
#ifdef _PDB_DEBUG_DATA_
  std::cout<<" CPDB: positive m_xyz = "<<m_xyz;
#endif
  for(uint i=0;i<__total_atoms;i++){
    m_xyz[i] = m_xyz_input[i] - min_xyz;
  }
  if(__is_periodic){
    max_xyz[0]=m_xyz.col_max(0);
    max_xyz[1]=m_xyz.col_max(1);
    max_xyz[2]=m_xyz.col_max(2);
#ifdef _PDB_DEBUG_DATA_
    std::cout<<" CPDB: is_periodic"<<std::endl;
    std::cout<<" CPDB: max_xyz = "<<max_xyz;
#endif
    for(uint i=0;i<__total_atoms;i++){
      //m_xyz[i]-=min_xyz;
      for(uint j=0; j<3; j++){
        m_xyz[i][j]+=(0.5*(scl_uvwTxyz[j][j]-max_xyz[j]));
        //m_xyz[i]+=0.25*max_xyz;
        //m_xyz[i]+=0.25*scl_uvwTxyz[j];
      }
    }

  }else{
    real _vacuum = 2.2*__vacuum_space;
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

    unit_uvwTxyz[0]=uvwTxyz[0]/uvwTxyz[0].magnitude();
    unit_uvwTxyz[1]=uvwTxyz[1]/uvwTxyz[1].magnitude();
    unit_uvwTxyz[2]=uvwTxyz[2]/uvwTxyz[2].magnitude();
  }
#ifdef _PDB_DEBUG_DATA_
    std::cout<<" CPDB: boxed m_xyz = "<<m_xyz;
    std::cout<<" CPDB: min_xyz = "<<min_xyz;
    std::cout<<" CPDB: max_xyz = "<<max_xyz;
#endif
  //std::cout<<" xyz = "<<m_xyz;
  //uvwTxyz[0][0]=1.1*max_xyz[0];
  //uvwTxyz[1][1]=1.1*max_xyz[1];
  //uvwTxyz[2][2]=1.1*max_xyz[2];
}

// END

