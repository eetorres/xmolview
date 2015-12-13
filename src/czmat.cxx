//========================================================================
// FILE: CZmat.cxx
//
// This an utility program to manipulate and generate structure files
//
// Copyright 2012-2015 by Edmanuel Torres
// email:   eetorres@gmail.com
//
//========================================================================
//  This file is part of xmolview                                       //
//                                                                      //
//  Foobar is free software: you can redistribute it and/or modify      //
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
#include<czmat.h>
#include<ctype.h>
#include<algorithm>

//#define _ZMAT_INFO_MESSAGES_
//#define _ZMAT_DEBUG_MESSAGES_
//#define _ZMAT_SHOW_DATA_

// formats
// 1 zmat no fragments
// 2 zmat including fragments

CZmat::CZmat(){
  clear();
}

void CZmat::clear(void){
  //__xyz_filename = "xyz";
  __system_title = "System Title";
  __export_format=0;
  __file_format=1; // standard format; not fragments
  __header_lines=0;
  //__is_selected_dynamics=false;
  __is_direct=false;
  __is_periodic=false;
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
  //scl_uvwTxyz.resize(3,3);
}

bool CZmat::get_zmat_format(void){
  int res1, res2, var1, var2, var3;
  int l, h, j, k;
  uint i;
  uint config;
  float a, b, c;//, f;
  char str[20], sa[18], sb[18], sc[18];
  std::string text_line;
  std::ifstream zmat;
  try{
    zmat.open(__filename.c_str());
    if(!zmat.is_open()){
#ifdef _ZMAT_DEBUG_MESSAGES_
      std::cout<<"CZmat: The zmat file: "<<__filename<<std::endl;
      std::cout<<"CZmat: was not found\n"<<std::endl;
#endif
      return false;
    }//else{
    res1=0;
    __header_lines=0;
    v_variable_table.resize(0);
    v_variable_config.resize(0);
    while(res1!=2 && !zmat.eof()){
      //zmat.ignore(256,'\n');
      std::getline(zmat,text_line);
      res1 = std::sscanf((const char*)text_line.c_str(),"%i %i %*s",&k,&l);
      __header_lines++;
#ifdef _ZMAT_DEBUG_MESSAGES_
     std::cout<<"CZmat:  res1: "<<res1<<std::endl;
#endif
    }
#ifdef _ZMAT_DEBUG_MESSAGES_
    std::cout<<"CZmat:  res1: "<<res1<<std::endl;
#endif

#ifdef _ZMAT_DEBUG_MESSAGES_
    std::cout<<"CZmat: file format: "<<__file_format<<std::endl;
#endif
    // atoms need to be counted
    //zmat.open(__filename.c_str());
    __total_atoms=0;
    while((res1 >= 1 && res1 <= 8) && !zmat.eof()){
      std::getline(zmat,text_line);
#ifdef _ZMAT_DEBUG_MESSAGES_
      std::cout<<"CZmat: ----- new atom -----"<<std::endl;
      std::cout<<"CZmat: text "<<text_line<<std::endl;
#endif
      res1 = std::sscanf((const char*)text_line.c_str(),"%s %i %f %i %f %i %f %i",str,&h,&a,&i,&b,&j,&c,&k);
      config = 0;
#ifdef _ZMAT_DEBUG_MESSAGES_
      std::cout<<"CZmat: res1="<<res1<<std::endl;
#endif
      //
      res2 = std::sscanf((const char*)text_line.c_str(),"%s %*s %f %*s",str,&a);
      var1 = std::sscanf((const char*)text_line.c_str(),"%s %*s %s %*s",str,sa);
      if((var1!=res2) && (var1==2)){
        //std::cout<<" --- sa:"<<sa<<std::endl;
        add_variable(sa);
        config += 1;
      }
      res2 = std::sscanf((const char*)text_line.c_str(),"%s %*s %*s %*s %f %*s",str,&b);
      var2 = std::sscanf((const char*)text_line.c_str(),"%s %*s %*s %*s %s %*s",str,sb);
      if((var2!=res2) && (var2==2)){
        //std::cout<<" --- sb:"<<sb<<std::endl;
        add_variable(sb);
        config += 2;
      }
      res2 = std::sscanf((const char*)text_line.c_str(),"%s %*s %*s %*s %*s %*s %f %*s",str,&c);
      var3 = std::sscanf((const char*)text_line.c_str(),"%s %*s %*s %*s %*s %*s %s %*s",str,sc);
      if((var3!=res2) && (var3==2)){
        //std::cout<<" --- sc:"<<sc<<std::endl;
        add_variable(sc);
        config += 4;
      }
      if(res1 >= 1 && res1 <= 8 ){
        __total_atoms++;
        v_variable_config.push_back(config);
      }
#ifdef _ZMAT_DEBUG_MESSAGES_
      std::cout<<"CZmat: total atoms "<<__total_atoms<<std::endl;
      std::cout<<"CZmat: var1 "<<var1<<" var2 "<<var2<<" var3 "<<var3<<std::endl;
      std::cout<<"CZmat: config = "<<config<<std::endl;
#endif
      if((config==0) && (res1!=-1)){
        if((__total_atoms==1) && (res1!=1))
          return false;
        else if((__total_atoms==2) && (res1!=3))
          return false;
        else if((__total_atoms==3) && (res1!=5))
          return false;
        else if((__total_atoms>3) && ((res1!=7) && (res1!=8)) )
          return false;
      }
    }
    v_variable_values.resize(v_variable_table.size());

    uint nvars=0;
    while((!zmat.eof() && (nvars<v_variable_table.size())) && !zmat.eof()){
      std::getline(zmat,text_line);
      res1 = std::sscanf((const char*)text_line.c_str(),"%s %f",sa,&a);
#ifdef _ZMAT_DEBUG_MESSAGES_
      std::cout<<"CZmat: ----- res1 "<<res1<<" -----"<<std::endl;
      std::cout<<"CZmat: text "<<text_line<<std::endl;
      std::cout<<"CZmat: sa="<<sa<<" a="<<a<<std::endl;
#endif
      if(res1==2){
        for(i=0; i<v_variable_table.size(); i++)
          if(sa==v_variable_table[i]){
            v_variable_values[i]=a;
            nvars++;
          }
      }
      //std::cout<<"CZmat: ---------------------------"<<std::endl;
    }
#ifdef _ZMAT_DEBUG_MESSAGES_
    std::cout<<"CZmat: v_variable_config = "<<v_variable_config<<std::endl;
    std::cout<<"CZmat: v_variable_values = "<<v_variable_values<<std::endl;
    //std::cout<<"CZmat: ---------------------------"<<std::endl;
#endif
    if(zmat.is_open()){
      zmat.close();
    }
#ifdef _ZMAT_DEBUG_MESSAGES_
    std::cout<<"CZmat: number of atoms found "<<__total_atoms<<std::endl;
    std::cout<<"CZmat: The zmat file was successful closed"<<std::endl;
#endif
    if((__total_atoms>0) && (v_variable_config.size()==__total_atoms))
      return true;
    else
      return false;
  }catch(...){
    return false;
  }
}

bool CZmat::add_variable(std::string s){
  bool found = false;
  std::string::iterator it;
  if(s[0]=='-'){
    it=s.begin();
    s.erase(it);
  }
  for(uint i=0; i<v_variable_table.size(); i++){
    if(v_variable_table[i]==s)
      found = true;
  }
  if(!found)
    v_variable_table.push_back(s);
  //std::cout<<" v_variable_table="<<v_variable_table;
  return found;
}

real CZmat::get_variable(std::string s){
  bool found = false;
  real value = 0;
  std::string::iterator it;
  if(s[0]=='-'){
    it=s.begin();
    s.erase(it);
    found = true;
  }
  for(uint i=0; i<v_variable_table.size(); i++){
    if(v_variable_table[i]==s)
      value = v_variable_values[i];
  }
  if(found)
    value *= -1.0;
  return value;
}


real CZmat::get_distance(uint i, uint j){
  TVector<real> v1, v2, v3;
  v1 = get_cartesian(i);
  v2 = get_cartesian(j);
  v3 = v2-v1;
  return v3.magnitude();
}

real CZmat::get_angle(uint i, uint j, uint k){
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

real CZmat::get_dihedral(uint i, uint j, uint k, uint l){
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

// load the zmat file
bool CZmat::read_file(std::string d, std::string f){
  //real data, a, b, c;
  real a, b, c;
  //uint unum;
#ifdef _ZMAT_DEBUG_MESSAGES_
  int res1;
#endif
  int h, i, j, k, config;
  TVector<real> v1(3), v2(3), v3(3), vij, vjk;
  std::string symbol, text_line;
  char str[2], sa[16], sb[16], sc[16];
  __filename = d+"/"+f;
  std::ifstream zmat;
  //
  try{
    if(!get_zmat_format()){
      return false;
    }
    //
    zmat.open(__filename.c_str());
#ifdef _ZMAT_DEBUG_MESSAGES_
    std::cout<<"CZmat: Total atoms: "<<__total_atoms<<std::endl;
#endif
    v_atomic_symbols.resize(__total_atoms);
    v_fragment_table.resize(__total_atoms);
    v_atom_table.resize(__total_atoms);
    m_xyz_input.resize(__total_atoms,3);
    m_xyz.resize(__total_atoms,3);
    m_uvw.resize(__total_atoms,3);
    //while(xyz.peek()=='\n')
      //xyz.ignore(0,'\n');
    for(uint m=0; m<__header_lines; m++){
      zmat.ignore(256,'\n');
      //std::getline(zmat,text_line);
      //res1 = std::sscanf((const char*)text_line.c_str(),"%i %i %*s",str,&k,&l);
      //__header_lines++;
    }
    //zmat.ignore(0,'\n');
    //zmat.ignore(0,'\n');
#ifdef _ZMAT_INFO_MESSAGES_
    std::cout<<"CZmat: loading atomic coordinates"<<std::endl;
#endif
    for(uint f=0;f<__total_atoms;f++){
      std::getline(zmat,text_line);
      config=v_variable_config[f];
#ifdef _ZMAT_DEBUG_MESSAGES_
      std::cout<<" text line: "<<text_line<<" config="<<config<<std::endl;
#endif
      if(config==0){
#ifdef _ZMAT_DEBUG_MESSAGES_
        res1 = std::sscanf((const char*)text_line.c_str(),"%s %i %lf %i %lf %i %lf %i",str,&i,&a,&j,&b,&k,&c,&h);
#else
        std::sscanf((const char*)text_line.c_str(),"%s %i %lf %i %lf %i %lf %i",str,&i,&a,&j,&b,&k,&c,&h);
#endif
      }else if(config==1){
#ifdef _ZMAT_DEBUG_MESSAGES_
        res1 = std::sscanf((const char*)text_line.c_str(),"%s %i %s %i %lf %i %lf %i",str,&i,sa,&j,&b,&k,&c,&h);
#else
        std::sscanf((const char*)text_line.c_str(),"%s %i %s %i %lf %i %lf %i",str,&i,sa,&j,&b,&k,&c,&h);
#endif
        a=get_variable(sa);
      }else if(config==2){
#ifdef _ZMAT_DEBUG_MESSAGES_
        res1 = std::sscanf((const char*)text_line.c_str(),"%s %i %lf %i %s %i %lf %i",str,&i,&a,&j,sb,&k,&c,&h);
#else
        std::sscanf((const char*)text_line.c_str(),"%s %i %lf %i %s %i %lf %i",str,&i,&a,&j,sb,&k,&c,&h);
#endif
        b=get_variable(sb);
      }else if(config==4){
#ifdef _ZMAT_DEBUG_MESSAGES_
        res1 = std::sscanf((const char*)text_line.c_str(),"%s %i %lf %i %lf %i %s %i",str,&i,&a,&j,&b,&k,sc,&h);
#else
        std::sscanf((const char*)text_line.c_str(),"%s %i %lf %i %lf %i %s %i",str,&i,&a,&j,&b,&k,sc,&h);
#endif
        c=get_variable(sc);
      }else if(config==3){
#ifdef _ZMAT_DEBUG_MESSAGES_
        res1 = std::sscanf((const char*)text_line.c_str(),"%s %i %s %i %s %i %lf %i",str,&i,sa,&j,sb,&k,&c,&h);
#else
        std::sscanf((const char*)text_line.c_str(),"%s %i %s %i %s %i %lf %i",str,&i,sa,&j,sb,&k,&c,&h);
#endif
        a=get_variable(sa);
        b=get_variable(sb);
      }else if(config==5){
#ifdef _ZMAT_DEBUG_MESSAGES_
        res1 = std::sscanf((const char*)text_line.c_str(),"%s %i %s %i %lf %i %s %i",str,&i,sa,&j,&b,&k,sc,&h);
#else
        std::sscanf((const char*)text_line.c_str(),"%s %i %s %i %lf %i %s %i",str,&i,sa,&j,&b,&k,sc,&h);
#endif
        a=get_variable(sa);
        c=get_variable(sc);
      }else if(config==6){
#ifdef _ZMAT_DEBUG_MESSAGES_
        res1 = std::sscanf((const char*)text_line.c_str(),"%s %i %lf %i %s %i %s %i",str,&i,&a,&j,sb,&k,sc,&h);
#else
        std::sscanf((const char*)text_line.c_str(),"%s %i %lf %i %s %i %s %i",str,&i,&a,&j,sb,&k,sc,&h);
#endif
        b=get_variable(sb);
        c=get_variable(sc);
      }else if(config==7){
#ifdef _ZMAT_DEBUG_MESSAGES_
        res1 = std::sscanf((const char*)text_line.c_str(),"%s %i %s %i %s %i %s %i",str,&i,sa,&j,sb,&k,sc,&h);
#else
        std::sscanf((const char*)text_line.c_str(),"%s %i %s %i %s %i %s %i",str,&i,sa,&j,sb,&k,sc,&h);
#endif
        a=get_variable(sa);
        b=get_variable(sb);
        c=get_variable(sc);
      }
#ifdef _ZMAT_DEBUG_MESSAGES_
      std::cout<<"CZmat: ["<<f<<"] res1= "<<res1;
      if( res1 > 0 ){
        std::cout<<" : str= "<<str;
      }
      if( res1 > 2 ){
        std::cout<<" : i= "<<i;
        std::cout<<" : a= "<<a;
      }
      if( res1 > 4 ){
        std::cout<<" : j= "<<j;
        std::cout<<" : b= "<<b;
      }
      if( res1 > 6 ){
        std::cout<<" : k= "<<k;
        std::cout<<" : c= "<<c;
        std::cout<<" : h= "<<h;
      }
      std::cout<<std::endl;
#endif
//    zmat>>symbol;
      symbol=str;
      // make sure it uses the standar atom symbols
      transform(symbol.begin(), symbol.begin()+1, symbol.begin(), toupper);
      transform(symbol.begin()+1, symbol.end(), symbol.begin()+1, tolower);
      v_atomic_symbols[f]=symbol;
      if(f==1){
        m_xyz_input[f][0]=a;
      }else if(f==2){
        if(i!=j){
          vij = m_xyz_input[j-1]-m_xyz_input[i-1];
          real norm = vij.magnitude();
          //norm = sqrt(x1^2 + y1^2 + z1^2);
          vij = a*(vij/norm);
          //x1 /= norm; y1 /= norm; z1 /= norm;
        }else{
          vij[0]=0;
          vij[1]=0;
          vij[2]=a;
        }
        //v1 = m_xyz[f];
        m_xyz_input[f]=vrot_y(vij,b*DEG_RAD);
        m_xyz_input[f]+=m_xyz_input[i-1];
        //vij = vij + m_xyz_input[i-1];
      }else if(f>2){
        // First reference vector
        vij = m_xyz_input[j-1]-m_xyz_input[i-1];
        //x1 = x[k2] - x[k1]; y1 = y[k2] - y[k1]; z1 = z[k2] - z[k1];
        real norm = vij.magnitude();
        //norm = sqrt(x1^2 + y1^2 + z1^2);
        vij = vij/norm;
        //x1 /= norm; y1 /= norm; z1 /= norm;
        // Second reference vector
        vjk = m_xyz_input[k-1]-m_xyz_input[j-1];
        //x2 = x[k3] - x[k2]; y2 = y[k3] - y[k2]; z2 = z[k3] - z[k2];
        norm = vij*vjk;
        //norm = x1 * x2 + y1 * y2 + z1 * z2; // Project into perp plane
        vjk -= norm*vij;
        //x2 -= norm*x1; y2 -= norm*y1; z2 -= norm*z1;
        norm = vjk.magnitude();
        //norm = sqrt(x2^2 + y2^2 + z2^2); // Normalize if possible.
        if (norm > 0) { // Can skip if sin(a) == 0.
          //x2 /= norm; y2 /= norm; z2 /= norm;
          vjk = vjk/norm;
        }
        // Third reference vector
        //x3 = y1 * z2 - y2 * z1;
        //y3 = z1 * x2 - z2 * x1;
        //z3 = x1 * y2 - x2 * y1;
        v3[0] = vij[1] * vjk[2] - vjk[1] * vij[2];
        v3[1] = vij[2] * vjk[0] - vjk[2] * vij[0];
        v3[2] = vij[0] * vjk[1] - vjk[0] * vij[1];

        // Compute final position
        //x[i] = x[k1]; y[i] = y[k1]; z[i] = z[k1];
        m_xyz_input[f] = m_xyz_input[i-1];
        // x[i] += b * (cos(a) * x1 + sin(a) * (cos(t) * x2 - sin(t) * x3));
        // y[i] += b * (cos(a) * y1 + sin(a) * (cos(t) * y2 - sin(t) * y3));
        // z[i] += b * (cos(a) * z1 + sin(a) * (cos(t) * z2 - sin(t) * z3));
        b *= DEG_RAD;
        c *= DEG_RAD;
        m_xyz_input[f][0] += a * (cos(b) * vij[0] + sin(b) * (cos(c) * vjk[0] - sin(c) * v3[0]));
        m_xyz_input[f][1] += a * (cos(b) * vij[1] + sin(b) * (cos(c) * vjk[1] - sin(c) * v3[1]));
        m_xyz_input[f][2] += a * (cos(b) * vij[2] + sin(b) * (cos(c) * vjk[2] - sin(c) * v3[2]));
      }
    }
    eval_atomic_composition();
    center_coordinates();
    eval_atomic_number_table();
    //m_uvw=m_uvw-0.5;
    //////compute_general();
#ifdef _ZMAT_DEBUG_MESSAGES_
    std::cout<<"--- Summary of the zmat file ---"<<std::endl;
    std::cout<<"atom table ="<<v_atom_table;
    std::cout<<"atomic symbol table ="<<v_atomic_symbol_table;
    std::cout<<"atomic number table ="<<v_atomic_number_table;
    std::cout<<"atomic composition table ="<<v_atomic_composition_table;
    std::cout<<"atomic symbols ="<<v_atomic_symbols;
    std::cout<<"atomic numbers ="<<v_atomic_numbers;
    std::cout<<"uvwTxyz'="<<uvwTxyz;
    std::cout<<"--- Summary of the zmat file ---"<<std::endl;
#endif

#ifdef _ZMAT_SHOW_DATA_
    std::cout<<"Cartesian'="<<m_xyz;
#endif
    if(zmat.is_open()){
      zmat.close();
    }
#ifdef _ZMAT_DEBUG_MESSAGES_
    std::cout<<"CZmat: The zmat was successful closed!"<<std::endl;
#endif
    return true;
  }catch(...){
    return false;
  }
}

void CZmat::eval_atomic_composition(void){
  std::string _tstr;
  uint _an;
  v_atomic_symbol_table.resize(0);
  TVector<real> _tv;
  __vacuum_space=0;
  if(!__is_symbol_table){
    for(uint j=0; j<periodic_table_atoms; j++){
      for(uint f=0;f<__total_atoms;f++){
        if(!strcmp(v_atomic_symbols[f].c_str(),symbol[j].c_str())){
          v_atomic_symbol_table.push_back(symbol[j]);
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
      _an = get_atomic_number(v_atomic_symbols[j]);
      v_atomic_numbers[j]=_an;
      __vacuum_space=maxi(__vacuum_space,atom_rrgb[_an][0]);
      if(v_atomic_symbol_table[i]==v_atomic_symbols[j]){
        v_atom_table[j]=i;
        /*if(__is_sort && (j!=_f)){
           _tstr=v_atomic_symbols[_f];
           v_atomic_symbols[_f]=v_atomic_symbols[j];
           v_atomic_symbols[j]=_tstr;
           _tv = m_xyz[_f];
           m_xyz[_f]=m_xyz[j];
           m_xyz[j]=_tv;
        }*/
        //v_atomic_numbers[_f]=get_atomic_number(v_atomic_symbols[_f]);
        //_f++;
        _c++;
      }

    }
    v_atomic_composition_table[i]=_c;
  }
}

void CZmat::write_file(void){
  write_file(__filename);
}

void CZmat::write_file(std::string d, std::string f){
  std::string _f=d+'/'+f;
  write_file(_f);
}

//  save the zmat file
void CZmat::write_file(std::string _fn){
  std::string symbol;
  real data;
  TVector<real> v_scale(3);
  std::ofstream ozmat;
  ozmat.open(_fn.c_str());
  //ozmat.width(17);
  if(ozmat.bad())
    std::cout<<" CZmat: The file was not created"<<std::endl;
#ifdef _ZMAT_DEBUG_MESSAGES_
  std::cout<<" CZmat: file format: "<<__file_format<<std::endl;
#endif
  if(__export_format==OUTPUT_FORMAT_ATM_NFR || __export_format==OUTPUT_FORMAT_ATM_FRG){
    ozmat<<__total_atoms;
    ozmat<<std::endl;
    ozmat<<std::endl;
  }
  std::string stime = get_date();
  ozmat<<"z-matrix file generated by XMolView (www.xmol.org) ["<<stime<<"]"<<std::endl;
  ozmat<<"0 1"<<std::endl;
  ozmat.precision(6);
  //for(uint cell=0; cell<__total_cells; cell++){ // repetition cells
    for(uint f=0;f<__total_atoms;f++){
      symbol=v_atomic_symbols[f];
      if(strcmp(symbol.c_str(),"X") || __is_dummy){
        ozmat<<symbol;
        if(f>2){
		  ozmat.width(5);
		  ozmat<<std::fixed<<std::right<<"1";
		  data=get_distance(0,f);
		  ozmat.width(12);
          ozmat<<std::fixed<<std::right<<data;
          ozmat.width(5);
		  ozmat<<std::fixed<<std::right<<"2";
          data=get_angle(f,0,1);
		  ozmat.width(12);
          ozmat<<std::fixed<<std::right<<data;
          ozmat.width(5);
		  ozmat<<std::fixed<<std::right<<"3";
          data=get_dihedral(f,2,1,0);
		  ozmat.width(12);
          ozmat<<std::fixed<<std::right<<data;
		}else if(f>1){
		  ozmat.width(5);
		  ozmat<<std::fixed<<std::right<<"1";
		  data=get_distance(0,f);
		  ozmat.width(12);
          ozmat<<std::fixed<<std::right<<data;
          ozmat.width(5);
		  ozmat<<std::fixed<<std::right<<"2";
          data=get_angle(f,0,1);
		  ozmat.width(12);
          ozmat<<std::fixed<<std::right<<data;
		}else if(f>0){
		  ozmat.width(5);
		  ozmat<<std::fixed<<std::right<<"1";
		  data=get_distance(0,f);
		  ozmat.width(12);
          ozmat<<std::fixed<<std::right<<data;
          //for(uint c=0;c<3;c++){
            //data=m_xyz[f+cell*__total_atoms][c];
            //data=m_xyz[f][c];
            //ozmat.width(22);
            //ozmat<<std::fixed<<std::right<<data;
          //}
        }
        // export format
        // 1: zmatssian
        if(__export_format==OUTPUT_FORMAT_ATM_FRG || __export_format==OUTPUT_FORMAT_NAT_FRG){
          ozmat.width(5);
          ozmat<<std::fixed<<std::right<<v_fragment_table[f];
        }
        ozmat<<std::endl;
      }
    }
  //}
  if(__export_format==0)
    write_copyright(ozmat);
  ozmat.close();
}

// SET FUNCTIONS

void CZmat::set_file_format(const uint u){
  __file_format=u;
}

void CZmat::set_export_format(const uint u){
  __export_format=u;
}

void CZmat::set_total_cells(const uint x, const uint y, const uint z, const uint t){
  __x_cells=x;
  __y_cells=y;
  __z_cells=z;
  __total_cells=t;
}

void CZmat::set_fragments(const TVector<uint>& _v){
  v_fragment_table=_v;
#ifdef _ZMAT_DEBUG_MESSAGES_
  std::cout<<" CZmat: fragments: "<<v_fragment_table;
#endif
}

void CZmat::set_atomic_symbol_table(const TVector<std::string>& _v){
  v_atomic_symbol_table=_v;
  __atomic_species=v_atomic_symbol_table.size();
  v_atomic_composition_table.resize(__atomic_species);
  __is_symbol_table=true;
}

void CZmat::eval_atomic_number_table(void){
  uint z;
  v_atomic_number_table.resize(0);
  for( uint i=0; i<v_atomic_symbol_table.size(); i++){
    z=get_atomic_number(v_atomic_symbol_table[i]);
    v_atomic_number_table.push_back(z);
  }
}

void CZmat::set_cartesian(const TMatrix<real>& _m){
   m_xyz=_m;
#ifdef _ZMAT_DEBUG_MESSAGES_
  std::cout<<" CZmat: coordinates: "<<m_xyz;
#endif
}

void CZmat::set_atomic_symbols(const TVector<std::string>& v){
  v_atomic_symbols=v;
#ifdef _ZMAT_DEBUG_MESSAGES_
  std::cout<<" CZmat: symbols: "<<v_atomic_symbols;
#endif
}

//void CZmat::set_direct(const TMatrix<real>& _m){
//  m_uvw=_m;
//}

//////////////////////////////
//   Coordinate functions  //
//////////////////////////////

uint CZmat::get_atomic_number(std::string s){
  //transform(symbol.begin(), symbol.end(), symbol.begin(), toupper);
  //transform(s.begin(), s.begin()+1, s.begin(), toupper);
  for(uint j=0; j<periodic_table_atoms; j++){
    if(!strcmp(s.c_str(),symbol[j].c_str())){
//#ifdef _ZMAT_DEBUG_MESSAGES_
//      std::cout<<"CZmat: atomic number "<<j+1<<std::endl;
//#endif
      return j;
    }
  }
  return 118;
}

uint CZmat::get_atomic_composition(uint i){
  return  v_atomic_composition_table[i];
}

bool CZmat::is_direct(){
  return __is_direct;
}

bool CZmat::is_periodic(){
  return __is_periodic;
}

uint CZmat::get_total_species(void){
  return __atomic_species;
}

uint CZmat::get_total_atoms(void){
  return m_xyz.rows();
}

uint CZmat::get_format(void){
  return __file_format;
}

real CZmat::get_cartesian(uint x, uint y){
  return  m_xyz[x][y];
}

//real CZmat::get_direct(uint x, uint y){
//  return  m_uvw[x][y];
//}

TVector<std::string> CZmat::get_atomic_symbols(void){
  return  v_atomic_symbols;
}

TVector<std::string> CZmat::get_atomic_symbol_table(void){
  return  v_atomic_symbol_table;
}

TVector<uint> CZmat::get_atom_table(void){
  return  v_atom_table;
}

TVector<uint> CZmat::get_atomic_composition_table(void){
  return  v_atomic_composition_table;
}

TVector<uint> CZmat::get_atomic_number_table(void){
  return  v_atomic_number_table;
}

TVector<uint> CZmat::get_atomic_numbers(void){
  return  v_atomic_numbers;
}

TVector<uint> CZmat::get_fragment_table(void){
  return  v_fragment_table;
}

//TVector<real> CZmat::get_direct(uint x){
//  return  m_uvw[x];
//}

TVector<real> CZmat::get_cartesian(uint x){
  return  m_xyz[x];
}

TVector<real> CZmat::get_direct_basis(void){
  return v_xyz;
}

TMatrix<real> CZmat::get_xyz_input(void){
  return  m_xyz_input;
}

TMatrix<real> CZmat::get_cartesian(void){
  return  m_xyz;
}

TMatrix<real> CZmat::get_direct(void){
  return  m_uvw;
}

// Extra functions
void CZmat::write_copyright(std::ofstream& _o){
  time_t rawtime;
  char * _st, tmp_buffer[256];
  std::string stime;
  time (&rawtime);
  sprintf(tmp_buffer,"%s",ctime(&rawtime));
  _st = strtok(tmp_buffer,"\n");
  stime = _st;
  _o<<std::endl<<std::endl<<std::endl;
  _o<<"# ------------------------------------------------------------#"<<std::endl;
  _o<<"# zmatrix input file generated by XMolView (www.xmol.org)    #"<<std::endl;
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

void CZmat::center_coordinates(void){
  TVector<real> min_xyz(3);
  TVector<real> max_xyz(3);
  min_xyz[0]=m_xyz_input.col_min(0);
  min_xyz[1]=m_xyz_input.col_min(1);
  min_xyz[2]=m_xyz_input.col_min(2);
  real _vacuum=2.2*__vacuum_space;
#ifdef _debugging_data_
  std::cout<<" CXYZ: initial m_xyz = "<<m_xyz;
#endif
  for(uint i=0;i<__total_atoms;i++){
    m_xyz[i] = m_xyz_input[i] - min_xyz;
  }
#ifdef _debugging_data_
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
#ifdef _debugging_data_
  std::cout<<" CXYZ: boxed m_xyz = "<<m_xyz;
#endif
  unit_uvwTxyz[0]=uvwTxyz[0]/uvwTxyz[0].magnitude();
  unit_uvwTxyz[1]=uvwTxyz[1]/uvwTxyz[1].magnitude();
  unit_uvwTxyz[2]=uvwTxyz[2]/uvwTxyz[2].magnitude();
  //
#ifdef _debugging_data_
  std::cout<<" CXYZ: min_xyz = "<<min_xyz;
  std::cout<<" CXYZ: max_xyz = "<<max_xyz;
#endif

  //std::cout<<" xyz = "<<m_xyz; 
  //uvwTxyz[0][0]=1.1*max_xyz[0];
  //uvwTxyz[1][1]=1.1*max_xyz[1];
  //uvwTxyz[2][2]=1.1*max_xyz[2];
}

// END

