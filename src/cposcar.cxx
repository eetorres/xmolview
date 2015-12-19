//========================================================================
// FILE: cposcar.cxx
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
#include <config_debug.h>
#include <cposcar.h>

//#define _POSCAR_INFO_MESSAGES_  1
//#define _POSCAR_DEBUG_MESSAGES_    1
//#define _POSCAR_SHOW_DATA_          1

CPoscar::CPoscar(){
  clear();
}

void CPoscar::clear(void){
  __filename = "POSCAR";
  __system_title = "System Title";
  __is_atom_comments = false;
  __is_velocity_reset=false;
  __is_acceleration_reset=false;
  __is_selected_dynamics=false;
  __is_direct=false;
  __is_periodic=true;
  __is_sorted=false;
  __is_centered=false;
  __is_labels=false;
  __is_read=false;
  __fix_all_atoms=false;
  __is_potcar = false;
  __is_mit_format = false;
  //
  __lattice_constant=1.0;
  __head_lines=8;
  __atomic_species=0;
  //__potcar_species=0;
  __total_atoms=0;
  a1=0;

  m_coordinates_uvw.resize(0,3);
  m_dynamic.resize(0,3);
  m_xyz.resize(3,3);

  v_atomic_composition_table.resize(0);
  v_atomic_symbols.resize(0);
  //v_atomic_symbols.resize(0);
  v_xyz.resize(3);
  v_cells.resize(3);

  v_cells[0]=1;
  v_cells[1]=1;
  v_cells[2]=1;
  __total_cells=1;

  unit_uvwTxyz.resize(3,3);
  scl_uvwTxyz.resize(3,3);
  uvwTxyz.resize(3,3);
}

bool CPoscar::get_poscar_format(void){
  //int res2;
  int res1;
  //int k, l;
  //float a, b, c, d, e;
  float a, b, c;
  //float alpha, beta, gamma;
  char tmp_buffer[256], str1[16], str2[16], str3[16];
  std::string text_line;
  std::ifstream poscar;
  try{
    poscar.open(__poscar_file.c_str());
    if(!poscar.is_open()){
#ifdef _POSCAR_DEBUG_MESSAGES_
      setvbuf(stdout, NULL, _IONBF, 0);
      std::cout<<"CPOSCAR: The POSCAR file: "<<__poscar_file<<std::endl;
      std::cout<<"CPOSCAR: was not found\n"<<std::endl;
      std::cout<< std::flush;
#endif
      return false;
    }//else{
#ifdef _POSCAR_DEBUG_MESSAGES_
    printf("<-- The POSCAR was successful opened -->\n");
    printf("<--        Analysing the POSCAR      -->\n");
#endif
    res1=0;
    //res2=0;
    __head_lines=0;
    while((res1<3 || __head_lines<8) && !poscar.eof()){
      //poscar.ignore(256,'\n');
      std::getline(poscar,text_line,'\n');
      strcpy (poscar_head_buffer[__head_lines],text_line.c_str());
#ifdef _POSCAR_SHOW_DATA_
      std::cout<<"CPOSCAR: [header "<<__head_lines<<"] "<<text_line<<std::endl;
#endif
      res1 = std::sscanf((const char*)text_line.c_str(),"%f %f %f %s %s %s",&a,&b,&c,str1,str2,str3);
#ifdef _POSCAR_DEBUG_MESSAGES_
      std::cout <<" res1 ="<<res1<<std::endl;
#endif
      __head_lines++;
    }
    //
    // check for dynamic type
    // Read dynamic labels
    strcpy(tmp_buffer,poscar_head_buffer[__head_lines-3]);
    //if(strpbrk(tmp_buffer,"Ss")!=NULL){
#ifdef _POSCAR_DEBUG_MESSAGES_
    std::cout <<" tmp_buffer ="<<tmp_buffer<<std::endl;
#endif
    if((tmp_buffer[0]=='S') || (tmp_buffer[0]=='s')){
#ifdef _POSCAR_DEBUG_MESSAGES_
      std::cout<<" Selective Dynamics..."<<std::endl;10.5599956512451172
#endif
      __is_selected_dynamics = true;
      //poscar.getline(poscar_head_buffer[__head_lines-1],256,'\n');
      //strcpy(tmp_buffer,poscar_head_buffer[__head_lines-1]);
//#ifdef _POSCAR_DEBUG_MESSAGES_
      //printf("%s\n",poscar_head_buffer[__head_lines-1+res1]);
//#endif
      //__head_lines++;
    }else{
#ifdef _POSCAR_DEBUG_MESSAGES_
      std::cout<<" No Selective Dynamics...!!!"<<std::endl;
#endif
      __is_selected_dynamics = false;
      //strcpy(tmp_buffer,poscar_head_buffer[__head_lines-2+res1]);
    }
    //strcpy(tmp_buffer,poscar_head_buffer[__head_lines-1]);
#ifdef _POSCAR_DEBUG_MESSAGES_
    std::cout<<" POSCAR: buffer = "<<tmp_buffer<<std::endl;
#endif
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Check for coordinate type
    strcpy(tmp_buffer,poscar_head_buffer[__head_lines-2]);
    res1 = std::sscanf((const char*)tmp_buffer,"%s",str1);
#ifdef _POSCAR_DEBUG_MESSAGES_
    std::cout <<" tmp_buffer ="<<tmp_buffer<<std::endl;
    std::cout <<" res1 ="<<res1<<std::endl;
#endif
    if((tmp_buffer[0]=='D') || (tmp_buffer[0]=='d')){
      __is_direct = true;
#ifdef _POSCAR_DEBUG_MESSAGES_
      std::cout<<" Direct coordinates...!!!"<<std::endl;
#endif
    }else{
      __is_direct = false;
#ifdef _POSCAR_DEBUG_MESSAGES_
      std::cout<<" Cartesian coordinates...!!!"<<std::endl;
#endif
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////
    strcpy(tmp_buffer,poscar_head_buffer[__head_lines-1]);
    res1 = std::sscanf((const char*)tmp_buffer,"%f %f %f %s %s %s",&a,&b,&c,str1,str2,str3);
#ifdef _POSCAR_DEBUG_MESSAGES_
    std::cout <<" tmp_buffer ="<<tmp_buffer<<std::endl;
    std::cout <<" res1 ="<<res1<<std::endl;
#endif
    if(res1>=3){
      __total_atoms=1;
    }else{
      __total_atoms=0;
    }
    int last_res=3;
    while(res1>=3 && !poscar.eof()){
      //poscar.ignore(256,'\n');
      std::getline(poscar,text_line);
#ifdef _POSCAR_SHOW_DATA_
      std::cout<<"CPOSCAR: "<<text_line<<std::endl;
#endif
      last_res = res1;
      res1 = std::sscanf((const char*)text_line.c_str(),"%f %f %f %s %s %s",&a,&b,&c,str1,str2,str3);
#ifdef _POSCAR_DEBUG_MESSAGES_
      std::cout <<" res1 ="<<res1<<std::endl;
#endif
     if(res1 >= 3)
       __total_atoms++;
    }
    if(last_res>3 && __is_selected_dynamics)
      __is_atom_comments = false;
    else
      __is_atom_comments = true;
#ifdef _POSCAR_DEBUG_MESSAGES_
    std::cout<<"CPOSCAR: total header lines "<<__head_lines<<std::endl;
#endif
    //
    if(poscar.is_open())
      poscar.close();
#ifdef _POSCAR_DEBUG_MESSAGES_
    std::cout<<"CPOSCAR: number of atoms found "<<__total_atoms<<std::endl;
#endif
#ifdef _POSCAR_DEBUG_MESSAGES_
    std::cout<<"CPOSCAR: !!!The POSCAR file was successful closed!!!"<<std::endl;
#endif
    if(__total_atoms>0)
      return true;
    else
      return false;
  }catch(...){
    return false;
  }
  return false;
}

// load the POSCAR file
bool CPoscar::read_file(std::string d, std::string f){
  char tmp_buffer[256],*pch, dyn;
  std::string str;
  int res1=0, chem_pos;
  //string atom_symbol;
  uint ui;
  real data;
  TVector<real> v1, v2, v3;
  // Windows work only with double
  double x=0, y=0, z=0;
  //long double x=0, y=0, z=0;
  __dir = d;
  __filename = f;
  __poscar_file = __dir+"/"+__filename;
  std::ifstream poscar;
  try{
    if(!get_poscar_format()){
      return false;
    }
    poscar.open(__poscar_file.c_str());
    if(!poscar.is_open()){
#ifdef _POSCAR_DEBUG_MESSAGES_
      std::cout<<"The POSCAR: "<<__poscar_file<<std::endl;
      std::cout<<"was not found\n"<<std::endl;
#endif
      return false;
    }else{
#ifdef _POSCAR_DEBUG_MESSAGES_
      printf("<-- The POSCAR to read was successful opened -->\n");
#endif
    }
#ifdef _POSCAR_DEBUG_MESSAGES_
    std::cout<<" [Head file start] "<<std::endl;
#endif
    // read the header
    for(uint _i=0;_i<__head_lines-1;_i++){
      poscar.getline(poscar_head_buffer[_i],256,'\n');
#ifdef _POSCAR_DEBUG_MESSAGES_
      printf(" header: %s\n",poscar_head_buffer[_i]);
#endif
    }
#ifdef _POSCAR_DEBUG_MESSAGES_
    std::cout<<" [Head file end] "<<std::endl;
#endif
    __system_title = poscar_head_buffer[0];
    // Read structure constant
    //sscanf(poscar_head_buffer[1],"%lf",&__lattice_constant);
    char * cp = poscar_head_buffer[1];
    __lattice_constant = (real)strtold(cp,NULL);
#ifdef _POSCAR_DEBUG_MESSAGES_
    std::cout.precision(16);
    std::cout<<std::fixed<<"Structure scaling factor "<<__lattice_constant<<std::endl;
    std::cout<<std::fixed<<"Read lattice vectors"<<std::endl;
#endif
    // Read lattice vectors
    for(uint i=2; i<5; i++){
      std::cout<<std::fixed<<"Line ("<<i<<"]:"<<poscar_head_buffer[i]<<std::endl;
      // Windows use floats
      sscanf(poscar_head_buffer[i],"%lf %lf %lf",&x,&y,&z);
      ///sscanf(poscar_head_buffer[i],"%Lf %Lf %Lf",&x,&y,&z);
      //printf("V%i %Lf %Lf %Lf\n",i,x,y,z);
      scl_uvwTxyz[i-2][0] = x;
      scl_uvwTxyz[i-2][1] = y;
      scl_uvwTxyz[i-2][2] = z;
      //uvwTxyz[i-2]=m_uvw[i-2];
      v_xyz[i-2] = __lattice_constant*scl_uvwTxyz[i-2].magnitude();
      unit_uvwTxyz[i-2]=scl_uvwTxyz[i-2]/scl_uvwTxyz[i-2].magnitude();
    }
#ifdef _POSCAR_DEBUG_MESSAGES_
    //cout<<"uvwTxyz="<<uvwTxyz;
    std::cout<<"v_xyz="<<v_xyz;
    std::cout<<"Basis vectors = "<<scl_uvwTxyz;
#endif
    uvwTxyz = __lattice_constant*scl_uvwTxyz;
    uvwTxyz_t = uvwTxyz;
    uvwTxyz_t.transpose();
    //m_xyz = __lattice_constant*scl_uvwTxyz;
#ifdef _POSCAR_DEBUG_MESSAGES_
    std::cout<<"uvwTxyz = "<<uvwTxyz;
    //cout<<"Basis vectors XYZ = "<<uvwTxyz;
#endif
    //
    //v3=uvwTxyz[2]/uvwTxyz[2].magnitude();
    //unit_uvwTxyz[2]=_v3;
    //
    strcpy(tmp_buffer,poscar_head_buffer[5]);
    res1 = std::sscanf((const char*)tmp_buffer,"%u",&ui);
#ifdef _POSCAR_DEBUG_MESSAGES_
    std::cout <<" tmp_buffer ="<<tmp_buffer<<std::endl;
    std::cout <<" res ="<<res1<<std::endl;
#endif
    if(res1==1) chem_pos = 5;
    else{
      __is_potcar = true;
      pch = strtok(tmp_buffer," ");
      while (pch != NULL){
#ifdef _POSCAR_DEBUG_MESSAGES_
        printf ("%s\n",pch);
#endif
        str = pch;
        str.erase(std::remove(str.begin(), str.end(), '\r'), str.end());
        //__total_atoms +=  a1;
        v_atomic_symbols.push_back(str);
        //v_atomic_composition_table
        pch = strtok (NULL," ");
        //__atomic_species++;
        //c++;
      }
#ifdef _POSCAR_DEBUG_MESSAGES_
      std::cout<<" v_atomic_symbols="<<v_atomic_symbols;
#endif
      __is_mit_format = true;
      chem_pos = 6;
    }
    //if(strpbrk(tmp_buffer,"Dd")!=NULL){
/*    if((tmp_buffer[0]=='D') || (tmp_buffer[0]=='d')){
      __is_direct = true;
#ifdef _POSCAR_DEBUG_MESSAGES_
      std::cout<<" Direct coordinates...!!!"<<std::endl;
#endif
    }else{
      __is_direct = false;
#ifdef _POSCAR_DEBUG_MESSAGES_
      std::cout<<" Cartesian coordinates...!!!"<<std::endl;
#endif
    }10.5599956512451172
*/
    // Read Atomic Composition
    strcpy(tmp_buffer,poscar_head_buffer[chem_pos]);
#ifdef _POSCAR_DEBUG_MESSAGES_
    std::cout <<" tmp_buffer ="<<tmp_buffer<<std::endl;
#endif
    pch = strtok(tmp_buffer," ");
    while (pch != NULL){
#ifdef _POSCAR_DEBUG_MESSAGES_
      printf ("%s\n",pch);
#endif
      sscanf(pch,"%i",&a1);
      //__total_atoms +=  a1;
      v_atomic_composition_table.push_back(a1);
      pch = strtok (NULL," ");
      __atomic_species++;
      //c++;
    }
#ifdef _POSCAR_DEBUG_MESSAGES_
    printf ("Atomic species: %u\n",__atomic_species);
    std::cout<<"Number of atoms per specie: "<<v_atomic_composition_table;
    std::cout<<"Total atoms: "<<__total_atoms<<std::endl;
#endif
    if((__lattice_constant==0) || (v_xyz[0]==0)  || (v_xyz[1]==0) || (v_xyz[2]==0))
      return false;
    if(__is_selected_dynamics)
      m_dynamic.resize(__total_atoms,3);
    m_xyz_input.resize(__total_atoms,3);
    m_coordinates_uvw.resize(__total_atoms,3);
    m_coordinates_xyz.resize(__total_atoms,3);
    ////mv_coordinates_xyz.resize(__total_atoms,3);

#ifdef _POSCAR_INFO_MESSAGES_
    std::cout<<"loading atomic coordinates"<<std::endl;
#endif
    for(uint f=0;f<__total_atoms;f++){
      for(uint c=0;c<3;c++){
        poscar>>data;
#ifdef _POSCAR_SHOW_DATA_
        std::cout<<" "<<data;
#endif
        // centering the coordinates arround (x,y,z) because VASP
        // POSCAR use coordinates inside the positive unit cell.
        // WARNING: This may be a user choice outside the cposcar class
        // is implementesd here temporary
        //data-=0.5;
        //if(data<0) data = (data+1.0);
        //printf("%f-",data);
        //m_coordinates_uvw[f][c]=data*v_xyz[c];
        m_xyz_input[f][c]=data;
        //if(__is_direct)
          //m_coordinates_uvw[f][c]=data;
        //else
          //m_coordinates_xyz[f][c]=data;
        //mv_coordinates_xyz[f][c]=compute_general(data,c);
      }
      if(__is_selected_dynamics){
        for(uint c=0;c<3;c++){
          poscar>>dyn;
          //std::cout<<" "<<dyn;
          m_dynamic[f][c]=dyn;
        }
        if(poscar.peek()!='\n')
          poscar.ignore(1024,'\n');
      }else{
      //if(__is_mit_format){
        //poscar>>dyn;
        poscar.ignore(1024,'\n');
      }
#ifdef _POSCAR_SHOW_DATA_
      std::cout<<std::endl;
#endif
    }
    if(__is_direct){
      m_coordinates_uvw=m_xyz_input;
      //m_xyz_input = m_coordinates_uvw*uvwTxyz_t;
      m_coordinates_xyz=m_xyz_input*uvwTxyz;
    }else{
      m_coordinates_xyz=m_xyz_input;
    }
    //m_coordinates_uvw=m_coordinates_uvw-0.5;
    //////compute_general();
#ifdef _POSCAR_SHOW_DATA_
    //std::cout<<"m_xyz_input ="<<m_xyz_input;
    //std::cout<<"Coordinates XYZ ="<<m_coordinates_xyz;
    //std::cout<<"Coordinates UVW ="<<m_coordinates_uvw;
    //std::cout<<"Dynamic ="<<m_dynamic;
#endif
    if(poscar.is_open()){
      poscar.close();
    }
#ifdef _POSCAR_DEBUG_MESSAGES_
    printf("The POSCAR was successful loaded..!\n");
    //cout<<"uvwTxyz'="<<uvwTxyz;
#endif
    //
    if((__lattice_constant>0) && (v_xyz[0]>0)  && (v_xyz[1]>0) && (v_xyz[2]>0)){
      __is_read=true;
      return true;
    }else
      return false;
  }catch(...){
    return false;
  }
}

void CPoscar::save_file(void){
  save_file_as(__poscar_file);
}

void CPoscar::save_file_as(std::string d, std::string f){
  std::string _df=d+'/'+f;
  save_file_as(_df);
  // if POTCAR does not exist then create a basic one
  std::string _potcar=d+"/POTCAR";
  std::ifstream pf(_potcar.c_str());
  if(pf.is_open()){
    pf.close();
#ifdef _POSCAR_DEBUG_MESSAGES_
    std::cout<<"The POTCAR was found"<<std::endl;
#endif
  }else{
    write_potcar(_potcar);
#ifdef _POSCAR_DEBUG_MESSAGES_
    std::cout<<"The POTCAR was NOT found"<<std::endl;
#endif
  }
}

void CPoscar::save_file_as(std::string _fn){
//  save the POSCAR file
  TVector<real> v_scale(3);
  std::ofstream oposcar;
  real val;
  //string _pos_file = _pos_pth + _pos_nam;
#ifdef _POSCAR_DEBUG_MESSAGES_
  std::cout<<" saving "<<_fn<<std::endl;
#endif
  oposcar.open(_fn.c_str());
  if(oposcar.bad())
    std::cout<<"The file was not created"<<std::endl;
#ifdef _POSCAR_DEBUG_MESSAGES_
  else
    std::cout<<"The file was created"<<std::endl;
#endif
  // set scaling factors
  for(uint i=0; i<3; i++)
    v_scale[i]=1.0/uvwTxyz[i].magnitude();
  oposcar.precision(16);
  //poscar.precision(9);
  //poscar<<"Name of the system"<<endl;
  //oposcar<<__system_title<<std::endl;
  std::string stime = get_date();
  if(__is_read){
    //for(uint _i=0;_i<__head_lines-1;_i++){
    //poscar.getline(poscar_head_buffer[_i],256,'\n');
    oposcar<<poscar_head_buffer[0]<<std::endl;
    //}
  }else{
    oposcar<<"POSCAR structure file generated by XMolView (www.xmol.org) ["<<stime<<"]"<<std::endl;
  }
#ifdef _POSCAR_DEBUG_MESSAGES_
  std::cout<<"Structure constant: "<<__lattice_constant<<std::endl;
#endif
  //set_format(true);
  oposcar<<std::fixed<<std::right<<"    "<<(real)__lattice_constant<<std::endl;
  oposcar<<" ";
  // write lattice vectors
  for(uint f=0; f<3; f++){
    for(uint c=0; c<3; c++){ oposcar<<" ";
      if(__lattice_constant!=1.0) val=scl_uvwTxyz[f][c];
      else val=uvwTxyz[f][c];
      if(val>=0) oposcar<<" ";
      oposcar.width(20);
      oposcar<<std::fixed<<std::right<<val*v_cells[f];
      //oposcar<<" ";
    }
    if(f<2) oposcar<<"\n ";
  }

  oposcar<<std::endl<<" ";
  for(uint i=0; i<v_atomic_composition_table.size(); i++){
    oposcar<<v_atomic_composition_table[i]*__total_cells<<" ";
  }
  oposcar<<std::endl;
  if(__is_selected_dynamics) oposcar<<"Selective dynamics"<<std::endl;
  if(__is_direct)
    oposcar<<"Direct coordinates"<<std::endl;
  else
    oposcar<<"Cartesian coordinates"<<std::endl;
  if(!__is_sorted)
    eval_sorting();
  if(!__is_centered)
    eval_centering();
  // write the atoms coordinates
  for(uint f=0; f<__total_atoms; f++){
    for(uint cell=0; cell<__total_cells; cell++){ // repetition cells
      oposcar<<" ";
      for(uint c=0; c<3; c++){
        if(__is_direct) val=m_coordinates_uvw[f][c];
        else val=m_coordinates_xyz[v_sorted_position[f]+cell*__total_atoms][c];
        //
        if(val>=0) oposcar<<" ";
        //oposcar<<fixed<<right<<m_coordinates_uvw[_f][_c]*v_scale[_c];
        oposcar.width(20);
        oposcar<<std::fixed<<std::right<<val;
        oposcar<<" ";
      }
      if(__is_selected_dynamics){
        for(uint c=0; c<3; c++){
          if(__fix_all_atoms) oposcar<<"F";
          else oposcar<<m_dynamic[f][c];
          oposcar<<" ";
        }
      }
      oposcar<<std::endl;
    }
  }
  //write_copyright(oposcar);
  oposcar.close();
}

void CPoscar::write_potcar(std::string _fn){
  std::ofstream opotcar;
  opotcar.open(_fn.c_str());
  if(opotcar.bad())
    std::cout<<"The POTCAR: "<<_fn<<" was not created"<<std::endl;
#ifdef _POSCAR_DEBUG_MESSAGES_
  else
    std::cout<<"The POTCAR: "<<_fn<<" was generated"<<std::endl;
#endif
  for(uint i=0; i<v_atomic_symbol_table.size(); i++)
    opotcar<<"VRHFIN="<<v_atomic_symbol_table[i]<<std::endl;
  opotcar.close();
}
// set functions

void CPoscar::set_direct(const TMatrix<real>& _m){
  //std::cout<<m_coordinates_uvw;
  m_coordinates_uvw=_m;
  //std::cout<<m_coordinates_uvw;
}

void CPoscar::set_cartesian(const TMatrix<real>& _m){
  m_coordinates_xyz=_m;
}

void CPoscar::set_velocities_zero(void){
  __is_velocity_reset=true;
}

void CPoscar::set_accelerations_zero(void){
  __is_acceleration_reset=true;
};

void CPoscar::set_total_cells(const uint x, const uint y, const uint z, const uint t){
  v_cells[0]=x;
  v_cells[1]=y;
  v_cells[2]=z;
  __total_cells=t;
}

//////////////////////////////
//   Coordinate functions  //
//////////////////////////////

TVector<uint> CPoscar::get_composition(void){
  return v_atomic_composition_table;
}

uint CPoscar::get_composition(uint i){
  return  v_atomic_composition_table[i];
}

bool CPoscar::get_format(void){
  return __is_direct;
}

bool CPoscar::is_periodic(void){
  return __is_periodic;
}

uint CPoscar::get_total_atoms(void){
  return m_coordinates_xyz.rows();
}

uint CPoscar::get_total_species(void){
  return  __atomic_species;
}

real CPoscar::get_cartesian(uint x, uint y){
  return  m_coordinates_xyz[x][y];
}

real CPoscar::get_direct(uint x, uint y){
  return  m_coordinates_uvw[x][y];
}

TVector<real> CPoscar::get_direct(uint x){
  return  m_coordinates_uvw[x];
}

TVector<real> CPoscar::get_cartesian(uint x){
  return  m_coordinates_xyz[x];
}

TVector<real> CPoscar::get_direct_basis(void){
  return v_xyz;
}

TMatrix<real> CPoscar::get_xyz_input(void){
  return  m_xyz_input;
}

TMatrix<real> CPoscar::get_cartesian(void){
  return  m_coordinates_xyz;
}

TMatrix<real> CPoscar::get_direct(void){
  return  m_coordinates_uvw;
}

void CPoscar::eval_sorting(void){
  TVector<uint> z_table = v_atomic_composition_table;
  TVector<real> vtmp;
  //bool found;
  //uint utmp;
  uint acomp = v_atomic_composition_table.size();
  uint cont=0;
  //, pos=0;
  v_sorted_position.resize(__total_atoms);
  for(uint i=0; i<acomp; i++){
    //for(uint j=0; j< v_atomic_composition_table[i]; j++){
      for(uint k=0; k< __total_atoms; k++){
        if(v_atomic_numbers[k]==v_atomic_number_table[i]){
          v_sorted_position[cont]=k;
          cont++;
        }
      }
    //}
  }
}

void CPoscar::eval_centering(void){
  //TVector<real> min_xyz(3);
  //TVector<real> max_xyz(3);
  //
  //min_xyz[0]=m_coordinates_xyz.col_min(0);
  //min_xyz[1]=m_coordinates_xyz.col_min(1);
  //min_xyz[2]=m_coordinates_xyz.col_min(2);
  //
  //max_xyz[0]=m_coordinates_xyz.col_max(0);
  //max_xyz[1]=m_coordinates_xyz.col_max(1);
  //max_xyz[2]=m_coordinates_xyz.col_max(2);
//#ifdef _debugging_data_
//  std::cout<<" CXYZ: initial m_coordinates_xyz = "<<m_coordinates_xyz;
//#endif
  for(uint i=0;i<__total_atoms*__total_cells;i++){
    //m_coordinates_xyz[i]-=(min_xyz);
    //m_coordinates_xyz[i]+=(0.5*max_xyz);
    for(uint j=0; j<3; j++){
      m_coordinates_xyz[i]+=(0.5*uvwTxyz[j]);
    }
  }
  //unit_uvwTxyz[0]=uvwTxyz[0]/uvwTxyz[0].magnitude();
  //unit_uvwTxyz[1]=uvwTxyz[1]/uvwTxyz[1].magnitude();
  //unit_uvwTxyz[2]=uvwTxyz[2]/uvwTxyz[2].magnitude();
  //
//#ifdef _debugging_data_
//  std::cout<<" CXYZ: min_xyz = "<<min_xyz;
//  std::cout<<" CXYZ: max_xyz = "<<max_xyz;
//#endif
}

// Extra functions
void CPoscar::write_copyright(std::ofstream& _o){
  time_t rawtime;
  char * _st, tmp_buffer[256];
  std::string stime;
  time (&rawtime);
  sprintf(tmp_buffer,"%s",ctime(&rawtime));
  _st = strtok(tmp_buffer,"\n");
  stime = _st;
  _o<<std::endl<<std::endl<<std::endl;
  _o<<"# ------------------------------------------------------------#"<<std::endl;
  _o<<"# POSCAR input file generated by the XMolView (www.xmol.org)  #"<<std::endl;
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



// END

