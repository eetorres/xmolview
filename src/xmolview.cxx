//========================================================================
// XMOLview - www.molecular-explorer.com
//
// This program can manipulate and genarate structure files
//
//
// Features:
//   + visualization and edition
//   + Topology definition of fragments
//   + Configuration structure scanning
//
//
// Copyrigth 2002-2015 by Edmanuel Torres
//
// email:   eetorres@gmail.com
//
//========================================================================

#include<fl_xmol_view.h>


int main(int argc, char* argv[]){
  //std::string dir_file;
  //char full_dir_file[256];
  fl_xmol_view * mv = new fl_xmol_view();
  /*
  if ( argc > 1){
    //std::cout<<" argument :"<<argv[1]<<std::endl;
    dir_file = argv[1];
    fl_filename_absolute(full_dir_file, sizeof(full_dir_file), dir_file.c_str());
    mv->open_file(full_dir_file);
    std::cout<<" directory :"<<full_dir_file<<std::endl;
  }*/
  mv->show();
  Fl::run();
  return 0;
}

