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

#include<fl_xmol_view.h>

int main(int argc, char* argv[]){
#if defined (BUILD_FOR_LINUX) || defined (BUILD_FOR_MACOS)
  std::string dir_file;
  char full_dir_file[256];
#endif
  fl_xmol_view * mv = new fl_xmol_view();
#if defined (BUILD_FOR_LINUX) || defined (BUILD_FOR_MACOS)
  //if ( argc > 1){
    //std::cout<<" argument :"<<argv[1]<<std::endl;
    //dir_file = argv[1];
    dir_file = "";
    fl_filename_absolute(full_dir_file, sizeof(full_dir_file), dir_file.c_str());
    std::cout<<" directory :"<<full_dir_file<<std::endl;
    mv->open_file(full_dir_file);
  //}
#endif
  mv->show();
  //mv->show_view(full_dir_file);
  Fl::run();
  return 0;
}

