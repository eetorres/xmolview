//========================================================================
// FILE - Fl_Gl_Mol_View.cxx                                            //
// For the Fast Light Tool Kit (FLTK) - www.fltk.org                    //
//========================================================================
//                                                                      //
// OpenGL 3D visualization widget to display atomic structures.         //
//                                                                      //
// Copyright 2002-2016 by Edmanuel Torres                               //
// email: eetorres@gmail.com                                            //
//                                                                      //
// Lastest update: 01/03/2016                                           //
//======================================================================//
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

#include <fl_gl_mol_view.h>
#include <config_xmv.h>

// Selection Buffer
#define SelBufferSize 512

GLuint selectBuf[BUFSIZE];
GLint hits;
int render_mode = MODE_RENDER;
//int menu_mode = MODE_RENDER;
int cursorX, cursorY;

GLint viewport[4];
GLdouble modelview[16];
GLdouble projection[16];

//const int font_size_symbol=12;
//const int font_size_panel_label=12;
//const int font_size_slider_label=8;
//const int font_size_pie_label=12;

#define ALPHA 0.5
#define HAVE_GL 1

#if HAVE_GL
Fl_Gl_Mol_View::Fl_Gl_Mol_View(int x,int y,int w,int h,const char *l) : Fl_Gl_Atom(x,y,w,h,l)
#else
Fl_Gl_Mol_View::Fl_Gl_Mol_View(int x,int y,int w,int h,const char *l) : Fl_Box(x,y,w,h,l)
#endif // HAVE_GL
{
  clear();
  base_view = 10.0;
  shift_factor = 0.05;
  // initialize fonts
  font_size_symbol=12;
  font_size_pie_label=12;
  font_size_panel_label=12;
  font_size_slider_label=8;
  // View direction
  __x_ang = 0.0;
  __y_ang = 0.0;
  __z_ang = 0.0;
  // View position and scale
  size    = 10.0;
  x_shift = 0.0;
  y_shift = 0.0;
  z_shift = 0.0;
  y_off   = 1.0;
  zoom    = 1.0;
  // Foreground colors
  fgred   = 1.0; //0.5;
  fggreen = 1.0;
  fgblue  = 1.0; //0.5;
  // Background colors
  bgred   = 0.0;
  bggreen = 0.0;
  bgblue  = 0.1;
  // TrackBall
  tb_button   = -1;
  tb_tracking = GL_FALSE;
  tb_animate  = GL_TRUE;
  tb_angle    = 0.0;
  // Atom Appareance
  f_atom_brightness          = 0.8;
  f_atom_radius_scale        = 0.5;
  f_bond_radius_scale        = 1.0;
  f_atom_bond_delta          = 0.1;
  f_atom_bond_inv_delta      = 10.0;
  // Bond Appareance
  f_bond_brightness          = 0.8;
  f_background_brightness    = 0.0;
  f_highlight_brightness     = 1.0;
  f_select_brightness        = 1.0;
  f_highlight_brightness_max = 1.5;
  f_select_brightness_max    = 1.5;
  f_atom_brightness_max      = 1.5;
  //
  //eval_sphere(2);
  // for the inial world axes
  //eval_cylinder(20);
  // Locatios
  view_left   =  base_view;
  view_right  =  base_view;
  view_bottom =  base_view;
  view_top    =  base_view;
  view_near   = -1000;
  view_far    =  1000;
  view_axis_x = -base_view+8.0;
  view_axis_y = -base_view+2.0;
  // Mouse
  is_left_click=false;
  is_right_click=false;
  // 3D View area
  is_draw_world_axes_     = true;
  is_draw_labels_         = false;
  is_draw_symbols_        = false;
  is_draw_numbers_        = false;
  is_draw_tools_          = false;
  is_draw_molecular_axes_ = false;
  is_draw_molecular_axis_ = false;
  // GUI
  is_draw_processing      = false;
  is_draw_menu            = false;
  is_lock_controls        = false;
  // Properties
  is_eval_bonds           = true;
  is_update_bonds         = false;
  update_bonds_color      = false;
  is_draw_bonds_          = false;
  is_initialize_rot       = true;
  // Behaviour
  is_control_left_on      = false;
  last_action = 0;
  // Arrays
  v_bbox.resize(3);
  v_distance1.resize(3);
  v_distance2.resize(3);
  v_distance3.resize(3);
  for(uint i=0; i<NUMBER_OF_RADIOS; i++){
    is_radio_active[i]    = false;
  }
  is_draw_bbox_           = true;
  //init_font();
  //
#if !HAVE_GL
  label("OpenGL is required for this demo to operate.");
  align(FL_ALIGN_WRAP | FL_ALIGN_INSIDE);
#endif /* !HAVE_GL */
}

////////////////////////////OPENGL////////////////////////////////
#if HAVE_GL
Fl_Gl_Mol_View::~Fl_Gl_Mol_View(){
  //m_sphere_e1.clear();
  //m_sphere_e2.clear();
  //m_texture1.clear();
  //m_texture2.clear();
  clear();
}

void Fl_Gl_Mol_View::clear(void){
  m_atom_position.clear();
  m_atom_coordinates.clear();
  //m_text_position.clear();
  m_radius_color.clear();
  v_atom_symbols.clear();
  //m_atoms_strip.clear();
  //
  is_eval_bonds           = true;
  update_normal_color     = true;
  update_bonds_color      = false;
  update_dark_mask        = false;
  update_highlight_fragment = false;
  update_highlight_atom   = false;
  update_coordinates      = false;
  u_active_menu           = NOT_MENU;
  //
  is_graphics_on          = false;
  is_first_structure_     = true;
  //is_draw_bbox_           = true;
  is_draw_pie_menu        = false;
  is_draw_pie_submenu     = false;
  is_draw_line            = false;
  is_draw_point           = false;
  is_draw_controls        = true;
  //
  is_atom_picked          = false;
  is_menu_picked          = false;
  is_background_picked    = false;
  is_menu_position        = false;
  is_lock_dragging        = false;
  is_menu_pie_picked      = false;
  //
  is_mode_atom            = true;
  //
  //is_draw_molecular_axes_ = false;
  //is_draw_molecular_axis_ = false;
  is_draw_atoms_          = false;
  is_update_mask_rcolor   = false;
  is_update_radius        = false;
  //
  is_update_atomic_properties = false;
  is_normal_color_        = true;
  is_highlight_atom_      = true;
  //is_normal_atom_         = false;
  is_highlight_fragment_  = false;
  is_dark_mask_           = false;
  is_pbc                  = true;
  is_linked_cell          = true;
  //
  __fragment_active       = 1;
  __highlight_atom        = 0;
  __last_highlight_atom   = 0;
  __highlight_atom_a      = 0;
  __highlight_atom_b      = 0;
  __select_begin          = -1;
  __select_end            = -1;

  i_number_of_bonds       = 0;
  i_number_of_bonds_pbc   = 0;
  __axis_precession       = 0;
  __axis_tilt             = 0;
  //r_axis_bias = 0.0;
  v_axes_position.resize(3);
  //
  u_selected_index        = 0;
  // initialize sliders
  u_slider_index          = 0;
  for(uint i=0; i<NUMBER_OF_SLIDERS; i++){
    is_slider_active[i]   = false;
  }
  // default slider
  is_slider_active[0]     = true;
  // initialize radios
  u_radio_index           = 0;
  //for(uint i=0; i<NUMBER_OF_RADIOS; i++){
    //is_radio_active[i]    = false;
  //}
  // default radio
  //is_radio_active[0]      =true;
  //is_slider_active[1]     =false;
  //is_slider_active[2]     =false;
  //is_slider_active[3]     =false;
  //
  label="";
}

bool Fl_Gl_Mol_View::initialize(void){
  bool res=false;
  res = read_files_view();
  if(res){
    // initialize_fragments();
    // initialize data to visualize
    initialize_view();
    set_view_active_fragment(0);
    // atoms are counted from 1 in the scene
    // Beging GUI functions
    set_active_fragment_index(1);
    set_update_active_fragment(); // the same as zero above.
    // set the data to visualize
    set_bounding_box(get_view_bounding_box());
    is_draw_bbox(is_view_periodic());
    set_atomic_labels(get_view_atomic_labels());
    set_atomic_symbols(get_view_atomic_symbols());
    set_atomic_numbers(get_view_atomic_numbers());
    set_atom_table(get_view_atom_table());
    set_atomic_number_table(get_view_atomic_number_table());
    set_atomic_symbol_table(get_view_atomic_symbol_table());
    set_fragment_total(get_view_total_fragments());
    set_fragment_table(get_view_fragment_table());
    initialize_atomic_coordinates(get_view_cartesian());
    // End GUI functions
#ifdef _ATOM_DATA_MESSAGES_
    std::cout<<" GLMOL: atom labels = "<<v_atom_labels;
    std::cout<<" GLMOL: atom symbols = "<<v_atom_symbols;
#endif
    // allow to display the OpenGL scene
    is_graphics(true);
    set_update_coordinates(true);
    if(is_draw_bonds_){
      update_bonds_color=true;
      is_update_mask_rcolor=true;
    }
    //update_data();
  }
  return res;
}

void Fl_Gl_Mol_View::set_active_slider(uint u){
  for(uint i=0; i<NUMBER_OF_SLIDERS; i++){
    if(i==u){
      u_slider_index=i;
      is_slider_active[i]=true;
    }else is_slider_active[i]=false;
  }
}

void Fl_Gl_Mol_View::set_active_radio(uint u, bool b){
  is_radio_active[u]=b;
}

// delete the lists if they exist
void Fl_Gl_Mol_View::delete_sphere_dl(void){
  if(v_sphere_list.size() > 0){
    for(uint i=0; i<__total_species; i++){
      glDeleteLists(v_sphere_list[i],1);
    }
  }
}

void Fl_Gl_Mol_View::create_sphere_dl(void){
#ifdef _GLMOL_DEBUG_MESSAGES_
  std::cout<<" create_sphere_dl (0)"<<std::endl;
#endif
  real radius;
  // delete the lists only if it exists
  delete_sphere_dl();
  __total_species=v_atomic_number_table_gl.size();
#ifdef _GLMOL_DEBUG_MESSAGES_
  std::cout<<" total species = "<<__total_species<<std::endl;
#endif
#ifdef _GLMOL_DEBUG_MESSAGES_
  std::cout<<" create_sphere_dl (1)"<<std::endl;
#endif
  // Create the id for each sphere list
  v_sphere_list.resize(__total_species);
  for(uint i=0; i<__total_species; i++){
    radius = atom_rrgb[v_atomic_number_table_gl[i]][0];
    v_sphere_list[i]=glGenLists(1);
    // start list
    glNewList(v_sphere_list[i],GL_COMPILE);
    // call the function that contains the rendering commands
    initialize_sphere(radius);
    // end list
    glEndList();
  }
  //return(sphere_dl);
#ifdef _GLMOL_DEBUG_MESSAGES_
  std::cout<<" create_sphere_dl (2)"<<std::endl;
#endif
}

void Fl_Gl_Mol_View::create_cylinder_dl(void){
  real length;
  // Create the id for the cylinder list
  v_cylinder_list.resize(u_bond_types);
  for(uint i=0; i<u_bond_types; i++){
#ifdef _SHOW_DEBUG_CYLINDERS_
    //std::cout<<" bond index: "<<v_bond_number[i]<<std::endl;
    std::cout<<" bond length: "<<v_bond_table[i]<<std::endl;
#endif
    length = 0.5*f_atom_bond_delta*(real)v_bond_table[i];
    v_cylinder_list[i]=glGenLists(2);
    // start list
    glNewList(v_cylinder_list[i],GL_COMPILE);
    // call the function that contains the rendering commands
    initialize_cylinder(-length);
    //initialize_cylinder(-length);
    // end list
    glEndList();
    glNewList(v_cylinder_list[i]+1,GL_COMPILE);
    // call the function that contains the rendering commands
    //initialize_cylinder(length);
    initialize_cylinder(length);
    // end list
    glEndList();
  }
}

void Fl_Gl_Mol_View::eval_mask_rcolor(void){
  gm_rgb color;
  //TVector<real> vcolor;
#ifdef _GLMOL_DEBUG_BOND_COLORS_
  std::cout<<" EVAL MASK COLOR: eval_mask_rcolor (0)"<<std::endl;
#endif
  //uint _fhl;
  if(update_normal_color){
    for(int i=0; i<__number_of_atoms; i++){
      //if(is_highlight_fragment_ && v_fragment_table_gl[i]==__fragment_active){
        //color = palette.get_color(__fragment_active);
        //m_atom_rcolor[i][1] = f_atom_brightness*color.r;
        //m_atom_rcolor[i][2] = f_atom_brightness*color.g;
        //m_atom_rcolor[i][3] = f_atom_brightness*color.b;
      //}else
      if(is_mode_atom){
        //m_atom_rcolor[i] = f_atom_brightness*m_radius_color[i];
        m_atom_rcolor[i][1] = f_atom_brightness*m_radius_color[i][1];
        m_atom_rcolor[i][2] = f_atom_brightness*m_radius_color[i][2];
        m_atom_rcolor[i][3] = f_atom_brightness*m_radius_color[i][3];
      }else{
        //std::cout<<" fragment ="<<i<<std::endl;
        color = palette.get_color(v_fragment_table_gl[i]);
        //vcolor = palette.get_vcolor(v_fragment_table_gl[i]);
        //std::cout<<" vcolor ="<<vcolor<<std::endl;
        //m_atom_rcolor[i] = (f_atom_brightness*vcolor);
        m_atom_rcolor[i][1] = f_atom_brightness*color.r;
        m_atom_rcolor[i][2] = f_atom_brightness*color.g;
        m_atom_rcolor[i][3] = f_atom_brightness*color.b;
      }
    }
    update_normal_color=false;
  }
#ifdef _GLMOL_DEBUG_BOND_COLORS_
  std::cout<<" EVAL MASK COLOR: eval_mask_rcolor (1)"<<std::endl;
#endif
  update_bonds_color=true;
  if(update_bonds_color){
    for(uint i=0; i<i_number_of_bonds; i++){
      int idx0 = m_bond_indices[i][0];
      int idx1 = m_bond_indices[i][1];
      if(is_mode_atom){
        //m_bond_rcolor_0[i] = f_atom_brightness*m_radius_color[idx0];
        m_bond_rcolor_0[i][1] = f_atom_brightness*m_radius_color[idx0][1];
        m_bond_rcolor_0[i][2] = f_atom_brightness*m_radius_color[idx0][2];
        m_bond_rcolor_0[i][3] = f_atom_brightness*m_radius_color[idx0][3];
        //m_bond_rcolor_1[i] = f_atom_brightness*m_radius_color[idx1];
        m_bond_rcolor_1[i][1] = f_atom_brightness*m_radius_color[idx1][1];
        m_bond_rcolor_1[i][2] = f_atom_brightness*m_radius_color[idx1][2];
        m_bond_rcolor_1[i][3] = f_atom_brightness*m_radius_color[idx1][3];
      }else{
        //std::cout<<" bond("<<i<<")=["<<idx0<<","<<idx1<<"]="<<v_fragment_table_gl[idx0]<<std::endl;
        color = palette.get_color(v_fragment_table_gl[idx0]);
        //vcolor = palette.get_vcolor(v_fragment_table_gl[idx0]);
        //m_bond_rcolor_0[i] = f_atom_brightness*vcolor;
        m_bond_rcolor_0[i][1] = f_atom_brightness*color.r;
        m_bond_rcolor_0[i][2] = f_atom_brightness*color.g;
        m_bond_rcolor_0[i][3] = f_atom_brightness*color.b;
        //m_bond_rcolor_1[i] = f_atom_brightness*vcolor;
        m_bond_rcolor_1[i][1] = f_atom_brightness*color.r;
        m_bond_rcolor_1[i][2] = f_atom_brightness*color.g;
        m_bond_rcolor_1[i][3] = f_atom_brightness*color.b;
      }
    }
#ifdef _GLMOL_DEBUG_BOND_COLORS_
  std::cout<<" EVAL MASK COLOR: eval_mask_rcolor (2)"<<std::endl;
  std::cout<<" PBC bonds: "<<i_number_of_bonds_pbc<<std::endl;
#endif
    for(uint i=0; i<i_number_of_bonds_pbc; i++){
#ifdef _ATOM_DEBUG_BOND_COLORS_
        std::cout<<" bond number = "<<i<<std::endl;
#endif
      int idx0 = m_bond_indices_pbc[i][0];
      int idx1 = m_bond_indices_pbc[i][1];
#ifdef _ATOM_DEBUG_BOND_COLORS_
        std::cout<<" (idx0,idx1) = ("<<idx0<<","<<idx1<<")"<<std::endl;
#endif
      if(is_mode_atom){
        //m_bond_rcolor_pbc_0[i] = f_atom_brightness*m_radius_color[idx0];
        m_bond_rcolor_pbc_0[i][1] = f_atom_brightness*m_radius_color[idx0][1];
        m_bond_rcolor_pbc_0[i][2] = f_atom_brightness*m_radius_color[idx0][2];
        m_bond_rcolor_pbc_0[i][3] = f_atom_brightness*m_radius_color[idx0][3];
        //m_bond_rcolor_pbc_1[i] = f_atom_brightness*m_radius_color[idx1];
        m_bond_rcolor_pbc_1[i][1] = f_atom_brightness*m_radius_color[idx1][1];
        m_bond_rcolor_pbc_1[i][2] = f_atom_brightness*m_radius_color[idx1][2];
        m_bond_rcolor_pbc_1[i][3] = f_atom_brightness*m_radius_color[idx1][3];
      }else{
#ifdef _ATOM_DEBUG_BOND_COLORS_
        std::cout<<" bond("<<i<<")=["<<idx0<<","<<idx1<<"]="<<v_fragment_table_gl[idx0]<<std::endl;
#endif
        color = palette.get_color(v_fragment_table_gl[idx0]);
        //vcolor = palette.get_vcolor(v_fragment_table_gl[idx0]);
        //m_bond_rcolor_pbc_0[i] = f_atom_brightness*vcolor;
        m_bond_rcolor_pbc_0[i][1] = f_atom_brightness*color.r;
        m_bond_rcolor_pbc_0[i][2] = f_atom_brightness*color.g;
        m_bond_rcolor_pbc_0[i][3] = f_atom_brightness*color.b;
        //m_bond_rcolor_pbc_1[i] = f_atom_brightness*vcolor;
        m_bond_rcolor_pbc_1[i][1] = f_atom_brightness*color.r;
        m_bond_rcolor_pbc_1[i][2] = f_atom_brightness*color.g;
        m_bond_rcolor_pbc_1[i][3] = f_atom_brightness*color.b;
      }
    }
    update_bonds_color=false;
#ifdef _ATOM_DEBUG_BOND_COLORS_
  std::cout<<" EVAL MASK COLOR: eval_mask_rcolor (3)"<<std::endl;
#endif
  }
#ifdef _ATOM_DEBUG_BOND_COLORS_
  std::cout<<" f_bond_brightness = "<<f_bond_brightness<<std::endl;
  std::cout<<" Bond colors = "<<m_bond_rcolor_0;
  std::cout<<" Fragment table = "<<v_fragment_table_gl;
  //std::cout<<" Radius color = "<<m_radius_color;
  std::cout<<" Bond index = "<<m_bond_indices;
#endif
  // testing
  /*
  if(is_highlight_fragment_ && update_highlight_fragment){
    for(int i=0; i<__number_of_atoms; i++){
      if(v_fragment_table_gl[i]==__fragment_active){
        //glColor3f(0.2*atom_rrgb[i][1],0.2*atom_rrgb[i][2],0.2*atom_rrgb[i][3]);
        //m_atom_rcolor[i][1] = f_select_brightness*m_radius_color[i][1];
        //m_atom_rcolor[i][2] = f_select_brightness*m_radius_color[i][2];
        //m_atom_rcolor[i][3] = f_select_brightness*m_radius_color[i][3];
        m_atom_rcolor[i][1] = f_atom_brightness;
        m_atom_rcolor[i][2] = 0;
        m_atom_rcolor[i][3] = 0;
      }else{
        //glColor3f(0.8*atom_rrgb[i][1],0.8*atom_rrgb[i][2],0.8*atom_rrgb[i][3]);
        m_atom_rcolor[i][1] = f_atom_brightness*m_radius_color[i][1];
        m_atom_rcolor[i][2] = f_atom_brightness*m_radius_color[i][2];
        m_atom_rcolor[i][3] = f_atom_brightness*m_radius_color[i][3];
      }
    }
    update_highlight_fragment=false;
  }*/

/*
  if(is_dark_mask_ && update_dark_mask){
    for(int i=0; i<__number_of_atoms; i++){
      if((i>=__select_begin) && (i<=__select_end)){
        //glColor3f(0.2*atom_rrgb[i][1],0.2*atom_rrgb[i][2],0.2*atom_rrgb[i][3]);
        m_atom_rcolor[i][1] = f_select_brightness*m_radius_color[i][1];
        m_atom_rcolor[i][2] = f_select_brightness*m_radius_color[i][2];
        m_atom_rcolor[i][3] = f_select_brightness*m_radius_color[i][3];
      }else{
        //glColor3f(0.8*atom_rrgb[i][1],0.8*atom_rrgb[i][2],0.8*atom_rrgb[i][3]);
        m_atom_rcolor[i][1] = f_atom_brightness*m_radius_color[i][1];
        m_atom_rcolor[i][2] = f_atom_brightness*m_radius_color[i][2];
        m_atom_rcolor[i][3] = f_atom_brightness*m_radius_color[i][3];
      }
    }
    update_dark_mask=false;
  }
*/
  if((is_highlight_atom_ || is_highlight_fragment_) && update_highlight_atom && !is_draw_tools_){
    //glColor3f(1.5*atom_rrgb[i][1],1.5*atom_rrgb[i][2],1.5*atom_rrgb[i][3]);
    //if(is_highlight_fragment_){// && v_fragment_table_gl[__last_highlight_atom]==__fragment_active){
      //m_atom_rcolor[__last_highlight_atom][1] = f_atom_brightness;
      //m_atom_rcolor[__last_highlight_atom][2] = 0;
      //m_atom_rcolor[__last_highlight_atom][3] = 0;
    //}
    if(is_mode_atom){
      //m_atom_rcolor[__last_highlight_atom] = f_atom_brightness*m_radius_color[__last_highlight_atom];
      //m_atom_rcolor[__highlight_atom] = f_select_brightness*m_radius_color[__highlight_atom];
      m_atom_rcolor[__last_highlight_atom][1] = f_atom_brightness*m_radius_color[__last_highlight_atom][1];
      m_atom_rcolor[__last_highlight_atom][2] = f_atom_brightness*m_radius_color[__last_highlight_atom][2];
      m_atom_rcolor[__last_highlight_atom][3] = f_atom_brightness*m_radius_color[__last_highlight_atom][3];
      m_atom_rcolor[__highlight_atom][1] = f_select_brightness*m_radius_color[__highlight_atom][1];
      m_atom_rcolor[__highlight_atom][2] = f_select_brightness*m_radius_color[__highlight_atom][2];
      m_atom_rcolor[__highlight_atom][3] = f_select_brightness*m_radius_color[__highlight_atom][3];
    }else{
      uint __last_fragment = v_fragment_table_gl[__last_highlight_atom];
      color = palette.get_color(__last_fragment);
      //vcolor = palette.get_vcolor(__last_fragment);
      for(int i=0; i<__number_of_atoms; i++){
        if(v_fragment_table_gl[i]==__last_fragment){
          //m_atom_rcolor[i] = f_atom_brightness*vcolor;
          m_atom_rcolor[i][1] = f_atom_brightness*color.r;
          m_atom_rcolor[i][2] = f_atom_brightness*color.g;
          m_atom_rcolor[i][3] = f_atom_brightness*color.b;
        }
      }
      uint __new_fragment  = v_fragment_table_gl[__highlight_atom];
      color = palette.get_color(__new_fragment);
      //vcolor = palette.get_vcolor(__new_fragment);
      for(int i=0; i<__number_of_atoms; i++){
        if(v_fragment_table_gl[i]==__new_fragment){
          //m_atom_rcolor[i] = f_select_brightness*vcolor;
          m_atom_rcolor[i][1] = f_select_brightness*color.r;
          m_atom_rcolor[i][2] = f_select_brightness*color.g;
          m_atom_rcolor[i][3] = f_select_brightness*color.b;
        }
      }
      //m_atom_rcolor[__highlight_atom][1] = f_select_brightness*color.r;
      //m_atom_rcolor[__highlight_atom][2] = f_select_brightness*color.g;
      //m_atom_rcolor[__highlight_atom][3] = f_select_brightness*color.b;
    }
    update_highlight_atom=false;
  }/*else if(is_highlight_atom_ && (hits == 0)){
    if(is_highlight_fragment_ && v_fragment_table_gl[__highlight_atom]==__fragment_active){
      gm_rgb color = palette.get_color(__highlight_atom);
      m_atom_rcolor[__highlight_atom][1] = f_atom_brightness*color.r;
      m_atom_rcolor[__highlight_atom][2] = f_atom_brightness*color.g;
      m_atom_rcolor[__highlight_atom][3] = f_atom_brightness*color.b;
      //m_atom_rcolor[__highlight_atom][1] = f_atom_brightness;
      //m_atom_rcolor[__highlight_atom][2] = 0;
      //m_atom_rcolor[__highlight_atom][3] = 0;
    }else{
      m_atom_rcolor[__highlight_atom][1] = f_atom_brightness*m_radius_color[__highlight_atom][1];
      m_atom_rcolor[__highlight_atom][2] = f_atom_brightness*m_radius_color[__highlight_atom][2];
      m_atom_rcolor[__highlight_atom][3] = f_atom_brightness*m_radius_color[__highlight_atom][3];
    }
  }*/else if(is_draw_tools_ && update_selected_atoms){
    for(uint i=0; i<u_selected_index; i++){
      m_atom_rcolor[v_selected_atoms[i]][1] = f_select_brightness;
      m_atom_rcolor[v_selected_atoms[i]][2] = 0;
      m_atom_rcolor[v_selected_atoms[i]][3] = 0;
    }
    if(is_unselected_atom){
      // unselect if the atom was selected
      if(is_mode_atom){
        //m_atom_rcolor[u_unselected_atom] = f_atom_brightness*m_radius_color[u_unselected_atom];
        m_atom_rcolor[u_unselected_atom][1] = f_atom_brightness*m_radius_color[u_unselected_atom][1];
        m_atom_rcolor[u_unselected_atom][2] = f_atom_brightness*m_radius_color[u_unselected_atom][2];
        m_atom_rcolor[u_unselected_atom][3] = f_atom_brightness*m_radius_color[u_unselected_atom][3];
      }else{
        color = palette.get_color(v_fragment_table_gl[u_unselected_atom]);
        //vcolor = palette.get_vcolor(v_fragment_table_gl[u_unselected_atom]);
        //m_atom_rcolor[u_unselected_atom] = f_atom_brightness*vcolor;
        m_atom_rcolor[u_unselected_atom][1] = f_atom_brightness*color.r;
        m_atom_rcolor[u_unselected_atom][2] = f_atom_brightness*color.g;
        m_atom_rcolor[u_unselected_atom][3] = f_atom_brightness*color.b;
      }
      is_unselected_atom=false;
    }
  }
  //m_atom_rcolor = m_radius_color;
}

void Fl_Gl_Mol_View::draw_atoms(void){
  TVector<real> _x(3);
  TVector<real> _xyz;
  TVector<real> e(3),p(3);
  if(is_draw_atoms_){
    for(int c=0; c<get_total_cells(); c++){ // repetition in z
      for(int i=0; i<__number_of_atoms; i++){
        if(render_mode==MODE_SELECT){
          //glPushName(i+100); // the first 100 names are reserved for the menues
          ui_rgb color;
          color = index_palette.get_index(i+MENU_RESERVED_IDS);
          //std::cout<<" color["<<i<<"]: "<<color.r<<" "<<color.g<<" "<<color.b<<std::endl;
          glColor3ub(color.r,color.g,color.b);
        }else{
          glColor3f(m_atom_rcolor[i][1],m_atom_rcolor[i][2],m_atom_rcolor[i][3]);
          // testing the selection palette ==================================
          //ui_rgb color;
          //color = index_palette.get_index(i);
          //glColor3ub(color.r,color.g,color.b);
          ///std::cout<<" color["<<i<<"]: "<<color.r<<" "<<color.g<<" "<<color.b<<std::endl;
          // testing the selection palette ==================================
        }
        _xyz=m_atom_position[i+c*__number_of_atoms];
        glPushMatrix();
        glTranslatef(_xyz[0],_xyz[1],_xyz[2]);
        //glColor3f(m_atom_rcolor[i][1],m_atom_rcolor[i][2],m_atom_rcolor[i][3]);
        glCallList(v_sphere_list[v_atom_table[i]]);
        glPopMatrix();
        // to be removed
        //if(render_mode==MODE_SELECT)
          //glPopName();
        // to be removed
      }
    }
  }
}

// draw bonds between atoms closer than the sum of their van der Waals radius
void Fl_Gl_Mol_View::draw_bonds(void){
  //uint _bn;
#ifdef _SHOW_DEBUG_BONDS_
  std::cout<<" GLMOL: draw_bonds"<<std::endl;
#endif
  TVector<real> _xyz(3), _ang(2);
  //glColor3f(f_bond_brightness*0.8, f_bond_brightness*0.8, f_bond_brightness*0.8);
  // unit cell bonds
#ifdef _SHOW_DEBUG_BONDS_
  std::cout<<" GLMOL: ----- internal bonds -----"<<std::endl;
  //std::cout<<" f_bond_brightness = "<<f_bond_brightness<<std::endl;
  //std::cout<<" Used bond colors = "<<m_bond_rcolor_0;
#endif
  for(int x=neg_x_cells; x<pos_x_cells+1; x++){ // repetition in x
    for(int y=neg_y_cells; y<pos_y_cells+1; y++){ // repetition in y
      for(int z=neg_z_cells; z<pos_z_cells+1; z++){ // repetition in z
        for(uint i=0; i<i_number_of_bonds; i++){
          _ang = m_bond_angles[i];
          _xyz = m_bond_position[i];
          _xyz=_xyz+2.0*(x*_vu+y*_vv+z*_vw);
          glPushMatrix();
          glTranslatef(_xyz[0],_xyz[1],_xyz[2]);
          glRotatef(_ang[0],0,0,1); // rotation arround the z axis
          glRotatef(_ang[1],0,1,0); // rotation arround the y asis
          //glCallList(cylinder_dl);
          //glCallList(v_cylinder_list[0]);
          glColor3f(m_bond_rcolor_0[i][1],m_bond_rcolor_0[i][2],m_bond_rcolor_0[i][3]);
          glCallList(v_cylinder_list[v_bond_number[i]]);
          //
          glColor3f(m_bond_rcolor_1[i][1],m_bond_rcolor_1[i][2],m_bond_rcolor_1[i][3]);
          glCallList(v_cylinder_list[v_bond_number[i]]+1);
          glPopMatrix();
        }
      }
    }
  }
#ifdef _SHOW_DEBUG_PERIODIC_BONDS_
  std::cout<<" m_bond_position_pbc="<<m_bond_position_pbc;
#endif
  int xu, yu, zu;
  // across boundary bonds
#ifdef _SHOW_DEBUG_BONDS_
  std::cout<<" GLMOL: ----- periodic bonds -----"<<std::endl;
#endif
  for(int x=neg_x_cells; x<pos_x_cells+1; x++){ // repetition in x
    if( pos_x_cells!=0 ){
      if( x==neg_x_cells ) xu = -1;
      else xu = int(x/pos_x_cells);
    }else xu=2;
    //std::cout<<" xu = "<<xu<<std::endl;
    for(int y=neg_y_cells; y<pos_y_cells+1; y++){ // repetition in y
      if( pos_y_cells!=0 ){
        if( y==neg_y_cells ) yu = -1;
        else yu = int(y/pos_y_cells);
      }else yu=2;
      //std::cout<<" yu = "<<yu<<std::endl;
      for(int z=neg_z_cells; z<pos_z_cells+1; z++){ // repetition in z
        if( pos_z_cells!=0 ){
          if( z==neg_z_cells ) zu = -1;
          else zu = int(z/pos_z_cells);
        }else zu=2;
        //std::cout<<" zu = "<<zu<<std::endl;
#ifdef _SHOW_DEBUG_PERIODIC_BONDS_
        std::cout<<" ("<<xu<<","<<yu<<","<<zu<<")"<<std::endl;
#endif
        if(!(pos_x_cells==0 && pos_y_cells==0 && pos_z_cells==0)){  // avoid incomplete bonds in the unit cell
          for(uint i=0; i<i_number_of_bonds_pbc; i++){
            if((m_bond_boundary_pbc[i][0]==0 || (m_bond_boundary_pbc[i][0]!=xu && xu!=2)) \
            && (m_bond_boundary_pbc[i][1]==0 || (m_bond_boundary_pbc[i][1]!=yu && yu!=2)) \
            && (m_bond_boundary_pbc[i][2]==0 || (m_bond_boundary_pbc[i][2]!=zu && zu!=2))){
#ifdef _SHOW_DEBUG_PERIODIC_BONDS_
              std::cout<<" bond = "<<m_bond_boundary_pbc[i];
              std::cout<<" ("<<x<<","<<y<<","<<z<<")"<<std::endl;
              //std::cout<<" ("<<_vu<<","<<_vv<<","<<_vw<<")"<<std::endl;
#endif
              _ang = m_bond_angles_pbc[i];
              _xyz = m_bond_position_pbc[i];
#ifdef _SHOW_DEBUG_PERIODIC_BONDS_
              std::cout<<" _xyz = "<<_xyz<<std::endl;
              std::cout<<" bond cylinder = "<<v_bond_number_pbc[i];
#endif
              _xyz=(_xyz+2.0*((real)x*_vu+(real)y*_vv+(real)z*_vw));
              glPushMatrix();
              glTranslatef(_xyz[0],_xyz[1],_xyz[2]);
              glRotatef(_ang[0],0,0,1); // rotation arround the z axis
              glRotatef(_ang[1],0,1,0); // rotation arround the y asis
              //glCallList(cylinder_dl);
              //glCallList(v_cylinder_list[0]);
              //_bn = v_bond_number_pbc[i];
              //glColor3f(f_bond_brightness*0.8, f_bond_brightness*0.8, f_bond_brightness*0.8);
              glColor3f(m_bond_rcolor_pbc_0[i][1],m_bond_rcolor_pbc_0[i][2],m_bond_rcolor_pbc_0[i][3]);
              glCallList(v_cylinder_list[v_bond_number_pbc[i]]);
              //
              glColor3f(m_bond_rcolor_pbc_1[i][1],m_bond_rcolor_pbc_1[i][2],m_bond_rcolor_pbc_1[i][3]);
              glCallList(v_cylinder_list[v_bond_number_pbc[i]]+1);
              glPopMatrix();
#ifdef _SHOW_DEBUG_PERIODIC_BONDS_
              std::cout<<" bond cylinder = "<<v_bond_number_pbc[i];
#endif
              //
              //glCallList(v_cylinder_list[v_bond_number_pbc[i]]);
              //glPopMatrix();
            }
          }
        }
      }
    }
  }
  //update_bonds_color = false;
  is_update_bonds = false;
#ifdef _SHOW_DEBUG_PERIODIC_BONDS_
  std::cout<<" GLMOL: draw_bonds end"<<std::endl;
#endif
}

//Thu Dec 29 20:49:43 MST 2011
// this code have been optimized, but still in very alpha state
// dont remove the commented code!
// use a label vector to optimize this function
void Fl_Gl_Mol_View::draw_symbols(void){
  TVector<real> _xyz;
  glColor3f(0.0F,1.0F,0.0F); // text color
  //gl_font(FL_COURIER|FL_BOLD,font_size_symbol); // text font
  //gl_font(1,GLint(f_atom_radius_scale*24)); // text font
#if defined (BUILD_FOR_MACOS)
    gl_font(FL_HELVETICA,12); // text font
#elif defined (BUILD_FOR_WINDOWS)
    //gl_font(FL_HELVETICA,14); // text font
#else
  gl_font(FL_COURIER,12); // text font
#endif
  //gl_font(FL_COURIER,font_size_symbol); // text font
  glPushMatrix();
  GLboolean boolval;
  char buff[10];
  //GLfloat str_width;
  for(int i=0; i<__number_of_atoms; i++){
  //for(int i=0; i<2; i++){
    _xyz=m_atom_coordinates[i];
    if(is_draw_symbols_ && is_draw_numbers_){
      if(is_mode_atom)
        sprintf(buff,"%s-%i",v_atom_symbols[i].c_str(),i+1);
      else
        sprintf(buff,"%i-%s",v_fragment_table_gl[i],v_atom_symbols[i].c_str());
    }else if(is_draw_labels_ && is_draw_numbers_){
      if(is_mode_atom)
        sprintf(buff,"%s-%i",v_atom_labels[i].c_str(),i+1);
      else
        sprintf(buff,"%s-%i",v_atom_labels[i].c_str(),v_fragment_table_gl[i]);
    }else if(is_draw_symbols_){
      sprintf(buff,"%s",v_atom_symbols[i].c_str());
    }else if(is_draw_labels_){
      sprintf(buff,"%s",v_atom_labels[i].c_str());
    }else{
      if(is_mode_atom)
        sprintf(buff,"%i",i+1);
      else
        sprintf(buff,"%i-%i",v_fragment_table_gl[i],i+1);
    }
    glLoadIdentity();
    glTranslatef(x_shift, y_shift, z_shift);
    // zoom in and zoom out
    glScalef(zoom,zoom,zoom);
    // set label text offset
    glTranslatef(0.5*f_atom_radius_scale*m_radius_color[i][0],0.5*f_atom_radius_scale*m_radius_color[i][0],4.1*f_atom_radius_scale*m_radius_color[i][0]);
    glMultMatrixf((GLfloat*)rot_matrix);
    //glNormal3f(0,0,1);
    // text position
    glRasterPos3f(_xyz[0],_xyz[1],_xyz[2]);
      //glRasterPos3f(0,0,0);
      // check raster position validity
      glGetBooleanv(GL_CURRENT_RASTER_POSITION_VALID, &boolval);
      if(boolval == GL_TRUE) {
          //printf("raster pos valid\n");
            gl_draw(buff, strlen(buff));
      }
  }
  glPopMatrix();
}

void Fl_Gl_Mol_View::draw_selected_numbers(void){
  TVector<real> _xyz;
  glColor3f(0.0F,0.0F,0.0F); // text color
  gl_font(FL_COURIER,12);
  //gl_font(FL_COURIER,font_size_symbol); // text font
  //gl_font(1,GLint(f_atom_radius_scale*24)); // text font
  glPushMatrix();
  char buff[10];
  for(uint i=0; i<u_selected_index; i++){
    _xyz=m_atom_position[v_selected_atoms[i]];
    sprintf(buff,"%i",i+1);
    glLoadIdentity();
    glTranslatef(x_shift, y_shift, z_shift);
    // zoom in and zoom out
    glScalef(zoom,zoom,zoom);
    // set label text offset
    glTranslatef(0,0,2.0*f_atom_radius_scale*m_radius_color[v_selected_atoms[i]][0]);
    glMultMatrixf((GLfloat*)rot_matrix);
    glNormal3f(0,0,1);
    // text position
    glRasterPos3f(_xyz[0],_xyz[1],_xyz[2]);
    // render the text
    gl_draw(buff, strlen(buff));
  }
  glPopMatrix();
}

// draw world axes
void Fl_Gl_Mol_View::draw_axes(void){
  real length = base_view/10.0;
  real width = base_view/80.0;
  //std::cout<<" world axis = "<<v_axes_position;
  //std::cout<<" v_axes_position = "<<v_axes_position;
  //std::cout<<" length = "<<length<<std::endl;
  //std::cout<<" width = "<<width<<std::endl;
  // X red axis
  glColor3f(1.0, 0.0, 0.0);
  add_axis(v_axes_position,length,width,0.0,90.0*DEG_RAD);
  // Y green axis
  glColor3f(0.0, 1.0, 0.0);
  add_axis(v_axes_position,length,width,90.0*DEG_RAD,90.0*DEG_RAD);
  // Z blue axis
  glColor3f(0.2, 0.2, 1.0);
  add_axis(v_axes_position,length,width,0.0,0.0);
}

// draw the bounding box unsing the lattice vectors
void Fl_Gl_Mol_View::draw_box(void){
  TVector<real> _v1, _v2, _v3, _v4;
  TVector<real> _v5, _v6, _v7, _v8;
  TVector<real> _p1, _p2, _p3, _p4;
  TVector<real> _p5, _p6, _p7, _p8;
  //
  _p1 =  _vu+_vv-_vw;
  _p2 =  _vu-_vv-_vw;
  _p3 = -_vu-_vv-_vw;
  _p4 = -_vu+_vv-_vw;
  //
  _p5 =  _vu+_vv+_vw;
  _p6 =  _vu-_vv+_vw;
  _p7 = -_vu-_vv+_vw;
  _p8 = -_vu+_vv+_vw;
  //
#ifdef _GLMOL_DEBUG_BBOX_
  std::cout<<"p1="<<_p1;FL_HELVETICA
  std::cout<<"p2="<<_p2;
  std::cout<<"p3="<<_p3;
  std::cout<<"p4="<<_p4;
  std::cout<<"p5="<<_p5;
  std::cout<<"p6="<<_p6;
  std::cout<<"p7="<<_p7;
  std::cout<<"p8="<<_p8;
#endif
  //
  for(int x=neg_x_cells; x<pos_x_cells+1; x++){ // repetition in x
    for(int y=neg_y_cells; y<pos_y_cells+1; y++){ // repetition in y
      for(int z=neg_z_cells; z<pos_z_cells+1; z++){ // repetition in z
        _v1 = _p1+2.0*(x*_vu+y*_vv+z*_vw);
        _v2 = _p2+2.0*(x*_vu+y*_vv+z*_vw);
        _v3 = _p3+2.0*(x*_vu+y*_vv+z*_vw);
        _v4 = _p4+2.0*(x*_vu+y*_vv+z*_vw);
        //
        _v5 = _p5+2.0*(x*_vu+y*_vv+z*_vw);
        _v6 = _p6+2.0*(x*_vu+y*_vv+z*_vw);
        _v7 = _p7+2.0*(x*_vu+y*_vv+z*_vw);
        _v8 = _p8+2.0*(x*_vu+y*_vv+z*_vw);
        //
#ifdef _GLMOL_DEBUG_BBOX_
        std::cout<<"v1="<<_v1;
        std::cout<<"v2="<<_v2;
        std::cout<<"v3="<<_v3;
        std::cout<<"v4="<<_v4;
        std::cout<<"v5="<<_v5;
        std::cout<<"v6="<<_v6;
        std::cout<<"v7="<<_v7;
        std::cout<<"v8="<<_v8;
#endif
        //
        glColor3f(1.0,0.0,0.0);
        glBegin(GL_LINES);// X lines
        glVertex3f(_v1[0],_v1[1],_v1[2]); glVertex3f(_v4[0],_v4[1],_v4[2]);
        glVertex3f(_v2[0],_v2[1],_v2[2]); glVertex3f(_v3[0],_v3[1],_v3[2]);
        glVertex3f(_v7[0],_v7[1],_v7[2]); glVertex3f(_v6[0],_v6[1],_v6[2]);
        glVertex3f(_v8[0],_v8[1],_v8[2]); glVertex3f(_v5[0],_v5[1],_v5[2]);
        glEnd();
        //
        glColor3f(0.0,1.0,0.0);
        glBegin(GL_LINES);// Y lines
        glVertex3f(_v1[0],_v1[1],_v1[2]); glVertex3f(_v2[0],_v2[1],_v2[2]);
        glVertex3f(_v3[0],_v3[1],_v3[2]); glVertex3f(_v4[0],_v4[1],_v4[2]);
        glVertex3f(_v5[0],_v5[1],_v5[2]); glVertex3f(_v6[0],_v6[1],_v6[2]);
        glVertex3f(_v7[0],_v7[1],_v7[2]); glVertex3f(_v8[0],_v8[1],_v8[2]);
        glEnd();
        //
        //glColor3f(0.8,0.8,0.8);
        glColor3f(0.0,0.0,1.0);
        glBegin(GL_LINES);// Z lines
        glVertex3f(_v1[0],_v1[1],_v1[2]); glVertex3f(_v5[0],_v5[1],_v5[2]);
        glVertex3f(_v2[0],_v2[1],_v2[2]); glVertex3f(_v6[0],_v6[1],_v6[2]);
        glVertex3f(_v3[0],_v3[1],_v3[2]); glVertex3f(_v7[0],_v7[1],_v7[2]);
        glVertex3f(_v4[0],_v4[1],_v4[2]); glVertex3f(_v8[0],_v8[1],_v8[2]);
        glEnd();
      }
    }
  }
}

void Fl_Gl_Mol_View::draw_scene(void){
  //
  glPointSize(2.0);
  glPushMatrix();
  glLoadIdentity();
  glRotatef(tb_angle, tb_axis[0], tb_axis[1], tb_axis[2]);
  // What the line below does?
  ////>glMultMatrixf((real*)tb_transform);  //<---------//?
  glMultMatrixf((GLfloat*)rot_matrix);
  glGetFloatv(GL_MODELVIEW_MATRIX, (GLfloat*)rot_matrix);
  glPopMatrix();
  glTranslatef(x_shift, y_shift, z_shift);
  glMultMatrixf((GLfloat*)rot_matrix);
  // zoom in and zoom out
  glScalef(zoom,zoom,zoom);
  //glTranslatef(x_shift,y_shift,z_shift);
  if(is_draw_atoms_){
    eval_system_properties();
    draw_atoms();
    // draw atomics bonds using the van der Waals radius
    // bond selection not implemented yet
    // using render_mode!=MODE_SELECT to avoid picking
    if(is_draw_bonds_ && render_mode!=MODE_SELECT)
      draw_bonds();
    glNormal3f(0,0,1);
    // the same as above
    if(is_draw_molecular_axis_ && render_mode!=MODE_SELECT){
      glColor3f(0.0,0.0,1.0);  // axis
      add_axis(v_axis_position, 5.0, 0.05, __axis_precession, __axis_tilt);
      glColor3f(1.0,0.0,0.0);  // backbone plane
      add_axis(v_axis_position, 5.0, 0.05, __backbone_precession, __backbone_tilt);
    }
  }
  if((is_draw_symbols_ || is_draw_labels_ || is_draw_numbers_) && render_mode!=MODE_SELECT){
    draw_symbols();
  }
  // draw the bounding box using the lattice vectors
  if(is_draw_bbox_ && render_mode!=MODE_SELECT)
    draw_box();
  // Molecular axisFL_HELVETICA
  // 1) rotate y axis, second angle, tilt
  // 2) rotate z axis, first angle, precession
  if(is_draw_world_axes_ && render_mode!=MODE_SELECT){
    // world coordinate axes
    glPushMatrix();
    glLoadIdentity();
    glTranslatef(0.85*view_left,0.85*view_bottom,700); //<----------------
    glMultMatrixf((GLfloat*)rot_matrix); //<----------------
    draw_axes(); //<----------------
    // translate the world axes to the lower left corner
    glPopMatrix();
  }//
  // draw the pie menu
  glPushMatrix();
  glLoadIdentity();
  // draw the pie menu
  //if(is_draw_pie_menu && render_mode!=MODE_SELECT){
  if(is_draw_pie_menu){
    draw_pie_menu(cursorX,cursorY, 792,base_view/4.0,100);
  }
  // draw the subpie
  //if(is_draw_pie_submenu && render_mode!=MODE_SELECT){
  if(is_draw_pie_submenu){
    draw_pie_submenu(792,base_view/4.0,100);
  }
  if(is_draw_tools_ && render_mode!=MODE_SELECT){
    draw_tools(790);
    draw_selected_numbers();
  }
  // draw the slide controls
  if(is_draw_controls && render_mode!=MODE_SELECT){
    draw_settings(790);
    draw_information(790);
    draw_controls(790);
  }
  // draw the processing message
  if(is_draw_processing)
    draw_message(794);
  glPopMatrix();
}

void Fl_Gl_Mol_View::draw(){
  initialize_opengl();
  // restoring the original projection matrix
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  //glMatrixMode(GL_MODELVIEW);
  //glFlush();
  // Reset the coordinate system before modifying
  //glMatrixMode(GL_PROJECTION);
  //glPushMatrix();
  glLoadIdentity();
  if(render_mode==MODE_SELECT){
    // process start picking
    // turn off texturing, lighting and fog during picking
    glDisable(GL_TEXTURE_2D);
    glDisable(GL_FOG);
    glDisable(GL_LIGHTING);
    glDisable (GL_BLEND);
    //glDisable (GL_DITHER);
    glDisable (GL_LIGHTING);
    glShadeModel(GL_FLAT);
    glClearColor(1,1,1,1);   // Background Color
  }
  // Set the viewport to be the entire window
  view_reshape(w(),h());
  glViewport(0,0,w(),h());
  //std::cout<<" W="<<w()<<" H="<<h()<<std::endl;
  if(w() <= h()){
        y_factor=(real)h()/(real)w();
        view_left   = -base_view;
        view_right  =  base_view;
        view_bottom = -base_view*y_factor;
        view_top    =  base_view*y_factor;
        view_axis_x =  view_left+6.0;
        view_axis_y =  view_bottom+6.0;
        glOrtho(view_left, view_right, view_bottom, view_top, view_near, view_far);
  }else{
        x_factor=(real)w()/(real)h();
        view_left   = -base_view*x_factor;
        view_right  =  base_view*x_factor;
        view_bottom = -base_view;
        view_top    =  base_view;
        view_axis_x =  view_left+6.0;
        view_axis_y =  view_bottom+6.0;
        glOrtho(view_left, view_right, view_bottom, view_top, view_near, view_far);
  }
  //set_font_size();
  scaled_light_position[0]=view_left*light_position[0];
  scaled_light_position[1]=view_top*light_position[1];
  scaled_light_position[2]=view_far*light_position[2];
  scaled_light_position[3]=light_position[3];
  glLightfv(GL_LIGHT0, GL_POSITION, scaled_light_position);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glMatrixMode(GL_MODELVIEW);
  //
  if(is_graphics_on){
    glPushMatrix();
    glLoadIdentity();
    draw_scene();
    glPopMatrix();
    glFlush();
  }
  //glMultMatrixf(tb_transform);
  // new picking process using color index
  if(render_mode==MODE_SELECT){
    //process_stop_picking();
    // get color information from frame buffer
    glGetIntegerv(GL_VIEWPORT, viewport);
    glReadPixels(cursorX, viewport[3]-cursorY, 1, 1, GL_RGB, GL_UNSIGNED_BYTE, pixel);
    process_picking(pixel);
    // reenable OpenGL
    glEnable(GL_LIGHTING);
    glEnable (GL_BLEND);
    glEnable (GL_LIGHTING);
    glShadeModel(GL_SMOOTH);           // Use smooth shading
    render_mode=MODE_RENDER;
    glClearColor(bgred,bggreen,bgblue,1.0);   // Background Color
    draw();
    if(is_handle_atom_){
      handle_atom_menu();
      is_handle_atom_=false;
    }
    if(is_handle_main_){
      handle_main_menu();
      is_handle_main_=false;
    }
    draw();
  }else{
    redraw();
  }
  //redraw();
}

void Fl_Gl_Mol_View::resize(int X,int Y,int W,int H) {
  //std::cout<<"X="<<X<<" Y="<<Y<<" W="<<W<<" H="<<H<<std::endl;
  Fl_Gl_Window::resize(X,Y,W,H);
  //is_draw_controls=true;
  //is_draw_tools_= false;
  //is_menu_position=true;
  redraw();
  redraw();
}

#endif /* HAVE_GL */

////////////////////////////HANDLE EVENTS////////////////////////////////
int Fl_Gl_Mol_View::handle(int event){
  static int last_x;
  static int last_y;
  int delta_x, delta_y;
  int key;
  //int key = Fl::event_key();
  //std::cout<<"Event: "<<key<<std::endl;
  //... position in Fl::event_x() and Fl::event_y()
  // get current mouse position and process event
  int x = Fl::event_x();
  int y = Fl::event_y();
  //int ret = Fl_Gl_Atoms::handle(event);
  int ret = Fl_Gl_Window::handle(event);
  switch(event){
  case FL_FOCUS:
    //if(Fl::focus() != this)
    //    Fl::focus(this);
    return 1;
  case FL_UNFOCUS:
    //... Return 1 if you want keyboard events, 0 otherwise
    //redraw();FL_HELVETICA
    return 1;
  case FL_PUSH:
    //... mouse down event ...
    //std::cout<<" push"<<std::endl;
    if(Fl::focus() != this)
        Fl::focus(this);
    cursorX = x;
    cursorY = y;
    render_mode=MODE_SELECT;
    if(Fl::event_button()==FL_RIGHT_MOUSE){
      //menu_mode=MODE_MENU;
      if(!is_draw_pie_menu)
        is_draw_pie_menu=true;
      is_menu_position=true;
      is_draw_line=false;
      is_draw_point=false;
      is_left_click=false;
      is_right_click=true;
      //is_menu_picked=false;
    }else if(Fl::event_button()==FL_LEFT_MOUSE){
    //else
      //is_menu_picked=true;
      is_menu_position=false;
      is_left_click=true;
      is_right_click=false;
        //if(is_draw_menu)
          //is_draw_menu=false;
    }
    // save mouse position to track drag events
    last_x = x;
    last_y = y;
    point_to_vector(last_x, last_y, tb_width, tb_height, tb_lastposition);
    ////>initialize_transform_matrix();
    redraw();
    return 1;
  case FL_DRAG:
    if((Fl::event_button()==FL_LEFT_MOUSE) && !is_lock_dragging){
      delta_x = x - last_x;
      delta_y = y - last_y;
      last_x = x;
      last_y = y;
      if(Fl::event_shift()){
        x_shift += shift_factor*delta_x*cos(DEG_RAD*__x_ang);
        y_shift -= shift_factor*(delta_y*cos(DEG_RAD*__y_ang)-delta_x*sin(DEG_RAD*__x_ang)*sin(DEG_RAD*__y_ang));
        z_shift += shift_factor*(delta_x*sin(DEG_RAD*__x_ang)*cos(DEG_RAD*__y_ang)+delta_y*sin(DEG_RAD*__y_ang));
      }else{
        set_mouse_motion(last_x,last_y);
      }
      if(is_draw_pie_menu && !Fl::event_is_click()){
        is_draw_pie_menu=false;
        is_draw_line=false;
        is_draw_point=false;
        u_active_menu=NOT_MENU;
      }
      redraw();
    }
    if(is_draw_controls==true && is_lock_controls!=true){
       is_draw_controls=false;
    }
    redraw();
    return 1;
  case FL_RELEASE:
    initialize_transform_matrix();
    render_mode=MODE_RENDER;
    //render_mode=MODE_SELECT;
    //redraw();
    return 1;
  case FL_MOUSEWHEEL:
    if(Fl::event_shift()){
      set_zoom(zoom-0.01*Fl::event_dy());
    }else{
      set_zoom(zoom-0.1*Fl::event_dy());
    }
    is_draw_pie_menu=false;
    is_draw_line=false;
    is_draw_point=false;
    u_active_menu=NOT_MENU;
    //if(Fl::event_ctrl()){
    //}else{
    //}
    redraw();
    return 1;
  //case FL_KEYBOARD:
  case FL_KEYDOWN:
    if(Fl::focus() != this)
      Fl::focus(this);
    //... keypress, key is in Fl::event_key(), ascii in Fl::event_text()
    key = Fl::event_key();
#ifdef _SHOW_KEY_MESSAGES_
    std::cout<<" pressed key: "<<key<<std::endl;
#endif
    switch(key){
      case FL_Page_Up:
        //std::cout<<" pressed key: Page_Up"<<std::endl;
        switch(u_slider_index){
          case 0:
            set_atom_radius_scale(f_atom_radius_scale+0.05);
            break;
          case 1:
            set_bond_radius_scale(f_bond_radius_scale+0.05);
            break;
          case 2:
            //set_highlihght_brightness(f_highlight_brightness+0.05);
            set_atom_brightness(f_atom_brightness+0.05);
            break;
          case 3:
            set_select_brightness(f_select_brightness+0.05);
            break;
          default:
            break;
        }
        is_draw_controls=true;
        redraw();
        return 1;
        break;
      case FL_Page_Down:
        //std::cout<<" pressed key: Page_Down"<<std::endl;
        switch(u_slider_index){
          case 0:
            set_atom_radius_scale(f_atom_radius_scale-0.05);
            break;
          case 1:
            set_bond_radius_scale(f_bond_radius_scale-0.05);
            break;
          case 2:
            //set_highlihght_brightness(f_highlight_brightness-0.1);
            set_atom_brightness(f_atom_brightness-0.05);
            break;
          case 3:
            set_select_brightness(f_select_brightness-0.05);
            break;
          default:
            break;
        }
        is_draw_controls=true;
        redraw();
        return 1;
        break;
      case 97:  // a) world axis
        is_draw_world_axes(!is_draw_world_axes_);
        is_draw_controls=true;
        redraw();
        return 1;
        break;
      case 98:  // b) bonds
        is_draw_bonds(!is_draw_bonds_);
        is_draw_controls=true;
        redraw();
        return 1;
        break;
      case 106:  // j) molecular axis
        is_draw_molecular_axis(!is_draw_molecular_axis_);
        //is_draw_labels(!is_draw_labels_);
        is_draw_controls=true;
        redraw();
        return 1;
        break;
      case 108:  // l) Labels
        is_draw_labels(!is_draw_labels_);
        is_draw_controls=true;
        redraw();
        return 1;
        break;
      case 107:  // k) lock controls
        set_lock_controls(!is_lock_controls);
        if(is_lock_controls)
          is_draw_controls=true;
        //is_lock_controls=false;
        redraw();
        return 1;
        break;
      case 110:  // n) numbers
        is_draw_numbers(!is_draw_numbers_);
        is_draw_controls=true;
        redraw();
        return 1;
        break;
      case 115:  // s) symbols
        is_draw_symbols(!is_draw_symbols_);
        is_draw_controls=true;
        redraw();
        return 1;
        break;
      case 116:  // t) tools
        is_draw_tools(!is_draw_tools_);
        is_draw_controls=true;
        redraw();
        return 1;
        break;
      case 118:  // v) volumen
        is_draw_bbox(!is_draw_bbox_);
        is_draw_controls=true;
        redraw();
        return 1;
        break;
      case 120:  // x) yz views
        if(Fl::event_shift()){
          set_view_xy_front();
        }else{
          set_view_xy_back();
        }
        is_draw_pie_menu=false;
        is_draw_line=false;
        is_draw_point=false;
        u_active_menu=NOT_MENU;
        redraw();
        break;
      case 121: // y) xz views
        if(Fl::event_shift()){
          set_view_yz_front();
        }else{
          set_view_yz_back();
        }
        is_draw_pie_menu=false;
        is_draw_line=false;
        is_draw_point=false;
        u_active_menu=NOT_MENU;
        redraw();
        break;
      case 122: // z) xy views
        if(Fl::event_shift()){
          set_view_zx_front();
        }else{
          set_view_zx_back();
        }
        is_draw_pie_menu=false;
        is_draw_line=false;
        is_draw_point=false;
        u_active_menu=NOT_MENU;
        redraw();
        break;
      // use the navigation keys to move an change the control menu
      case FL_Escape:
        clear_scene();
        redraw();
        break;
      case FL_Up:
        return 1;
      case FL_Down:
        return 1;
      case FL_Left:
        return 1;
      case FL_Right:
        return 1;
      case FL_Tab:
        is_slider_active[u_slider_index]=false;
        if(Fl::event_ctrl()) u_slider_index--;
        else u_slider_index++;
        if(u_slider_index > 3) u_slider_index = 0;
        else if(u_slider_index < 0) u_slider_index = 3;
        is_slider_active[u_slider_index]=true;
        is_draw_controls=true;
        redraw();
        return 1;
      default:
        return 1;
    }
    redraw();
    //... Return 1 if you understand/use the keyboard event, 0 otherwise...
    return 1;
  case FL_ENTER:
  case FL_LEAVE:
    return 1;
  case FL_SHORTCUT:
    //... shortcut, key is in Fl::event_key(), ascii in Fl::event_text()
    //... Return 1 if you understand/use the shortcut event, 0 otherwise...
    return 1;
  default:
    // pass other events to the base class...
    redraw();
    return ret;
  }
}

void Fl_Gl_Mol_View::set_mouse_motion(int x, int y){
  GLfloat current_position[3], dx, dy, dz;
  //assert(tb_button != -1);
  //if (tb_tracking == GL_FALSE)
    //return;
  point_to_vector(x, y, tb_width, tb_height, current_position);
  // calculate the angle to rotate by (directly proportional to the
  // length of the mouse movement.
  dx = current_position[0]-tb_lastposition[0];
  dy = current_position[1]-tb_lastposition[1];
  dz = current_position[2]-tb_lastposition[2];
  tb_angle=(90.0*sqrt(dx*dx+dy*dy+dz*dz));
  // calculate the axis of rotation (cross product)
  tb_axis[0]=tb_lastposition[1]*current_position[2]-tb_lastposition[2]*current_position[1];
  tb_axis[1]=tb_lastposition[2]*current_position[0]-tb_lastposition[0]*current_position[2];
  tb_axis[2]=tb_lastposition[0]*current_position[1]-tb_lastposition[1]*current_position[0];
  // reset for next time
  // tb_lasttime = glutGet(GLUT_ELAPSED_TIME);
  tb_lastposition[0] = current_position[0];
  tb_lastposition[1] = current_position[1];
  tb_lastposition[2] = current_position[2];
  // remember to draw new position
}

void Fl_Gl_Mol_View::point_to_vector(int x, int y, int width, int height, GLfloat v[3]){
  real d, a;
  // project x, y onto a hemi-sphere centered within width, height.
  v[0] = (2.0 * x - width) / width;
  v[1] = (height - 2.0 * y) / height;
  d = sqrt(v[0] * v[0] + v[1] * v[1]);
  v[2] = cos((3.14159265 / 2.0) * ((d < 1.0) ? d : 1.0));
  a = 1.0 / sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  v[0] *= a;
  v[1] *= a;
  v[2] *= a;
}

void Fl_Gl_Mol_View::set_view(real x,real y,real z){
  //init_rot_matrix();
  initialize_transform_matrix();
  //glMatrixMode(GL_PROJECTION);
  glMatrixMode(GL_MODELVIEW);
  // put the identity in the trackball transform
  glPushMatrix();
  glLoadIdentity();
  //glRotatef(__x_ang,1,0,0);
  glRotatef(x,1,0,0);
  glRotatef(y,0,1,0);
  glRotatef(z,0,0,1);
  glGetFloatv(GL_MODELVIEW_MATRIX,(GLfloat*)rot_matrix);
  glPopMatrix();
}

void Fl_Gl_Mol_View::view_reshape(int width, int height){
  //assert(tb_button != -1);
  tb_width  = width;
  tb_height = height;
}

void Fl_Gl_Mol_View::set_font_size(void){
  //GLint scale = 0;
  //scale = (GLint)(0.1*mini(w(), h()));
  //std::cout<<" Font scale: "<<scale<<std::endl;
  //font_size_symbol=scale+4;
  int offset = 0;
  font_size_symbol=14+offset;
  font_size_pie_label=14+offset;
  //font_size_panel_label=scale+4+offset;
  //font_size_panel_label=scale+2;
  font_size_panel_label=14+offset;
  font_size_slider_label=14+offset;
  /*
  font_size_symbol=12;
  font_size_pie_label=12;
  font_size_panel_label=10;
  font_size_slider_label=8;
  */
  //gl_font(FL_COURIER,12);  // text font
}
// Initialize GL
void Fl_Gl_Mol_View::initialize_opengl(void){
  if(!valid()) {
    glClearColor(bgred,bggreen,bgblue,1.0);   // Background Color
    glClearDepth(1.0f);
    //if(gm_depth)
    glEnable(GL_DEPTH_TEST);           // Enable Depth testing
    glDepthFunc(GL_LEQUAL);
    glShadeModel(GL_SMOOTH);           // Use smooth shading
    // Set the smooth shaiding to the best we can have
    glHint(GL_SHADE_MODEL, GL_NICEST);
    // Line options
    glEnable(GL_LINE_SMOOTH);               // Use smooth lines
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST); // Use smooth lines //
    glEnable(GL_BLEND); //<----TEST
    //}
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_NORMALIZE);
    //glEnable(GL_CULL_FACE);
    scaled_light_position[0]=base_view*light_position[0];
    scaled_light_position[1]=base_view*light_position[1];
    scaled_light_position[2]=base_view*light_position[2];
    scaled_light_position[3]=light_position[3];
    //  Enable Lighting
    glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
    //glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    glLightfv(GL_LIGHT0, GL_POSITION, scaled_light_position);
    // enable color tracking
    // set material properties which will be assigned by glColor
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    // material properties
    glMaterialf(GL_FRONT, GL_SHININESS, 50);
    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    glNormal3f(0,0,1);
    glPointSize(1.1);
    glLineWidth(1.1);
    //initialize_transform_matrix();
  }
}

// PICKING HANDLING //////////////////////////////////////////////////////////////

// Fri Sep  5 11:43:22 EDT 2014
// alpha version
void Fl_Gl_Mol_View::process_picking(unsigned char pc[3]){
  ui_rgb color;
#ifdef _GLMOL_DEBUG_PICKING_
  std::cout<<" pixel ---> "<<(uint)pc[0]<<"-"<<(uint)pc[1]<<"-"<<(uint)pc[2]<<std::endl;
#endif
  // to keep it safe while the menu picking is not implemented
  // it does not take care of the axes
  // all the intesive computation must be out of this function
  // here only process the picking by the mouse
  if((uint)pc[2] == 0){
      color.r = pc[0];
      color.g = pc[1];
      color.b = pc[2];
      uint idx = index_palette.get_index_rgb(color);
#ifdef _GLMOL_DEBUG_PICKING_
      std::cout<<" color picked ---> "<<color.r<<"-"<<color.g<<"-"<<color.b<<std::endl;
      std::cout<<" atom picked ---> "<<idx<<std::endl;
      std::cout<<" !You picked an atom!"<<std::endl;
#endif
      if(!is_lock_controls)
        is_draw_controls=false;
      set_highlight_atom(idx);
      if(is_draw_tools_)
        set_selected_atom(idx);
      is_atom_picked=true;
      if(is_draw_pie_menu){
        u_active_menu=ATOM_MENU;
        is_lock_dragging=true;
      }else{
        is_lock_dragging=false;
      }
      is_draw_point=true;
      is_draw_pie_submenu=false;
    }else if((uint)pc[2] == 255){ // Pie Menu or Controls
      is_atom_picked=false;
      is_lock_dragging=true;
      //if((uint)pc[0] >= 0 && (uint)pc[0] < 6 && !is_draw_pie_submenu){
      if((uint)pc[0] < 6 && !is_draw_pie_submenu){
        u_menu_index=(uint)pc[0];
#ifdef _GLMOL_DEBUG_PICKING_
        std::cout<<" Active menu: "<<u_active_menu<<std::endl;
        std::cout<<" You picked a menu: "<<u_menu_index<<std::endl;
        std::cout<<" Label picked: "<<legends[u_menu_index]<<std::endl;
#endif
        is_menu_pie_picked=true;
        if(!is_lock_controls)
          is_draw_controls=false;
        switch(u_menu_index){
          case CLOSE_MENU:
            is_menu_pie_picked=false;
            u_active_menu=NOT_MENU;
            is_draw_pie_menu=false;
            is_lock_dragging=false;
            break;
          default:
            if(u_active_menu==MAIN_MENU || u_active_menu==CONTROL_MENU){
              is_draw_pie_submenu=true;
#ifdef _GLMOL_DEBUG_PICKING_
              std::cout<<" draw submenu"<<std::endl;
#endif
            }else if(u_active_menu==ATOM_MENU){
              if(u_menu_index!=1){
                is_draw_pie_submenu=false;
                is_draw_pie_menu=false;
                //---->handle_atom_menu();
                is_handle_atom_ = true;
              }else{
                is_draw_pie_submenu=true;
              }
            }
#ifdef _GLMOL_DEBUG_PICKING_
            if(is_draw_pie_submenu){
              std::cout<<" [draw submenu]"<<std::endl;
            }
#endif
            break;
        }
      }else if((uint)pc[0] >= 6 && (uint)pc[0] < 12){ // Pie Menu
        u_submenu_index=(uint)pc[0]-6;
#ifdef _GLMOL_DEBUG_PICKING_
        std::cout<<"You picked a submenu: "<<u_submenu_index<<std::endl;
        std::cout<<" label picked: "<<sub_legends[u_submenu_index]<<std::endl;
#endif
        //is_menu_pie_picked=true;
        is_draw_pie_submenu=false;
        if(!is_lock_controls)
          is_draw_controls=false;
        switch(u_active_menu){
          case MAIN_MENU:
            //---->handle_main_menu();
            is_handle_main_ = true;
          break;
          case ATOM_MENU:
            is_handle_atom_ = true;
            //----> handle_atom_menu();
          break;
        }
      }else{
#ifdef _GLMOL_DEBUG_PICKING_
        std::cout<<" [background picked]"<<std::endl;
#endif
        if(is_draw_tools_){
          u_selected_index=0;
          r_distance1=0.0;
          r_distance2=0.0;
          r_distance3=0.0;
          v_distance1.zero();
          v_distance2.zero();
          v_distance3.zero();
          r_angle1=0.0;
          r_angle2=0.0;
          r_dihedral=0.0;
          update_normal_color=true;
          is_update_mask_rcolor=true;
        }
        is_background_picked=true;
        if(is_right_click){
          u_active_menu=MAIN_MENU;
        }else if(is_left_click){
          u_active_menu=NOT_MENU;
          is_menu_pie_picked=false;
          is_draw_pie_menu=false;
          is_lock_dragging=false;
        }
        if(!is_lock_controls)
          is_draw_controls=false;
        is_atom_picked=false;
        is_lock_dragging=false;
        is_draw_pie_submenu=false;
        is_draw_line=false;
        is_draw_point=false;
      }
  }
}

// handle funtions on the main menu
void Fl_Gl_Mol_View::handle_main_menu(void){
  //uint _af;
  switch(u_menu_index){
    case MAIN_MODE_SUBMENU:
        switch(u_submenu_index){
          case MODE_PBC_BUTTON:
            redraw();
            break;
          case MODE_SYMMETRY_BUTTON:
            redraw();
            break;
          case MODE_FRAGMENTS_BUTTON:
            is_mode_atom=false;
            is_highlight_fragment(true);
            redraw();
            break;
          case MODE_ATOMS_BUTTON:
            is_mode_atom=true;
            is_highlight_fragment(false);
            redraw();
            break;
          case MODE_BLANK_BUTTON:
            //is_highlight_fragment(true);
            //redraw();
            //break;
          default:
            break;
    }
    break;
    case MAIN_VIEW_SUBMENU:
        switch(u_submenu_index){
          case VIEW_XY_BUTTON:
            set_view_xy_front();
            break;
          case VIEW_YZ_BUTTON:
            set_view_yz_front();
            break;
          case VIEW_ZX_BUTTON:
            set_view_zx_front();
            break;
          case VIEW_YX_BUTTON:
            set_view_xy_back();
            break;
          case VIEW_ZY_BUTTON:
            set_view_yz_back();
            break;
          case VIEW_XZ_BUTTON:
            set_view_zx_back();
            break;
          default:
            break;
    }
    redraw();
    break;
    // view submenu
    case MAIN_SHOW_SUBMENU:
      switch(u_submenu_index){
          case SHOW_LABELS_BUTTON:
            is_draw_labels(!is_draw_labels_);
            redraw();
            break;
          case SHOW_BOX_BUTTON:
            is_draw_bbox(!is_draw_bbox_);
            redraw();
            break;
          case SHOW_SYMBOLS_BUTTON:
            is_draw_symbols(!is_draw_symbols_);
            redraw();
            break;
          case SHOW_BONDS_BUTTON:
            is_draw_bonds(!is_draw_bonds_);
            redraw();
            break;
          case SHOW_NUMBERS_BUTTON:
            is_draw_numbers(!is_draw_numbers_);
            break;
          default:
            //is_lock_dragging=false;
            break;
      }
      break;
    case MAIN_TOOLS_SUBMENU:
      switch(u_submenu_index){
          case TOOLS_MESURE_BUTTON:
            is_draw_tools(!is_draw_tools_);
            redraw();
            break;
          case TOOLS_BOX_BUTTON:
            is_draw_bbox(!is_draw_bbox_);
            redraw();
            break;
          case TOOLS_AXES_BUTTON:
            is_draw_world_axes(!is_draw_world_axes_);
            redraw();
            break;
          case TOOLS_FRAGMENT_BUTTON:
            //update_coordinates=true;
            //update_data();
            is_draw_processing=true;
            render_mode=MODE_RENDER;
            redraw();
            //is_draw_processing=true;
            //redraw();
            //Fl::wait(0.1);
            compute_vdw_fragments();
            //set_update_coordinates(true);
            //set_active_fragment_index(1);
            set_update_active_fragment(); // the same as zero above.
            //update_data();
            is_draw_processing=false;
            // switch to fragment mode
            if(is_mode_atom){
              is_mode_atom=false;
              is_highlight_fragment(true);
            }
            // switch to fragment mode
            render_mode=MODE_RENDER;
            redraw();
            break;
          case TOOLS_ATOMS_BUTTON:
            is_draw_processing=true;
            render_mode=MODE_RENDER;
            redraw();
            compute_atom_fragments();
            set_update_active_fragment(); // the same as zero above.
            is_draw_processing=false;
            // switch to fragment mode
            if(is_mode_atom){
              is_mode_atom=false;
              is_highlight_fragment(true);
            }
            render_mode=MODE_RENDER;
            redraw();
            break;
          default:
            //is_lock_dragging=false;
            break;
      }
      break;
    default:
      break;
  }
}

// handle operation on the picked atom
void Fl_Gl_Mol_View::handle_atom_menu(void){
  //uint _af;
  switch(u_menu_index){
    /*case ATOM_ELEMENT_SUBMENU:
      set_active_fragment(__highlight_atom);
      //_af= v_fragment_table_gl[__highlight_atom];
      //set_map_active_fragment(_af-1);
      //set_update_active_fragment(_af);
      //update_coordinates=true;
      //
      compute_vdw_fragment(__highlight_atom);
      //set_atom_fragment(__highlight_atom);
      //set_update_coordinates(true);
      //set_map_active_fragment(_af-1);
      //set_update_active_fragment(_af); // the same as zero above.
      //
      set_update_active_fragment();
      //render_mode=MODE_RENDER;
      //is_draw_processing=true;
      //redraw();
      //Fl::wait(0.1);
      //set_update_coordinates(true);
      //set_map_active_fragment(_af-1);
      //set_update_active_fragment(_af); // the same as zero above.
      //is_draw_processing=false;
      // switch to fragment mode
      if(is_mode_atom){
        is_mode_atom=false;
        is_highlight_fragment(true);
      }
      // switch to fragment mode
      render_mode=MODE_RENDER;
      redraw();
      break;*/
    case ATOM_AXIS_SUBMENU:
      switch(u_submenu_index){
          case AXIS_SHOW_BUTTON:
            is_draw_molecular_axis(!is_draw_molecular_axis_);
            redraw();
            break;
          case AXIS_HEAD_BUTTON:
            //set_fragment_axis(__highlight_atom);
            set_active_fragment(__highlight_atom);
            set_update_active_fragment();
            //redraw();
            supercell.set_fragment_axis(0,__highlight_atom);
            supercell.set_fragment_axis(2,__highlight_atom);
            set_active_fragment(__highlight_atom);
            set_update_active_fragment();
            render_mode=MODE_RENDER;
            redraw();
            //is_draw_processing=false;
            break;
          case AXIS_TAIL_BUTTON:
            set_active_fragment(__highlight_atom);
            set_update_active_fragment();
            supercell.set_fragment_axis(1,__highlight_atom);
            set_active_fragment(__highlight_atom);
            set_update_active_fragment();
            render_mode=MODE_RENDER;
            redraw();
            //is_draw_processing=false;
            break;
          case AXIS_PLANE_BUTTON:
            set_active_fragment(__highlight_atom);
            set_update_active_fragment();
            supercell.set_fragment_axis(3,__highlight_atom);
            set_active_fragment(__highlight_atom);
            set_update_active_fragment();
            render_mode=MODE_RENDER;
            redraw();
            //is_draw_processing=false;
            break;
          default:
            //is_lock_dragging=false;
            break;
      }
      break;
    //case ATOM_ACTIVE_SUBMENU:
      //switch(u_submenu_index){
    case ATOM_ACTIVE_SUBMENU:
      set_active_fragment(__highlight_atom);
      //_af= v_fragment_table_gl[__highlight_atom];
      //set_map_active_fragment(_af-1);
      //set_update_active_fragment(_af);
      //update_coordinates=true;
      //update_data();
      set_update_active_fragment();
      render_mode=MODE_RENDER;
      //is_mode_atom=false;
      redraw();
      break;
    case ATOM_DELETE_SUBMENU:
      //is_draw_pie_submenu=false;
      break;
    case ATOM_FRAGMENT_SUBMENU:
      switch(u_submenu_index){
        case ATOM_FRAGMENT_VDW_BUTTON:
          set_active_fragment(__highlight_atom);
          compute_vdw_fragment(__highlight_atom);
          set_update_active_fragment();
          // switch to fragment mode
          if(is_mode_atom){
            is_mode_atom=false;
            is_highlight_fragment(true);
          }
          // switch to fragment mode
          render_mode=MODE_RENDER;
          redraw();
          break;
        case ATOM_FRAGMENT_ATOM_BUTTON:
          set_active_fragment(__highlight_atom);
          compute_atom_fragment(__highlight_atom);
          set_update_active_fragment();
          //is_draw_processing=true;
          //is_draw_processing=false;
          // switch to fragment mode
          if(is_mode_atom){
            is_mode_atom=false;
            is_highlight_fragment(true);
          }
          // switch to fragment mode
          render_mode=MODE_RENDER;
          redraw();
          break;
        case ATOM_FRAGMENT_RADIAL_BUTTON:
          if ( r_distance1 > 0.0 ){
            set_active_fragment(__highlight_atom);
            compute_radial_fragment(__highlight_atom,1.1*r_distance1);
            set_update_active_fragment();
            //is_draw_processing=true;
            //is_draw_processing=false;
            // switch to fragment mode
            if(is_mode_atom){
              is_mode_atom=false;
              is_highlight_fragment(true);
            }
            // switch to fragment mode
            render_mode=MODE_RENDER;
            redraw();
          }
          break;
        case ATOM_FRAGMENT_MERGE_BUTTON:
          //cell.eval_merge_fragment(__highlight_atom,false);
          compute_merge_fragments(__highlight_atom);
          set_update_active_fragment();
          // switch to fragment mode
          if(is_mode_atom){
            is_mode_atom=false;
            is_highlight_fragment(true);
          }
          render_mode=MODE_RENDER;
          redraw();
          break;
        default:
          //is_lock_dragging=false;
          break;
      }
      break;
    case ATOM_ELEMENT_SUBMENU:
    default:
      break;
  }
}
// GL_Toolkit /////////////////////////////////////////////////////////////////

// control box
void Fl_Gl_Mol_View::draw_controls(GLfloat z){
  char buff[10];
  GLfloat y_hight = 0.33*view_bottom;
  GLfloat x_width = 0.8*base_view;
  GLfloat x_pos = -0.5*x_width;//0.8*view_left;
  GLfloat y_pos = 0.65*view_bottom;
  GLfloat z_pos = z;
  GLfloat x_pos_end = x_pos+x_width;
  GLfloat y_pos_end = y_pos+y_hight;
  //std::cout<<" y_hight="<<y_hight<<" x_width="<<x_width<<std::endl;
  //glDisable(GL_DEPTH_TEST);           // Enable Depth testing
  //glEnable(GL_BLEND);
  //glEnable(GL_LINE_SMOOTH);
  //glEnable(GL_POLYGON_SMOOTH);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glNormal3f(0,0,1);
  glColor4f(0.3,0.3,0.3,0.7);
  // background for the control menu
  glBegin(GL_POLYGON); // draw in triangle strips
  glVertex3f( x_pos,     y_pos, z);
  glVertex3f( x_pos_end, y_pos, z);
  glVertex3f( x_pos_end, (y_pos+y_hight), z);
  glVertex3f( x_pos,     (y_pos+y_hight), z);
  glEnd();
  // draw the control box
  // solid border for the control menu
  //
  //gl_font(FL_COURIER,12); // text font
  //gl_font(FL_SYMBOL,12);
#if defined (BUILD_FOR_MACOS)
    gl_font(FL_HELVETICA,12); // text font
#elif defined (BUILD_FOR_WINDOWS)
    gl_font(FL_COURIER,12); // text font
    //gl_font(FL_HELVETICA,14); // text font
#else
    gl_font(FL_COURIER,14); // text font
#endif
  //gl_font(FL_COURIER,font_size_panel_label); // text font
  glColor4f(0.0,1.0,0.0,1.0); // text color
  sprintf(buff,"%s","Controls");
  glRasterPos3f(x_pos, 0.638*view_bottom,z);
  gl_draw(buff, strlen(buff));
  glColor4f(0.8,0.8,0.8,1.0);
  glBegin(GL_LINE_LOOP);
  glVertex3f( x_pos,     y_pos, z_pos);
  glVertex3f( x_pos_end, y_pos, z_pos);
  glVertex3f( x_pos_end, (y_pos+y_hight), z);
  glVertex3f( x_pos,     (y_pos+y_hight), z);
  glEnd();
  //glColor4f(0.5,0.45,0.27,0.7);
  //glColor4f(0.98,0.87,0.34,0.8);
  //draw_slider(GLfloat x1, GLfloat y1, GLfloat y2, GLfloat z1, GLfloat val, char* l)
  // slider bar 1 /////////////////////////////////////////////////////////////////
  draw_slider(x_pos, y_pos, y_pos_end, z, f_atom_radius_scale, is_slider_active[0], (char*)" ");
  // slider bar 2 /////////////////////////////////////////////////////////////////
  draw_slider(x_pos+0.10*base_view, y_pos, y_pos_end, z, f_bond_radius_scale, is_slider_active[1], (char*)" ");
  // slider bar 3 /////////////////////////////////////////////////////////////////
  draw_slider(x_pos+0.20*base_view, y_pos, y_pos_end, z, f_atom_brightness/f_atom_brightness_max, is_slider_active[2], (char*)" ");
  // slider bar 4 /////////////////////////////////////////////////////////////////
  draw_slider(x_pos+0.30*base_view, y_pos, y_pos_end, z, f_select_brightness/f_select_brightness_max, is_slider_active[3], (char*)" ");
  /////////////////////////////////////////////////////////////////////////////////
  // Draw the radios for keyboard shortcuts
  // Radio 1 (axis)    [a]
  draw_radio_button(x_pos+0.40*base_view, y_pos, y_pos_end, z, 1.0, is_radio_active[0], (char*)"a");
  // Radio 2 (bonds)   [b]
  draw_radio_button(x_pos+0.50*base_view, y_pos, y_pos_end, z, 1.0, is_radio_active[1], (char*)"b");
  // Radio 3 (symbols) [s]
  draw_radio_button(x_pos+0.40*base_view, y_pos+0.05*view_bottom, y_pos_end, z, 1.0, is_radio_active[2], (char*)"s");
  // Radio 4 (volume)  [v]
  draw_radio_button(x_pos+0.50*base_view, y_pos+0.05*view_bottom, y_pos_end, z, 1.0, is_radio_active[3], (char*)"v");
  // Radio 5 (numbers) [n]
  draw_radio_button(x_pos+0.40*base_view, y_pos+0.10*view_bottom, y_pos_end, z, 1.0, is_radio_active[4], (char*)"n");
  // Radio 6 (tools)   [t]
  draw_radio_button(x_pos+0.50*base_view, y_pos+0.10*view_bottom, y_pos_end, z, 1.0, is_radio_active[5], (char*)"t");
  // Radio 7 (labels)  [l]
  draw_radio_button(x_pos+0.40*base_view, y_pos+0.15*view_bottom, y_pos_end, z, 1.0, is_radio_active[6], (char*)"l");
  // Radio 8 (tools)   [t]
  draw_radio_button(x_pos+0.50*base_view, y_pos+0.15*view_bottom, y_pos_end, z, 1.0, is_radio_active[7], (char*)"j");
  // Radio 9 (labels)  [m]
  draw_radio_button(x_pos+0.40*base_view, y_pos+0.20*view_bottom, y_pos_end, z, 1.0, is_radio_active[8], (char*)"k");
  //glDisable(GL_BLEND);
  //glEnable(GL_DEPTH_TEST);           // Enable Depth testing
  //glDisable(GL_POLYGON_SMOOTH);
}

void Fl_Gl_Mol_View::draw_message(GLfloat z){
  char buff[50];
  GLfloat y_hight = 0.1*view_bottom;
  GLfloat x_width = 0.5*base_view;
  GLfloat x_pos = -0.5*x_width;//0.8*view_left;
  GLfloat y_pos = 0.06*view_top;
  GLfloat z_pos = z;
  GLfloat x_pos_end = x_pos+x_width;
  //GLfloat y_pos_end = y_pos+y_hight;
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glNormal3f(0,0,1);
  glColor4f(0.3,0.3,0.3,0.7);
  // background for the control menu
  glBegin(GL_POLYGON); // draw in triangle strips
  glVertex3f( x_pos,     y_pos, z);
  glVertex3f( x_pos_end, y_pos, z);
  glVertex3f( x_pos_end, (y_pos+y_hight), z);
  glVertex3f( x_pos,     (y_pos+y_hight), z);
  glEnd();
  // draw the control box
  // solid border for the control menu
  //gl_font(FL_COURIER,font_size_panel_label); // text font
  gl_font(FL_COURIER,14); // text font
  glColor4f(0.0,1.0,0.0,1.0); // text color
  sprintf(buff,"%s","Processing...");
  glRasterPos3f(x_pos+0.3*x_width, 0.0,z);
  gl_draw(buff, strlen(buff));
  glColor4f(0.8,0.8,0.8,1.0);
  glBegin(GL_LINE_LOOP);
  glVertex3f( x_pos,     y_pos, z_pos);
  glVertex3f( x_pos_end, y_pos, z_pos);
  glVertex3f( x_pos_end, (y_pos+y_hight), z);
  glVertex3f( x_pos,     (y_pos+y_hight), z);
  glEnd();
}

void Fl_Gl_Mol_View::draw_tools(GLfloat z){
  char buff[10];
  //GLdouble winX=200, winY=0, winZ=791;
  //GLdouble posX, posY, posZ;
  //gluUnProject( winX, winY, winZ, modelview, projection, viewport, &posX, &posY, &posZ);
  GLfloat y_hight = -0.26*base_view;
  //GLfloat x_width = 0.45*base_view;
  //GLfloat x_width = -0.45*view_left;
  //GLfloat x_width = 5*x_factor;
  //GLfloat x_width = posX;
  GLfloat x_pos = 0.95*view_left;//0.8*view_left;
#ifdef PLATFORM_MAC
  GLfloat y_pos = 0.94*view_top;
#else
  GLfloat y_pos = 0.955*view_top;
#endif
  //GLfloat x_pos_end = x_pos+x_width;
  GLfloat y_pos_end = y_pos+y_hight;
  //glDisable(GL_DEPTH_TEST);           // Enable Depth testing
  //glEnable(GL_BLEND);
  //glEnable(GL_LINE_SMOOTH);
  //glEnable(GL_POLYGON_SMOOTH);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glNormal3f(0,0,1);
  glColor4f(0.3,0.3,0.3,1.0);
  // background for the tool box
  //glBegin(GL_POLYGON); /draw in triangle strips
  //glVertex3f( x_pos,     y_pos, z);
  //glVertex3f( x_pos_end, y_pos, z);
  //glVertex3f( x_pos_end, (y_pos+y_hight), z);
  //glVertex3f( x_pos,     (y_pos+y_hight), z);
  //glEnd();
  // draw the tool box
  // solid border for the control menu
#if defined (BUILD_FOR_MACOS)
    gl_font(FL_HELVETICA,12); // text font
#elif defined (BUILD_FOR_WINDOWS)
    //gl_font(FL_HELVETICA,14); // text font
#else
    gl_font(FL_COURIER,12); // text font
#endif
  //gl_font(FL_COURIER,font_size_panel_label); // text font
  glColor4f(0.0,1.0,0.0,1.0); // text color
  sprintf(buff,"%s","Tools");
  glRasterPos3f(x_pos, y_pos,z);
  gl_draw(buff, strlen(buff));
  //glColor4f(0.8,0.8,0.8,1.0);
  //glBegin(GL_LINE_LOOP);
  //glVertex3f( x_pos,     y_pos, z_pos);
  //glVertex3f( x_pos_end, y_pos, z_pos);
  //glVertex3f( x_pos_end, (y_pos+y_hight), z);
  //glVertex3f( x_pos,     (y_pos+y_hight), z);
  //glEnd();
  // widgets
  // text output 1
  GLfloat y_delta = 0.04*base_view;
  y_pos-=y_delta;
  x_pos+=(0.03*base_view);
  int a = v_selected_atoms[0] + 1;
  int b = v_selected_atoms[1] + 1;
  widget_float_output_xy(x_pos, y_pos, y_pos_end, z, a, b, r_distance1, (char*)"Distance (1,2):");
  y_pos-=y_delta;
  widget_vector_output(x_pos, y_pos, y_pos_end, z, v_distance1, (char*)"Cartesian (1,2):");
  y_pos-=y_delta;
  widget_float_output(x_pos, y_pos, y_pos_end, z, r_distance2, (char*)"Distance (2,3):");
  y_pos-=y_delta;
  widget_vector_output(x_pos, y_pos, y_pos_end, z, v_distance2, (char*)"Cartesian (2,3):");
  y_pos-=y_delta;
  widget_float_output(x_pos, y_pos, y_pos_end, z, r_distance3, (char*)"Distance (3,4):");
  y_pos-=y_delta;
  widget_vector_output(x_pos, y_pos, y_pos_end, z, v_distance3, (char*)"Cartesian (3,4):");
  y_pos-=y_delta;
  widget_float_output(x_pos, y_pos, y_pos_end, z, r_angle1, (char*)"Angle  (1,2,3):");
  y_pos-=y_delta;
  widget_float_output(x_pos, y_pos, y_pos_end, z, r_angle2, (char*)"Angle  (2,3,4):");
  y_pos-=y_delta;
  widget_float_output(x_pos, y_pos, y_pos_end, z, r_dihedral, (char*)"Dihedral:      ");
}

void Fl_Gl_Mol_View::draw_settings(GLfloat z){
  char buff[10];
  GLfloat y_hight = -0.26*base_view;
  GLfloat x_pos = 0.95*view_left;//0.8*view_left;
#ifdef PLATFORM_MAC
  GLfloat y_pos = 0.52*view_top;
#else
  GLfloat y_pos = 0.55*view_top;
#endif
  //GLfloat z_pos = z;
  //GLfloat x_pos_end = x_pos+x_width;
  GLfloat y_pos_end = y_pos+y_hight;
  //glDisable(GL_DEPTH_TEST);           // Enable Depth testing
  //glEnable(GL_BLEND);
  //glEnable(GL_LINE_SMOOTH);
  //glEnable(GL_POLYGON_SMOOTH);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glNormal3f(0,0,1);
  glColor4f(0.3,0.3,0.3,1.0);
  // background for the tool box
  //glBegin(GL_POLYGON); /draw in triangle strips
  //glVertex3f( x_pos,     y_pos, z);
  //glVertex3f( x_pos_end, y_pos, z);
  //glVertex3f( x_pos_end, (y_pos+y_hight), z);
  //glVertex3f( x_pos,     (y_pos+y_hight), z);
  //glEnd();
  // draw the tool box
  // solid border for the control menu
  //gl_font(FL_COURIER,font_size_panel_label); // text font
#if defined (BUILD_FOR_MACOS)
    gl_font(FL_HELVETICA,12); // text font
#elif defined (BUILD_FOR_WINDOWS)
    gl_font(FL_HELVETICA,14); // text font
#else
  gl_font(FL_COURIER,12); // text font
#endif
  //gl_font(FL_HELVETICA,font_size_pie_label); // text font
  glColor4f(0.5,1.0,0.5,1.0); // text color
  sprintf(buff,"%s","Settings");
  glRasterPos3f(x_pos, y_pos,z);
  gl_draw(buff, strlen(buff));
  //glColor4f(0.8,0.8,0.8,1.0);
  //glBegin(GL_LINE_LOOP);
  //glVertex3f( x_pos,     y_pos, z_pos);
  //glVertex3f( x_pos_end, y_pos, z_pos);
  //glVertex3f( x_pos_end, (y_pos+y_hight), z);
  //glVertex3f( x_pos,     (y_pos+y_hight), z);
  //glEnd();
  // widgets
  // text output 1
  GLfloat y_delta = 0.04*base_view;
  x_pos+=(0.03*base_view);
  // slider bar 1 /////////////////////////////////////////////////////////////////
  y_pos-=y_delta;
  widget_float_output(x_pos, y_pos, y_pos_end, z, f_atom_radius_scale, (char*)"Atom Radius:");
  // slider bar 2 /////////////////////////////////////////////////////////////////
  y_pos-=y_delta;
  widget_float_output(x_pos, y_pos, y_pos_end, z, f_bond_radius_scale, (char*)"Bond Radius:");
  // slider bar 3 /////////////////////////////////////////////////////////////////
  y_pos-=y_delta;
  widget_float_output(x_pos, y_pos, y_pos_end, z, f_atom_brightness/f_atom_brightness_max, (char*)"Brightness:");
  // slider bar 4 /////////////////////////////////////////////////////////////////
  y_pos-=y_delta;
  widget_float_output(x_pos, y_pos, y_pos_end, z, f_select_brightness/f_select_brightness_max, (char*)"Selection:");
  y_pos-=y_delta;
  draw_switch_output(x_pos, y_pos, y_pos_end, z,  is_radio_active[0],(char*)"Axes:    <A>");
  y_pos-=y_delta;
  draw_switch_output(x_pos, y_pos, y_pos_end, z, is_radio_active[2], (char*)"Symbols: <S>");
  y_pos-=y_delta;
  draw_switch_output(x_pos, y_pos, y_pos_end, z, is_radio_active[4], (char*)"Numbers: <N>");
  y_pos-=y_delta;
  draw_switch_output(x_pos, y_pos, y_pos_end, z, is_radio_active[6], (char*)"Labels:  <L>");
  y_pos-=y_delta;
  draw_switch_output(x_pos, y_pos, y_pos_end, z, is_radio_active[1], (char*)"Bonds:   <B>");
  y_pos-=y_delta;
  draw_switch_output(x_pos, y_pos, y_pos_end, z, is_radio_active[3], (char*)"Volume:  <V>");
  y_pos-=y_delta;
  draw_switch_output(x_pos, y_pos, y_pos_end, z, is_radio_active[5], (char*)"Tools:   <T>");
  y_pos-=y_delta;
  draw_switch_output(x_pos, y_pos, y_pos_end, z, is_radio_active[7], (char*)"Axis:    <J>");
}

void Fl_Gl_Mol_View::draw_information(GLfloat z){
  char buff[20];
  GLfloat y_hight = -0.26*base_view;
  GLfloat x_pos = 0.95*view_left;//0.8*view_left;
#ifdef PLATFORM_MAC
  GLfloat y_pos = -0.03*view_top;
#else
  GLfloat y_pos = 0.03*view_top;
#endif
  //GLfloat z_pos = z;
  //GLfloat x_pos_end = x_pos+x_width;
  GLfloat y_pos_end = y_pos+y_hight;
  //glDisable(GL_DEPTH_TEST);           // Enable Depth testing
  //glEnable(GL_BLEND);
  //glEnable(GL_LINE_SMOOTH);
  //glEnable(GL_POLYGON_SMOOTH);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glNormal3f(0,0,1);
  glColor4f(0.3,0.3,0.3,0.7);
  // background for the tool box
  //glBegin(GL_POLYGON); /draw in triangle strips
  //glVertex3f( x_pos,     y_pos, z);
  //glVertex3f( x_pos_end, y_pos, z);
  //glVertex3f( x_pos_end, (y_pos+y_hight), z);
  //glVertex3f( x_pos,     (y_pos+y_hight), z);
  //glEnd();
  // draw the tool box
  // solid border for the control menu
#if defined (BUILD_FOR_MACOS)
    gl_font(FL_HELVETICA,12); // text font
#elif defined (BUILD_FOR_WINDOWS)
    gl_font(FL_HELVETICA,14); // text font
#else
    gl_font(FL_COURIER,12); // text font
#endif
  //gl_font(FL_COURIER,font_size_panel_label); // text font
  glColor4f(0.0,1.0,0.0,1.0); // text color
  sprintf(buff,"%s","Information");
  glRasterPos3f(x_pos, y_pos,z);
  gl_draw(buff, strlen(buff));
  //glColor4f(0.8,0.8,0.8,1.0);
  //glBegin(GL_LINE_LOOP);
  //glVertex3f( x_pos,     y_pos, z_pos);
  //glVertex3f( x_pos_end, y_pos, z_pos);
  //glVertex3f( x_pos_end, (y_pos+y_hight), z);
  //glVertex3f( x_pos,     (y_pos+y_hight), z);
  //glEnd();
  // widgets
  // text output 1
  GLfloat y_delta = 0.04*base_view;
  x_pos+=y_delta;
  y_pos-=1.5*y_delta;
  widget_int_output(x_pos, y_pos, y_pos_end, z, __number_of_atoms, (char*)"Atoms:");
  y_pos-=y_delta;
  widget_int_output(x_pos, y_pos, y_pos_end, z, __fragment_total, (char*)"Fragments:");
  y_pos-=y_delta;
  widget_vector_output(x_pos, y_pos, y_pos_end, z, v_bbox, (char*)"Box:");
  //widget_float_output(x_pos, y_pos, y_pos_end, z, f_bond_radius_scale, (char*)"Bond Radius:");
  //y_pos-=y_delta;
  //widget_float_output(x_pos, y_pos, y_pos_end, z, f_atom_brightness/f_atom_brightness_max, (char*)"Brightness:");
  //y_pos-=y_delta;
  //widget_float_output(x_pos, y_pos, y_pos_end, z, f_select_brightness/f_select_brightness_max, (char*)"Selection:");
  //y_pos-=y_delta;
  //draw_switch_output(x_pos, y_pos, y_pos_end, z,  is_radio_active[0],(char*)"Axes:    <A>");
  //y_pos-=y_delta;
  //draw_switch_output(x_pos, y_pos, y_pos_end, z, is_radio_active[2], (char*)"Symbols: <S>");
  //y_pos-=y_delta;
  //draw_switch_output(x_pos, y_pos, y_pos_end, z, is_radio_active[4], (char*)"Numbers: <N>");
  //y_pos-=y_delta;
  //draw_switch_output(x_pos, y_pos, y_pos_end, z, is_radio_active[6], (char*)"Labels:  <L>");
  //y_pos-=y_delta;
  //draw_switch_output(x_pos, y_pos, y_pos_end, z, is_radio_active[1], (char*)"Bonds:   <B>");
  //y_pos-=y_delta;
  //draw_switch_output(x_pos, y_pos, y_pos_end, z, is_radio_active[3], (char*)"Volume:  <V>");
  //y_pos-=y_delta;
  //draw_switch_output(x_pos, y_pos, y_pos_end, z, is_radio_active[5], (char*)"Tools:   <T>");
  //y_pos-=y_delta;
  //draw_switch_output(x_pos, y_pos, y_pos_end, z, is_radio_active[7], (char*)"Axis:    <J>");
}

// slider widget
void Fl_Gl_Mol_View::draw_slider(GLfloat x1, GLfloat y1, GLfloat y2, GLfloat z, GLfloat val, bool active, char* l){
  GLfloat side_x1, side_x2;
  GLfloat side_y1, side_y2;
  glNormal3f(0,0,1);
  if(active)
    glColor4f(0.78,0.68,0.0,0.8);
  else
    glColor4f(0.98,0.87,0.64,0.7);
  // slider bar
  side_x1 = x1+0.11*base_view;
  side_x2 = x1+0.14*base_view;
  side_y1 = y1+0.30*view_bottom;
  side_y2 = (side_y1-val*0.24*view_bottom);
  glBegin(GL_POLYGON);
  glVertex3f(side_x1, side_y2, z);
  glVertex3f(side_x2, side_y2, z);
  glVertex3f(side_x2, side_y1, z);
  glVertex3f(side_x1, side_y1, z);
  glEnd();
  // slider box
  side_x1 = x1+0.10*base_view;
  side_x2 = x1+0.15*base_view;
  side_y1 = y1+0.31*view_bottom;
  side_y2 = (side_y1-0.26*view_bottom);
  glRasterPos3f(side_x1, 0.69*view_bottom,z);
  //gl_font(FL_COURIER,font_size_slider_label); // text font
  gl_draw(l, strlen(l));
  glBegin(GL_LINE_LOOP);
  glVertex3f(side_x1, side_y2, z); // left corner of the roof
  glVertex3f(side_x2, side_y2, z); // right corner of the roof
  glVertex3f(side_x2, side_y1, z); // bottom right corner of the house
  glVertex3f(side_x1, side_y1, z); // bottom left corner of the house
  glEnd();
}

// radio button widget
void Fl_Gl_Mol_View::draw_radio_button(GLfloat x1, GLfloat y1, GLfloat y2, GLfloat z, GLfloat val, bool active, char* l){
  GLfloat side_x1, side_x2;
  GLfloat side_y1, side_y2;
  glNormal3f(0,0,1);
  // lower left corner
  side_x1 = x1+0.11*base_view;
  side_x2 = x1+0.13*base_view;
  // uper right corner
  side_y1 = y1+0.09*view_bottom;
  side_y2 = (side_y1+0.02*base_view);
  if(active)
    glColor4f(0.08,0.08,0.08,0.8);
  else
    glColor4f(0.98,0.87,0.64,0.7);
  // radio solid mark
  glBegin(GL_POLYGON);
  glVertex3f(side_x1, side_y2, z);
  glVertex3f(side_x2, side_y2, z);
  glVertex3f(side_x2, side_y1, z);
  glVertex3f(side_x1, side_y1, z);
  glEnd();
  // radio box
  glColor4f(0.98,0.87,0.64,0.7);
  glBegin(GL_LINE_LOOP);
  glVertex3f(side_x1, side_y2, z); // left corner of the roof
  glVertex3f(side_x2, side_y2, z); // right corner of the roof
  glVertex3f(side_x2, side_y1, z); // bottom right corner of the house
  glVertex3f(side_x1, side_y1, z); // bottom left corner of the house
  glEnd();
  //
  // radio label
  glRasterPos3f(side_x2+0.01*base_view, side_y1/*+0.3*(side_y2-side_y1)*/,z);
  //gl_font(FL_COURIER,font_size_slider_label); // text font
  gl_draw(l, strlen(l));
}

// switch output widget
void Fl_Gl_Mol_View::draw_switch_output(GLfloat x1, GLfloat y1, GLfloat y2, GLfloat z, bool val, char* l){
  GLfloat side_x1;//, side_x2;
  GLfloat side_y1;//, side_y2;
  char buff[60];
  glNormal3f(0,0,1);
  // lower left corner
  side_x1 = x1;
  // uper right corner
  side_y1 = y1;
  //
  //if(active)
    glColor4f(0.78,0.68,0.0,1.0);
    //glColor4f(0.3,0.3,0.3,0.7);
  //else
    //glColor4f(0.98,0.87,0.64,0.7);
  // radio label
  //font_size_pie_label
  //gl_font(FL_COURIER,font_size_slider_label); // text font
  //gl_font(FL_COURIER,font_size_panel_label); // text font
  //gl_font(FL_COURIER,font_size_panel_label); // text font
  if(val)
    sprintf(buff,"%-15s ON",l);
  else
    sprintf(buff,"%-15s OFF",l);
  glRasterPos3f(side_x1, side_y1,z);
  gl_draw(buff, strlen(buff));
}

// float output widget
void Fl_Gl_Mol_View::widget_float_output(GLfloat x1, GLfloat y1, GLfloat y2, GLfloat z, float val, char* l){
  GLfloat side_x1;//, side_x2;
  GLfloat side_y1;//, side_y2;
  char buff[60];
  glNormal3f(0,0,1);
  // radio box
  // lower left corner
  side_x1 = x1;
  //side_x2 = x1+0.14*base_view;
  // uper right corner
  side_y1 = y1;
  //side_y2 = (side_y1+0.04*base_view);
  //
  //glColor4f(0.98,0.87,0.64,0.7);
  //glBegin(GL_LINE_LOOP);
  //glVertex3f(side_x2-1, side_y2, z); // left corner of the roof
  //glVertex3f(side_x2, side_y2, z); // right corner of the roof
  //glVertex3f(side_x2, side_y1, z); // bottom right corner of the house
  //glVertex3f(side_x2-1, side_y1, z); // bottom left corner of the house
  //glEnd();
  //
  //if(active)
  glColor4f(0.78,0.68,0.5,1.0);
    //glColor4f(0.3,0.3,0.3,0.7);
  //else
    //glColor4f(0.98,0.87,0.64,0.7);
  // radio label
  //font_size_pie_label
  //gl_font(FL_COURIER,font_size_slider_label); // text font
  //gl_font(FL_COURIER,font_size_panel_label); // text font
  //gl_font(FL_COURIER,font_size_panel_label); // text font
  sprintf(buff,"%-15s %f",l,val);
  glRasterPos3f(side_x1, side_y1,z);
  gl_draw(buff, strlen(buff));
  // radio mark
  //side_x1 = x1+0.11*base_view;
  //side_x2 = x1+0.13*base_view;
  //side_y1 = y1+0.09*view_bottom;
  //side_y2 = (side_y1+val*0.02*base_view);
  //glBegin(GL_POLYGON);
  //glVertex3f(side_x1, side_y2, z);
  //glVertex3f(side_x2, side_y2, z);
  //glVertex3f(side_x2, side_y1, z);
  //glVertex3f(side_x1, side_y1, z);
  //glEnd();
}

// float output widget
void Fl_Gl_Mol_View::widget_float_output_xy(GLfloat x1, GLfloat y1, GLfloat y2, GLfloat z, int x, int y, float val, char* l){
  GLfloat side_x1;//, side_x2;
  GLfloat side_y1;//, side_y2;
  char buff[512];
  glNormal3f(0,0,1);
  // radio box
  // lower left corner
  side_x1 = x1;
  //side_x2 = x1+0.14*base_view;
  // uper right corner
  side_y1 = y1;
  //side_y2 = (side_y1+0.04*base_view);
  //
  //glColor4f(0.98,0.87,0.64,0.7);
  //glBegin(GL_LINE_LOOP);
  //glVertex3f(side_x2-1, side_y2, z); // left corner of the roof
  //glVertex3f(side_x2, side_y2, z); // right corner of the roof
  //glVertex3f(side_x2, side_y1, z); // bottom right corner of the house
  //glVertex3f(side_x2-1, side_y1, z); // bottom left corner of the house
  //glEnd();
  //
  //if(active)
  glColor4f(0.78,0.68,0.5,1.0);
    //glColor4f(0.3,0.3,0.3,0.7);
  //else
    //glColor4f(0.98,0.87,0.64,0.7);
  // radio label
  //font_size_pie_label
  //gl_font(FL_COURIER,font_size_slider_label); // text font
  //gl_font(FL_COURIER,font_size_panel_label); // text font
  //gl_font(FL_COURIER,font_size_panel_label); // text font
  sprintf(buff,"%-15s %i %i %f",l,x,y,val);
  glRasterPos3f(side_x1, side_y1,z);
  gl_draw(buff, strlen(buff));
  // radio mark
  //side_x1 = x1+0.11*base_view;
  //side_x2 = x1+0.13*base_view;
  //side_y1 = y1+0.09*view_bottom;
  //side_y2 = (side_y1+val*0.02*base_view);
  //glBegin(GL_POLYGON);
  //glVertex3f(side_x1, side_y2, z);
  //glVertex3f(side_x2, side_y2, z);
  //glVertex3f(side_x2, side_y1, z);
  //glVertex3f(side_x1, side_y1, z);
  //glEnd();
}

void Fl_Gl_Mol_View::widget_text_output(GLfloat x1, GLfloat y1, GLfloat y2, GLfloat z, char* text, char* l){
  GLfloat side_x1;//, side_x2;
  GLfloat side_y1;//, side_y2;
  char buff[60];
  glNormal3f(0,0,1);
  // radio box
  // lower left corner
  side_x1 = x1;
  //side_x2 = x1+0.14*base_view;
  // uper right corner
  side_y1 = y1;
  //side_y2 = (side_y1+0.04*base_view);
  //
  //glColor4f(0.98,0.87,0.64,0.7);
  //glBegin(GL_LINE_LOOP);
  //glVertex3f(side_x2-1, side_y2, z); // left corner of the roof
  //glVertex3f(side_x2, side_y2, z); // right corner of the roof
  //glVertex3f(side_x2, side_y1, z); // bottom right corner of the house
  //glVertex3f(side_x2-1, side_y1, z); // bottom left corner of the house
  //glEnd();
  //
  //if(active)
    glColor4f(0.88,0.98,0.5,1.0);
    //glColor4f(0.3,0.3,0.3,0.7);
  //else
    //glColor4f(0.98,0.87,0.64,0.7);
  // radio label
  //font_size_pie_label
  //gl_font(FL_COURIER,font_size_slider_label); // text font
  //gl_font(FL_COURIER,font_size_panel_label); // text font
  //gl_font(FL_COURIER,font_size_panel_label); // text font
  sprintf(buff,"%s %s",l,text);
  glRasterPos3f(side_x1, side_y1,z);
  gl_draw(buff, strlen(buff));
  // radio mark
  //side_x1 = x1+0.11*base_view;
  //side_x2 = x1+0.13*base_view;
  //side_y1 = y1+0.09*view_bottom;
  //side_y2 = (side_y1+val*0.02*base_view);
  //glBegin(GL_POLYGON);
  //glVertex3f(side_x1, side_y2, z);
  //glVertex3f(side_x2, side_y2, z);
  //glVertex3f(side_x2, side_y1, z);
  //glVertex3f(side_x1, side_y1, z);
  //glEnd();
}

void Fl_Gl_Mol_View::widget_int_output(GLfloat x1, GLfloat y1, GLfloat y2, GLfloat z, int val, char* l){
  char buff[20];
  sprintf(buff,"%i",val);
  widget_text_output(x1, y1, y2, z, buff, l);
}

void Fl_Gl_Mol_View::widget_vector_output(GLfloat x1, GLfloat y1, GLfloat y2, GLfloat z, TVector<real> v, char* l){
  char buff[60];
  sprintf(buff,"%10.6f %10.6f %10.6f",v[0],v[1],v[2]);
  widget_text_output(x1, y1, y2, z, buff, l);
}

// pie menu interaction/handle
void Fl_Gl_Mol_View::draw_pie_menu(GLfloat cx, GLfloat cy, GLfloat z, GLfloat r, GLint num_segments){
  GLdouble winX, winY; //, winZ;
  //GLdouble posX=0, posY=0, posZ=0;
  GLdouble posZ=0;
  float delta_ang, x1, y1;
  float theta;
  char buff[10];
  delta_ang = 2.0f * 3.1415926f / float(num_segments); //get the angle step
  //
  glGetDoublev( GL_MODELVIEW_MATRIX, modelview );
  glGetDoublev( GL_PROJECTION_MATRIX, projection );
  glGetIntegerv( GL_VIEWPORT, viewport );
  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glNormal3f(0,0,1);
  //
  // Did you account for the inverted Y coordinate? Most window systems (Microsoft Windows, X Windows,
  // others?) usually return mouse coordinates to your program with Y=0 at the top of the window, while
  // OpenGL assumes Y=0 is at the bottom of the window. Assuming you're using a default viewport,
  // transform the Y value from window system coordinates to OpenGL coordinates as (windowHeight-y).
  if(is_menu_position){
    menu_pos_cx=cx;
    menu_pos_cy=cy;
    winX = menu_pos_cx;
    winY = viewport[3]-menu_pos_cy;
    //glReadPixels( winX, int(winY), 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ );
    //GLbyte color[4];
    //GLfloat depth;
    //GLdouble pX, pY, pZ;
    //GLuint index;
    //glReadPixels(winX, int(winY), 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, color);
    //glReadPixels(winX, int(winY), 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &depth);
    //glReadPixels(winX, int(winY), 1, 1, GL_STENCIL_INDEX, GL_UNSIGNED_INT, &index);
    //printf("Clicked on pixel %d, %d, color %02hhx%02hhx%02hhx%02hhx, depth %f, stencil index %u\n", winX, winY, color[0], color[1], color[2], color[3], depth, index);
    //printf("Clicked on pixel %d, %d, color %x %x %x %x, depth %f, stencil index %u\n", x, y, color[0], color[1], color[2], color[3], depth, index);
    //printf("Clicked on pixel %d, %d, depth %f, stencil index %u\n", x, y, depth, index);
    //
    //gluUnProject( winX, winY, (int)depth, modelview, projection, viewport, &pX, &pY, &pZ);
    //printf("Window coordinates %f, %f, %f\n", pX, pY, pZ);
    // How do I use gluProject() and gluUnProject()?
    // Both routines take a ModelView matrix, Projection matrix, and OpenGL Viewport as parameters.
    // gluProject() also takes an XYZ-object space coordinate. It returns the transformed XYZ window
    // (or device) coordinate equivalent.
    // gluUnProject() does the opposite. It takes an XYZ window coordinate and returns the
    // back-transformed XYZ object coordinate equivalent.
    // The concept of window space Z is often confusing. It's the depth buffer value expressed as a
    // GLdouble in the range 0.0 to 1.0. Assuming a default glDepthRange(), a window coordinate with
    // a Z value of 0.0 corresponds to an eye coordinate located on the zNear clipping plane.
    // Similarly, a window space Z value of 1.0 corresponds to an eye space coordinate located on the
    // zFar plane. You can obtain any window space Z value by reading the depth buffer with glReadPixels().
    //
    //gluUnProject( winX, winY, winZ, modelview, projection, viewport, &menu_pos_x, &menu_pos_y, &menu_pos_z);
#if defined (BUILD_FOR_MACOS)
    //std::cout<<"Mac OS X code here"<<std::endl;
    gluUnProject( winX, winY, 1, modelview, projection, viewport, &menu_pos_x, &menu_pos_y, &menu_pos_z);
#else
    gluUnProject( winX, winY, 1, modelview, projection, viewport, &menu_pos_x, &menu_pos_y, &menu_pos_z);
#endif
    //is_menu_position=false;
  }
  //

  if(is_atom_picked){
    // clicked point position
    winX = cx;
    winY = viewport[3]-cy;
    //glReadPixels( cx, int(winY), 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ );
    //gluProject( winX, winY, 800, modelview, projection, viewport, &tmpX, &tmpY, &tmpZ);
    //gluUnProject( winX, winY, winZ, modelview, projection, viewport, &click_pos_x, &click_pos_y, &posZ);

#if defined (BUILD_FOR_MACOS)
    //std::cout<<"Mac OS X code here"<<std::endl;
    gluUnProject( winX, winY, 1, modelview, projection, viewport, &click_pos_x, &click_pos_y, &posZ);
#else
    gluUnProject( winX, winY, 1, modelview, projection, viewport, &click_pos_x, &click_pos_y, &posZ);
#endif
    GLfloat dist=(click_pos_x-menu_pos_x)*(click_pos_x-menu_pos_x)+(click_pos_y-menu_pos_y)*(click_pos_y-menu_pos_y);
    if(dist > r){
      //std::cout<<" distance :"<<dist<<std::endl;
      //glColor3f(0.0,0.0,1.0);
      GLdouble delta_x = menu_pos_x-click_pos_x;
      GLdouble delta_y = menu_pos_y-click_pos_y;
      GLdouble ang = atan2(delta_y,delta_x)+ C_PI;
      //
      side_pos_x = r * cosf(ang);          //calculate the x component
      side_pos_y = r * sinf(ang);          //calculate the y component
      is_draw_line=true;
    }else{
      is_draw_line=false;
    }
  }
  // draw the connection line
  //if((is_menu_pie_picked || is_atom_picked) && is_draw_line){
  if(is_draw_line){
    glColor3f(0.0,1.0,0.0);
    glBegin(GL_LINES);
    glVertex3f(menu_pos_x+side_pos_x,menu_pos_y+side_pos_y,z);
    glVertex3f(click_pos_x,click_pos_y,z);
    glEnd();
  }
  // current atom point
  //if((is_menu_pie_picked || is_atom_picked) && is_draw_point){
  if(is_draw_point){
    glColor3f(0.0F,1.0F,0.0F);
    glBegin(GL_POLYGON);
    for(int ii = 0; ii < num_segments+1; ii++){
      theta = delta_ang * float(ii);        //get the current angle
      x1 = 0.05*r * cosf(theta);            //calculate the x component
      y1 = 0.05*r * sinf(theta);            //calculate the y component
      glVertex3f(x1 + click_pos_x, y1 + click_pos_y,z); //output vertex
    }
    glEnd();
  }
  // May 32, 2014: Becuase a bug on the Mac this lines are required
  //for(int i=0; i<6; i++)
    //legends[i] = l0[i];
  // May 32, 2014: Becuase a bug on the Mac this lines are required
  if(!is_menu_pie_picked){
    switch(u_active_menu){
      case MAIN_MENU:
#ifdef _SHOW_INFO_
        std::cout<<" backgroud legends"<<std::endl;
#endif
        label = " "; // blank label
        for(int i=0; i<6; i++)
          legends[i] = l0[i];
        break;
      case CONTROL_MENU:
#ifdef _SHOW_INFO_
        std::cout<<" menu legends"<<std::endl;
#endif
        label = " "; // blank label
        for(int i=0; i<6; i++)
          legends[i] = l1[i];
      break;
      case ATOM_MENU:
//#ifdef _SHOW_INFO_
//        std::cout<<" atom legends"<<std::endl;
//#endif
        // set the label to the picked atom
        sprintf(buff,"%i-%s",__highlight_atom+1,v_atom_symbols[__highlight_atom].c_str());
        label=buff;
        for(int i=0; i<6; i++)
          legends[i] = l2[i];
      break;
    }
  }
  draw_pie(menu_pos_cx,menu_pos_cy,z,legends,6,r,num_segments);
  is_atom_picked=false;
  is_menu_picked=false;
  is_menu_pie_picked=false;
  is_menu_position=false;
  is_background_picked=false;
  //u_active_menu=NOT_MENU;
}

void Fl_Gl_Mol_View::draw_pie_submenu(GLfloat z, GLfloat r, GLint num_segments){
  // blank label
  switch(u_active_menu){
    case MAIN_MENU:
      switch(u_menu_index){
        case MAIN_MODE_SUBMENU:
          set_pie_labels(sl_mode,"mode");
          break;
        case MAIN_VIEW_SUBMENU:
          set_pie_labels(sl_view,"view");
          break;
        case MAIN_SHOW_SUBMENU:
          set_pie_labels(sl_show,"show");
          break;
        case MAIN_TOOLS_SUBMENU:
          set_pie_labels(sl_tools,"show");
          break;
        default:
          set_pie_labels(lblank,"");
          break;
      }
    break;
    case CONTROL_MENU:
      switch(u_menu_index){
        default:
          set_pie_labels(lblank,"");
        break;
      }
    break;
    case ATOM_MENU:
      switch(u_menu_index){
        //case FRAGMENT_AXIS_BUTTON:
          //set_pie_labels(l_frag,"frag.");
          //break;
        case ATOM_AXIS_SUBMENU:
          set_pie_labels(l_axis,"axis");
          break;
        case ATOM_FRAGMENT_SUBMENU:
          set_pie_labels(l_frag,"Frags");
          break;
        default:
          set_pie_labels(lblank,"");
        break;
      }
    break;
  }
  draw_sub_pie(submenu_pos_cx,submenu_pos_cy,z,legends,6,r,num_segments);
}

void Fl_Gl_Mol_View::set_pie_labels(const std::string s[], std::string l){
  sub_label = l;
  for(int i=0; i<6; i++)
    sub_legends[i] = s[i];
}

// pie menu interaction/handle
void Fl_Gl_Mol_View::draw_sub_pie(GLfloat cx, GLfloat cy, GLfloat z, std::string l[], GLint nl, GLfloat r, GLint n){
  float delta_ang, x1, x2, y1, y2;
  float theta;
  // draw the buttons only
  if(render_mode==MODE_SELECT){
    glNormal3f(0,0,1);
    glColor4f(0.3,0.2,0.6,0.5);
    delta_ang = 2.0f * 3.1415926f / float(nl);//get the current angle
    for(int ii = 0; ii < 6; ii++){
      //glPushName(ii+7); // the first 7 names are reseverd for the main menu
      ui_rgb color;
      color = index_palette.get_index(ii+6);
      glColor3ub(color.r,color.g,color.b);
      glBegin(GL_QUADS);
      theta = delta_ang * float(ii);                   //get the current angle
      x1 = 0.2*r * cosf(theta);                        //calculate the x component
      y1 = 0.2*r * sinf(theta);                        //calculate the y component
      x2 = 1.2*r * cosf(theta);                        //calculate the x component
      y2 = 1.2*r * sinf(theta);                        //calculate the y component
      glVertex3f(x1 + submenu_pos_x, y1 + submenu_pos_y, z+4);//output vertex
      glVertex3f(x2 + submenu_pos_x, y2 + submenu_pos_y, z+4);//output vertex
      theta = delta_ang * float(ii+1);                 //get the current angle
      x1 = 0.2*r * cosf(theta);                        //calculate the x component
      y1 = 0.2*r * sinf(theta);                        //calculate the y component
      x2 = 1.2*r * cosf(theta);                        //calculate the x component
      y2 = 1.2*r * sinf(theta);                        //calculate the y component
      glVertex3f(x2 + submenu_pos_x, y2 + submenu_pos_y, z+4);//output vertex
      glVertex3f(x1 + submenu_pos_x, y1 + submenu_pos_y, z+4);//output vertex
      glEnd();
      //glPopName();
    }
  }else{
    // draw the pie disk
    draw_pie_disk(submenu_pos_x,submenu_pos_y,z+2,r,n);
    // draw the labels
    draw_pie_labels(submenu_pos_cx,submenu_pos_cy,z+2,r,sub_legends,sub_label,nl);
  }
  //is_menu_pie_picked=false;
}

// pie menu widget handle
void Fl_Gl_Mol_View::draw_pie(GLfloat cx, GLfloat cy, GLfloat z, std::string l[], GLint nl, GLfloat r, GLint n){
  float delta_ang, x1, x2, y1, y2;
  float theta;
  // draw the buttons
  if(render_mode==MODE_SELECT){
    glNormal3f(0,0,1);
    //glColor4f(0.5,0.3,0.1,0.5);
    delta_ang = 2.0f * 3.1415926f / float(nl);//get the current angle
    for(int ii = 0; ii < 6; ii++){
      ui_rgb color;
      color = index_palette.get_index(ii);
      glColor3ub(color.r,color.g,color.b);
      //glPushName(ii+1); // the first 100 names are reseverd for the menues
      glBegin(GL_QUADS);
      theta = delta_ang * float(ii);                   //get the current angle
      x1 = 0.2*r * cosf(theta);                        //calculate the x component
      y1 = 0.2*r * sinf(theta);                        //calculate the y component
      x2 = 1.2*r * cosf(theta);                        //calculate the x component
      y2 = 1.2*r * sinf(theta);                        //calculate the y component
      glVertex3f(x1 + menu_pos_x, y1 + menu_pos_y, z+2);//output vertex
      glVertex3f(x2 + menu_pos_x, y2 + menu_pos_y, z+2);//output vertex
      theta = delta_ang * float(ii+1);                 //get the current angle
      x1 = 0.2*r * cosf(theta);                        //calculate the x component
      y1 = 0.2*r * sinf(theta);                        //calculate the y component
      x2 = 1.2*r * cosf(theta);                        //calculate the x component
      y2 = 1.2*r * sinf(theta);                        //calculate the y component
      glVertex3f(x2 + menu_pos_x, y2 + menu_pos_y, z+2);//output vertex
      glVertex3f(x1 + menu_pos_x, y1 + menu_pos_y, z+2);//output vertex
      glEnd();
      //glPopName();
    }
  }else{
    // draw the pie disk
    draw_pie_disk(menu_pos_x,menu_pos_y,z,r,n);
    // draw the labels
    draw_pie_labels(menu_pos_cx,menu_pos_cy,z+1,r,legends,label,nl);
  }
}

// pie disk
void Fl_Gl_Mol_View::draw_pie_disk(GLfloat x, GLfloat y, GLfloat z, GLfloat r, GLint n){
  GLdouble winX, winY, winZ;
  float delta_ang, x1, x2, y1, y2;
  float theta;
  delta_ang = 2.0f * 3.1415926f / float(n);//get the current angle
  glNormal3f(0,0,1);
  glColor4f(0.0,0.0,0.0,0.9);
  glBegin(GL_QUAD_STRIP);
  //glBegin(GL_POLYGON);
  for(int ii=0; ii < n+1; ii++){
    theta = delta_ang * float(ii); //get the current angle
    x1 = 0.2*r * cosf(theta);      //calculate the x component
    y1 = 0.2*r * sinf(theta);      //calculate the y component
    x2 = r * cosf(theta);          //calculate the x component
    y2 = r * sinf(theta);          //calculate the y component
    glVertex3f(x1 + x, y1 + y, z); //output vertex
    glVertex3f(x2 + x, y2 + y, z); //output vertex
  }
  glEnd();
  //glDisable (GL_BLEND);
  // draw the iner circle
  //glColor3f(0.3,0.3,0.3);
  glColor3f(0.8,0.8,0.8);
  glBegin(GL_LINE_LOOP);
  for(int ii = 0; ii < n+1; ii++){
    theta = delta_ang * float(ii); //get the current angle
    x1 = 0.2*r * cosf(theta);      //calculate the x component
    y1 = 0.2*r * sinf(theta);      //calculate the y component
    glVertex3f(x1 + x, y1 + y, z); //output vertex
  }
  glEnd();
  // draw the outer circle
  glBegin(GL_LINE_LOOP);
  for(int ii = 0; ii < n+1; ii++){
    theta = delta_ang * float(ii); //get the current angle
    x2 = r * cosf(theta);          //calculate the x component
    y2 = r * sinf(theta);          //calculate the y component
    glVertex3f(x2 + x, y2 + y, z); //output vertex
  }
  glEnd();
  // menu divisions
  //glColor3f(0.7,0.7,0.7);
  glColor3f(0.8,0.8,0.8);
  delta_ang = 2.0f * 3.1415926f / float(6);          //get the angle step
  glBegin(GL_LINES);
  for(int ii = 0; ii < 6; ii++){
    theta = delta_ang * float(ii);                   //get the current angle
    x1 = 0.3*r * cosf(theta);                        //calculate the x component
    y1 = 0.3*r * sinf(theta);                        //calculate the y component
    x2 = 0.9*r * cosf(theta);                        //calculate the x component
    y2 = 0.9*r * sinf(theta);                        //calculate the y component
    glVertex3f(x1 + x, y1 + y, z);//output vertex
    glVertex3f(x2 + x, y2 + y, z);//output vertex
  }
  glEnd();
  /*
  if(render_mode==MODE_SELECT){
    //glNormal3f(0,0,1);
    glColor4f(0.0,0.0,0.0,0.0);
    //delta_ang = 2.0f * 3.1415926f / float(n);//get the current angle
    for(int ii = 0; ii < 6; ii++){
      glPushName(ii+1); // the first 100 names are reseverd for the menues
      glBegin(GL_QUADS);
      theta = delta_ang * float(ii);                   //get the current angle
      x1 = 0.2*r * cosf(theta);                        //calculate the x component
      y1 = 0.2*r * sinf(theta);                        //calculate the y component
      x2 = 1.2*r * cosf(theta);                        //calculate the x component
      y2 = 1.2*r * sinf(theta);                        //calculate the y component
      glVertex3f(x1 + x, y1 + y, z+2);//output vertex
      glVertex3f(x2 + x, y2 + y, z+2);//output vertex
      theta = delta_ang * float(ii+1);                 //get the current angle
      x1 = 0.2*r * cosf(theta);                        //calculate the x component
      y1 = 0.2*r * sinf(theta);                        //calculate the y component
      x2 = 1.2*r * cosf(theta);                        //calculate the x component
      y2 = 1.2*r * sinf(theta);                        //calculate the y component
      glVertex3f(x2 + x, y2 + y, z+2);//output vertex
      glVertex3f(x1 + x, y1 + y, z+2);//output vertex
      glEnd();
      glPopName();
    }
  }
  */
  // compute the coordinates for the pie subumenu
  if(is_menu_pie_picked){
    //
    glGetDoublev( GL_MODELVIEW_MATRIX, modelview );
    glGetDoublev( GL_PROJECTION_MATRIX, projection );
    glGetIntegerv( GL_VIEWPORT, viewport );
    //
    delta_ang = 2.0f * 3.1415926f / float(360);
    theta = delta_ang * float(u_menu_index*60+30); //get the current angle
    x1 = 0.5*r * cosf(theta);                             //calculate the x component
    y1 = 0.5*r * sinf(theta);                             //calculate the y component
    submenu_pos_x=x+x1;
    submenu_pos_y=y+y1;
    winX = submenu_pos_x;
    //winY = viewport[3]-submenu_pos_cy;
    winY = submenu_pos_y;
    //glReadPixels( int(winX), int(winY), 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ );
#if defined (BUILD_FOR_MACOS)
    gluProject( winX, winY, z+2, modelview, projection, viewport, &submenu_pos_cx, &submenu_pos_cy, &winZ);
    //std::cout<<"Mac OS X code here"<<std::endl;
#else
    gluProject( winX, winY, z+2, modelview, projection, viewport, &submenu_pos_cx, &submenu_pos_cy, &winZ);
#endif
    submenu_pos_cy = viewport[3]-submenu_pos_cy;
  }
}

void Fl_Gl_Mol_View::draw_pie_labels(GLfloat cx, GLfloat cy, GLfloat z, GLfloat r, std::string l[6], std::string m, GLint nl){
  GLdouble winX, winY, winZ=791;
  GLdouble posX, posY, posZ;
  float delta_ang, x1, y1;
  float theta, str_width;
  char buff[10];
  int str_length;
  winX = cx;
  winY = viewport[3]-cy;
  delta_ang = 2.0f * 3.1415926f / float(360);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glNormal3f(0,0,1);
  glColor3f(1.0F,1.0F,1.0F); // text color
  //gl_font(FL_COURIER,14); // text font
  //gl_font(FL_COURIER,14); // text font
  //glColor4f(0.0,1.0,0.0,1.0); // text color
  //gl_font(FL_ZAPF_DINGBATS,12); // text font
#ifdef PLATFORM_MAC
  gl_font(FL_HELVETICA,font_size_pie_label); // text font
#endif
  for(int i=0; i<nl; i++){
    //theta = delta_ang * float(90+i*60); //get the current angle
    theta = delta_ang * float(i*60+30); //get the current angle
    x1 = 0.7*r * cosf(theta);           //calculate the x component
    y1 = 0.7*r * sinf(theta);           //calculate the y component
    sprintf(buff,"%s",l[i].c_str());
    str_length = strlen(buff);
    str_width = gl_width(buff);
    winX = cx-0.5*str_width;
#if defined (BUILD_FOR_MACOS)
    gluUnProject( winX, winY, winZ, modelview, projection, viewport, &posX, &posY, &posZ);
    //std::cout<<"Mac OS X code here"<<std::endl;
#else
    gluUnProject( winX, winY, winZ, modelview, projection, viewport, &posX, &posY, &posZ);
#endif
    glRasterPos3f(posX+x1,posY+y1,z);
    gl_draw(buff, str_length);
  }
  // set the label
  //if(is_atom_picked){
  //if(render_mode!=MODE_SELECT){
    glColor3f(1.0F,1.0F,0.0F); // text color
    //gl_font(FL_COURIER,12);  // text font
    sprintf(buff,"%s",m.c_str());
    str_length = strlen(buff);
    str_width = gl_width(buff);
    //winX = menu_pos_cx-0.5*str_width;
    winX = cx-0.5*str_width;
#if defined (BUILD_FOR_MACOS)
    gluUnProject( winX, winY, winZ, modelview, projection, viewport, &posX, &posY, &posZ);
    //std::cout<<"Mac OS X code here"<<std::endl;
#else
    gluUnProject( winX, winY, winZ, modelview, projection, viewport, &posX, &posY, &posZ);
#endif
    glRasterPos3f(posX,posY+0.4*r,z);
    gl_draw(buff, str_length);
  //}
}

//////////////////////////////UTILS///////////////////////////////
/*
// Tue Feb 26 12:25:32 MST 2013
// Migrated from fl_gl_mol_view.h
void Fl_Gl_Mol_View::update_data(void){
  map_update_active_fragment();
  if(update_coordinates){
    update_atomic_coordinates(get_view_cartesian());
    set_axis_position(get_view_centered_position_cartesian());
    set_axis_precession(get_view_axis_precession());
    set_axis_tilt(get_view_axis_tilt());
    set_backbone_precession(get_view_backbone_precession());
    set_backbone_tilt(get_view_backbone_tilt());
    update_coordinates=false;
  }
}*/

void Fl_Gl_Mol_View::eval_system_properties(void){
#ifdef _GLMOL_DEBUG_MESSAGES_
    std::cout<<" GLMOL: eval_system_properties (0)"<<std::endl;
#endif
  if(is_update_atomic_properties){
    ////eval_initial_properties();
    //eval_atomic_positions(); // <---------
    is_update_atomic_properties=false;
    if(is_draw_tools_){
      eval_tool_parameters();
    }
#ifdef _GLMOL_DEBUG_MESSAGES_
    std::cout<<" GLMOL: update atomic positions"<<std::endl;
#endif
  }
#ifdef _GLMOL_DEBUG_MESSAGES_
    std::cout<<" GLMOL: eval_system_properties (1)"<<std::endl;
#endif
  if(is_update_radius){
    //eval_atomic_radius(); // <-----------
    is_update_radius=false;
#ifdef _GLMOL_DEBUG_MESSAGES_
    std::cout<<" GLMOL: update radius"<<std::endl;
#endif
  }
#ifdef _GLMOL_DEBUG_MESSAGES_
    std::cout<<" GLMOL: eval_system_properties (2)"<<std::endl;
#endif
  // Sun Feb 24 16:11:36 MST 2013
  // it evaluates the bonds only once then
  // use efficient update the bonds
  if(is_draw_bonds_ && is_update_bonds){
    if(is_eval_bonds){
      eval_atomic_bonds();  // <----------
      is_eval_bonds=false;
    }else{
      update_atomic_bonds();  // <----------
    }
    create_cylinder_dl();  // <----------
    is_update_bonds = false;
#ifdef _SHOW_DEBUG_BONDS_
    std::cout<<" GLMOL: update bonds"<<std::endl;
#endif
  }
#ifdef _GLMOL_DEBUG_MESSAGES_
    std::cout<<" GLMOL: eval_system_properties (3)"<<std::endl;
#endif
  if(is_update_mask_rcolor){
    eval_mask_rcolor();
    is_update_mask_rcolor=false;
#ifdef _GLMOL_DEBUG_MESSAGES_
    std::cout<<" GLMOL: update mask rcolor"<<std::endl;
#endif
  }
#ifdef _GLMOL_DEBUG_MESSAGES_
    std::cout<<" GLMOL: eval_system_properties (4)"<<std::endl;
#endif
  //if(is_update_mask_rcolor || is_update_radius || is_update_atomic_properties){
    create_sphere_dl();    // <----------
  //}
  //
#ifdef _GLMOL_DEBUG_MESSAGES_
    std::cout<<" GLMOL: eval_system_properties (5)"<<std::endl;
#endif
}

void Fl_Gl_Mol_View::set_highlight_atom(int i){
#ifdef _GLMOL_DEBUG_MESSAGES_
  std::cout<<" GLMOL: Highlight atom: "<<i<<std::endl;
#endif
  __last_highlight_atom=__highlight_atom;
  __highlight_atom = i;
  update_highlight_atom=true;
  is_update_mask_rcolor=true;
}

void Fl_Gl_Mol_View::set_selected_atom(uint u){
  is_unselected_atom=false;
  //uint pos;
  for(uint i=0; i<u_selected_index; i++){
      if(v_selected_atoms[i]==u){
        is_unselected_atom=true;
        // unselect if the atom was selected
        u_unselected_atom=v_selected_atoms[i];
        for(uint j=i; j<u_selected_index-1; j++){
          v_selected_atoms[j]=v_selected_atoms[j+1];
        }
        u_selected_index--;
        //update_normal_color=true;
        update_selected_atoms=true;
      }
  }
  if(u_selected_index<4){
    // add a new selected atom
    if(!is_unselected_atom){
      v_selected_atoms[u_selected_index]=u;
      //std::cout<<" selected: "<<u<<std::endl;
      u_selected_index++;
      if(u_selected_index==2){
        v_distance1=supercell.get_vector_diff(v_selected_atoms[0],v_selected_atoms[1]);
        r_distance1=v_distance1.magnitude();  //get_distance(v_selected_atoms[0],v_selected_atoms[1]);
        //std::cout<<" distance = "<<r_distance1<<std::endl;
      }else if(u_selected_index==3){
        v_distance2=supercell.get_vector_diff(v_selected_atoms[1],v_selected_atoms[2]);
        r_distance2=v_distance2.magnitude();  //get_distance(v_selected_atoms[1],v_selected_atoms[2]);
        r_angle1=supercell.get_angle(v_selected_atoms[0],v_selected_atoms[1],v_selected_atoms[2]);
        //std::cout<<" angle = "<<r_angle1<<std::endl;
      }else if(u_selected_index==4){
        v_distance3=supercell.get_vector_diff(v_selected_atoms[2],v_selected_atoms[3]);
        r_distance3=v_distance3.magnitude();  //get_distance(v_selected_atoms[2],v_selected_atoms[3]);
        r_angle2=supercell.get_angle(v_selected_atoms[1],v_selected_atoms[2],v_selected_atoms[3]);
        r_dihedral=supercell.get_dihedral(v_selected_atoms[0],v_selected_atoms[1],v_selected_atoms[2],v_selected_atoms[3]);
        //std::cout<<" dihedral = "<<r_dihedral<<std::endl;
      }
      update_selected_atoms=true;
    }
  }
}

void Fl_Gl_Mol_View::eval_tool_parameters(void){
  for(uint uindex=0; uindex<=u_selected_index; uindex++){
    if(uindex==2){
        v_distance1=supercell.get_vector_diff(v_selected_atoms[0],v_selected_atoms[1]);
        r_distance1=v_distance1.magnitude();  //get_distance(v_selected_atoms[0],v_selected_atoms[1]);
        //std::cout<<" distance = "<<r_distance1<<std::endl;
    }else if(uindex==3){
        v_distance2=supercell.get_vector_diff(v_selected_atoms[1],v_selected_atoms[2]);
        r_distance2=v_distance2.magnitude();  //get_distance(v_selected_atoms[1],v_selected_atoms[2]);
        r_angle1=supercell.get_angle(v_selected_atoms[0],v_selected_atoms[1],v_selected_atoms[2]);
        //std::cout<<" angle = "<<r_angle1<<std::endl;
    }else if(uindex==4){
        v_distance3=supercell.get_vector_diff(v_selected_atoms[2],v_selected_atoms[3]);
        r_distance3=v_distance3.magnitude();  //get_distance(v_selected_atoms[2],v_selected_atoms[3]);
        r_angle2=supercell.get_angle(v_selected_atoms[1],v_selected_atoms[2],v_selected_atoms[3]);
        r_dihedral=supercell.get_dihedral(v_selected_atoms[0],v_selected_atoms[1],v_selected_atoms[2],v_selected_atoms[3]);
        //std::cout<<" dihedral = "<<r_dihedral<<std::endl;
    }
  }
  //update_selected_atoms=true;
}

uint Fl_Gl_Mol_View::get_highlight_atom(void){
  return __highlight_atom;
}

void Fl_Gl_Mol_View::set_update_active_fragment(void){
  //__fragment_active=u;
  if(is_highlight_fragment_){
    update_normal_color=true;
    is_update_mask_rcolor=true;
  }
}

void Fl_Gl_Mol_View::set_highlight_atom_a(int i){
  __highlight_atom_a = i;
}

void Fl_Gl_Mol_View::set_highlight_atom_b(int i){
  __highlight_atom_b = i;
}

void Fl_Gl_Mol_View::set_select_begin(int i){
  __select_begin=i;
}

void Fl_Gl_Mol_View::set_select_end(int i){
  __select_end=i;
}

/*
void Fl_Gl_Mol_View::set_fragment_total(uint u){
  __fragment_total=u;
}*/

void Fl_Gl_Mol_View::set_bond_brightness(real f){
  f_bond_brightness = f;
  update_bonds_color=true;
}

void Fl_Gl_Mol_View::set_background_brightness(real f){
  f_background_brightness = f;
}

void Fl_Gl_Mol_View::set_highlihght_brightness(real f){
  if(f >= 0.0 && f<=f_highlight_brightness_max)
    f_highlight_brightness = f;
  else if(f>1.0)
    f_highlight_brightness = f_highlight_brightness_max;
  else
    f_highlight_brightness = 0.0;
}

void Fl_Gl_Mol_View::set_atom_brightness(real f){
  f_atom_brightness = f;
  if(f >= 0.0 && f<=f_atom_brightness_max)
    f_atom_brightness = f;
  else if(f>1.0)
    f_atom_brightness = f_atom_brightness_max;
  else
    f_atom_brightness = 0.0;
  update_normal_color=true;
  is_update_mask_rcolor=true;
}

void Fl_Gl_Mol_View::set_select_brightness(real f){
  real f_select_brightness_max = 1.5;
  if(f >= 0.0 && f<=f_select_brightness_max)
    f_select_brightness = f;
  else if(f>f_select_brightness_max)
    f_select_brightness = f_select_brightness_max;
  else
    f_select_brightness = 0.0;
  update_highlight_atom=true;
  is_update_mask_rcolor=true;
}

void Fl_Gl_Mol_View::set_atom_radius_scale(real f){
  if(f >= 0.0 && f<=1.0)
    f_atom_radius_scale = f;
  else if(f>1.0)
    f_atom_radius_scale = 1.0;
  else
    f_atom_radius_scale = 0.0;
  //std::cout<<" f_atom_radius_scale = "<<f_atom_radius_scale<<std::endl;
  is_update_radius=true;
}

void Fl_Gl_Mol_View::set_bond_radius_scale(real f){
  if(f >= 0.0 && f<=1.0)
    f_bond_radius_scale = f;
  else if(f>1.0)
    f_bond_radius_scale = 1.0; // maximun bound
  else
    f_bond_radius_scale = 0.0;
  create_cylinder_dl();  // <----------
  //is_update_bonds=true;
}

void Fl_Gl_Mol_View::is_graphics(bool b){
  is_graphics_on=b;
  redraw();
}

void Fl_Gl_Mol_View::is_highlight_atom(bool b){
  is_highlight_atom_ = b;
  //if(is_highlight_atom_)
    //update_highlight_atom=true;
  //else
    //update_normal_color=true;
  //is_update_mask_rcolor=true;
}

void Fl_Gl_Mol_View::is_highlight_fragment(bool b){
  is_highlight_fragment_ = b;
  if(is_highlight_fragment_){
    update_highlight_fragment=true;
  }
  update_bonds_color=true;
  update_normal_color=true;
  is_update_mask_rcolor=true;
}

void Fl_Gl_Mol_View::set_pcb(bool b){
  is_pbc = b;
}

void Fl_Gl_Mol_View::is_mask_atoms(bool b){
  is_dark_mask_ = b;
}

void Fl_Gl_Mol_View::is_draw_bbox(bool b){
  is_draw_bbox_ = b;
  set_active_radio(3,is_draw_bbox_);
}

void Fl_Gl_Mol_View::is_draw_world_axes(bool b){
  is_draw_world_axes_ = b;
  set_active_radio(0,is_draw_world_axes_);
}

void Fl_Gl_Mol_View::is_draw_molecular_axes(bool b){
  is_draw_molecular_axes_ = b;
}

void Fl_Gl_Mol_View::is_draw_molecular_axis(bool b){
  is_draw_molecular_axis_ = b;
  set_active_radio(7,is_draw_molecular_axis_);
}

void Fl_Gl_Mol_View::is_draw_bonds(bool b){
  is_draw_bonds_ = b;
  if(is_draw_bonds_){
    update_bonds_color=true;
    is_update_mask_rcolor=true;
    if(f_atom_radius_scale > 0.40)
      set_atom_radius_scale(0.25);
    set_active_slider(1);
  }else{
    set_atom_radius_scale(0.50);
    set_active_slider(0);
  }
  set_active_radio(1,is_draw_bonds_);
}

void Fl_Gl_Mol_View::is_draw_labels(bool b){
  is_draw_labels_ = b;
  set_active_radio(6,is_draw_labels_);
}

void Fl_Gl_Mol_View::is_draw_symbols(bool b){
  is_draw_symbols_ = b;
  set_active_radio(2,is_draw_symbols_);
}

void Fl_Gl_Mol_View::is_draw_tools(bool b){
  is_draw_tools_ = b;
  if(!is_draw_tools_){
    u_selected_index=0;
    r_distance1=0.0;
    r_distance2=0.0;
    r_distance3=0.0;
    r_angle1=0.0;
    r_angle2=0.0;
    r_dihedral=0.0;
    update_normal_color=true;
    is_update_mask_rcolor=true;
  }
  set_active_radio(5,is_draw_tools_);
  update_normal_color=true;
}

void Fl_Gl_Mol_View::is_draw_numbers(bool b){
  is_draw_numbers_ = !is_draw_numbers_;
  set_active_radio(4,is_draw_numbers_);
}

void Fl_Gl_Mol_View::set_lock_controls(bool b){
  is_lock_controls = b;
  set_active_radio(8,is_lock_controls);
}

void Fl_Gl_Mol_View::set_atomic_labels(const TVector<std::string>& v){
  v_atom_labels = v;
}

void Fl_Gl_Mol_View::set_atomic_symbols(const TVector<std::string>& v){
  v_atom_symbols = v;
}

void Fl_Gl_Mol_View::set_atomic_symbol_table(const TVector<std::string>& v){
  v_atomic_symbol_table_gl = v;
}

void Fl_Gl_Mol_View::set_atomic_numbers(const TVector<uint>& v){
  v_atom_numbers = v;
}

void Fl_Gl_Mol_View::set_atom_table(const TVector<uint>& v){
  v_atom_table = v;
}

void Fl_Gl_Mol_View::set_atomic_number_table(const TVector<uint>& v){
  v_atomic_number_table_gl = v;
  r_cut_radius=0;
  //std::cout<<" v_atomic_number_table_gl = "<<v_atomic_number_table_gl;
  for(uint i=0; i<v_atomic_number_table_gl.size(); i++){
    r_cut_radius=maxi(r_cut_radius,atom_rrgb[v_atomic_number_table_gl[i]][0]);
  }
  r_cut_radius*=2.0;
  r_cut_radius_2 = r_cut_radius*r_cut_radius;
#ifdef _SHOW_DEBUG_LINKED_
  std::cout<<" Cut Radius = "<<r_cut_radius<<std::endl;
#endif
}

/*
void Fl_Gl_Mol_View::set_fragment_table(const TVector<uint>& v){
  v_fragment_table_gl = v;
}*/

void Fl_Gl_Mol_View::set_background_color(real r,real g,real b){
  bgred=r; bggreen=g; bgblue=b;
}

void Fl_Gl_Mol_View::set_foreground_color(real r,real g,real b){
  fgred=r; fggreen=g; fgblue=b;
}

void Fl_Gl_Mol_View::clear_all(void){
  //clear_cl3d();
  //clear_cn3d();
}

void Fl_Gl_Mol_View::clear_scene(void){
  is_draw_pie_menu=false;
  is_draw_line=false;
  is_draw_point=false;
  x_shift=0;
  y_shift=0;
  z_shift=0;
  u_active_menu=NOT_MENU;
  if(!is_lock_controls)
    is_draw_controls=false;
  is_atom_picked=false;
  is_lock_dragging=false;
  is_draw_pie_submenu=false;
  is_draw_line=false;
  is_draw_point=false;
  if(is_draw_tools_){
    u_selected_index=0;
    r_distance1=0.0;
    r_distance2=0.0;
    r_distance3=0.0;
    r_angle1=0.0;
    r_angle2=0.0;
    r_dihedral=0.0;
    //update_normal_color=true;
    //update_bonds_color=true;
    //is_update_mask_rcolor=true;
  }
  update_normal_color=true;
  is_update_mask_rcolor=true;
}

void Fl_Gl_Mol_View::set_view_xy_front(void){
  set_view(0,0,0);
}

void Fl_Gl_Mol_View::set_view_yz_front(void){
  set_view(-90,0,-90);
}

void Fl_Gl_Mol_View::set_view_zx_front(void){
  set_view(0,90,90);
}

void Fl_Gl_Mol_View::set_view_xy_back(void){
  set_view(0,180,180);
}

void Fl_Gl_Mol_View::set_view_yz_back(void){
  set_view(-90,0,90);
}

void Fl_Gl_Mol_View::set_view_zx_back(void){
  set_view(0,90,-90);
}

// images functions //////////////////////////////////////////////////////////////

/*
   Write the current view to a file
   The multiple fputc()s can be replaced with
      fwrite(image,width*height*3,1,fptr);
   If the memory pixel order is the same as the destination file format.
*/
int Fl_Gl_Mol_View::WindowDump(void){
   int i,j;
   FILE *fptr;
   static int counter = 0; /* This supports animation sequences */
   char fname[32];
   unsigned char *image;
   unsigned int width =  w();
   unsigned int height = h();
   //unsigned int pos_x =  x();
   //unsigned int pos_y =  y();
#ifdef _SHOW_INFO_
   std::cout<<" w = "<<(uint)width<<"  h = "<<(uint)height<<std::endl;
#endif

   /* Allocate our buffer for the image */
   if ((image = (unsigned char*) malloc(3*width*height*sizeof(char))) == NULL) {
   //if ((image = malloc(3*sizeof(char))) == NULL) {
      fprintf(stderr,"Failed to allocate memory for image\n");
      return(false);
   }
   glPixelStorei(GL_PACK_ALIGNMENT,1);
   // Open the file
   //if (stereo)
      sprintf(fname,"image_%04d.ppm",counter);
   //else
      //sprintf(fname,"C_%04d.raw",counter);
   if ((fptr = fopen(fname,"w")) == NULL) {
      fprintf(stderr,"Failed to open file for window dump\n");
      return(false);
   }
   // Copy the image into our buffer
   glReadBuffer(GL_BACK_LEFT);
   glReadPixels(0,0,width,height,GL_RGB,GL_UNSIGNED_BYTE,image);

   // Write the raw file
   fprintf(fptr,"P6\n%d %d\n255\n",width,height); //for ppm
   for (j=(int)height-1;j>=0;j--) {
      for (i=0;i<(int)width;i++) {
         fputc(image[3*j*width+3*i+0],fptr);
         fputc(image[3*j*width+3*i+1],fptr);
         fputc(image[3*j*width+3*i+2],fptr);
      }
   }
   fclose(fptr);

   //
   // Clean up
   counter++;
   free(image);
   //return(TRUE);
  return 1;
}


// END

