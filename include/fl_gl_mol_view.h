//========================================================================
// FILE - Fl_Atom_view.cxx                                              //
// For the Fast Light Tool Kit (FLTK) - www.fltk.org                    //
//========================================================================
//                                                                      //
// OpenGL 3D visualization widget to display atomic structures.         //
//                                                                      //
// Copyright 2002-2015 by Edmanuel Torres                               //
// email: eetorres@gmail.com                                            //
//                                                                      //
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

#ifndef _FL_GLMOL_VIEW_H_
#define _FL_GLMOL_VIEW_H_

#include <stdlib.h>
#include <stdio.h>

#include <atom_color.h>
#include <atom_symbol.h>
#include <atom_name.h>
#include <msmvtl/const.h>
#include <msmvtl/linalg.h>
#include <assert.h>

#include <FL/Fl.H>


#define HAVE_GL 1
#if HAVE_GL
    #include <fl_gl_atom.h>
//    #include <stroke.h>
#else
    #include <FL/Fl_Box.H>
#endif // HAVE_GL

#define MODE_RENDER    1
#define MODE_SELECT    2
#define MODE_MENU      3

#define BUFSIZE     1024

// Main menues
#define MAIN_MENU    0
#define CONTROL_MENU 1
#define ATOM_MENU    2
#define NOT_MENU     3

const uint NUMBER_OF_SLIDERS  = 4;
const uint NUMBER_OF_RADIOS   = 9;

// Default menu
const std::string lblank[6] = { "blank", "blank", "blank", "blank",   "close", "blank" }; // 0
//                                 0        1        2        3          4        5
// Main menu
const std::string l0[6] = { "mode", "view",  "show",  "tools", "close", "colors" }; // 0
//                             0       1        2        3        4        5
// Controls menu
const std::string l1[6] = { "menu", "color", "zoom",  "move", "close", "copy" };   // 1
//                            0        1        2        3      4        5
// Atom menu
const std::string l2[6] = { "edit", "axis", "select", "frag", "close", "fix" };   // 2
//                            0      1        2        3        4       5
//
// Main menu buttons
#define MAIN_MODE_SUBMENU     0
#define MAIN_VIEW_SUBMENU     1
#define MAIN_SHOW_SUBMENU     2
#define MAIN_TOOLS_SUBMENU    3
#define MAIN_COLOR_SUBMENU    5
// Atom menu buttons
#define ATOM_ELEMENT_SUBMENU  0
#define ATOM_AXIS_SUBMENU     1
#define ATOM_ACTIVE_SUBMENU   2
#define ATOM_FRAGMENT_SUBMENU 3
#define ATOM_DELETE_SUBMENU   5
//
// Submenues                      0       1         2         3        4         5
// Main menu
// Mode submenu                     PBC     symetry
const std::string sl_mode[6] = { "--",    "--",    "frags",   "atoms", "cancel", "--" };
// View submenu
const std::string sl_view[6] = { "XY (z)", "YZ (x)", "ZX (y)", "YX (Z)", "ZY (X)", "XZ (Y)" };
// Show submenu
const std::string sl_show[6] = { "labels", "box",  "symbols", "bonds", "cancel", "numbers" };
// Tools submenu
const std::string sl_tools[6] = { "params", "box",  "axes", "split", "cancel", "atoms" };


// Submenues                       0       1         2         3        4         5
// Atom menu
// Fragment submenu
const std::string l_frag[6]  = { "vdW",   "atom", "delete", "merge", "cancel", "radial" };
// Axis submenu
const std::string l_axis[6]  = { "show",   "head", "tail", "plane", "cancel", "zero" };
// Edit submenu
const std::string l_edit[6]  = { "delete",   "copy", "symbol", "color", "cancel", "size" };

#define CLOSE_MENU             4
#define CLOSE_SUBMENU          4

// Mode submenu buttons
#define MODE_PBC_BUTTON        0
#define MODE_SYMMETRY_BUTTON   1
#define MODE_FRAGMENTS_BUTTON  2
#define MODE_ATOMS_BUTTON      3
#define MODE_BLANK_BUTTON      5

// View submenu buttons
#define VIEW_XY_BUTTON         0
#define VIEW_YZ_BUTTON         1
#define VIEW_ZX_BUTTON         2
#define VIEW_YX_BUTTON         3
#define VIEW_ZY_BUTTON         4
#define VIEW_XZ_BUTTON         5

// Show submenu buttons
#define SHOW_LABELS_BUTTON     0
#define SHOW_BOX_BUTTON        1
#define SHOW_SYMBOLS_BUTTON    2
#define SHOW_BONDS_BUTTON      3
#define SHOW_NUMBERS_BUTTON    5

// Tools submenu buttons
#define TOOLS_MESURE_BUTTON    0
#define TOOLS_BOX_BUTTON       1
#define TOOLS_AXES_BUTTON      2
#define TOOLS_FRAGMENT_BUTTON  3
#define TOOLS_ATOMS_BUTTON     5

// Fragment submenu buttons
#define ATOM_FRAGMENT_VDW_BUTTON    0
#define ATOM_FRAGMENT_ATOM_BUTTON   1
#define ATOM_FRAGMENT_DELETE_BUTTON 2
#define ATOM_FRAGMENT_MERGE_BUTTON  3
#define ATOM_FRAGMENT_RADIAL_BUTTON 5

// Axis submenu buttons
#define AXIS_SHOW_BUTTON       0
#define AXIS_HEAD_BUTTON       1
#define AXIS_TAIL_BUTTON       2
#define AXIS_PLANE_BUTTON      3
#define AXIS_ZERO_BUTTON       5

//#define FRAGMENT_MENU 3
#define AUTO_SUBMENU           11


class Fl_Gl_Mol_View : public Fl_Gl_Atom{

public:

  real size;
  //gm_rgb oldrgb, rgb;
  Fl_Gl_Mol_View();
  Fl_Gl_Mol_View(int,int,int,int,const char* l=0);
  ~Fl_Gl_Mol_View();
  //void initData(void);
  void clear(void);
  bool initialize(void);
  void graph_cb(void);
  //  void redraw(void){draw();};
  void set_background_color(real,real,real);
  void set_foreground_color(real,real,real);
  // Scale function
  void set_zoom(real f){ if(f>0.0001) zoom = f;};
  void set_zoom_step(real f){ if(f>0.0001) zoom_step = f;};
  // the rotation about the vertical (x) axis
  void x_angle(real f){ __x_ang = f;redraw();};
  real x_angle(){return __x_ang;};
  // the rotation about the horizontal (y) axis
  void y_angle(real f){ __y_ang = f;redraw();};
  real y_angle(){return __y_ang ;};
  // the rotation about the (z) axis
  void z_angle(real f){ __z_ang = f;redraw();};
  real z_angle(){return __z_ang ;};
  // Begin GUI controls
  void set_position_x(real f){x_shift=f;};
  void set_position_y(real f){y_shift=f;};
  void set_position_z(real f){z_shift=f;};
  //
  void set_view_xy_front(void);
  void set_view_yz_front(void);
  void set_view_zx_front(void);
  void set_view_xy_back(void);
  void set_view_yz_back(void);
  void set_view_zx_back(void);
  void set_view(real,real,real);
  void set_highlight_atom(int);
  void set_highlight_atom_a(int);
  void set_highlight_atom_b(int);
  void set_select_begin(int);
  void set_select_end(int);
  //
  void set_selected_atom(uint);
  void set_active_slider(uint);
  //void set_active_fragment(uint);
  //void set_new_fragment(uint);
  void set_active_radio(uint,bool);
  //
  //void set_fragment_total(uint);
  void set_update_active_fragment(void);
  //
  void set_atom_brightness(real);
  void set_bond_brightness(real);
  void set_background_brightness(real);
  void set_highlihght_brightness(real);
  void set_select_brightness(real);

  // Set radius scaling factor
  void set_atom_radius_scale(real);
  void set_bond_radius_scale(real);
  //void set_axis_precession(real f){ __axis_precession=f;};
  //void set_axis_tilt(real f){ __axis_tilt=f;};
  //void set_backbone_precession(real f){ __backbone_precession=f;};
  //void set_backbone_tilt(real f){ __backbone_tilt=f;};
  //
  //void set_basis_vectors(const TVector<real>&,const TVector<real>&,const TVector<real>&);
  // End GUI controls
  void set_atomic_labels(const TVector<std::string>&);
  void set_atomic_symbols(const TVector<std::string>&);
  void set_atomic_symbol_table(const TVector<std::string>&);
  void set_atomic_numbers(const TVector<uint>&);
  void set_atom_table(const TVector<uint>&);
  void set_atomic_number_table(const TVector<uint>&);
  //void set_fragment_table(const TVector<uint>&);
  //
  uint get_highlight_atom(void);
  uint get_action(void){ return last_action;};
  //
  void is_graphics(bool);
  void is_highlight_atom(bool);
  void is_mask_atoms(bool);
  void is_draw_bbox(bool);
  void is_draw_world_axes(bool);
  void is_draw_molecular_axes(bool);
  void is_draw_molecular_axis(bool);
  void is_draw_bonds(bool);
  void is_draw_labels(bool);
  void is_draw_symbols(bool);
  void is_draw_tools(bool);
  void is_draw_numbers(bool);
  void set_lock_controls(bool);
  void is_highlight_fragment(bool);
  void set_pcb(bool);
  //inline int sign_of(int i){ return (i==0)?0:(i<0?-1:1); }
  int WindowDump(void);
  int handle(int);
  void handle_main_menu(void);
  void handle_atom_menu(void);
  void handle_controls_menu(void);
  //
  // void set_update_coordinates(bool b){ update_coordinates=b;};
  virtual void view_redraw(void){ redraw();};
  //

#if HAVE_GL
#else
#endif // HAVE_GL
#if HAVE_GL
    void draw();
#endif // HAVE_GL

private:

  // graphic quality control
  int  __sphere_strip_size, _draw_mpts, _draw_data;
  int  i_sphere_resolution;
  int  __highlight_atom_a, __highlight_atom_b;
  int  __select_begin, __select_end;
  //
  int font_size_symbol;
  int font_size_panel_label;
  int font_size_slider_label;
  int font_size_pie_label;
  //
  int u_slider_index;
  int u_radio_index;
  //
  uint last_action;
  uint __highlight_atom;
  uint __last_highlight_atom;
  //
  //uint __fragment_total;
  //uint __fragment_active;
  uint __total_species;
  //
  uint u_menu_index;
  uint u_submenu_index;
  uint u_selected_index;
  uint u_unselected_atom;
  //
  uint u_active_menu;
  // Gl window parameters
  real view_left;
  real view_right;
  real view_bottom;
  real view_top;
  real view_near;
  real view_far;
  real x_factor, y_factor;
  real view_axis_x, view_axis_y;
  //
  real fgred, fggreen, fgblue;
  real bgred, bggreen, bgblue;
  real f_atom_brightness;
  real f_bond_brightness;
  real f_background_brightness;
  real f_highlight_brightness;
  real f_select_brightness;
  real f_atom_brightness_max;
  real f_highlight_brightness_max;
  real f_select_brightness_max;
  //
  real x_shift, y_shift, z_shift;
  real __x_ang, __y_ang, __z_ang;
  real y_off, _scl, zoom_step;
  real shift_factor;
  //real __axis_precession;
  //real __axis_tilt;
  //real __backbone_precession;
  //real __backbone_tilt;
  //
  real r_distance1;
  real r_distance2;
  real r_distance3;
  real r_angle1;
  real r_angle2;
  real r_dihedral;
  //
  TVector<real> v_distance1;
  TVector<real> v_distance2;
  TVector<real> v_distance3;
  //
  GLdouble menu_pos_cx, menu_pos_cy;
  GLdouble menu_pos_x, menu_pos_y, menu_pos_z;
  GLdouble click_pos_x, click_pos_y;
  GLdouble side_pos_x, side_pos_y;
  GLdouble submenu_pos_cx, submenu_pos_cy;
  GLdouble submenu_pos_x, submenu_pos_y, submenu_pos_z;
  GLdouble label_atom_pos_x, label_atom_y, label_atom_z;
  GLdouble label_symbol_pos_x, label_symbol_y, label_symbol_z;
  GLdouble label_menu_pos_x, label_menu_y, label_menu_z;
  //
  uint v_selected_atoms[4];
  unsigned char pixel[3];
  //
  std::string label;
  std::string legends[6];
  std::string sub_label;
  std::string sub_legends[6];
  //
  //GLint viewport[4];
  //GLdouble modelview[16];
  //GLdouble projection[16];
  // atomic spheres
  GLuint    sphere_dl;
  GLuint    cylinder_dl;
  // what is displayed
  bool is_graphics_on, is_draw_world_axes_;
  bool is_draw_molecular_axes_, is_draw_molecular_axis_;
  bool is_normal_color_, is_highlight_atom_;
  bool is_dark_mask_, is_highlight_fragment_;
  bool is_pbc;
  bool is_draw_processing;
  bool is_draw_menu;
  bool is_draw_controls;
  bool is_lock_controls;
  bool is_draw_tools_;
  bool is_draw_pie_menu;
  bool is_draw_pie_submenu;
  bool is_draw_line;
  bool is_draw_point;
  bool is_atom_picked;
  bool is_background_picked;
  bool is_menu_picked;
  bool is_menu_pie_picked;
  bool is_submenu_pie_picked;
  bool is_menu_position;
  bool is_lock_dragging;
  //
  bool is_mode_atom;
  bool is_mode_fragment;
  //
  bool is_left_click;
  bool is_right_click;
  //
  //bool is_update_mask_rcolor;
  bool is_update_position;
  //
  //bool update_coordinates;
  bool update_bonds_color;
  bool update_normal_color;
  bool update_dark_mask;
  bool update_highlight_fragment;
  bool update_highlight_atom;
  bool update_selected_atoms;
  //
  bool is_unselected_atom;
  //
  bool is_handle_atom_;
  bool is_handle_main_;
  //
  bool is_control_left_on;
  bool is_slider_active[NUMBER_OF_SLIDERS];
  bool is_radio_active[NUMBER_OF_RADIOS];
  //
  void eval_system_properties(void);
  void eval_mask_rcolor(void);
  void eval_tool_parameters(void);
  // drawing functions
  void clear_all(void);
  void clear_scene(void);
  void draw_scene(void);
  void draw_atoms(void);
  void draw_axes(void);
  void draw_bonds(void);
  void draw_symbols(void);
  void draw_selected_numbers(void);
  void draw_box(void);
  // GL-GUI /////////////////////////////////
  void resize(int X,int Y,int W,int H);
  void draw_pie_menu(GLfloat cx, GLfloat cy, GLfloat z, GLfloat r, GLint num_segments);
  void draw_pie_submenu( GLfloat z, GLfloat r, GLint num_segments);
  void set_pie_labels(const std::string l[], std::string);
  void draw_controls(GLfloat z);
  void draw_message(GLfloat z);
  void draw_tools(GLfloat z);
  void draw_settings(GLfloat z);
  void draw_information(GLfloat z);
  void draw_slider(GLfloat x1, GLfloat y1, GLfloat y2, GLfloat z, GLfloat val, bool active, char* l);
  void draw_radio_button(GLfloat x1, GLfloat y1, GLfloat y2, GLfloat z, GLfloat val, bool active, char* l);
  void draw_switch_output(GLfloat x1, GLfloat y1, GLfloat y2, GLfloat z, bool val, char* l);
  void widget_float_output(GLfloat x1, GLfloat y1, GLfloat y2, GLfloat z, float val, char* l);
  void widget_float_output_xy(GLfloat x1, GLfloat y1, GLfloat y2, GLfloat z, int x, int y, float val, char* l);
  void widget_text_output(GLfloat x1, GLfloat y1, GLfloat y2, GLfloat z, char* text, char* l);
  void widget_int_output(GLfloat x1, GLfloat y1, GLfloat y2, GLfloat z, int val, char* l);
  void widget_vector_output(GLfloat x1, GLfloat y1, GLfloat y2, GLfloat z, TVector<real> v, char* l);
  //void draw_pie_submenu(GLfloat cx, GLfloat cy, GLfloat z, GLfloat r, GLint num_segments);
  //void draw_pie_sub_menu(GLfloat,GLfloat,GLfloat,GLint);
  void draw_pie(GLfloat x1, GLfloat y1, GLfloat z, std::string l[], GLint nl, GLfloat r, GLint n);
  void draw_sub_pie(GLfloat x1, GLfloat y1, GLfloat z, std::string l[], GLint nl, GLfloat r, GLint n);
  void draw_pie_disk(GLfloat x, GLfloat y, GLfloat z, GLfloat r, GLint n);
  void draw_pie_labels(GLfloat cx, GLfloat cy, GLfloat z, GLfloat r, std::string l[], std::string m, GLint nl);
  //
  // GUI
  // begin sphere
  void displace_scale(point*,float,point*);
  // end sphere
  inline void view_reshape(int width, int height);
  void set_font_size(void);
  void initialize_opengl(void);
  void process_picking(unsigned char pc[3]);
  // deprecated begin
  void process_start_picking(void);
  void process_stop_picking(void);
  void process_mouse_hits(GLint, GLuint*);
  // deprecated end
  //
  void create_sphere_dl(void);
  void delete_sphere_dl(void);
  void create_cylinder_dl(void);
  //
  inline int sign_of(int i){
    return (i<0 ?-1:1);
  };
  //
  inline void set_mouse_motion(int,int);
  inline void point_to_vector(int, int, int, int,GLfloat v[3]);
  // font functions
  //  drawLetter() interprets the instructions from the array
  //  for that letter and renders the letter with line segments.
  //void drawLetter(CPfont *l);
  // Create a display list for each of 6 characters
  //void init_font (void);
  //void printStrokedString(char *s);
};

#endif //

