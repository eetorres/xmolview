//========================================================================
// FILE - Fl_Gl_Atom.cxx                                                //
// Low level primitives for atomic structure visualization using OpenGL //
//========================================================================
//                                                                      //
// OpenGL 3D visualization widget to display atomic structures.         //
//                                                                      //
// Copyright 2002-2015 by Edmanuel Torres                               //
// email: eetorres@gmail.com                                            //
//                                                                      //
// Lastest update: 13/07/2015                                           //
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

#ifndef _FL_GL_ATOM_H_
#define _FL_GL_ATOM_H_

#include <config_debug.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctimer.h>
#include <atom_color.h>
#include <atom_symbol.h>
#include <atom_name.h>
#include <msmvtl/const.h>
#include <msmvtl/linalg.h>
#include <assert.h>
#include <FL/Fl.H>

#define HAVE_GL 1
#if HAVE_GL
    #include <FL/Fl_Gl_Window.H>
    #include <FL/gl.h>
    #include <FL/glu.h>
#else
    #include <FL/Fl_Box.H>
#endif // HAVE_GL

#include <cpalette.h>
#include <cviewmol.h>
#include <sphere.h>

// Neighbouring cells
const int neighbor_cells[27][3] = {
{  0, 0, 0}, // 0
{  1, 0, 0}, // 1
{  1, 1, 0}, // 2
{  0, 1, 0}, // 3
{ -1, 1, 0}, // 4
{  0, 0, 1}, // 5
{  1, 0, 1}, // 6
{ -1, 0, 1}, // 7
{  1, 1, 1}, // 8
{  0, 1, 1}, // 9
{ -1, 1, 1}, // 10
{ -1,-1, 1}, // 11
{  0,-1, 1}, // 12
{  1,-1, 1}, // 13
{  0, 0,-1}, // 14
{  1, 0,-1}, // 15
{ -1, 0,-1}, // 16
{  1, 1,-1}, // 17
{  0, 1,-1}, // 18
{ -1, 1,-1}, // 19
{ -1,-1,-1}, // 20
{  0,-1,-1}, // 21
{  1,-1,-1}, // 22
{ -1, 0, 0}, // 23
{ -1,-1, 0}, // 24
{  0,-1, 0}, // 25
{  1,-1, 0}  // 26
};

static const GLfloat light_ambient[]  = { 0.3, 0.3, 0.3, 1.0};
static const GLfloat light_diffuse[]  = { 0.5, 0.5, 0.5, 1.0};
static const GLfloat light_specular[] = { 1.0, 1.0, 1.0, 1.0};
static const GLfloat light_position[] = { 1.0, 1.0, 1.5, 1.0};

static const GLfloat mat_specular[] = { 0.2, 0.2, 0.2, 0.8};
static const GLfloat mat_diffuse[] = {0.5, 0.5, 0.5, 1};

//static GLfloat scaled_light_ambient[4];
//static GLfloat scaled_light_diffuse[4];
//static GLfloat scaled_light_specular[4];
//static GLfloat scaled_light_position[4];

#define GLV(v)    { glNormal3f(v[0],v[1],v[2]); \
                    glVertex3f(scale*v[0],scale*v[1],scale*v[2]); \
                  }

#define NORMV(v,n)  { point x; x = v; \
                    normalize_point(&x); \
                    vt[0]=x.x; \
                    vt[1]=x.y; \
                    vt[2]=x.z; \
                    m_sphere[n]=vt; \
                    }

const uint MENU_RESERVED_IDS  = 100;

class Fl_Gl_Atom: public Fl_Gl_Window, public CViewmol{

public:

  Fl_Gl_Atom();
  Fl_Gl_Atom(int,int,int,int,const char* l=0);
  ~Fl_Gl_Atom(){};
  int handle(int);
  //
  void initialize_transform_matrix(void);
  void initialize_rotation_matrix(void);
  void initialize_atomic_coordinates(const TMatrix<real>&);
  // cell list
  void set_cells(void);
  void set_inverse_cell(void);
  void set_cell_list(void);
  void set_bounding_box(const TMatrix<real>&);
  //
  void set_x_cells(int);
  void set_y_cells(int);
  void set_z_cells(int);
  void set_sphere_resolution(uint);
  void set_xyz_cells(void);
  //
  void set_fragment_total(const uint);
  void set_fragment_active(uint);
  void set_active_fragment_index(const uint u){ __fragment_active=u;};
  void set_active_fragment(const uint);
  //void set_new_fragment(const uint);
  //void set_atom_fragment(const uint);
  //
  void set_fragment_table(const TVector<uint>&);
  // Tue Feb 26 12:25:32 MST 2013
  // Migrated from fl_gl_mol_view.h
  void set_axis_position(const TVector<real>&);
  void set_axis_precession(real f){ __axis_precession=f;};
  void set_axis_tilt(real f){ __axis_tilt=f;};
  void set_backbone_precession(real f){ __backbone_precession=f;};
  void set_backbone_tilt(real f){ __backbone_tilt=f;};
  uint get_bond_index(const real f){ return (uint)(f/f_atom_bond_delta);};
  //
  ////////////////////////////////////////////////////////////////////////////////////////
  // Evaluation functions
  ////////////////////////////////////////////////////////////////////////////////////////
  void eval_initial_properties(void);
  void eval_linked_list(void);
  void eval_atomic_bonds(void);
  void eval_sphere(uint);
  void eval_cylinder(uint);

  ////////////////////////////////////////////////////////////////////////////////////////
  // Compute functions
  ////////////////////////////////////////////////////////////////////////////////////////
  void compute_vdw_fragment(uint);
  void compute_atom_fragment(uint);
  void compute_radial_fragment(uint,real);
  void compute_vdw_fragments(void);
  void compute_atom_fragments(void);
  void compute_merge_fragments(const uint);

  ////////////////////////////////////////////////////////////////////////////////////////
  // Update functions
  ////////////////////////////////////////////////////////////////////////////////////////
  void update_atomic_coordinates(const TMatrix<real>&);
  void update_atomic_bonds(void);
  void update_fragments(uint,bool);
  //
  ////////////////////////////////////////////////////////////////////////////////////////
  // OpenGL functions
  ////////////////////////////////////////////////////////////////////////////////////////
  void initialize_sphere(real);
  void initialize_cylinder(real);
  void add_stick(const TVector<real>&,real,real,real,real);
  void add_axis(const TVector<real>&,real,real,real,real);
  //
  inline uint check_bond(uint);
  inline void linearly_interpolate(point*,point*,float,point*);
  inline void normalize_point(point*);
  //
  int get_total_cells(void){ return total_cells;};
  //
  int  pos_x_cells;
  int  pos_y_cells;
  int  pos_z_cells;
  int  neg_x_cells;
  int  neg_y_cells;
  int  neg_z_cells;
  int  total_cells;
  int  x_cells;
  int  y_cells;
  int  z_cells;
  //
  // virtual functions
  virtual void update_data(void);
  virtual void set_update_coordinates(bool b){ update_coordinates=b;};
  virtual void save_wysiwyg_as(std::string);
  virtual void save_wysiwyg_extension(std::string);
  virtual void save_wysiwyg_as(std::string,std::string);
  //
protected:
  //
  CPalette palette;
  CPalette index_palette;
  //
  bool is_first_structure_;
  bool is_initialize_rot;
  bool is_eval_bonds;
  bool is_eval_sphere;
  //
  bool is_draw_labels_;
  bool is_draw_symbols_;
  bool is_draw_numbers_;
  bool is_draw_atoms_;
  bool is_draw_bonds_;
  bool is_draw_bbox_;
  //
  bool is_update_atomic_properties;
  bool is_update_radius;
  bool is_update_bonds;
  bool is_linked_cell;
  bool is_update_mask_rcolor;
  bool update_coordinates;
  //
  int  __number_of_atoms;
  uint  i_number_of_bonds;
  uint i_number_of_bonds_pbc;
  uint u_sphere_resolution;
  int  u_sphere_rows;
  int  __sphere_strip_length;
  int  __cylinder_strip_length;
  int  u_cylinder_resolution;
  int  i_neighbor_cells;
  //
  uint u_bond_types;
  uint __fragment_total;
  uint __fragment_active;
  //
  real f_atom_bond_delta;
  real f_atom_bond_inv_delta;
  real f_atom_radius_scale;
  real f_bond_radius_scale;
  real base_view;
  real r_axes_position;
  real zoom;
  //
  real __axis_precession;
  real __axis_tilt;
  real __backbone_precession;
  real __backbone_tilt;
  // begin trackball
  GLboolean tb_tracking;
  GLboolean tb_animate;
  GLuint    tb_width;
  GLuint    tb_height;
  GLint     tb_button;
  GLuint    tb_lasttime;
  GLfloat   tb_angle;
  GLfloat   tb_axis[3];
  GLfloat   tb_transform[4][4];
  GLfloat   rot_matrix[4][4];
  GLfloat   tb_lastposition[3];
  GLfloat   scaled_light_position[4];
  // end trackball
  // axis
  TVector<real> v_axes_position;
  TVector<real> v_axis_position;
  TVector<real> v_bond_length;
  TVector<uint> v_atom_numbers; // <-----------------
  TVector<uint> v_atom_table;
  TVector<uint> v_bond_number;
  TVector<uint> v_bond_number_pbc;
  TVector<uint> v_bond_table;
  TVector<uint> v_fragment_table_gl;
  TVector<uint> v_atomic_number_table_gl;
  TVector<GLuint> v_sphere_list;
  TVector<GLuint> v_cylinder_list;
  TVector<GLuint> v_cylinder_list_pbc;
  //
  TVector<real> v_bbox, _vu, _vv, _vw;
  //
  // linked cell variables
  uint u_cell_number, u_particle_number;
  real r_cut_radius, r_cut_radius_2;
  TVector<real> v_cell_frac, v_box_size, v_box_middle;
  TVector<uint> v_cell_number;
  TVector<int>  v_cell_side, v_cell_list, v_cell_head;
  TVector<int>  v_neighbor_cell, v_md_pbc;
  //
  TVector<std::string> v_atom_labels; // <--------------
  TVector<std::string> v_atom_symbols; // <--------------
  TVector<std::string> v_atomic_symbol_table_gl;
  //
  TMatrix<uint> m_bond_indices;
  TMatrix<uint> m_bond_indices_pbc;
  // bonds with PBC
  TMatrix<int>  m_bond_boundary_pbc;
  TMatrix<int>  neighbor_cells_xyz;
  //
  TMatrix<real> m_sphere;
  TMatrix<real> m_cylinder;
  // cylinder
  // this must be removed, still used in "add_stick"
  TMatrix<real> m_cylinder_e1;
  // this may be useful?
  TMatrix<real> m_cylinder_texture1;
  // atoms
  TMatrix<real> m_atom_coordinates;
  TMatrix<real> m_atom_position;
  TMatrix<real> m_radius_color;
  TMatrix<real> m_atom_rcolor;
  //
  TMatrix<real> m_bond_rcolor_0;
  TMatrix<real> m_bond_rcolor_1;
  TMatrix<real> m_bond_rcolor_pbc_0;
  TMatrix<real> m_bond_rcolor_pbc_1;
  // bonds no PBC
  TMatrix<real> m_bond_position;
  TMatrix<real> m_bond_angles;
  // bonds with PBC
  TMatrix<real> m_bond_position_pbc;
  TMatrix<real> m_bond_angles_pbc;
  // bounding box
  TMatrix<real> m_bbox;
  TMatrix<real> u_bbox;
  TMatrix<real> u_inv_bbox;
  //
  private:
    CTimer gl_atom_clock;
};

#endif //

