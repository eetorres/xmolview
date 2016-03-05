//========================================================================
// FILE - fl_gl_atom.cxx                                              //
// For the Fast Light Tool Kit (FLTK) - www.fltk.org                    //
//========================================================================
//                                                                      //
// OpenGL 3D visualization widget to display atomic structures.         //
//                                                                      //
// Copyright 2002-2015 by Edmanuel Torres                               //
// email: eetorres@gmail.com                                            //
//                                                                      //
// Lastest update: 13/07/2015
//======================================================================//

#include <config_debug.h>
#include <fl_gl_atom.h>


// Selection Buffer
//#define SelBufferSize 512

// Picking Stuff //
//#define MODE_RENDER    1
//#define MODE_SELECT    2
//#define BUFSIZE     1024

//GLuint selectBuf[BUFSIZE];
//GLint hits;
//int render_mode = MODE_RENDER;
//int cursorX, cursorY;

//#define ALPHA 0.5
//#define HAVE_GL 1

#if HAVE_GL
Fl_Gl_Atom::Fl_Gl_Atom(int x,int y,int w,int h,const char *l) : Fl_Gl_Window(x,y,w,h,l)
#else
Fl_Gl_Atom::Fl_Gl_Atom(int x,int y,int w,int h,const char *l) : Fl_Box(x,y,w,h,l)
#endif // HAVE_GL
{
  pos_x_cells = 0;
  pos_y_cells = 0;
  pos_z_cells = 0;
  neg_x_cells = 0;
  neg_y_cells = 0;
  neg_z_cells = 0;
  u_sphere_resolution = 0;
  u_cylinder_resolution = 10;
  i_number_of_bonds = 0;
  i_number_of_bonds_pbc = 0;
  is_eval_sphere=true;

  total_cells = 1;
  x_cells = 1;
  y_cells = 1;
  z_cells = 1;
#if !HAVE_GL
  label("OpenGL is required for this demo to operate.");
  align(FL_ALIGN_WRAP | FL_ALIGN_INSIDE);
#endif // !HAVE_GL
}

int Fl_Gl_Atom::handle(int event){
  int ret = Fl_Gl_Window::handle(event);
  /*switch(event){
  case FL_FOCUS:
    //std::cout<<" got the focus"<<std::endl;
    //if(Fl::focus() != this)
    //    Fl::focus(this);
    //return 1;
  case FL_UNFOCUS:
    //if(Fl::focus() != this)
    //    Fl::focus(this);
    //... Return 1 if you want keyboard events, 0 otherwise
    //redraw();
    return 1;
  }*/
  return ret;
}

void Fl_Gl_Atom::eval_initial_properties(void){
  real _radius;
  TVector<real> _atom_xyz;
  TVector<real> _rcolor(4);
  m_radius_color.resize(__number_of_atoms,4);
  //index_palette.set(__number_of_atoms+MENU_RESERVED_IDS);
  //index_palette.set(__number_of_atoms);
  index_palette.set_color(0);
  index_palette.initialize(0,__number_of_atoms+MENU_RESERVED_IDS,__number_of_atoms+MENU_RESERVED_IDS);
  index_palette.update_palette_index();
  int i_z=0;
#ifdef _ATOM_DEBUG_MESSAGES_
  std::cout<<" ATOM: Number of atoms: "<<__number_of_atoms<<std::endl;
#endif
  if(is_eval_sphere){
    set_sphere_resolution(2);
    //eval_sphere(u_sphere_resolution);
    is_eval_sphere=false;
  }
  if(is_draw_atoms_){
    for(int i=0; i<__number_of_atoms; i++){
      i_z=v_atom_numbers[i];
      _radius=atom_rrgb[i_z][0];
      _rcolor[0]=_radius;
      _atom_xyz = m_atom_coordinates[i];
      //m_text_position.add_row(_atom_xyz); // text position
      _rcolor[1] = atom_rrgb[i_z][1];
      _rcolor[2] = atom_rrgb[i_z][2];
      _rcolor[3] = atom_rrgb[i_z][3];
      //m_radius_color.add_row(_rcolor);
      m_radius_color[i]=_rcolor;
    }
  }
  m_atom_rcolor = m_radius_color;
//#ifdef _ATOM_DEBUG_MESSAGES_
//  std::cout<<" m_atom_rcolor = "<<m_atom_rcolor;
//#endif
  //eval_atom_spheres(); //<--------------
}

void Fl_Gl_Atom::set_axis_position(const TVector<real>& v){
  v_axis_position=v;
}

void Fl_Gl_Atom::update_atomic_coordinates(const TMatrix<real>& m){
  m_atom_coordinates = m;
  is_update_bonds = true;
  is_update_atomic_properties = true;
  is_update_mask_rcolor = true;
  set_xyz_cells();
}

void Fl_Gl_Atom::initialize_atomic_coordinates(const TMatrix<real>& m){
  //clear();
#ifdef _ATOM_DEBUG_MESSAGES_
  std::cout<<" ATOM: Initialize"<<std::endl;
#endif
  //m_atom_coordinates = m;
  update_atomic_coordinates(m);
  __number_of_atoms = m_atom_coordinates.rows();
  if(__number_of_atoms<100){
    is_linked_cell=false;
  }
  //cout<<" ATOM: New number of atoms: "<<__number_of_atoms<<std::endl;
  //std::cout<<" atoms: "<<m_atom_coordinates;
  if(is_first_structure_){
    //std::cout<<" ATOM: Initialize 1-1"<<std::endl;
#ifdef _ATOM_DEBUG_MESSAGES_
    std::cout<<" ATOM: First structures"<<std::endl;
#endif
    //std::cout<<"Sphere level: "<<__sphere_strip_size<<std::endl;
    //cout<<"m_atoms_strip.rows = "<<m_atoms_strip.rows()<<std::endl;
    if(is_initialize_rot)
      initialize_rotation_matrix();
    initialize_transform_matrix();
    is_draw_atoms_=true;
    is_first_structure_=false;
    //std::cout<<" ATOM: Initialize 1-2"<<std::endl;
  }
  //m_atoms_strip.resize(__number_of_atoms*__sphere_strip_length,3);
  //std::cout<<" ATOM: Initialize 1-3"<<std::endl;
  is_update_atomic_properties = true;
  if(is_draw_bonds_){
    is_update_bonds=true;
    //update_bonds_color=true;
  }
  is_update_radius=true;
  eval_initial_properties(); // <---------------------!!!@
  //set_xyz_cells();
  //eval_system_properties();
#ifdef _ATOM_DEBUG_MESSAGES_
  std::cout<<" ATOM: Initialized"<<std::endl;
#endif
}

// this the function in charge of the scene actualization
void Fl_Gl_Atom::update_data(void){
  update_view();
  if(update_coordinates){
    update_atomic_coordinates(get_view_cartesian());
    set_axis_position(get_view_centered_position_cartesian());
    set_axis_precession(get_view_axis_precession());
    set_axis_tilt(get_view_axis_tilt());
    set_backbone_precession(get_view_backbone_precession());
    set_backbone_tilt(get_view_backbone_tilt());
    update_coordinates=false;
  }
  Fl::wait(0.1);
}

void Fl_Gl_Atom::update_atomic_bonds(void){
#ifdef _ATOM_DEBUG_MESSAGES_
  std::cout<<" ATOM: start - update_atomic_bonds(void)"<<std::endl;
#endif
  //bool use_pbc;
  uint i, j;//,  nbonds;
  //real rl, r, r2, rr, ri, rj, rlz, rlxy;
  //real r, r2;
  //real ri, rj, rlz, rlxy;
  //real ri, rlz, rlxy;
  real rlz, rlxy;
  //TVector<uint> vidx(2);
  TVector<real> vi, vj, vij, vang(2);
  TVector<real> vi_uvw, vj_uvw, vij_uvw;
  //
  TVector<int>  v_pbc;
  //vidx[0]=i;
  //vidx[1]=j;
  //r = sqrt(r2);
  //vij = (vj-vi);
  //rlz = vij[2];
  //vij[2] = 0;
  //rlxy = vij.magnitude();
  //vang[0] = RAD_DEG*atan2(vij[1],vij[0]);       // bond precession
  //vang[1] = RAD_DEG*fabs(atan2(rlxy,rlz));      // bond tilt
  //vij = 0.5*(vj+vi);
  //k = (uint)(r*f_atom_bond_inv_delta);
  //m_bond_indices.add_row(vidx);

  // Inside the box
  //m_bond_angles.add_row(vang);
  //m_bond_position.add_row(vij);
  //v_bond_number.push_back(check_bond(k));
  //m_bond_indices.add_row(vidx);
  //nbonds = m_bond_indices.rows();
  //for(uint n=0; n<nbonds; n++){
  for(uint n=0; n<i_number_of_bonds; n++){
    i=m_bond_indices[n][0];
    j=m_bond_indices[n][1];
    //
    vi = m_atom_coordinates[i];
    //vi_uvw = (vi*u_inv_bbox);
    //ri = m_radius_color[i][0];
    vj = m_atom_coordinates[j];
    //vj_uvw = (vj*u_inv_bbox);
    //rj = m_radius_color[i][1];
    // Sun Feb 24 16:35:37 MST 2013
    // should check boundaries only if the atoms were moved
    // right now is checking all, low efficiency
    //use_pbc = false;
    //vij_uvw = (vj_uvw-vi_uvw);
    /* bonds inside the volume box dont need PBC
    for(uint coord=0; coord<3; coord++){
      //v_pbc[coord] = 0;
      if(vij_uvw[coord] <= -v_bbox[coord]){
        vj += 2.0*m_bbox[coord];       // PBC
        //use_pbc = true;
        //v_pbc[coord] = 1;
      }else if(vij_uvw[coord] > v_bbox[coord]){
        vj -= 2.0*m_bbox[coord];       // PBC
        //use_pbc = true;
        //v_pbc[coord] = -1;
      }
    }*/
    vij = (vj-vi);
    // check if the bond is still within the range
    //r2 = vij.magnitude();
    //r = sqrt(r2);
    //uint k = get_bond_index(r);
    // working here...!!!!
    rlz = vij[2];
    vij[2] = 0;
    rlxy = vij.magnitude();
    vang[0] = RAD_DEG*atan2(vij[1],vij[0]);       // bond precession
    vang[1] = RAD_DEG*fabs(atan2(rlxy,rlz));      // bond tilt
    vij = 0.5*(vj+vi);
    //
    //k = (uint)(r*f_atom_bond_inv_delta);
    m_bond_angles[n]=vang;
    m_bond_position[n]=vij;
  }

  // PBC
  //m_bond_boundary_pbc.add_row(v_pbc);
  //m_bond_angles_pbc.add_row(vang);
  //m_bond_position_pbc.add_row(vij);
  //v_bond_number_pbc.push_back(check_bond(k));
  //m_bond_indices_pbc.add_row(vidx);

  // THERE IS A BIG BUG HERE
  // PBC BONDS DONT SHOW CORRECT
  //nbonds = m_bond_indices_pbc.rows();
  for(uint n=0; n<i_number_of_bonds_pbc; n++){
    i=m_bond_indices_pbc[n][0];
    j=m_bond_indices_pbc[n][1];
    //
    vi = m_atom_coordinates[i];
    vi_uvw = (vi*u_inv_bbox);
    //ri = m_radius_color[i][0];
    vj = m_atom_coordinates[j];
    for(uint coord=0; coord<3; coord++){
      //vj[coord]+=(neighbor_cells[_m][coord]*v_box_size[coord]);                    // PBC
      //vj[coord]+=(neighbor_cells[_m][coord]*2.0*(m_bbox[0][coord]+m_bbox[1][coord]+m_bbox[2][coord]));                    // PBC
      if ( m_bond_boundary_pbc[n][coord] != 0 )
        vj += m_bond_boundary_pbc[n][coord]*2.0*m_bbox[coord];       // PBC
      //v_pbc[coord] = neighbor_cells[_m][coord];
    }
    vj_uvw = (vj*u_inv_bbox);
    //rj = m_radius_color[i][1];
    //
    //r = sqrt(r2);
    vij = (vj-vi);
    rlz = vij[2];
    vij[2] = 0;
    rlxy = vij.magnitude();
    vang[0] = RAD_DEG*atan2(vij[1],vij[0]);       // bond precession
    vang[1] = RAD_DEG*fabs(atan2(rlxy,rlz));      // bond tilt
    vij = 0.5*(vj+vi);
    //k = (uint)(r*f_atom_bond_inv_delta);
    //
    m_bond_angles_pbc[n]=vang;
    m_bond_position_pbc[n]=vij;
  }
#ifdef _ATOM_DEBUG_MESSAGES_
  std::cout<<" ATOM: end - update_atomic_bonds(void)"<<std::endl;
#endif
}

void Fl_Gl_Atom::set_x_cells(int i){
  if(i>=0){
    pos_x_cells = i;
    neg_x_cells = 0;
  }else{
    pos_x_cells = abs(i);
    neg_x_cells = i;
  }
  set_xyz_cells();
  supercell.set_gsf_modified(true);
}

void Fl_Gl_Atom::set_y_cells(int i){
  if(i>=0){
    pos_y_cells = i;
    neg_y_cells = 0;
  }else{
    pos_y_cells = abs(i);
    neg_y_cells = i;
  }
  set_xyz_cells();
  supercell.set_gsf_modified(true);
}

void Fl_Gl_Atom::set_z_cells(int i){
  if(i>=0){
    pos_z_cells = i;
    neg_z_cells = 0;
  }else{
    pos_z_cells = abs(i);
    neg_z_cells = i;
  }
  set_xyz_cells();
  supercell.set_gsf_modified(true);
}

void Fl_Gl_Atom::set_sphere_resolution(uint u){
  if(u_sphere_resolution != u){
    u_sphere_resolution = u;
    eval_sphere(u_sphere_resolution);
    u_cylinder_resolution=(5*u_sphere_resolution+10);
    eval_cylinder(u_cylinder_resolution);
  }
}

void Fl_Gl_Atom::set_xyz_cells(void){
  TVector<real> _x(3);
  TVector<real> _xyz;
  TVector<real> e(3),p(3);
  uint cont = 0;
  x_cells = (pos_x_cells-neg_x_cells+1);
  y_cells = (pos_y_cells-neg_y_cells+1);
  z_cells = (pos_z_cells-neg_z_cells+1);
  total_cells = x_cells*y_cells*z_cells;
#ifdef _ATOM_DEBUG_MESSAGES_
  std::cout<<" total cells = "<<total_cells<<std::endl;
#endif
  m_atom_position.resize(total_cells*__number_of_atoms,3);
  if(is_draw_atoms_){
    for(int x=neg_x_cells; x<pos_x_cells+1; x++){ // repetition in x
      for(int y=neg_y_cells; y<pos_y_cells+1; y++){ // repetition in y
        for(int z=neg_z_cells; z<pos_z_cells+1; z++){ // repetition in z
          for(int i=0; i<__number_of_atoms; i++){
            _xyz=m_atom_coordinates[i];
            _xyz=_xyz+2.0*(x*_vu+y*_vv+z*_vw);
            m_atom_position[cont]=_xyz;
            cont++;
          }
        }
      }
    }
  }
}

void Fl_Gl_Atom::save_wysiwyg_as(std::string _p, std::string _f){
  // in case the bonds have not been evaluated
  if(((supercell.get_output_file_format() == OUTPUT_FORMAT_ATM_NFR) || (supercell.get_output_file_format() == OUTPUT_FORMAT_ATM_FRG)) && (supercell.get_output_file_type() == OUTPUT_FILE_TYPE_DLP)){
    //std::cout<<"bonds have not been evaluated"<<std::endl;
    if(is_eval_bonds){
      eval_atomic_bonds();  // <----------
      is_eval_bonds=false;
    }
    supercell.eval_connections(m_bond_indices,i_number_of_bonds);
  }
#ifdef _ATOM_DEBUG_MESSAGES_
  std::cout<<" ATOM: do save_wysiwyg_as (1)"<<std::endl;
  std::cout<<" ATOM: m_atom_position = "<<m_atom_position;
#endif
  supercell.save_as_file(_p,_f,is_draw_symbols_,is_draw_numbers_,m_atom_position,x_cells,y_cells,z_cells,total_cells,is_draw_bbox_);
}

void Fl_Gl_Atom::save_wysiwyg_as(std::string _f){
  //if((get_output_file_format() == OUTPUT_FORMAT_ATM_NFR) || (get_output_file_format() == OUTPUT_FORMAT_ATM_FRG))
    //std::cout<<"hola";
#ifdef _ATOM_DEBUG_MESSAGES_
  std::cout<<" ATOM: do save_wysiwyg_as (2)"<<std::endl;
#endif
  supercell.save_as_file(_f,is_draw_symbols_,is_draw_numbers_,m_atom_position,x_cells,y_cells,z_cells,total_cells,is_draw_bbox_);
}

void Fl_Gl_Atom::save_wysiwyg_extension(std::string _f){
  std::string _fext;
  switch(supercell.get_output_file_type()){
    case OUTPUT_FILE_TYPE_VSP:
      _fext = _f + ".vsp";
    break;
    case OUTPUT_FILE_TYPE_XYZ:
      _fext = _f + ".xyz";
    break;
    case OUTPUT_FILE_TYPE_GAU:
      _fext = _f + ".gau";
    break;
    case OUTPUT_FILE_TYPE_PDB:  // reserved for xmatrix
      _fext = _f + ".pdb";
    break;
    case OUTPUT_FILE_TYPE_DLP:
      _fext = _f + ".dlp";
    break;
  }
  supercell.save_as_file(_fext,is_draw_symbols_,is_draw_numbers_,m_atom_position,x_cells,y_cells,z_cells,total_cells,is_draw_bbox_);
}

void Fl_Gl_Atom::set_bounding_box(const TMatrix<real>& m){
  m_bbox = m;
  u_bbox.resize(3,3);
  _vu = m_bbox[0];
  _vv = m_bbox[1];
  _vw = m_bbox[2];
  // set up the half box sides
  // v_bbox is the half of the box
  v_bbox[0]=_vu.magnitude();
  v_bbox[1]=_vv.magnitude();
  v_bbox[2]=_vw.magnitude();
  // the full box is used for the bond search method
  v_box_size = 2.0*v_bbox;
  //v_box_size = v_bbox;
  //v_bbox[0]=m_bbox[0][0];
  //v_bbox[1]=m_bbox[1][1];
  //v_bbox[2]=m_bbox[2][2];
  //u_bbox[0] = _vu/_vu.magnitude();
  //u_bbox[1] = _vv/_vv.magnitude();
  //u_bbox[2] = _vw/_vw.magnitude();
  u_bbox=supercell.get_unit_uvw_to_xyz();
  //m_bbox.transpose();
  u_inv_bbox=u_bbox.inverse();
  // set view size
  base_view = maxi(_vu.magnitude(), _vv.magnitude());
  base_view = maxi(base_view, _vw.magnitude());
  r_axes_position = 0.9*base_view;
  base_view *= 1.1;
#ifdef _ATOM_DEBUG_MESSAGES_
  std::cout<<" Unit cell = "<<m_bbox<<std::endl;
  std::cout<<" Half Box = "<<v_bbox<<std::endl;
  std::cout<<" Full Box = "<<v_box_size<<std::endl;
  std::cout<<" Unit Box = "<<u_bbox<<std::endl;
  std::cout<<" Unit Inv Box = "<<u_inv_bbox<<std::endl;
#endif
}

void Fl_Gl_Atom::initialize_sphere(real r){
    int s;
    real scale = 2.0*f_atom_radius_scale*r;
    // iterate over the 20 sides of the icosahedron
    int cont = 0;
    for(s = 0; s < 20; s++){
        for(int i = 0; i < u_sphere_rows; i++){
            // strip the ith trapezoid block
            glBegin(GL_TRIANGLE_STRIP);
            for(int j=0; j<(2*i+3); j++){
                // add one vertice at a time
                GLV(m_sphere[cont]);
                cont++;
            }
            glEnd();
        }
    }
}

void Fl_Gl_Atom::initialize_cylinder(real r){
  TVector<real> t(3);
  glBegin(GL_QUAD_STRIP);
  real scale = f_bond_radius_scale*0.2;
  for (int i=0;i<__cylinder_strip_length;i++){
    t = m_cylinder[i];
    t[2]=0; // -r
    glNormal3f(t[0],t[1],t[2]);
    glVertex3f(scale*t[0],scale*t[1],t[2]);
    t[2]=r;
    glNormal3f(t[0],t[1],t[2]);
    glVertex3f(scale*t[0],scale*t[1],t[2]);
  }
  glEnd();
}

void Fl_Gl_Atom::add_stick(const TVector<real>& c,real l, real r, real a1, real a2){
  TVector<real> e(3), p(3), t(3);
#ifdef _SHOW_DEBUG_ADD_STICK_
  std::cout<<" ATOM: add_stick"<<std::endl;
#endif
  glBegin(GL_QUAD_STRIP);
  for (int i=0;i<=u_cylinder_resolution;i++) {
    e = m_cylinder_e1[i];
    p[0] = r * e[0];
    p[1] = r * e[1];
    p[2] = 0;
    t = vrot_y(p,a2);
    t = vrot_z(t,a1);
    t = c+t;
    e = vrot_y(e,a2);
    e = vrot_z(e,a1);
    glNormal3f(e[0],e[1],e[2]);
    //glTexCoord2f(m_cylinder_texture1[i][0],m_cylinder_texture1[i][1]);
    glVertex3f(t[0],t[1],t[2]);
    //p[0] = r * e[0];
    //p[1] = r * e[1];
    p[2] = l;
    t = vrot_y(p,a2);
    t = vrot_z(t,a1);
    t = c+t;
    glNormal3f(e[0],e[1],e[2]);
    //glTexCoord2f(m_cylinder_texture1[i][0],m_cylinder_texture1[i][1]);
    glVertex3f(t[0],t[1],t[2]);
  }
#ifdef _SHOW_DEBUG_ADD_STICK_
  std::cout<<" ATOM: add_stick end"<<std::endl;
#endif
  glEnd();
}

void Fl_Gl_Atom::add_axis(const TVector<real>& c,real l, real r, real a1, real a2){
#ifdef _SHOW_DEBUG_ADD_AXIS_
  std::cout<<" ATOM: add_axis"<<std::endl;
#endif
  TVector<real> e(3),p(3), t1(2), t2(2);
#ifdef _SHOW_DEBUG_ADD_AXIS_
  std::cout<<" ATOM: check r"<<std::endl;
#endif
  if(r < 0) r = -r;
  if(r <= 0) {
    glBegin(GL_POINTS);
    glVertex3f(c[0],c[1],c[2]);
    glEnd();
    return;
  }
#ifdef _SHOW_DEBUG_ADD_AXIS_
  std::cout<<" ATOM: init"<<std::endl;
#endif
  //glColor3f(0.4,0.9,0.1);
  //add_stick(c,10,l,__axis_precession,__axis_tilt);
  add_stick(c,l,r,a1,a2);
#ifdef _SHOW_DEBUG_ADD_AXIS_
  std::cout<<" ATOM: strip"<<std::endl;
#endif
  glBegin(GL_QUAD_STRIP);
  real tip = l/3.0;
  for (int i=0;i<=u_cylinder_resolution;i++) {
    e = m_cylinder_e1[i];
    p[0] = 2.5*r* e[0];
    p[1] = 2.5*r* e[1];
    p[2] = l;
    p = vrot_y(p,a2);
    p = vrot_z(p,a1);
    p = c+p;
    glNormal3f(e[0],e[1],0);
    glTexCoord2f(m_cylinder_texture1[i][0],m_cylinder_texture1[i][1]);
    glVertex3f(p[0],p[1],p[2]);
    p[0] = 0;
    p[1] = 0;
    p[2] = l+tip;
    p = vrot_y(p,a2);
    p = vrot_z(p,a1);
    p = c+p;
    glNormal3f(e[0],e[1],0);
    glTexCoord2f(m_cylinder_texture1[i][0],m_cylinder_texture1[i][1]);
    glVertex3f(p[0],p[1],p[2]);
  }
  glEnd();
}

// Linked and shell cell configuration functions
void Fl_Gl_Atom::set_cells(void){
  TVector<int> v1(3);
  int xpcb=-2, ypcb=-2, zpcb=-2;
  v_cell_side = iVScale(v_box_size, (1.0/r_cut_radius));
#ifdef _SHOW_DEBUG_LINKED_
  std::cout<<" LINKED: Cell = "<<v_cell_side;
#endif
  // set the necesary neighbor cells
  // neighbor_cells_xyz
  if(v_cell_side[0]==2) xpcb=-1;
  if(v_cell_side[1]==2) ypcb=-1;
  if(v_cell_side[2]==2) zpcb=-1;
  neighbor_cells_xyz.resize(0,3);
  //int count=0;
  i_neighbor_cells=0;
  for(int k=0; k<27; k++){
	v1[0] = neighbor_cells[k][0];
	v1[1] = neighbor_cells[k][1];
	v1[2] = neighbor_cells[k][2];
	if((v1[0] > xpcb) && (v1[1] > ypcb) && (v1[2] > zpcb)){
	  //std::cout<<v1<<std::endl;
	  neighbor_cells_xyz.add_row(v1);
	  i_neighbor_cells++;
	}
  }
  //std::cout<<" Neighbouring cells = "<<i_neighbor_cells<<std::endl;
  //std::cout<<" Neighbouring matrix = "<<neighbor_cells_xyz<<std::endl;
  //
  is_linked_cell=true;
  for(int i=0; i<3; i++){
    if(v_cell_side[1]<=2){
      is_linked_cell=false;
#ifdef _ATOM_DEBUG_MESSAGES_
  std::cout<<" LINKED CELL too small [turned off]"<<std::endl;
#endif
    }
  }
}

void Fl_Gl_Atom::set_inverse_cell(void){
  v_cell_frac = fVDiv(v_cell_side,v_box_size);
#ifdef _SHOW_DEBUG_LINKED_
  std::cout<<" LINKED: Cell Frac = "<<v_cell_frac;
#endif
}

void Fl_Gl_Atom::set_cell_list(void){
  u_cell_number = (uint)vVol(v_cell_side);
#ifdef _SHOW_DEBUG_LINKED_
  std::cout<<" LINKED: Number of used cells = "<<u_cell_number<<std::endl;
#endif
  v_cell_head.resize(u_cell_number);
  v_cell_list.resize(__number_of_atoms);
}

void Fl_Gl_Atom::eval_linked_list(void){
  set_cells();
  set_inverse_cell();
  set_cell_list();
#ifdef _SHOW_DEBUG_LINKED_
  std::cout<<" LINKED: BEGIN: eval_cell_list"<<std::endl;
#endif
  uint _n, u_icell;
  TVector<real> v_positive_r;
  TVector<int>  v_integer_r;
  for(_n=0; _n<u_cell_number; _n++)
    v_cell_head[_n] = -1;
#ifdef _SHOW_DEBUG_LINKED_
  std::cout<<" LINKED: build the linked cell"<<std::endl;
#endif
  for(_n=0; _n<(uint)__number_of_atoms; _n++){
#ifdef _SHOW_DEBUG_LINKED_EVAL_
    std::cout<<" atom coordinates["<<_n<<"] = "<<m_atom_coordinates[_n];
#endif
    v_positive_r=m_atom_coordinates[_n];
    // uvw coordinates
    v_positive_r =  (v_positive_r*u_inv_bbox);
#ifdef _SHOW_DEBUG_LINKED_EVAL_
    std::cout<<" direct r = "<<v_positive_r;
#endif
    // v_positive_r = fVAdd(m_atom_coordinates[_n],v_box_middle);
    // use half of the box here
    v_positive_r = fVAdd(v_positive_r,v_bbox);
#ifdef _SHOW_DEBUG_LINKED_EVAL_
    std::cout<<" positive r = "<<v_positive_r;
#endif
    // apply PBC to place all the atoms inside the box
    /*
    for(uint coord=0; coord<3; coord++){
      if(v_positive_r[coord] <= 0){
        v_positive_r[coord]+= v_box_size[coord];
      }else if(v_positive_r[coord] > v_box_size[coord]){
        v_positive_r[coord]-=v_box_size[coord];
      }
    }*/
    v_integer_r = iVMul(v_positive_r,v_cell_frac);
#ifdef _SHOW_DEBUG_LINKED_EVAL_
    std::cout<<" v_integer_r = "<<v_integer_r;
#endif
    u_icell = iVLinear(v_integer_r,v_cell_side);
#ifdef _SHOW_DEBUG_LINKED_EVAL_
    std::cout<<" u_icell = "<<u_icell<<std::endl;
#endif
    v_cell_list[_n] = v_cell_head[u_icell];
    v_cell_head[u_icell] = _n;
  }
#ifdef _SHOW_DEBUG_LINKED_
  std::cout<<" v_cell_list="<<v_cell_list<<std::endl;
  std::cout<<" END: eval_cell_list"<<std::endl;
#endif
}

// Fri Jan 13 16:55:51 MST 2012
// alpha version
// find bonds between atoms closer than the sum of their van der Waals radius
// it needs linked list to be faster (ASAP). (DONE)
// !!!This Bug was fixed!!!
// there is a bug in this code. It can be reproduce when I load the sam_asphaltane.xyz
// and then I add the bonds.
// !!!This Bug was fixed!!!
//
// Tue May 29 21:40:47 MDT 2012
// !!!This Bug was fixed!!!
// bonds broken when periodic images of the unit cell is used.
// this can be solved using PBC when the bonds are evald.
// !!!This Bug was fixed!!!
//
// !!!This Bug was fixed!!!
// Sun Oct 21 15:04:24 MDT 2012
// bond to the last atom somtimes is missing
// !!!This Bug was fixed!!!
void Fl_Gl_Atom::eval_atomic_bonds(void){
#ifdef _ATOM_DEBUG_MESSAGES_
  std::cout<<" ATOM: eval_atomic_bonds "<<std::endl;
#endif
#ifdef _SHOW_TIME_
  gl_atom_clock.start();
#endif
  int j;
  bool use_pbc;
  //real dx, r2;
  real rl, r, r2, rr, ri, rj, rlz, rlxy;
  TVector<uint> vidx(2);
  TVector<real> vi, vj, vij, vang(2);
  TVector<real> vi_uvw, vj_uvw, vij_uvw;
  TVector<uint> v_ft;
  v_ft = supercell.get_fragmol_fragment_table();
  //
  int u_icell;
  uint k;
  TVector<real> v_positive_r;
  TVector<int>  v_integer_r;
  TVector<int>  v_neighbor_cell;
  TVector<int>  v_pbc;
  //
  uint max_bonds = 15*(__number_of_atoms);
#ifdef _ATOM_DEBUG_MESSAGES_
  std::cout<<" ATOM: estimated max bonds "<<max_bonds<<std::endl;
#endif
  //v_bond_number.push_back(check_bond(k));
  //m_bond_indices.add_row(vidx);
  //
  v_bond_table.resize(0);
  //v_bond_length.resize(0);
  i_number_of_bonds = 0;
  i_number_of_bonds_pbc = 0;
  ////////////////////////////////////////
  //v_bond_number.resize(0);
  //m_bond_indices.resize(0,2);
  v_bond_number.resize(max_bonds);
  m_bond_indices.resize(max_bonds,2);
  //v_bond_number_pbc.resize(0);
  //m_bond_indices_pbc.resize(0,2);
  v_bond_number_pbc.resize(uint(max_bonds/2));
  m_bond_indices_pbc.resize(uint(max_bonds/2),2);
  m_bond_boundary_pbc.resize(uint(max_bonds/2),2);
  ////////////////////////////////////////
  //m_bond_angles.resize(0,2);
  //m_bond_position.resize(0,3);
  //
  //m_bond_angles_pbc.resize(0,2);
  //m_bond_position_pbc.resize(0,3);
  //m_bond_boundary_pbc.resize(0,3);
  v_pbc.resize(3);
  //
#ifdef _ATOM_DEBUG_BONDS_
  std::cout<<" FL_GL_ATOM: LINKED: Box = "<<v_box_size;
#endif
  if(is_linked_cell){
#ifdef _ATOM_DEBUG_BONDS_
    std::cout<<" FL_GL_ATOM: LINKED CELL USED"<<std::endl;
#endif
    eval_linked_list();
#ifdef _ATOM_DEBUG_BONDS_
    std::cout<<" FL_GL_ATOM: LINKED CELL READY"<<std::endl;
#endif
  }
#ifdef _ATOM_DEBUG_BONDS_
  else{
    std::cout<<" FL_GL_ATOM: LINKED CELL NOT USED"<<std::endl;
  }
#endif
  for(int i=0; i<__number_of_atoms-1; i++){
#ifdef _ATOM_DEBUG_BONDS_
    std::cout<<" ATOM: i="<<i<<std::endl;
#endif
    if(strcmp(v_atomic_symbol_table_gl[v_atom_table[i]].c_str(),"X")){
      vi = m_atom_coordinates[i];
      vi_uvw = (vi*u_inv_bbox);
      ri = m_radius_color[i][0];
      if(is_linked_cell){
        v_positive_r = fVAdd(vi_uvw,v_bbox);
        v_integer_r  = iVMul(v_positive_r,v_cell_frac);
        // using the neighbour list
        // searching inside the neighbour cells
        //for (int _m=0; _m<27; _m++){
        for (int _m=0; _m<i_neighbor_cells; _m++){
          //use_pbc = true;
          // Neighbour cell
          //v_neighbor_cell = iVAdd(v_integer_r,neighbor_cells[_m]);
          //std::cout<<" Neighbouring cells = "<<i_neighbor_cells<<std::endl;
          //std::cout<<" Neighbouring matrix = "<<neighbor_cells_xyz<<std::endl;
          v_neighbor_cell = iVAdd(v_integer_r,neighbor_cells_xyz[_m]);
          // Used to apply PCB to cells in each dimension
          /////////////////////////////////////////////////////////////////////////////////
          for (uint coord=0; coord<3; coord++){
            if(v_neighbor_cell[coord] >= v_cell_side[coord]){ // check if  PBC is necessary
              v_neighbor_cell[coord]= 0;                      // apply PBC to each  cell
              //if(is_pbc) use_pbc = false;
            }else if(v_neighbor_cell[coord] < 0){             // check if  PBC is necessary
              v_neighbor_cell[coord]= v_cell_side[coord]-1;   // apply PBC to each cell
              //if(is_pbc) use_pbc = false;
            }
          }
          /////////////////////////////////////////////////////////////////////////////////
          //if(use_pbc){
          u_icell = iVLinear(v_neighbor_cell,v_cell_side);    // head atom index in the the cell
          if(u_icell>=0 && u_icell < (int)u_cell_number){          // inside of a cells
            j = v_cell_head[u_icell];                         // head atom in the actual cell
          }else{                                              // out of the box
            j = -1;                                           // outside of a cells
          }
          while(1){                                           // over all the particles in the cell
            if(j<0) break;                                    // stop searching in the cell
            //if(_m!=0 || j>i){                               // avoid self-interaction
            if((j>i) && (v_ft[j] == v_ft[i])){                // avoid self-interaction and double bond
            //if(j>i){                // avoid self-interaction and double bond
#ifdef _ATOM_DEBUG_BONDS_
              std::cout<<" FL_GL_ATOM: test bond for [i-j]="<<i<<","<<j<<std::endl;
#endif
              r2 = 0;                                         // set distance to cero
              vj = m_atom_coordinates[j];
              vj_uvw = vj*u_inv_bbox;
              rj = m_radius_color[j][0];
              r = (ri+rj);
              rr = 1.2*(r*r);
              use_pbc = false;
              vij_uvw = (vj_uvw-vi_uvw);
              for(uint coord=0; coord<3; coord++){
                v_pbc[coord] = 0;
                if(vij_uvw[coord] <= -v_bbox[coord]){
                  vj += 2.0*m_bbox[coord];       // PBC
                  use_pbc = true;
                  v_pbc[coord] = 1;
                }else if(vij_uvw[coord] > v_bbox[coord]){
                  vj -= 2.0*m_bbox[coord];       // PBC
                  use_pbc = true;
                  v_pbc[coord] = -1;
                }
              }
              vij = (vj-vi);
              r2 = vij.norm();
              if(r2 < rr){                                    // atoms inside de cut radius
#ifdef _ATOM_DEBUG_BONDS_
                std::cout<<" FL_GL_ATOM: Bond found [i-j]="<<i<<","<<j<<std::endl;
#endif
                vidx[0]=i;
                vidx[1]=j;
                r = sqrt(r2);
                vij = (vj-vi);
                rlz = vij[2];
                vij[2] = 0;
                rlxy = vij.magnitude();
                vang[0] = RAD_DEG*atan2(vij[1],vij[0]);       // bond precession
                vang[1] = RAD_DEG*fabs(atan2(rlxy,rlz));      // bond tilt
                //vij = 0.5*(vj+vi);
                //k = (uint)(r*f_atom_bond_inv_delta);
                k = get_bond_index(r);
                //m_bond_indices.add_row(vidx);
                if(use_pbc){
                  v_bond_number_pbc[i_number_of_bonds_pbc]=check_bond(k);
                  m_bond_indices_pbc[i_number_of_bonds_pbc]=vidx;
                  m_bond_boundary_pbc[i_number_of_bonds_pbc]=v_pbc;
                  i_number_of_bonds_pbc++;
                  //m_bond_angles_pbc.add_row(vang);
                  //m_bond_position_pbc.add_row(vij);
                  //v_bond_number_pbc.push_back(check_bond(k));
                  //m_bond_indices_pbc.add_row(vidx);
                  //m_bond_boundary_pbc.add_row(v_pbc);
                }else{
                  v_bond_number[i_number_of_bonds]=check_bond(k);
                  m_bond_indices[i_number_of_bonds]=vidx;
                  i_number_of_bonds++;
                  //v_bond_number.push_back(check_bond(k));
                  //m_bond_indices.add_row(vidx);
                  //count_bonds++;
                  //m_bond_angles.add_row(vang);
                  //m_bond_position.add_row(vij);
                }
                //i_number_of_bonds++;
              }
            }
            j = v_cell_list[j];
          }
        }
      }else{
#ifdef _ATOM_DEBUG_BONDS_
        std::cout<<" FL_GL_ATOM: NOT LINKED CELL USED"<<std::endl;
#endif
        // the code below can be used for small number of atoms.
        // searching inside the neighbour cells
        for (int _m=0; _m<27; _m++){
          for(int j=i+1; j<__number_of_atoms; j++){
            if(v_ft[j] == v_ft[i]){  // avoid bonds between fragments
#ifdef _ATOM_DEBUG_BONDS_
            std::cout<<" ATOM: j="<<j<<" i="<<i<<" m="<<_m<<std::endl;
#endif
            if(strcmp(v_atomic_symbol_table_gl[v_atom_table[j]].c_str(),"X")){
              vj = m_atom_coordinates[j];
              for(uint coord=0; coord<3; coord++){
                //vj[coord]+=(neighbor_cells[_m][coord]*v_box_size[coord]);                                        // PBC
                //vj[coord]+=(neighbor_cells[_m][coord]*2.0*(m_bbox[0][coord]+m_bbox[1][coord]+m_bbox[2][coord])); // PBC
                vj += neighbor_cells[_m][coord]*2.0*m_bbox[coord];       // PBC
                v_pbc[coord] = neighbor_cells[_m][coord];
              }
              rj = m_radius_color[j][0];
              vj_uvw = vj*u_inv_bbox;
              vij = (vj-vi);
              r = (ri+rj);
              rr = 1.2*(r*r);
              rl = vij.norm();
              if( (rl <= rr) &&  (rl > 0.1) ){
#ifdef _ATOM_DEBUG_BONDS_
                std::cout<<" ATOM: j="<<j<<" i="<<i<<" bonded"<<std::endl;
#endif
                vidx[0]=i;
                vidx[1]=j;
                rlz = vij[2];
                vij[2] = 0;
                rlxy = vij.magnitude();
                vang[0] = RAD_DEG*atan2(vij[1],vij[0]);  // precession
                vang[1] = RAD_DEG*fabs(atan2(rlxy,rlz)); // tilt
                //vij = 0.5*(vj+vi);
                r = sqrt(rl);
                //uint k = (uint)(r/f_atom_bond_delta);
                uint k = get_bond_index(r);
                //
                //m_bond_angles.add_row(vang);
                //m_bond_position.add_row(vij);
                //v_bond_number.push_back(check_bond(k));
                //m_bond_indices.add_row(vidx);
                if(_m>0){
                  v_bond_number_pbc[i_number_of_bonds_pbc]=check_bond(k);
                  m_bond_indices_pbc[i_number_of_bonds_pbc]=vidx;
                  m_bond_boundary_pbc[i_number_of_bonds_pbc]=v_pbc;
                  i_number_of_bonds_pbc++;
                  //count_bonds++;
                  //m_bond_boundary_pbc.add_row(v_pbc);
                  //v_bond_number_pbc.push_back(check_bond(k));
                  //m_bond_indices_pbc.add_row(vidx);
                  //m_bond_angles_pbc.add_row(vang);
                  //m_bond_position_pbc.add_row(vij);
#ifdef _ATOM_DEBUG_BONDS_
                  std::cout<<" ATOM: bond pcb type: "<<check_bond(k)<<std::endl;
#endif
                }else{
                  v_bond_number[i_number_of_bonds]=check_bond(k);
                  m_bond_indices[i_number_of_bonds]=vidx;
                  i_number_of_bonds++;
                  //v_bond_number.push_back(check_bond(k));
                  //m_bond_indices.add_row(vidx);
                  //count_bonds++;
                  //m_bond_angles.add_row(vang);
                  //m_bond_position.add_row(vij);
#ifdef _SHOW_DEBUG_NOPBC_BONDS_
                  std::cout<<" ATOM: bond type: "<<check_bond(k)<<std::endl;
#endif
                }
                //i_number_of_bonds++; //
              }
            }
            }
          }
        }
      }
    }
    //i_number_of_bonds+=v_bond_number.size();
  }
  //i_number_of_bonds=v_bond_number.size()+v_bond_number_pbc.size();
  //i_number_of_bonds=v_bond_number.size();
  //i_number_of_bonds--;
  //i_number_of_bonds_pbc=v_bond_number_pbc.size();
  m_bond_rcolor_0.resize(i_number_of_bonds,4);
  m_bond_rcolor_1.resize(i_number_of_bonds,4);
  m_bond_rcolor_pbc_0.resize(i_number_of_bonds_pbc,4);
  m_bond_rcolor_pbc_1.resize(i_number_of_bonds_pbc,4);
  //
  m_bond_angles.resize(i_number_of_bonds,2);
  m_bond_angles_pbc.resize(i_number_of_bonds_pbc,2);
  m_bond_position.resize(i_number_of_bonds,3);
  m_bond_position_pbc.resize(i_number_of_bonds_pbc,3);
  //std::cout<<" m_bond_indices="<<m_bond_indices;
#ifdef _ATOM_DEBUG_BONDS_
  std::cout<<" ATOM: i_number_of_bonds="<<i_number_of_bonds<<std::endl;
  std::cout<<" ATOM: i_number_of_bonds_pbc="<<i_number_of_bonds_pbc<<std::endl;
  std::cout<<" ATOM: v_bond_number="<<v_bond_number;
  std::cout<<" ATOM: v_bond_number_pbc="<<v_bond_number_pbc;
  std::cout<<" ATOM: v_bond_table="<<v_bond_table;
  std::cout<<" ATOM: m_bond_position_pbc="<<m_bond_position_pbc;
#endif
  u_bond_types=v_bond_table.size();
#ifdef _ATOM_DEBUG_MESSAGES_
  std::cout<<" ATOM: max number of bonds "<<max_bonds<<std::endl;
  std::cout<<" ATOM: number of bonds "<<i_number_of_bonds<<std::endl;
  std::cout<<" ATOM: number of PBC bonds "<<i_number_of_bonds_pbc<<std::endl;
  std::cout<<" ATOM: eval_atomic_bonds "<<std::endl;
#endif
  update_atomic_bonds();
#ifdef _SHOW_TIME_
  gl_atom_clock.stop();
  gl_atom_clock.show();
#endif
  // send update color flag
  //update_bonds_color=true;
}

void Fl_Gl_Atom::set_fragment_total(uint u){
  __fragment_total=u;
  palette.set(u);
  palette.set_color(4);
  palette.initialize(0,u,u);
  palette.update_palette_real();
}

void Fl_Gl_Atom::set_fragment_table(const TVector<uint>& v){
  v_fragment_table_gl = v;
#ifdef _SHOW_DEBUG_FRAGMENTS_
  std::cout<<" FL_GL_ATOM: fragment table: "<<v_fragment_table_gl;
#endif
}

void Fl_Gl_Atom::update_fragments(uint _u, bool _sw){
  set_fragment_total(get_view_total_fragments());
  set_fragment_table(get_view_fragment_table());
  //__fragment_active=v_fragment_table_gl[_u];
  if(_sw)
    set_active_fragment_index(v_fragment_table_gl[_u]);
  else
    set_active_fragment_index(_u);
#ifdef _ATOM_DEBUG_MESSAGES_
  std::cout<<" total fragments: "<<__fragment_total<<std::endl;
#endif
  set_update_coordinates(true);
  // fragments are counted from 1
  set_map_active_fragment(__fragment_active-1);
  //set_fragment_active(_af); // the same as zero above.
  is_eval_bonds=true;
  is_update_bonds=true;
  //update_bonds_color=true;
  update_data();
}

void Fl_Gl_Atom::compute_atom_fragment(const uint _u){
  //supercell.eval_atom_fragment(_u);
  supercell.eval_scaled_fragment(_u,0.1);
  // Use fragment table
  update_fragments(_u,true);
}

void Fl_Gl_Atom::compute_radial_fragment(const uint _u, const real _r){
  supercell.eval_radial_fragment(_u,true,_r);
  // Use fragment table
  update_fragments(_u,true);
}

void Fl_Gl_Atom::compute_vdw_fragment(const uint _u){
  //supercell.eval_vdw_fragment(_u);
  supercell.eval_scaled_fragment(_u,1.1);
  update_fragments(_u,true);
}

void Fl_Gl_Atom::compute_atom_fragments(void){
  //unsigned int _n;
  supercell.eval_scaled_fragments(0.1);
  // Use fragment number
  update_fragments(1,false);
  //set_fragment_total(get_view_total_fragments());
#ifdef _ATOM_DEBUG_MESSAGES_
  std::cout<<" total fragments: "<<__fragment_total<<std::endl;
#endif
  //set_fragment_table(get_view_fragment_table());
}

void Fl_Gl_Atom::compute_vdw_fragments(void){
  //unsigned int _n;
  supercell.eval_scaled_fragments(1.1);
  // Use fragment number
  update_fragments(1,false);
  //set_fragment_total(get_view_total_fragments());
#ifdef _ATOM_DEBUG_MESSAGES_
  std::cout<<" total fragments: "<<__fragment_total<<std::endl;
#endif
  //set_fragment_table(get_view_fragment_table());
}

void Fl_Gl_Atom::compute_merge_fragments(const uint _u){
  //supercell.eval_atom_fragment(_u);
  supercell.eval_merge_fragment(_u,false);
  // Use fragment table
  update_fragments(_u,true);
}

void Fl_Gl_Atom::set_active_fragment(const uint _a){
  uint _af;
  // atoms are counted from 1 in the scene
  _af= v_fragment_table_gl[_a];
  __fragment_active=_af;
  // fragments are counted from 1
  set_map_active_fragment(_af-1);
  //set_fragment_active(_af);
  //update_coordinates=true;
  set_update_coordinates(true);
  update_data(); //<-------------------------
}

/*
void Fl_Gl_Atom::set_atom_fragment(const uint _a){
  //uint _af;
  //eval_new_fragment(_a);
  //compute_new_fragment(_a);
  // atoms are counted from 1 in the scene
  //_af= v_fragment_table_gl[_a];
  v_fragment_table_gl[_a]=__fragment_total;
  __fragment_active=__fragment_total;
  set_fragment_total(__fragment_total+1);
  set_update_coordinates(true);
  // fragments are counted from 1
  //set_map_active_fragment(__fragment_total-1);
  //set_fragment_active(_af); // the same as zero above.
  //update_data();
}*/

uint Fl_Gl_Atom::check_bond(uint u){
  uint size = v_bond_table.size();
  for(uint i=0; i<size; i++){
    if(v_bond_table[i]==u){
      return i;
    }
  }
  v_bond_table.push_back(u);
  return size;
}

void Fl_Gl_Atom::eval_sphere(uint maxlevel){
    u_sphere_rows = 1 << maxlevel;
    int s, cont = 0;
    __sphere_strip_length=20*(u_sphere_rows*(u_sphere_rows-1)+(u_sphere_rows*3));
    m_sphere.resize(__sphere_strip_length,3);
    //m_sphere.resize(0,3);
    TVector<real> vt(3);
    /* iterate over the 20 sides of the icosahedron */
    for(s = 0; s < 20; s++) {
        int i;
        triangle *t = (triangle *)&icosahedron[s];
        for(i = 0; i < u_sphere_rows; i++) {
            /* create a tstrip for each row */
            /* number of triangles in this row is number in previous +2 */
            /* strip the ith trapezoid block */
            point v0, v1, v2, v3, va, vb;
            int j;
            linearly_interpolate(&t->pt[1], &t->pt[0], (float)(i+1)/u_sphere_rows, &v0);
            linearly_interpolate(&t->pt[1], &t->pt[0], (float)i/u_sphere_rows, &v1);
            linearly_interpolate(&t->pt[1], &t->pt[2], (float)(i+1)/u_sphere_rows, &v2);
            linearly_interpolate(&t->pt[1], &t->pt[2], (float)i/u_sphere_rows, &v3);
            //glBegin(GL_TRIANGLE_STRIP);

            NORMV(v0,cont);
            cont++;
            NORMV(v1,cont);
            cont++;
            for(j = 0; j < i; j++) {
                /* calculate 2 more vertices at a time */
                linearly_interpolate(&v0, &v2, (float)(j+1)/(i+1), &va);
                linearly_interpolate(&v1, &v3, (float)(j+1)/i, &vb);
                NORMV(va,cont);
                cont++;
                NORMV(vb,cont);
                cont++;
            }
            NORMV(v2,cont);
            cont++;
            //std::cout<<s<<") "<<__sphere_strip_length<<" + "<<(2*i+3)<<" = ";
            //__sphere_strip_length+=(2*i+3);
            //std::cout<<__sphere_strip_length<<std::endl;
        }
    }
    //m_atoms_strip.resize(__sphere_strip_length,3);
#ifdef _ATOM_DEBUG_MESSAGES_
  std::cout<<" FL_GL_ATOM: sphere rows = "<<u_sphere_rows<<std::endl;
  std::cout<<" FL_GL_ATOM: sphere_strip_length = "<<__sphere_strip_length<<std::endl;
#endif
}

// cylinder with unitary length
void Fl_Gl_Atom::eval_cylinder(uint n){
  //real theta1,theta2;
  real theta3;
  TVector<real> e(3),t1(2);
  m_cylinder_e1.resize(0,3);
  m_cylinder_texture1.resize(0,2);
  m_cylinder.resize(0,3);
  __cylinder_strip_length=0;
  // check resolution
  if(n < 4){
      n = 4;
  }
  //else if(n < 0){
  //    n = -n;
  //}
  __cylinder_strip_length=n;
  for (uint i=0;i<=n;i++){
    theta3 = i * C_2PI / n;
    e[0] = cos(theta3);
    e[1] = sin(theta3);
    e[2] = 0;
    m_cylinder_e1.add_row(e); // still used
    m_cylinder.add_row(e);
    t1[0] = i/(real)n;
    t1[1] = 1;
    m_cylinder_texture1.add_row(t1);
  }
  m_cylinder.add_row(m_cylinder[n-1]);
  __cylinder_strip_length++;
}

void Fl_Gl_Atom::initialize_transform_matrix(void){
  //tb_button = button;
  tb_angle = 0.0;
  glMatrixMode(GL_PROJECTION);
  // put the identity in the trackball transform
  glPushMatrix();
  glLoadIdentity();
  glGetFloatv(GL_MODELVIEW_MATRIX, (GLfloat*)tb_transform);
  glPopMatrix();
}

void Fl_Gl_Atom::initialize_rotation_matrix(void){
  //tb_button = button;
  tb_angle = 0.0;
  //glMatrixMode(GL_PROJECTION);
  glMatrixMode(GL_MODELVIEW);
  /* put the identity in the trackball transform */
  glPushMatrix();
  glLoadIdentity();
  glGetFloatv(GL_MODELVIEW_MATRIX,(GLfloat*)rot_matrix);
  glPopMatrix();
  is_initialize_rot=false;
}

// linearly interpolate between a & b, by fraction f
void Fl_Gl_Atom::linearly_interpolate(point *a, point *b, float f, point *r) {
    r->x=a->x+f*(b->x-a->x);
    r->y=a->y+f*(b->y-a->y);
    r->z=a->z+f*(b->z-a->z);
}

// normalize_point point r
void Fl_Gl_Atom::normalize_point(point *r) {
    float mag;
    mag = r->x * r->x + r->y * r->y + r->z * r->z;
    if (mag != 0.0f) {
        mag = 1.0f / sqrt(mag);
        r->x *= mag;
        r->y *= mag;
        r->z *= mag;
    }
}

// END
