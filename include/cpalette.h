// Latest modification: Wed Sep  3 14:13:54 EDT 2014
// version: 0.1
// name:    cpalette.h
// Copyrigth 2008 by Edmanuel Torres A. (eetorres@gmail.com)
//========================================================================
// FILE - cpalette.h                                                    //
// For the Fast Light Tool Kit (FLTK) - www.fltk.org                    //
//========================================================================
//                                                                      //
//                                                                      //
//                                                                      //
// Copyright 2008-2014 by Edmanuel Torres                               //
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

#ifndef _CPalette_H_
#define _CPalette_H_

#include <msmvtl/tmatrix.h>
#include <msmvtl/tmmath.h>

typedef struct {
  double r,g,b;
} gm_rgb;

typedef struct {
  unsigned int r,g,b;
} ui_rgb;

class CPalette{

public:
    //
    unsigned int xcells, ycells, _lvls, _ipalette;
    double xdelta1, ydelta1, xdelta2, ydelta2;
    double y0, y1, x0, x1, x2, x3, x4, x5;
    double w0, dz, zmin, zmax;
    gm_rgb rgb, cf1, cf2;
    ui_rgb id_rgb;

    std::vector<gm_rgb> _color_palette;
    std::vector<ui_rgb> _index_palette;
    TMatrix<real> m_color_palette;

    CPalette();
    ~CPalette(){};
    // set
    void set_color(unsigned int);
    // get
    gm_rgb get_color(unsigned int);
    ui_rgb get_index(unsigned int);
    unsigned int get_index_rgb(ui_rgb);
    unsigned int get_index_rg(ui_rgb);
    TVector<real> get_vcolor(double);
    //
    void initialize(double mn, double mx, unsigned int l){
      zmin=mn;
      zmax=mx;
      _lvls=l;
    };

//private:
    //
    gm_rgb  index_palette_(double val);
    gm_rgb  linear_palette_(double val);
    gm_rgb  hsv_palette_(double val);
    gm_rgb  rgb_palette_(double val);
    gm_rgb  earth_palette_(double val);
    gm_rgb  terrain_palette_(double val);
    //
    gm_rgb  palette_color_(double val);
    double color_interpolation_(double,double,double);
    //
    // Color functions
    void update_palette_real(void);
    void update_palette_index(void);
    void set(double w);
};

#endif

////


