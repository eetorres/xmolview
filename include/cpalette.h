// Last modification 21/03/2008
// ========================================================================
// version: 0.1
// name:    fl_pallete.h
//
// Copyrigth 2008-2015 by Edmanuel Torres A. (eetorres@gmail.com)
//
//========================================================================

#ifndef _CPalette_H_
#define _CPalette_H_

#include <msmvtl/tmatrix.h>
#include <msmvtl/tmmath.h>

typedef struct {
  double r,g,b;
} gm_rgb;

class CPalette{

public:
    //
    unsigned int xcells, ycells, _lvls, _ipalette;
    double xdelta1, ydelta1, xdelta2, ydelta2;
    double y0, y1, x0, x1, x2, x3, x4, x5;
    double w0, dz, zmin, zmax;
    gm_rgb rgb, cf1, cf2;

    std::vector<gm_rgb> _color_palette;
    TMatrix<real> m_color_palette;

    CPalette();
    ~CPalette(){};
    // set
    void set_color(unsigned int);
    // get
    gm_rgb get_color(double);
    TVector<real> get_vcolor(double);
	//
    void initialize(double mn, double mx, unsigned int l){
      zmin=mn;
      zmax=mx;
      _lvls=l;
    };

//private:
    //
    gm_rgb  hsv_palette_(double val);
    gm_rgb  rgb_palette_(double val);
    gm_rgb  linear_palette_(double val);
    gm_rgb  earth_palette_(double val);
    gm_rgb  terrain_palette_(double val);
    //
    gm_rgb  palette_selection_(double val);
    double color_interpolation_(double,double,double);
    //
    // Color functions
    void update_palette_(void);
    void set(double w);
};

#endif

////


