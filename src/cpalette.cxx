// Latest modification: Wed Sep  3 14:13:54 EDT 2014
// version: 0.1
// name:    cpalette.cxx
// Copyrigth 2008 by Edmanuel Torres A. (eetorres@gmail.com)
//========================================================================
// FILE - cpalette.cxx                                                  //
// For the Fast Light Tool Kit (FLTK) - www.fltk.org                    //
//========================================================================
//                                                                      //
//                                                                      //
//                                                                      //
// Copyright 2002-2014 by Edmanuel Torres                               //
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

#include <cpalette.h>

//#define _SHOW_MESSAGE_ 1

CPalette::CPalette(void){
    _ipalette = 0;
}

void CPalette::set(double w){
  xdelta1 = w/w0;
  xdelta2 = xdelta1+xdelta1;
  x0=xdelta1;
  x1=xdelta1+x0;
  x2=xdelta1+x1;
  x3=xdelta1+x2;
  x4=xdelta1+x3;
  x5=xdelta1+x4;
  ydelta1 = 255.0;
}

void CPalette::set_color(unsigned int u){
  _ipalette = u;
  switch(u){
  case 0:
    w0 = 6.0;//Linear
    break;
  case 3:
    w0 = 6.0;//HSV
    break;
  default:
    w0 = 6.0;//HSV
    break;
  }
};

void CPalette::update_palette_real(void){
  double lv;
  _color_palette.clear();
  //dz=zmax/_lvls;
  set(zmax);
  _color_palette.resize(_lvls+1);
  m_color_palette.resize(_lvls+1,3);
#ifdef _SHOW_MESSAGE_
  printf("_lvls = %i\n",_lvls);
  //printf("dz = %f\n",dz);
#endif
  for(unsigned int _p=0; _p<(_lvls+1); _p++){
    lv = (_p);
    rgb = palette_color_(lv);
    rgb.r/=255.0;
    rgb.g/=255.0;
    rgb.b/=255.0;
    //printf("r=%f, g=%f, b=%f\n",rgb.r,rgb.g,rgb.b);
    _color_palette[_p] = rgb;
    m_color_palette[_p][0] = rgb.r;
    m_color_palette[_p][1] = rgb.g;
    m_color_palette[_p][2] = rgb.b;
  }
}

void CPalette::update_palette_index(void){
  _index_palette.clear();
  set(zmax);
  _index_palette.resize(_lvls+1);
#ifdef _SHOW_MESSAGE_
  printf("_lvls = %i\n",_lvls);
  printf("dz = %f\n",dz);
#endif
  id_rgb.r=0;
  id_rgb.g=0;
  id_rgb.b=255;
  for(unsigned int _p=0; _p<(_lvls+1); _p++){
    if ( _p == 100 ){
      id_rgb.r=0;
      id_rgb.b=0;
    }
    _index_palette[_p] = id_rgb;
    id_rgb.r++;
    if(id_rgb.r > 255)
    {
        id_rgb.r = 0;
        id_rgb.g++;
        if(id_rgb.g > 255){
          id_rgb.g = 0;
          id_rgb.b++;
        }
    }
  }
}

gm_rgb CPalette::palette_color_(double _x){
  switch(_ipalette){
  case 0:
    return index_palette_(_x);
    break;
  case 1:
    return terrain_palette_(_x);
    break;
  case 2:
    return earth_palette_(_x);
    break;
  case 3:
    return hsv_palette_(_x);
    break;
  default:
    return rgb_palette_(_x);
    break;
  }

}

gm_rgb CPalette::get_color(unsigned int c){
  //unsigned int c = (unsigned int)(val/dz);
  //unsigned int c = (unsigned int)val;
  return _color_palette[c];
}

ui_rgb CPalette::get_index(unsigned int c){
  return _index_palette[c];
}

unsigned int CPalette::get_index_rgb(ui_rgb u){
  uint index_=0;
  if(u.b > 0){
    index_=u.b*256;
  }
  if(u.g > 0){
    index_+=u.g*256;
  }
  if(u.r > 0){
    index_+=u.r;
  }
  return index_;
}

unsigned int CPalette::get_index_rg(ui_rgb u){
  uint index_=0;
  if(u.g > 0){
    index_+=u.g*256;
  }
  if(u.r > 0){
    index_+=u.r;
  }
  return index_;
}

TVector<real> CPalette::get_vcolor(double val){
  unsigned int c = (unsigned int)(val);
  return m_color_palette[c];
}

double CPalette::color_interpolation_(double dx, double dy, double _cl){
  return (dy/dx)*_cl;
}

//// Palettes

// HSV palette
gm_rgb CPalette::hsv_palette_(double val){
    double rcf1=100.0;
    double bcf2=200.0;
    double gcf3=255.0;
    ydelta1 = 255.0;
    if(val >= 0 && val <= x0){ //1
      rgb.r=rcf1;
      rgb.g=0.0;
      rgb.b=color_interpolation_(xdelta1, bcf2, val);
    }else if(val > x0 && val <= x1){ //2
      rgb.r=rcf1-color_interpolation_(xdelta1, rcf1, val-x0);
      rgb.g=0.0;
      rgb.b=bcf2;
    }else if(val > x1 && val <= x2){ // 3
      rgb.r=0.0;
      rgb.g=color_interpolation_(xdelta1, gcf3, val-x1);
      rgb.b=bcf2;
    }else if(val > x2 && val <= x3){ // 4
      rgb.r=0.0;
      rgb.g=gcf3;
      rgb.b=bcf2-color_interpolation_(xdelta1, bcf2, val-x2);
    }else if(val > x3 && val <= x4){ // 5
      rgb.r=color_interpolation_(xdelta1, rcf1, val-x3);
      rgb.g=gcf3;
      rgb.b=0.0;
    }else if(val > x4 && val <= x5){ // 6
      rgb.r=rcf1;
      rgb.g=gcf3-color_interpolation_(xdelta1, gcf3, val-x4);
      rgb.b=0.0;
    }
    return rgb;
}

// Earth palette
gm_rgb CPalette::earth_palette_(double val){
    cf1.r=250.0; cf1.g=250.0;  cf1.b=255.0;
    ydelta1 = 255.0;
    if(val >= 0 && val <= x0){ // 1
      rgb.r=0.0;
      rgb.g=0.0;
      rgb.b=color_interpolation_(xdelta1, cf1.b, val);
    }else if(val >= x0 && val <= x1){ // 2
      rgb.r=0.0;
      rgb.g=color_interpolation_(xdelta1, cf1.g, val-x0);
      rgb.b=0.0;
    }else if(val > x1 && val <= x2){ //3
      rgb.r=0.0;
      rgb.g=cf1.g;
     rgb.b=0.0;
    }else if(val > x2 && val <= x3){ //4
      rgb.r=color_interpolation_(xdelta1, cf1.r, val-x2);
      rgb.g=cf1.g;
      rgb.b=0.0;
    }else if(val > x3 && val <= x4){ //5
      rgb.r=0.5*cf1.r+0.25*color_interpolation_(xdelta1, cf1.r, val-x3);
      rgb.g=0.5*cf1.r+0.25*color_interpolation_(xdelta1, cf1.g, val-x3);
      rgb.b=0.5*cf1.r+0.25*color_interpolation_(xdelta1, cf1.b, val-x3);
    }else if(val > x4){ //6
      rgb.r=0.75*cf1.r+0.25*color_interpolation_(xdelta1, cf1.r, val-x4);
      rgb.g=0.75*cf1.g+0.25*color_interpolation_(xdelta1, cf1.g, val-x4);
      rgb.b=0.75*cf1.b+0.25*color_interpolation_(xdelta1, cf1.b, val-x4);
    }
    return rgb;
}

// Terrain palette
gm_rgb CPalette::terrain_palette_(double val){
    cf1.r=250.0; cf1.g=250.0;  cf1.b=255.0;
    ydelta1 = 255.0;
    if(val >= 0 && val <= x0){ // 1
      rgb.r=color_interpolation_(xdelta1, 184.0, val);
      rgb.g=color_interpolation_(xdelta1, 134.0, val);
      rgb.b=color_interpolation_(xdelta1, 11.2, val);
    }else if(val >= x0 && val <= x1){ // 2
      rgb.r=0.0;
      rgb.g=0.6*color_interpolation_(xdelta1, cf1.g, val-x0);
      rgb.b=0.0;
    }else if(val > x1 && val <= x2){ //3
      rgb.r=0.0;
      rgb.g=0.6*cf1.g;
      rgb.b=0.0;
    }else if(val > x2 && val <= x3){ //4
      rgb.r=34.0+0.3*color_interpolation_(xdelta1, cf1.r, val-x2);
      rgb.g=0.6*cf1.g;
      rgb.b=34.0+0.3*color_interpolation_(xdelta1, cf1.b, val-x2);;
    }else if(val > x3 && val <= x4){ //5
      rgb.r=0.5*cf1.r+0.25*color_interpolation_(xdelta1, cf1.r, val-x3);
      rgb.g=0.5*cf1.r+0.25*color_interpolation_(xdelta1, cf1.g, val-x3);
      rgb.b=0.5*cf1.r+0.25*color_interpolation_(xdelta1, cf1.b, val-x3);
    }else if(val > x4){ //6
      rgb.r=0.75*cf1.r+0.25*color_interpolation_(xdelta1, cf1.r, val-x4);
      rgb.g=0.75*cf1.g+0.25*color_interpolation_(xdelta1, cf1.g, val-x4);
      rgb.b=0.75*cf1.b+0.25*color_interpolation_(xdelta1, cf1.b, val-x4);
    }
    return rgb;
}




// RGB palette
gm_rgb CPalette::rgb_palette_(double val){
    cf1.r=250.0; cf1.g=250.0;  cf1.b=255.0;
    ydelta1 = 255.0;
    if(val >= 0 && val <= x0){ // 1
      rgb.r=100.0;
      rgb.g=100.0;
      rgb.b=color_interpolation_(xdelta1, cf1.b, val);
    }else if(val >= x0 && val <= x1){ // 2
      rgb.r=100.0;
      rgb.g=color_interpolation_(xdelta1, cf1.g, val-x0);
      rgb.b=cf1.b;
    }else if(val > x1 && val <= x2){ //3
      rgb.r=100.0;
      rgb.g=cf1.g;
      rgb.b=cf1.b-color_interpolation_(xdelta1, cf1.b, val-x1);
    }else if(val > x2 && val <= x3){ //4
      rgb.r=color_interpolation_(xdelta1, cf1.r, val-x2);
      rgb.g=cf1.g;
      rgb.b=100.0;
    }else if(val > x3 && val <= x4){ //5
      rgb.r=cf1.r;
      rgb.g=cf1.g-color_interpolation_(xdelta1, cf1.g, val-x3);
      rgb.b=100.0;
    }else if(val > x4){ //6
      rgb.r=cf1.r-color_interpolation_(xdelta1, cf1.r, val-x4);
      rgb.g=100.0;
      rgb.b=100.0;
    }
    return rgb;
}

// Linear palette
gm_rgb CPalette::linear_palette_(double val){
  double rydelta1, gydelta1, bydelta1;
  cf1.r=20.0;  cf1.g=20.0;  cf1.b=20.0;
  cf2.r=255.0; cf2.g=255.0; cf2.b=255.0;
  ydelta1 = 255.0;
  rydelta1 = cf2.r-cf1.r;
  gydelta1 = cf2.g-cf1.g;
  bydelta1 = cf2.b-cf1.b;
  if(val >= 0 && val <= x0){ // 1
    rgb.r=cf1.r+color_interpolation_(xdelta1, rydelta1, val);
    rgb.g=cf1.g+color_interpolation_(xdelta1, gydelta1, val);
    rgb.b=cf1.b+color_interpolation_(xdelta1, bydelta1, val);
  }
  return rgb;
}

// Index palette
gm_rgb CPalette::index_palette_(double val){
    if(val >= 0 && val <= x2){ //1
      rgb.r=color_interpolation_(x5, ydelta1, val);
      rgb.g=color_interpolation_(x2, ydelta1, val);
      rgb.b=255.0-color_interpolation_(x5, ydelta1, val);
    }else if(val > x2 && val <= x5){ // 4
      rgb.r=color_interpolation_(x5, ydelta1, val);
      rgb.g=255.0-color_interpolation_(x2, ydelta1, val-x2);
      rgb.b=255.0-color_interpolation_(x5, ydelta1, val);
    }
    return rgb;
}



/// RUBISH


