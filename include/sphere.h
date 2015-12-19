//========================================================================
// FILE - sphere.h the widget class file
// for the Fast Light Tool Kit (FLTK) - www.fltk.org
//
// Last modification: Sun Dec 18 18:04:34 MST 2011
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

#ifndef _SPHERE_H_
#define _SPHERE_H_

#include <msmvtl/libdef.h>

typedef struct {
    real x, y, z;
} point;

typedef struct {
    point pt[3];
} triangle;

// six equidistant points lying on the unit sphere //
#define UNIT_SPHERE_XPLUS {  1,  0,  0 }    //  X //
#define UNIT_SPHERE_XMIN  { -1,  0,  0 }    // -X //
#define UNIT_SPHERE_YPLUS {  0,  1,  0 }    //  Y //
#define UNIT_SPHERE_YMIN  {  0, -1,  0 }    // -Y //
#define UNIT_SPHERE_ZPLUS {  0,  0,  1 }    //  Z //
#define UNIT_SPHERE_ZMIN  {  0,  0, -1 }    // -Z //

// for icosahedron
#define IH_CZ (0.89442719099991)   //  2/sqrt(5) //
#define IH_SZ (0.44721359549995)   //  1/sqrt(5) //

#define IH_C1 (0.951056516)        // cos(18),  //
#define IH_S1 (0.309016994)        // sin(18) //
#define IH_C2 (0.587785252)        // cos(54),  //
#define IH_S2 (0.809016994)        // sin(54) //
#define IH_X1 (IH_C1*IH_CZ)
#define IH_Y1 (IH_S1*IH_CZ)
#define IH_X2 (IH_C2*IH_CZ)
#define IH_Y2 (IH_S2*IH_CZ)

#define IH_Ip0     {   0.0,    0.0,  1.0  }
#define IH_Ip1     {-IH_X2, -IH_Y2,  IH_SZ}
#define IH_Ip2     { IH_X2, -IH_Y2,  IH_SZ}
#define IH_Ip3     { IH_X1,  IH_Y1,  IH_SZ}
#define IH_Ip4     {   0.0,  IH_CZ,  IH_SZ}
#define IH_Ip5     {-IH_X1,  IH_Y1,  IH_SZ}

#define IH_Im0     {-IH_X1, -IH_Y1, -IH_SZ}
#define IH_Im1     {   0.0, -IH_CZ, -IH_SZ}
#define IH_Im2     { IH_X1, -IH_Y1, -IH_SZ}
#define IH_Im3     { IH_X2,  IH_Y2, -IH_SZ}
#define IH_Im4     {-IH_X2,  IH_Y2, -IH_SZ}
#define IH_Im5     {   0.0,    0.0, -1.0  }

// vertices of a unit icosahedron
static const triangle icosahedron[20]= {
        // front pole
        { {IH_Ip0, IH_Ip1, IH_Ip2}, },
        { {IH_Ip0, IH_Ip5, IH_Ip1}, },
        { {IH_Ip0, IH_Ip4, IH_Ip5}, },
        { {IH_Ip0, IH_Ip3, IH_Ip4}, },
        { {IH_Ip0, IH_Ip2, IH_Ip3}, },
        // mid
        { {IH_Ip1, IH_Im0, IH_Im1}, },
        { {IH_Im0, IH_Ip1, IH_Ip5}, },
        { {IH_Ip5, IH_Im4, IH_Im0}, },
        { {IH_Im4, IH_Ip5, IH_Ip4}, },
        { {IH_Ip4, IH_Im3, IH_Im4}, },
        { {IH_Im3, IH_Ip4, IH_Ip3}, },
        { {IH_Ip3, IH_Im2, IH_Im3}, },
        { {IH_Im2, IH_Ip3, IH_Ip2}, },
        { {IH_Ip2, IH_Im1, IH_Im2}, },
        { {IH_Im1, IH_Ip2, IH_Ip1}, },
        // back pole
        { {IH_Im3, IH_Im2, IH_Im5}, },
        { {IH_Im4, IH_Im3, IH_Im5}, },
        { {IH_Im0, IH_Im4, IH_Im5}, },
        { {IH_Im1, IH_Im0, IH_Im5}, },
        { {IH_Im2, IH_Im1, IH_Im5}, },
};

#endif

// END SPHERE