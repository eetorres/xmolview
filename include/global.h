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
#ifndef _GLOBAL_H_
#define _GLOBAL_H_

#include <msmvtl/libdef.h>
#include "config.h"

//

const bool COORD_FORMAT_CARTES = 0;
const bool COORD_FORMAT_DIRECT = 1;
const bool COORD_FORMAT_ZMATRX = 2;

const uint INPUT_FILE_TYPE_UNKNOWN = 100;
const uint INPUT_FILE_TYPE_VSP = 0;
const uint INPUT_FILE_TYPE_XYZ = 1;
const uint INPUT_FILE_TYPE_GAU = 2;
const uint INPUT_FILE_TYPE_PDB = 3;
const uint INPUT_FILE_TYPE_DLP = 4;
const uint INPUT_FILE_TYPE_ZMT = 5;

const uint INPUT_FORMAT_STD     = 0;
const uint INPUT_FORMAT_ATM_NFR = 1;
const uint INPUT_FORMAT_ATM_FRG = 2;
const uint INPUT_FORMAT_NAT_NFR = 3;
const uint INPUT_FORMAT_NAT_FRG = 4;
const uint INPUT_FORMAT_NAT_CHG = 5;
const uint INPUT_FORMAT_ATM_CHG = 6;

const bool VASP_FORMAT_CARTES = false;
const bool VASP_FORMAT_DIRECT = true;
//

const uint OUTPUT_FILE_TYPE_VSP = 0;
const uint OUTPUT_FILE_TYPE_XYZ = 1;
const uint OUTPUT_FILE_TYPE_GAU = 2;
const uint OUTPUT_FILE_TYPE_PDB = 3;
const uint OUTPUT_FILE_TYPE_DLP = 4;
const uint OUTPUT_FILE_TYPE_ZMT = 5;

const uint OUTPUT_FORMAT_STD     = 0;
const uint OUTPUT_FORMAT_ATM_NFR = 1;
const uint OUTPUT_FORMAT_ATM_FRG = 2;
const uint OUTPUT_FORMAT_NAT_NFR = 3;
const uint OUTPUT_FORMAT_NAT_FRG = 4;
const uint OUTPUT_FORMAT_NAT_CHG = 5;
const uint OUTPUT_FORMAT_ATM_CHG = 6;

#endif

// END

