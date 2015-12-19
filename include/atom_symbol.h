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
// Atomic Symbol.c

#ifndef _ATOM_SYMBOL_H_
#define _ATOM_SYMBOL_H_

#include<string>

const unsigned int periodic_table_atoms=110;

const std::string symbol[periodic_table_atoms] = {
"H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
"Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca",
"Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
"Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr",
"Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
"Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
"Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
"Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg",
"Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
"Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
"Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "X",
}; // last symbol is for a dummy atom

const std::string symbol_sorted[periodic_table_atoms] = {
"He", "Li", "Be", "Ne", "Na", "Mg", "Al", "Si", "Cl", "Ar",
"Ca", "Sc", "Ti", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
"Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Zr", "Nb",
"Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb",
"Te", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm",
"Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf",
"Ta", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi",
"Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "Np", "Pu",
"Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf",
"Db", "Sg", "Bh", "Hs", "Mt",
"H",  "B",  "C",  "N",  "O",  "F",  "P",  "S",  "K",  "V",
"Y",  "I",  "W",  "U",  "X",
}; // last symbol is for a dummy atom

const size_t symbol_t[periodic_table_atoms] = {
1, 2, 2, 2, 1, 1, 1, 1, 1, 2,
2, 2, 2, 2, 1, 1, 2, 2, 1, 2,
2, 2, 1, 2, 2, 2, 2, 2, 2, 2,
2, 2, 2, 2, 2, 2, 2, 2, 1, 2,
2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
2, 2, 1, 2, 2, 2, 2, 2, 2, 2,
2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
2, 2, 2, 1, 2, 2, 2, 2, 2, 2,
2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
2, 2, 2, 2, 2, 2, 2, 2, 2, 1,
}; // last symbol is for a dummy atom

enum atomic_number {
H=1,  He, Li, Be, B,  C,  N,  O,  F,  Ne,
Na, Mg, Al, Si, P,  S,  Cl, Ar, K,  Ca,
Sc, Ti, V,  Cr, Mn, Fe, Co, Ni, Cu, Zn,
Ga, Ge, As, Se, Br, Kr, Rb, Sr, Y,  Zr,
Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn,
Sb, Te, I,  Xe, Cs, Ba, La, Ce, Pr, Nd,
Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb,
Lu, Hf, Ta, W,  Re, Os, Ir, Pt, Au, Hg,
Tl, Pb, Bi, Po, At, Rn, Fr, Ra, Ac, Th,
Pa, U,  Np, Pu, Am, Cm, Bk, Cf, Es, Fm,
Md, No, Lr, Rf, Db, Sg, Bh, Hs, Mt, X,
};

#endif

