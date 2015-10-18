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

