//========================================================================
// FILE: cpotcar.h
//
// This an utility program to manipulate and genarate structure files
//
// Copyright 2011-2015 by Edmanuel Torres
// email:   eetorres@gmail.com
//
//========================================================================

#ifndef _POTCAR_H_
#define _POTCAR_H_

#include<config_debug.h>
#include<time.h>
#include<string.h>
#include<msmvtl/tmatrix.h>
#include<msmvtl/tvmath.h>
#include<atom_symbol.h>
//#include<const.h>

class CPotcar{

public:

  CPotcar();
  ~CPotcar(){};

  // read/write files
  void clear(void);
  // POSCAR VASP file
  bool read_potcar(std::string,std::string);
  //
  void eval_atomic_numbers(void);
  // set functions

  // get function
  bool format();
  //
  uint get_total_species(void);
  //
  void set_symbols(TVector<std::string> _v);
  //
  std::string get_symbol(uint);
  uint get_atomic_number(uint);
  TVector<uint> get_atomic_numbers(void);
  TVector<std::string> get_symbols(void);

private:
  //
  uint __potcar_species;
  TVector<std::string>  v_atomic_symbols;
  TVector<uint> v_atomic_numbers;
  //
};

#endif
