# what flags you want to pass to the C compiler & linker
#AM_CFLAGS = --pedantic -Wall -std=c99 -O2
#AM_LDFLAGS = 

# this lists the binaries to produce, the (non-PHONY, binary) targets in
# the previous manual Makefile
bin_PROGRAMS = xmolview
xmolview_SOURCES = xmolview.cxx fl_xmol_view.cxx fl_table.cxx fl_gl_mol_view.cxx fl_gl_atom.cxx \
      cviewmol.cxx cmapcell.cxx csupercell.cxx cfragment.cxx cfile.cxx ctopmol.cxx \
      cposcar.cxx cpotcar.cxx cxyz.cxx cgau.cxx czmat.cxx cpdb.cxx cdlp.cxx \
      sphere.cxx cpalette.cxx

#xmolview_LDFLAGS = $(shell fltk-config  --cxxflags --use-forms --use-glut --use-gl --ldstaticflags)