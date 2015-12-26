# xmolview
A graphical editor for molecular dynamics and ab-initio simulations

http://www.molecular-explorer.com

To compile xmolview the following two dependencies are needed:

Get the GPL open source MSMVTL template library from: https://github.com/eetorres/msmvtl, and add it to the include directory.

Get the GPL open source FLTK library from: http://www.fltk.org.

On Linux compile with:

$ ./configure

On Mac OS X

$ ./configure --enable-macos=yes

On Windows

$ ./configure  --enable-windows=yes

Then build the binary

$ make

run the executable

$ ./src/xmolview

You may also use the OS dependent Makefiles. Edit the makeinclude file according to your system.

=====

Molecular Explorer is a powerful WYSIWYG OpenGL graphical editor with several unique features to visualize structures for atomistic calculations. We are very close to a stable beta release.

Get the source code from:

https://github.com/eetorres/xmolview

Here a short list of features:

+ Load file from CLI.

+ Generate DL_FIELD molecular structure with DL_POLY format.

+ Periodic and non-periodic systems.

+ Uni/two-dimensional structure scans.

+ Fragment manipulation.

+ Based on the FLTK library therefore very fast and small.

+ Automatic fragmentation of non-bonded parts for counterpoise calculations.

+ Automatic molecular integrity recognition of split structures due periodic boundary conditions.

+ Read several standard files, such as: xyz, pdb, Gaussian, VASP, DL_POLY.

+ Convert structure files between the supported formats.

+ Build periodic systems using images of the unit cell.

+ Series of structures for Potential Energy Surfaces (PES) calculations.

+ Fragment definition to manipulate specific parts of the atomic structures.

+ Display labels for the atomic species (e.g., DL_POLY).

+ Distance and angle tools.

+ Keyboard shortcuts. 
