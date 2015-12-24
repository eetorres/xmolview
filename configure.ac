#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([xmolview], [15.12], [bugs@xmo.org])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/xmolview.cxx])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_MAKE_SET

# Checks for libraries.
# FLTK
unset fltkconfig
AC_CHECK_PROGS([fltkconfig], [fltk-config13 fltk-config])

if test "x$fltkconfig" = x; then
    AC_MSG_ERROR([fltk-config not found])
fi

VER=$($fltkconfig --version)
MINOR=$(echo $VER | cut -d. -f2)
PATCH=$(echo $VER | cut -d. -f3)

if test "$MINOR" -lt 3 || test "$PATCH" -lt 3; then
    AC_MSG_ERROR([FLTK >= 1.3.3 not found, got $VER])
fi

CXXFLAGS=["$CXXFLAGS $($fltkconfig --cxxflags | sed 's@-I@-isystem @g')"]
LIBS=["$LIBS $($fltkconfig --use-forms --use-glut --use-gl --ldstaticflags)"]

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([sqrt strcspn strdup strpbrk strrchr strstr strtol])

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
