#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([xmolview], [15.12], [bugs@xmol.org])
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

# Check for particular build options
AC_ARG_ENABLE([linux],
   AS_HELP_STRING([--enable-linux], [compile for linux [[yes]]]),
   [
      case "${enableval}" in
         yes) BUILD_FOR_LINUX=yes ;;
         no) BUILD_FOR_LINUX=no ;;
         *) AC_MSG_ERROR(bad value ${enableval} for --enable-linux) ;;
      esac
   ], [BUILD_FOR_LINUX=yes])  dnl Default value

AC_ARG_ENABLE([windows],
   AS_HELP_STRING([--enable-windows], [compile for windows [[no]]]),
   [
      case "${enableval}" in
         yes) BUILD_FOR_WINDOWS=yes ;;
         no) BUILD_FOR_WINDOWS=no ;;
         *) AC_MSG_ERROR(bad value ${enableval} for --enable-windows) ;;
      esac
   ], [BUILD_FOR_WINDOWS=no])  dnl Default value

AC_ARG_ENABLE([macos],
   AS_HELP_STRING([--enable-macos], [compile for MacOS [[no]]]),
   [
      case "${enableval}" in
         yes) BUILD_FOR_MACOS=yes ;;
         no) BUILD_FOR_MACOS=no ;;
         *) AC_MSG_ERROR(bad value ${enableval} for --enable-macos) ;;
      esac
   ], [BUILD_FOR_MACOS=no])  dnl Default value

# Conditional setup for Windows build
if test "$BUILD_FOR_WINDOWS" = "yes" ; then
   BUILD_FOR_LINUX="no"
   BUILD_FOR_MACOS="no"
   OS="mingw32"
   AC_DEFINE([BUILD_FOR_WINDOWS],[1],[use Windows OS components])
fi

# Conditional setup for MacOS build
if test "$BUILD_FOR_MACOS" = "yes" ; then
   BUILD_FOR_LINUX="no"
   BUILD_FOR_WINDOWS="no"
   OS="macos"
#  LIBS=["$LIBS -framework GLKit"]
   AC_DEFINE([BUILD_FOR_MACOS],[1],[use Mac OS X components])
fi

# Conditional setup for Linux build
if test "$BUILD_FOR_LINUX" = "yes" ; then
   OS="linux"
#  SYS=$(hostname)
   AC_DEFINE([BUILD_FOR_LINUX],[1],[use Linux components])
fi

# Make all the proper substitutions
AC_SUBST(OS)
#AC_SUBST(SYS)

# Generating files
AC_CONFIG_FILES([Makefile
                 src/Makefile])

AC_OUTPUT
