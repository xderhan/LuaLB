AC_PREREQ(2.53)

AC_INIT([lualb], [0.1-alpha1], [ngamsaad.waipot@gmail.com])

AC_CONFIG_SRCDIR([src/lua/lb-lua.c])
AM_CONFIG_HEADER(config.h)

AM_INIT_AUTOMAKE(AC_PACKAGE_NAME,AC_PACKAGE_VERSION)
AM_MAINTAINER_MODE

AC_PROG_CC
#AC_PROG_SWIG(1.3.26)
#AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_PROG_RANLIB

dnl
dnl  asserts
dnl

AC_ARG_ENABLE(asserts,
              [AC_HELP_STRING([--enable-asserts],
                              [turn on assert() and lb_assert() [default=yes]])],,
	      [enable_asserts="yes"])

GCC_WARNINGS="-Wall -Werror-implicit-function-declaration -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations"

changequote(,)dnl
if test "x$GCC" = "xyes"; then
  CFLAGS="$CFLAGS $GCC_WARNINGS"
fi

if test "x$enable_asserts" != "xyes"; then
  case " $CPPFLAGS " in
  *[\ \	]-DNDEBUG[\ \	]*) ;;
  *) CPPFLAGS="-DNDEBUG=yes $CPPFLAGS" ;;
  esac
  CPPFLAGS="-DLB_DISABLE_ASSERTS=yes $CPPFLAGS"
fi
changequote([,])dnl

dnl
dnl  useful headers && things
dnl

AC_HEADER_STDC

AC_C_CONST
AC_C_INLINE

dnl
dnl  MPI
dnl

AC_ARG_ENABLE(mpi,
              [AC_HELP_STRING([--enable-mpi],
                              [enable MPI parallelization [default=no]])],,
	      [enable_mpi="no"])

if test "x$enable_mpi" = "xyes"; then
  AC_CHECK_FUNC(MPI_Init,,
                AC_MSG_ERROR([failed to compile MPI program]))

  AC_DEFINE(LB_ENABLE_MPI, 1, [Define if you have & link an MPI library.])
fi

AM_CONDITIONAL(LB_ENABLE_MPI, [test "x$enable_mpi" = "xyes"])

dnl
dnl  HDF5
dnl

AC_ARG_WITH(hdf5,
            AC_HELP_STRING([--with-hdf5=DIRECTORY],
	                   [prefix, where HDF5 is installed]))

if test "x$with_hdf5" = "xyes" ; then
  AC_MSG_ERROR([
*** directory must be specified for --with-hdf5])
fi

if test "${with_hdf5+set}" = set ; then
  HDF5_CFLAGS="-I$with_hdf5/include"
  HDF5_LIBS="-L$with_hdf5/lib -lhdf5 -lz -lm"
else
  HDF5_CFLAGS=""
  HDF5_LIBS="-lhdf5 -lz -lm"
fi

saved_CPPFLAGS=$CPPFLAGS

CPPFLAGS="$CPPFLAGS $HDF5_CFLAGS"
AC_CHECK_HEADER(hdf5.h,,
  AC_MSG_ERROR([
************************************************************
**  HDF5 is required (see http://hdf.ncsa.uiuc.edu/HDF5)  **
************************************************************]))

saved_LIBS=$LIBS
LIBS="$LIBS $HDF5_LIBS"
AC_MSG_CHECKING([for HDF5 libraries])
AC_TRY_RUN([#include <hdf5.h>
 	    int main (int argc, char** argv) {
		hid_t id = H5Tcreate(H5T_COMPOUND, 8);
		return 0;
	    }],
	    AC_MSG_RESULT(ok),
	    [AC_MSG_RESULT(failed)
	    AC_MSG_ERROR([
************************************************************
**  failed to run HDF5 test program (check 'config.log')  **
************************************************************])])

if test "x$enable_mpi" = "xyes" ; then
  AC_CHECK_FUNCS(H5Pset_mpi H5Pset_fapl_mpio)
fi

CPPFLAGS=$saved_CPPFLAGS
LIBS=$saved_LIBS
 
AC_SUBST(HDF5_CFLAGS)
AC_SUBST(HDF5_LIBS)

dnl
dnl  Lua 5.0
dnl

AC_ARG_WITH(lua5,
            AC_HELP_STRING([--with-lua5=DIRECTORY],
	                   [prefix, where Lua 5 is installed]))

if test "x$with_lua5" = "xyes" ; then
  AC_MSG_ERROR([
*** directory must be specified for --with-lua5])
fi

if test "${with_lua5+set}" = set ; then
  LUA5_CFLAGS="-I$with_lua5/include"
  LUA5_LIBS="-L$with_lua5/lib -llualib -llua -lm"
else
  LUA5_CFLAGS=""
  LUA5_LIBS="-llualib -llua -lm"
fi

saved_CPPFLAGS=$CPPFLAGS

CPPFLAGS="$CPPFLAGS $LUA5_CFLAGS"
AC_CHECK_HEADER(lua.h,,
  AC_MSG_ERROR([
************************************************************
**        Lua 5 is required (see http://www.lua.org)      **
************************************************************]))

saved_LIBS=$LIBS
LIBS="$LIBS $LUA5_LIBS"
AC_MSG_CHECKING([for Lua 5 libraries])
AC_TRY_RUN([#include <lua.h>
 	    int main (int argc, char** argv) {
	        lua_State* l = lua_open();
		return 0;
	    }],
	    AC_MSG_RESULT(ok),
	    [AC_MSG_RESULT(failed)
	    AC_MSG_ERROR([
************************************************************
** failed to run Lua 5 test program (check 'config.log')  **
************************************************************])])

CPPFLAGS=$saved_CPPFLAGS
LIBS=$saved_LIBS
 
AC_SUBST(LUA5_CFLAGS)
AC_SUBST(LUA5_LIBS)

dnl
dnl  'lb.rgb' class
dnl

AC_ARG_ENABLE(rgb,
              [AC_HELP_STRING([--enable-rgb],
                              [build 'rgb' class [default=yes]])],,
	      [enable_rgb="yes"])

dnl
dnl  check for PNG library
dnl

PNG_LIBS=

if test "x$enable_rgb" = "xyes"; then
  AC_CHECK_LIB(png, png_read_info,
    [AC_CHECK_HEADER(png.h,PNG_LIBS="-lpng -lz -lm",enable_rgb=no)],
    AC_MSG_WARN(*** 'rgb' class will not be built (PNG library not found) ***), -lz -lm)
fi

AC_SUBST(PNG_LIBS)

if test "x$enable_rgb" = "xyes"; then
  AC_DEFINE(LB_ENABLE_RGB, 1, [Define if you want build 'rgb' class])
fi

AM_CONDITIONAL(LB_ENABLE_RGB, [test "x$enable_rgb" = "xyes"])

AC_CONFIG_FILES([
  Makefile
  doc/Makefile
  examples/Makefile
  src/Makefile
  src/core/Makefile
  src/graphics/Makefile
  src/BGK/Makefile
  src/Mix/Makefile
  src/Landau/Makefile
  src/lua/Makefile
])

AC_OUTPUT

dnl
dnl  amusements
dnl

echo "<=============================================================>
 Configuration summary :
<------------------------------------------------------------->

    COMPILER = ${CC}
    CFLAGS   = ${CFLAGS}
    CPPFLAGS = ${CPPFLAGS}
    LDFLAGS  = ${LDFLAGS}

         MPI = ${enable_mpi}

         RGB = ${enable_rgb}

 HDF5_CFLAGS = ${HDF5_CFLAGS}
   HDF5_LIBS = ${HDF5_LIBS}

 LUA5_CFLAGS = ${LUA5_CFLAGS}
   LUA5_LIBS = ${LUA5_LIBS}

     ASSERTS = ${enable_asserts}
      PREFIX = ${prefix}

<=============================================================>"
