AC_INIT([mvMAPIT],[0.1.1])

# Find the compiler and compiler flags used by R.
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CXX=`"${R_HOME}/bin/R" CMD config CXX`
CXX11=`"${R_HOME}/bin/R" CMD config CXX11`
BLAS_LIBS=`"${R_HOME}/bin/R" CMD config BLAS_LIBS`
FLIBS=`"${R_HOME}/bin/R" CMD config FLIBS`
LAPACK_LIBS=`"${R_HOME}/bin/R" CMD config LAPACK_LIBS`
PKG_CXX11STD=`"${R_HOME}/bin/R" CMD config CXX11STD`
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
AC_LANG(C++)
AC_PROG_CPP

# Search for a system library, in this case the zlib compression library (which
# contains the function 'deflate'). This sets the variable 'ac_cv_search_deflate'
# to what must be passed to ${PKG_LIBS}.
AC_SEARCH_LIBS(deflate, z, [], [AC_MSG_ERROR(The zlib library is required.)])
AC_CHECK_HEADERS(zlib.h, [], [AC_MSG_ERROR(The zlib library headers are required.)])

# Verify the correct src directory exists from where the file is being run.
AC_CONFIG_SRCDIR([src])

### Determine subdirectory source files
# Dynamically generate list of sources from subdirectories via shell
MAIN_SOURCES="$(cd src/ && ls *.cpp | tr '\n' ' ')"
SUBDIR_SOURCES="$(cd src/ && ls {mapit,mqs,gsm,tests,logging}/*.cpp | tr '\n' ' ')"
# Supply it as a variable
AC_SUBST(MAIN_SOURCES)
AC_SUBST(SUBDIR_SOURCES)

# Write the flags into the src/Makevars file.
AC_SUBST([BLAS_LIBS], ["${BLAS_LIBS}"])
AC_SUBST([FLIBS], ["${FLIBS}"])
AC_SUBST([LAPACK_LIBS], ["${LAPACK_LIBS}"])
AC_SUBST([PKG_CXX11STD], ["${PKG_CXX11STD}"])
AC_SUBST([PKG_CPPFLAGS], ["${PKG_CPPFLAGS}"])
AC_SUBST([PKG_LIBS], ["${LIBS} ${PKG_LIBS} ${ac_cv_search_deflate}"])
AC_SUBST([CXX11], ["${CXX11}"])
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

echo "
  --------------------------------------------------
  Configuration for ${PACKAGE_NAME} ${PACKAGE_VERSION}

    cppflags: ${CPPFLAGS} ${PKG_CPPFLAGS}
    cxxflags: ${PKG_CXX11STD} ${CXXFLAGS} ${BLAS_LIBS}
    libs:     ${PKG_LIBS} ${BLAS_LIBS} ${LAPACK_LIBS} ${FLIBS}

  --------------------------------------------------
"
