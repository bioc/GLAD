
AC_INIT([GLAD])

dnl Select include path for gsl
dnl AC_ARG_WITH([include-path], AC_HELP_STRING([--with-include-path], the location header files), [include_path=$withval])


dnl dnl Select library path for gsl
dnl AC_ARG_WITH([library-path], AC_HELP_STRING([--with-library-path], the location of library files), [library_path=$withval])

dnl AC_ARG_WITH([link-options], AC_HELP_STRING([--with-link-options], link options), [library_options=$withval])

dnl GLAD_LIBS_PATH=$library_path
dnl GLAD_INCLUDE_PATH=$include_path
dnl GLAD_LIBS_OPTIONS=$library_options

dnl echo Syntax example to pass compilation options:
dnl echo 'R CMD INSTALL --configure-args='--with-library-path="-L/bioinfo/local/gsl/lib/ -L/opt/SUNWspro/lib" --with-link-options="-lgsl -lRblas"' '
dnl echo ""

dnl if test [ -z "$GLAD_LIBS_OPTIONS" ] ; then
dnl echo "No link option has been given"
dnl echo Default will be used: -lgsl -lblas
dnl GLAD_LIBS_OPTIONS="-lgsl -lblas"
dnl fi

dnl if test [ -z "$GLAD_INCLUDE_PATH" ] ; then
dnl echo "No include path directory have been given"
dnl echo "default include path will be used"
dnl fi

dnl if test [ -z "$GLAD_LIBS_PATH" ] ; then
dnl echo "No library path directory have been given"
dnl echo "default library path will be used"
dnl fi


OS=$(uname -s)

case $OS in
     "Darwin")
	echo $OS
	CXXFLAGS="-DIS_MAC_OS"
	;;
      "Linux" | "SunOS")
	echo operating system: $OS
	# find gsl library
	: ${PKGCONFIG=pkg-config}
	if ${PKGCONFIG} --version >/dev/null 2>&1; then
		if ${PKGCONFIG} gsl; then
			CXXFLAGS="$(${PKGCONFIG} --cflags gsl)"
			LIBS="$(${PKGCONFIG} --libs gsl)"
		else
			pkg-config --cflags gsl

		fi
	fi
	;;
       *)
	"echo operating system: $OS"
	## BLAS.
	AC_ARG_WITH([blas],
	[AC_HELP_STRING([--with-blas],
	                [use system gsl BLAS library (if available)])],
			[BLAS_LIBS=$with_blas],
			[BLAS_LIBS="-lgslcblas"])

	PKG_CHECK_MODULES(GSL, [gsl >= 1.6], [HAVEGSL="yes"], [HAVEGSL="no"])
	if test "$HAVEGSL" = "yes"; then
		LIBS=$GSL_LIBS
	else 
		AM_PATH_GSL(1.6,[LIBS="${LIBS} -lgsl ${BLAS_LIBS} -lm"],[AC_MSG_ERROR([Cannot find Gnu Scientific Library >=1.6])])
	fi

	;;
esac


GLAD_LIBS_OPTIONS="$CXXFLAGS $LIBS"


dnl Substitution
dnl AC_SUBST(GLAD_LIBS_PATH)
AC_SUBST(GLAD_LIBS_OPTIONS)
dnl AC_SUBST(GLAD_INCLUDE_PATH)

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
