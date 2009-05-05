AC_INIT([locmod])

dnl Select include path for gsl
AC_ARG_WITH([include-path], AC_HELP_STRING([--with-include-path], the location header files), [include_path=$withval])


dnl if test [ -d $gsl_include_path ] ; then
dnl echo "gsl_include_path=$gsl_include_path"
dnl else
dnl echo $gsl_include_path does not exist
dnl fi

dnl Select library path for gsl
AC_ARG_WITH([library-path], AC_HELP_STRING([--with-library-path], the location of library files), [library_path=$withval])

AC_ARG_WITH([link-options], AC_HELP_STRING([--with-link-options], link options), [library_options=$withval])

GLAD_LIBS_PATH=$library_path
GLAD_INCLUDE_PATH=$include_path
GLAD_LIBS_OPTIONS=$library_options

if test [ -z "$GLAD_LIBS_OPTIONS" ] ; then
echo "No link option has been given"
echo Default will be used: -lgsl -lblas
GLAD_LIBS_OPTIONS="-lgsl -lblas"
fi

if test [ -z "$GLAD_INCLUDE_PATH" ] ; then
echo "No include path directoty have been given"
echo "default include path will be used"
fi

if test [ -z "$GLAD_LIBS_PATH" ] ; then
echo "No library path directoty have been given"
echo "default library path will be used"
fi

dnl Substitution
AC_SUBST(GLAD_LIBS_PATH)
AC_SUBST(GLAD_LIBS_OPTIONS)
AC_SUBST(GLAD_INCLUDE_PATH)

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
