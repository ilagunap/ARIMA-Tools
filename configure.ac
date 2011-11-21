# ARIMA Tools autoconf file  
AC_INIT([arimatools],[0.1],[ilaguna@purdue.edu])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

# Init Libtool	
LT_INIT

# Headers' configuration
AC_CONFIG_SRCDIR([src/arima_tools.C])
AC_CONFIG_HEADERS([config.h])

# Check for programs
AC_PROG_CXX
AC_PROG_INSTALL
AM_PROG_CC_C_O

# Define C++ language
AC_LANG([C++])

# Checks for header files.
AC_CHECK_HEADERS([stdio.h stdlib.h math.h])

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT

echo \
"-------------------------------------------------
 ${PACKAGE_NAME} Version ${PACKAGE_VERSION}

 Prefix: '${prefix}'
 Compiler: '${CXX}'
 CXXFLAGS: '${CXXFLAGS}'
 LIBS: '${LIBS}'
 LDFLAGS: '${LDFLAGS}'

 Now type 'make @<:@<target>@:>@'
   where the optional <target> is:
     all                - build all binaries
     install            - install everything
--------------------------------------------------"
