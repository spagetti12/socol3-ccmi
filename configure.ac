dnl Process this file with autoconf to produce a configure script.
dnl
dnl Luis Kornblueh August 31 1998
dnl Uwe Schulzweida November 4 2002
dnl  -add cross compiling for NEC (use --host=sx6)
dnl  -use autoconf version 2.52 
dnl  -add cross compiling for Earth Simulator, Yokohama (use --host=es)
dnl

AC_INIT(echam,5.3.02,schulzweida@dkrz.de,echam-5.3.02)
AC_PREREQ(2.59)

dnl Checks for package options.
dnl ######################################################################

AC_PREFIX_DEFAULT(.)

AC_CONFIG_AUX_DIR(config)

AC_CANONICAL_SYSTEM

dnl AM_INIT_AUTOMAKE

dnl ######################################################################

host_frag=""
case "${host}" in
    alpha-*-*)           host_frag=config/mh-alpha ;;
    5000-fujitsu-*)      host_frag=config/mh-fujitsu-5000 ;;
    *-fujitsu-*)         host_frag=config/mh-fujitsu ;;
    sr8000-hitachi-*)    host_frag=config/mh-hitachi-sr8000 ;;
    *90-cray-*)          host_frag=config/mh-cray-c90 ;;
    sv1-cray-*)          host_frag=config/mh-cray-c90 ;;
    nv1-cray-*)          host_frag=config/mh-cray-x1 ;;
    t3e-cray-*)          host_frag=config/mh-cray-t3e ;;
    xd1-cray-*)          host_frag=config/mh-cray-xd1 ;;
    xt3-cray-*)          host_frag=config/mh-cray-xt3 ;;
    mips*-sgi-irix6.*)   host_frag=config/mh-irix6 ;;
    *-*-solaris2*)       host_frag=config/mh-solaris ;;
    sx*-*-*)             host_frag=config/mh-nec ;;
    es*-*-*)             host_frag=config/mh-es;;
    *-ibm-aix*)          host_frag=config/mh-aix ;;
    ia64-*-linux-*)      host_frag=config/mh-linux-ia64 ;; 
    x86_64-*-linux-*)    host_frag=config/mh-linux-x64 ;;  
    i*86-*-linux-*)      host_frag=config/mh-linux ;;
    *-*-hpux*)           host_frag=config/mh-hpux ;;
    *) echo This configuration is not supported.
       echo Please set a valid config/mh-* file up.
       ;;	 	
esac

if test ! -f $host_frag ; then
   echo "machine dependent configuration file $host_frag does not exist!"
   exit 1
fi

if test "$build" != "$host" ; then
   CROSS=yes
else
   if test "$ARCH" = "CRAY_XT3" ; then
      CROSS=yes
   else
      CROSS=""
   fi
fi

if test -z "$CROSS" ; then
AC_CONFIG_HEADER(config/config.h)
fi

AC_GET_MH($host_frag)

dnl ######################################################################
dnl Checks for programs.
dnl ######################################################################

AC_PROGRAMS_CHECK(F90, f95 f90 ifc pgf90 frt xlf sunf95 lf95)
AC_PROG_RANLIB

AC_PATH_PROG(PERL, perl, /bin/false, \
                $PATH:$prefix/bin:/usr/bin:/usr/local/bin)
AC_PATH_PROG(GMAKE, gmake, /bin/false, \
                $PATH:$prefix/bin:/usr/bin:/usr/local/bin)

dnl ######################################################################

package_library=$LIBS
LIBS=""

if test -z "$CROSS" ; then
  AC_PROG_CC
fi

dnl ######################################################################
dnl Checks for NETCDF.
dnl ######################################################################

NETCDF_LIB=""
NETCDF_INCLUDE=""

AC_CHECK_FILE($NETCDFROOT/include/netcdf.inc,[],AC_MSG_ERROR([NETCDFROOT not properly defined]))
if test -n "$NETCDFROOT" ; then
  if test -r "$NETCDFROOT/lib/libnetcdf.a" ; then
     NETCDF_LIB='-L$(NETCDFROOT)/lib -lnetcdf'
  else
     echo "cannot read $NETCDFROOT/lib/libnetcdf.a"  
     exit 1    
  fi
  if test -r "$NETCDFROOT/include/netcdf.inc" ; then
     NETCDF_INCLUDE='-I$(NETCDFROOT)/include'
  else
     echo "cannot read $NETCDFROOT/include/netcdf.inc"  
     exit 1    
  fi
else
  echo "error: the variable NETCDFROOT is not set in $host_frag"  
  echo "       set NETCDFROOT to the netCDF root directory"
  exit 1
fi

dnl ######################################################################
dnl Checks for MPI.
dnl ######################################################################

AC_ARG_ENABLE(MPI,MPI)

MPI_INCLUDE=""

if test "$enable_MPI" = "no" ; then
  MPI="${DEFOPT}NOMPI"
  MPI_LIB=""
else
  if test -n "$MPIROOT" ; then
    AC_CHECK_FILE($MPIROOT/include/mpif.h,[],AC_MSG_ERROR([MPIROOT not properly defined]))
    MPI_LIB='-L$(MPIROOT)/lib '$MPI_LIB
    if test -r "$MPIROOT/include/mpif.h" ; then
       MPI_INCLUDE='-I$(MPIROOT)/include'
    else
       echo "cannot read $MPIROOT/include/mpif.h"  
       exit 1    
    fi
  else
    MPI="${DEFOPT}NOMPI"
    MPI_LIB=""
  fi
fi

dnl ######################################################################
dnl Checks for header files.

if test -z "$CROSS" ; then

   AC_HEADER_STDC
   AC_CHECK_HEADERS(fcntl.h fortran.h limits.h malloc.h netdb.h pwd.h unistd.h  sys/param.h sys/time.h sys/unistd.h sys/utsname.h )

dnl Checks for typedefs, structures, and compiler characteristics.
   AC_C_CONST
   AC_HEADER_TIME

dnl Checks for library functions.
   AC_CHECK_FUNCS(getrusage gettimeofday sysconf uname valloc)

dnl Check for byte ordering

   AC_C_BIGENDIAN

dnl check for C datatypes

   AC_CHECK_SIZEOF(char)
   AC_CHECK_SIZEOF(short) 
   AC_CHECK_SIZEOF(int) 
   AC_CHECK_SIZEOF(long) 
   AC_CHECK_SIZEOF(long long) 
   AC_CHECK_SIZEOF(float) 
   AC_CHECK_SIZEOF(double) 
   AC_CHECK_SIZEOF(long double)
   AC_CHECK_SIZEOF(int *)

dnl Check for Fortran types:

dnl Defines FORT_INT1=1, FORT_INT2=1, FORT_INT4=1, FORT_INT8=1, FORT_INT16=1
dnl         FORT_REAL4=1, FORT_REAL8=1, FORT_REAL16=1
dnl         FORT_COMPLEX8=1, FORT_COMPLEX16=1, FORT_COMPLEX32=1
dnl         FORT_DOUBLE_COMPLEX=1 if the types are available.
    
   AC_TEST_FORTTYPES()

dnl For each defined type, find the corresponding C datatype.
dnl its a bit complicated, but necessary.
dnl If none can be found, use 0.

   if test -z "$FORT_INTEGER_LEN" ; then
      AC_FORTRAN_GET_INTEGER_SIZE(fort_integer_len)
      if test -z "$fort_integer_len" ; then
   	 echo "Unable to determine the length of a FORTRAN integer"
	 echo "You can set the environment variable FORT_INTEGER_LEN"
	 echo "with the number of bytes in a FORTRAN integer and rerun"
	 echo "configure."
      fi
   else 
      fort_integer_len=$FORT_INTEGER_LEN
   fi

   if test -z "$FORT_REAL_LEN" ; then
      AC_FORTRAN_GET_REAL_SIZE(fort_real_len)
      if test -z "$fort_real_len" ; then
	 echo "Unable to determine the length of a FORTRAN real"
	 echo "You can set the environment variable FORT_REAL_LEN"
	 echo "with the number of bytes in a FORTRAN real and rerun"
	 echo "configure."
      fi
   else 
      fort_real_len=$FORT_REAL_LEN
   fi

dnl Set the corresponding datatypes

   char_def=1
   short_def=4
   int_def=6
   long_def=8
   long_long_def=13
   float_def=10
   double_def=11
   long_double_def=13

dnl Warning: float_cplx is Fortran COMPLEX, not (float,float)

   float_cplx_def=23

dnl Warning: double_cplx is Fortran DOUBLE COMPLEX, not (double,double)

   double_cplx_def=24

   for len in 1 2 4 8 16 ; do
      eval testval=\$"FORT_INT$len"
      AC_MSG_CHECKING(for Fortran INTEGER*$len)
      found=0
      if test "$testval" = 1 ; then
         for c_type in long_long long int short char ; do
	    eval typelen=\$"ac_cv_sizeof_${c_type}"
            if test $len = "$typelen" ; then
                # Fortran type FORT_INT$len corresponds to C type
                # $ctype
                eval FORT_INT$len=\$"${c_type}_def"
                AC_MSG_RESULT(found)
                found=1
                break
            fi 
         done
      fi
      if test $found = 0 ; then
         eval FORT_INT$len=0
         AC_MSG_RESULT(not supported)
      fi
   done

   for len in 4 8 16 ; do
      # Use len2 for Complex.
      len2=`expr $len + $len`
      lenby2=`expr $len / 2`
      eval testval=\$"FORT_REAL$len"
      AC_MSG_CHECKING(for Fortran REAL*$len and COMPLEX*$len2)
      found=0
      if test "$testval" = 1 ; then
         for c_type in long_double double float ; do
	    eval typelen=\$"ac_cv_sizeof_${c_type}"
            if test $len = "$typelen" ; then
                # Fortran type FORT_REAL$len corresponds to C type
                # $ctype
                eval FORT_REAL$len=\$"${c_type}_def"
                # Notes that we may not have been able to determine 
                # fort_real_len
                # The corresponding Complex type
                if test "$fort_real_len" = "$len" ; then 
                    eval FORT_COMPLEX$len2=$float_cplx_def
                elif test "$fort_real_len" = "$lenby2" ; then
                    eval FORT_COMPLEX$len2=$double_cplx_def
                else
                    eval FORT_COMPLEX$len2=0
                fi
                AC_MSG_RESULT(found)
                found=1
                break
            fi 
         done
      fi
      if test $found = 0 ; then
         eval FORT_REAL$len=0
         eval FORT_COMPLEX$len2=0
         AC_MSG_RESULT(not supported)
      fi
   done

   AC_MSG_CHECKING(for Fortran INTEGER)
   found=0
   if test "$fort_integer_len" -gt 1 ; then
      for c_type in long_long long int short char ; do
         eval typelen=\$"ac_cv_sizeof_${c_type}"
         if test $fort_integer_len = "$typelen" ; then
           # Fortran type fort_integer_len corresponds to C type $ctype
           eval FORT_INT=\$"${c_type}_def"
           AC_MSG_RESULT(found)
           found=1
           break
         fi 
      done
   fi
   if test $found = 0 ; then
      eval FORT_INT=0
      AC_MSG_RESULT(not supported)
   fi

   AC_MSG_CHECKING(for Fortran REAL)
   found=0
   if test "$fort_real_len" -gt 1 ; then
      for c_type in float double long_double  ; do
         eval typelen=\$"ac_cv_sizeof_${c_type}"
         if test $fort_real_len = "$typelen" ; then
           # Fortran type fort_real_len corresponds to C type $ctype
           eval FORT_REAL=\$"${c_type}_def"
           AC_MSG_RESULT(found)
           found=1
           break
         fi 
      done
   fi
   if test $found = 0 ; then
      eval FORT_REAL=0
      AC_MSG_RESULT(not supported)
   fi


   AC_DEFINE_UNQUOTED(FORT_INT, $FORT_INT)
   AC_DEFINE_UNQUOTED(FORT_INT1, $FORT_INT1)
   AC_DEFINE_UNQUOTED(FORT_INT2, $FORT_INT2)
   AC_DEFINE_UNQUOTED(FORT_INT4, $FORT_INT4)
   AC_DEFINE_UNQUOTED(FORT_INT8, $FORT_INT8)
   AC_DEFINE_UNQUOTED(FORT_INT16, $FORT_INT16)

   AC_DEFINE_UNQUOTED(FORT_REAL, $FORT_REAL)
   AC_DEFINE_UNQUOTED(FORT_REAL4, $FORT_REAL4)
   AC_DEFINE_UNQUOTED(FORT_REAL8, $FORT_REAL8)
   AC_DEFINE_UNQUOTED(FORT_REAL16, $FORT_REAL16)

   AC_DEFINE_UNQUOTED(FORT_COMPLEX8, $FORT_COMPLEX8)
   AC_DEFINE_UNQUOTED(FORT_COMPLEX16, $FORT_COMPLEX16)
   AC_DEFINE_UNQUOTED(FORT_COMPLEX32, $FORT_COMPLEX32)

   AC_DEFINE_UNQUOTED(FORT_INTEGER_LEN, $fort_integer_len)
   AC_DEFINE_UNQUOTED(FORT_REAL_LEN, $fort_real_len)

   AC_GET_FORTNAMES
   AC_DEFINE_UNQUOTED($FORTRANNAMES)

fi # ! CROSS


LIBS=$package_library

if test -z "$AR" ; then
  AR="ar"
fi  
if test -z "$AS" ; then
  AS="as"
fi  
if test -z "$ARFLAGS" ; then
  ARFLAGS="crv"
fi

AC_SUBST(ARCH,$ARCH)

AC_SUBST(NETCDFROOT,$NETCDFROOT)
AC_SUBST(NETCDF_LIB,$NETCDF_LIB)
AC_SUBST(NETCDF_INCLUDE,$NETCDF_INCLUDE)

AC_SUBST(MPIROOT,$MPIROOT)
AC_SUBST(MPI_LIB,$MPI_LIB)
AC_SUBST(MPI_INCLUDE,$MPI_INCLUDE)

AC_SUBST(AR,$AR)
AC_SUBST(AS,$AS)

AC_SUBST(CFLAGS,$CFLAGS)
AC_SUBST(FFLAGS,$FFLAGS)
AC_SUBST(F90FLAGS,$F90FLAGS)
AC_SUBST(ARFLAGS,$ARFLAGS)
AC_SUBST(LIBS,$LIBS)
AC_SUBST(MODOPT,$MODOPT)
LDFLAGS=$F90FLAGS
AC_SUBST(LDFLAGS,$LDFLAGS)
AC_SUBST(MPI,$MPI)
AC_SUBST(SRCDIRS,$SRCDIRS)

AC_OUTPUT(Makefile)

if test -n "$PERL" ; then
  AC_SUBST(PERL,$PERL)
   PERLS="util/f2html.pl util/F2f.pl util/f2F.pl util/createMakefiles.pl"
  AC_OUTPUT($PERLS)
  chmod 755 $PERLS
else
   echo "PERL is undefined!"
  exit 1
fi

if test -n "$CROSS" ; then
   case "${host}" in	
       es*-*-*)             ARCH=ES ;;
   esac
   cross_config=config/config.h.${ARCH}
   echo "config.status: creating config/config.h from $cross_config"
   if test ! -f $cross_config ; then
      echo "cross configuration file $cross_config does not exist!"
      exit 1
   fi   
   cp $cross_config config/config.h
fi

echo ""
echo "Configuration done. Run gmake to build ECHAM."
echo ""
