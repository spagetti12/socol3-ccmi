#!/usr/bin/perl
# -*-perl-*-
#
# Usage: createMakefiles.pl
#
# Generate a Makefile in the modules and src directory from ECHAM.
# To run createMakefiles, it will be necessary to modify the first line of 
# this script to point to the actual location of Perl on your system.
#
# Written by Uwe Schulzweida <schulzweida@dkrz.de> June, 1999
#
$EXECUTABLE="echam5";
#
$PROG=`basename "$0"`;
#
# Create Makefile in modules
#
$dir = modules; 
#
chdir $dir;
#
rename "Makefile", "Makefile.old";
open(MAKEFILE, "> Makefile");
print "create new Makefile in: $dir\n";
print MAKEFILE "# Generated automatically by $PROG \n";
#
# Source listing
#
print MAKEFILE "SRCS =\t";
@srcs = <*.f90 *.f *.F *.c>;
&PrintWords(8, 0, @srcs);
print MAKEFILE "\n\n";
#
# Object listing
#
#print MAKEFILE "OBJS =\t";
@objs = @srcs;
foreach (@objs) { s/\.[^.]+$/.o/ };
@mobjs = @objs;
#&PrintWords(8, 0, @objs);
print MAKEFILE 'OBJS := $(SRCS:.f90=.o)';
print MAKEFILE "\n\n";
#
# make
#
print MAKEFILE "all: \$(OBJS)\n\n";
#
# make clean
#
print MAKEFILE "clean:\n";
print MAKEFILE "\trm -f \$(OBJS) *.mod\n\n";
#
# Make .f90 a valid suffix
#
print MAKEFILE ".SUFFIXES: \$(SUFFIXES) .f90\n\n";
#
# .f90 -> .o
#
print MAKEFILE "%.o: %.f90\n";
print MAKEFILE "\t\$(F90) \$(F90FLAGS) -c \$<\n\n";
#
# Dependency listings
#
&MakeDependsf90($ARGV[1]);
#
#
# mo_spitfire: needs additional compiler options for inlining (SX bug)
# mo_tpcore:   needs additional compiler options for inlining (SX bug)
#
print MAKEFILE "\nifeq (\$(strip \$(ARCH)), SX)\n";
print MAKEFILE "mo_spitfire.o: mo_spitfire.f90\n";
print MAKEFILE "\t\$(F90) \$(F90FLAGS) -pi auto exp=minmod exp=medan exp=putyslice exp=cfdot1dp2 exp=cfint1x2 line=2000 -c mo_spitfire.f90\n";
print MAKEFILE "mo_tpcore.o: mo_tpcore.f90\n";
print MAKEFILE "#warning: Don't change line to 2000. This will give wrong code on the SX!!!\n";
print MAKEFILE "\t\$(F90) \$(F90FLAGS) -pi auto exp=xmist,fxppm,kmppm,lmppm,xtp noexp=map1_ppm_gp,ppm2m,steepz nest=3 line=1000 -c mo_tpcore.f90\n";
print MAKEFILE "mo_transpose.o: mo_transpose.f90\n";
print MAKEFILE "\t\$(F90) \$(F90FLAGS) -Npi -c mo_transpose.f90\n";
print MAKEFILE "endif\n\n";
print MAKEFILE "\nifeq (\$(strip \$(ARCH)), ES)\n";
print MAKEFILE "mo_spitfire.o: mo_spitfire.f90\n";
print MAKEFILE "\t\$(F90) \$(F90FLAGS) -pi auto exp=minmod exp=medan exp=putyslice exp=cfdot1dp2 exp=cfint1x2 line=2000 -c mo_spitfire.f90\n";
print MAKEFILE "mo_tpcore.o: mo_tpcore.f90\n";
print MAKEFILE "#warning: Don't change line to 2000. This will give wrong code on the SX!!!\n";
print MAKEFILE "\t\$(F90) \$(F90FLAGS) -pi auto exp=xmist,fxppm,kmppm,lmppm,xtp noexp=map1_ppm_gp,ppm2m,steepz nest=3 line=1000 -c mo_tpcore.f90\n";
print MAKEFILE "mo_transpose.o: mo_transpose.f90\n";
print MAKEFILE "\t\$(F90) \$(F90FLAGS) -Npi -c mo_transpose.f90\n";
print MAKEFILE "endif\n";
#
# mo_buffer_fft: compile only with -Ovector1 (CF90 bug, CRAY C90)
#
print MAKEFILE "\nifeq (\$(strip \$(ARCH)), CRAY_PVP)\n";
print MAKEFILE "mo_buffer_fft.o: mo_buffer_fft.f90\n";
print MAKEFILE "\t\$(F90) \$(F90FLAGS) -Ovector1 -c mo_buffer_fft.f90\n";
print MAKEFILE "mo_grib.o: mo_grib.f90\n";
print MAKEFILE "\t\$(F90) \$(F90FLAGS) -Ovector1 -c mo_grib.f90\n";
print MAKEFILE "endif\n";
#
# mo_io: compile only with -Ovector1 (CF90 bug, CRAY X1)
#
print MAKEFILE "\nifeq (\$(strip \$(ARCH)), CRAY_X1)\n";
print MAKEFILE "\t\$(F90) \$(F90FLAGS) -Ovector1 -Oscalar1 -Ostream1 -c mo_io.f90\n";
print MAKEFILE "endif\n";
#
# mo_transpose:
#
# Note: PGI compiler seems to interrogate $F90FLAGS from the environment
# and a -fastsse flag specified there will override the command line option
print MAKEFILE "\nifeq (\$(strip \$(ARCH)), CRAY_XT3)\n";
print MAKEFILE "mo_transpose.o: mo_transpose.f90\n";
print MAKEFILE "\t( F90FLAGS=\" \" ; \$(F90) -O2  -c mo_transpose.f90 )\n";
print MAKEFILE "endif\n\n";
print MAKEFILE "ifeq (\$(strip \$(ARCH)), CRAY_XD1)\n";
print MAKEFILE "mo_transpose.o: mo_transpose.f90\n";
print MAKEFILE "\t( F90FLAGS=\" \" ; \$(F90) -O2  -c mo_transpose.f90 )\n";
print MAKEFILE "endif\n";
#
close MAKEFILE;
#
chdir "..";
#
# Create Makefile in src
#
$dir = src; 
#
chdir $dir;
#
rename "Makefile", "Makefile.old";
open(MAKEFILE, "> Makefile");
print "create new Makefile in: $dir\n";
print MAKEFILE "# Generated automatically by $PROG \n";
#
print MAKEFILE "PROG =\t../bin/$EXECUTABLE\n\n";
#
# Source listing
#
print MAKEFILE "SRCS =\t";
@srcs = <*.f90 *.f *.F *.c>;
&PrintWords(8, 0, @srcs);
print MAKEFILE "\n\n";
#
# Object listing
#
#print MAKEFILE "OBJS =\t";
@objs = @srcs;
foreach (@objs) { s/\.[^.]+$/.o/ };
#&PrintWords(8, 0, @objs);
print MAKEFILE 'OBJS := $(SRCS:.f90=.o)';
print MAKEFILE "\n\n";
#
# Module object listing
#
print MAKEFILE "MOBJS =\t";
foreach (@mobjs) { s-^-\$(MODULES)/- };
&PrintWords(3, 0, @mobjs);
print MAKEFILE "\n\n";
#
# make
#
print MAKEFILE "all: \$(PROG)\n\n";
print MAKEFILE "\$(PROG): \$(OBJS) \$(MOBJS) ../lib/libsupport.a\n";
print MAKEFILE "\t\$(", &LanguageCompiler($ARGV[1], @srcs);
print MAKEFILE ") \$(LDFLAGS) -o \$@ \$(OBJS) \$(MOBJS) \$(LIBS)\n\n";
#
# make clean
#
print MAKEFILE "clean:\n";
print MAKEFILE "\trm -f \$(PROG) \$(OBJS) *.mod\n\n";
#
# Make .f90 a valid suffix
#
print MAKEFILE ".SUFFIXES: \$(SUFFIXES) .f90\n\n";
#
# .f90 -> .o
#
print MAKEFILE "%.o: %.f90\n";
print MAKEFILE "\t\$(F90) \$(F90FLAGS) -c \$<\n\n";
#
# Dependency listings
#
&MakeDependsf90mod($ARGV[1]);
#
print MAKEFILE "\nifeq (\$(strip \$(ARCH)), SX)\n";
print MAKEFILE "lti.o: lti.f90\n";
print MAKEFILE "\t\$(F90) \$(F90FLAGS) -Npi -c lti.f90\n";
print MAKEFILE "endif\n\n";
print MAKEFILE "\nifeq (\$(strip \$(ARCH)), ES)\n";
print MAKEFILE "lti.o: lti.f90\n";
print MAKEFILE "\t\$(F90) \$(F90FLAGS) -Npi -c lti.f90\n";
print MAKEFILE "endif\n\n";
#
close MAKEFILE;
#
exit;
#
#
# &PrintWords(current output column, extra tab?, word list); --- print words
#    nicely
#
sub PrintWords {
   local($columns) = 78 - shift(@_);
   local($extratab) = shift(@_);
   local($wordlength);
   #
   print MAKEFILE @_[0];
   $columns -= length(shift(@_));
   foreach $word (@_) {
      $wordlength = length($word);
      if ($wordlength + 1 < $columns) {
         print MAKEFILE " $word";
         $columns -= $wordlength + 1;
         }
      else {
         #
         # Continue onto a new line
         #
         if ($extratab) {
            print MAKEFILE " \\\n\t\t$word";
            $columns = 62 - $wordlength;
            }
         else {
            print MAKEFILE " \\\n\t$word";
            $columns = 70 - $wordlength;
            }
         }
      }
   }

#
# &LanguageCompiler(compiler, sources); --- determine the correct language
#    compiler
#
sub LanguageCompiler {
   local($compiler) = &toLower(shift(@_));
   local(@srcs) = @_;
   #
   if (length($compiler) > 0) {
      CASE: {
         grep(/^$compiler$/, ("fc", "f77")) &&
            do { $compiler = "FC"; last CASE; };
         grep(/^$compiler$/, ("cc", "c"))   &&
            do { $compiler = "CC"; last CASE; };
         $compiler = "F90";
         }
      }
   else {
      CASE: {
         grep(/\.f90$/, @srcs)   && do { $compiler = "F90"; last CASE; };
         grep(/\.(f|F)$/, @srcs) && do { $compiler = "FC";  last CASE; };
         grep(/\.c$/, @srcs)     && do { $compiler = "CC";  last CASE; };
         $compiler = "???";
         }
      }
   $compiler;
   }

#
# &toLower(string); --- convert string into lower case
#
sub toLower {
   local($string) = @_[0];
   $string =~ tr/A-Z/a-z/;
   $string;
   }

#
# &uniq(sorted word list); --- remove adjacent duplicate words
#
sub uniq {
   local(@words);
   foreach $word (@_) {
      if ($word ne $words[$#words]) {
         push(@words, $word);
         }
      }
   @words;
   }

#
# &MakeDepends(language pattern, include file sed pattern); --- dependency
#    maker
#
sub MakeDepends {
   local(@incs);
   local($lang) = @_[0];
   local($pattern) = @_[1];
   local(@modules);
   #
   foreach $file (<${lang}>) {
      open(FILE, $file) || warn "Cannot open $file: $!\n";
      while (<FILE>) {
         /$pattern/i && push(@incs, "\$(INCLUDE)/$1");
         /^\s*use\s+([^\s,!]+)/i && push(@modules, "\$(MODULES)/$1.o");
      }
      if (defined @incs || defined @modules) {
         $file =~ s/\.[^.]+$/.o/;
         print MAKEFILE "$file: ";
         &PrintWords(length($file) + 2, 0, @modules, @incs);
         print MAKEFILE "\n";
         undef @incs;
         undef @modules;
         }
      }
   }

#
# &MakeDependsf90(f90 compiler); --- FORTRAN 90 dependency maker
#
sub MakeDependsf90 {
   local(@dependencies);
   local(%filename);
   local(@incs);
   local(@modules);
   local($objfile);
   #
   # Associate each module with the name of the file that contains it
   #
   foreach $file (<*.f90>) {
      open(FILE, $file) || warn "Cannot open $file: $!\n";
      while (<FILE>) {
         /^\s*module\s+([^\s!]+)/i &&
            ($filename{&toLower($1)} = $file) =~ s/\.f90$/.o/;
         }
      }
   #
   # Print the dependencies of each file that has one or more include's or
   # references one or more modules
   #
   foreach $file (<*.f90>) {
      ($objfile = $file) =~ s/\.f90$/.o/;
      open(FILE, $file);
      while (<FILE>) {
         /^\s*#*include\s+["\']([^"\']+)["\']/i && push(@incs, $1);
         /^\s*use\s+([^\s,!]+)/i && push(@modules, &toLower($1));
      }

      if (defined @incs) {
          $jn = $#incs;
          while ($jn >= 0) {
             if (-e "../include/$incs[$jn]") {
                $incs[$jn] = "\$(INCLUDE)/$incs[$jn]";
             }
             else {
                print "  $file: include file $incs[$jn] not in \$(INCLUDE)\n";
                pop(@incs);
             }
             $jn--;
          }
      }
      if (defined @incs && $#incs < 0) {undef @incs;}

      if (defined @modules) {
          $jn = $#modules;
          while ($jn >= 0) {
             if ("$objfile" eq "$filename{$modules[$jn]}") {
#                print "  $objfile: $modules[$jn]\n";
                pop(@modules);
             }
             $jn--;
          }
      }
      if (defined @modules && $#modules < 0) {undef @modules;}

      if (defined @incs || defined @modules) {
         print MAKEFILE "$objfile: ";
         undef @dependencies;
         foreach $module (@modules) {
            push(@dependencies, $filename{$module});
            }
         @dependencies = &uniq(sort(@dependencies));
         &PrintWords(length($objfile) + 2, 0,
                     @dependencies, &uniq(sort(@incs)));
         print MAKEFILE "\n";
         undef @incs;
         undef @modules;
         }
      }
   }

#
# &MakeDependsf90(f90 compiler); --- FORTRAN 90 dependency maker
#
sub MakeDependsf90mod {
   local(@dependencies);
   local(%filename);
   local(@incs);
   local(@modules);
   local($objfile);
   #
   # Associate each module with the name of the file that contains it
   #
   foreach $file (<*.f90>) {
      open(FILE, $file) || warn "Cannot open $file: $!\n";
      while (<FILE>) {
         /^\s*module\s+([^\s!]+)/i &&
            ($filename{&toLower($1)} = $file) =~ s/\.f90$/.o/;
         }
      }
   chdir "../modules";
   foreach $file (<*.f90>) {
      open(FILE, $file) || warn "Cannot open $file: $!\n";
      while (<FILE>) {
         /^\s*module\s+([^\s!]+)/i &&
            ($filename{&toLower($1)} = "\$(MODULES)/$file") =~ s/\.f90$/.o/;
         }
      }
   chdir "../src";
   #
   # Print the dependencies of each file that has one or more include's or
   # references one or more modules
   #
   foreach $file (<*.f90>) {
      ($objfile = $file) =~ s/\.f90$/.o/;
      open(FILE, $file);
      while (<FILE>) {
         /^\s*#*include\s+["\']([^"\']+)["\']/i && push(@incs, $1);
         /^\s*use\s+([^\s,!]+)/i && push(@modules, &toLower($1));
         }

      if (defined @incs) {
          $jn = $#incs;
          while ($jn >= 0) {
             if (-e "../include/$incs[$jn]") {
                $incs[$jn] = "\$(INCLUDE)/$incs[$jn]";
             }
             else {
                print "  $file: include file $incs[$jn] not in \$(INCLUDE)\n";
                pop(@incs);
             }
             $jn--;
          }
      }
      if (defined @incs && $#incs < 0) {undef @incs;}

      if (defined @modules) {
          $jn = $#modules;
          while ($jn >= 0) {
             if ("$objfile" eq "$filename{$modules[$jn]}") {
#                print "  $objfile: $modules[$jn]\n";
                pop(@modules);
             }
             $jn--;
          }
      }
      if (defined @modules && $#modules < 0) {undef @modules;}

      if (defined @incs || defined @modules) {
         print MAKEFILE "$objfile: ";
         undef @dependencies;
         foreach $module (@modules) {
            if ($filename{$module})  {
               push(@dependencies, $filename{$module});
            }
            else {
               print "  $file: module $module not in \$(MODULES)\n";
            }
         }
         @dependencies = &uniq(sort(@dependencies));
         &PrintWords(length($objfile) + 2, 0,
                     @dependencies, &uniq(sort(@incs)));
         print MAKEFILE "\n";
         undef @incs;
         undef @modules;
         }
      }
   }
