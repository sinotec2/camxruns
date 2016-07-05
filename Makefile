#
#.........................................................................
# Version "$Id: Makefile.nocpl.sed 1 2014-03-14 20:22:54Z coats $"
# EDSS/Models-3 M3TOOLS
#    (C) 1992-2002 MCNC and Carlie J. Coats, Jr.,
#    (C) 2003-2004 by Baron Advanced Meteorological Systems, and
#    (C) 2005-2014 Carlie J. Coats, Jr.
# Distributed under the GNU GENERAL PUBLIC LICENSE version 2
# See file "GPL.txt" for conditions of use.
#.........................................................................
#  Environment Variables:
#       BIN     machine/OS/compiler/mode type. Shows up as suffix
#               for "Makeinclude.${BIN}" to determine compilation
#               flags, and in ${OBJDIR} and $(INSTALL) to determine
#               binary directories
#       INSTALL installation-directory root, used for "make install":
#               "libioapi.a" and the tool executables will be installed
#               in $(INSTALL)/${BIN}
#.........................................................................
#  Directories:
#       ${BASEDIR}  is the root directory for the I/O API library source,
#                   the M3Tools and M3Test source,the  HTML documentation,
#                   and the (machine/compiler/flag-specific) binary
#                   object/library/executable directories.
#       $(SRCDIR)   is the source directory for the M3TOOLS
#       $(IODIR)    is the source directory for the I/O API library
#       ${OBJDIR}   is the current machine/compiler/flag-specific
#                   build-directory
#       $(F90DIR)   is the current machine/compiler/flag-specific
#                   build-directory for F90-based programs (SGI & Sun)
#.........................................................................
#
#       ---------------     Definitions:   -------------------------

.SUFFIXES: .m4 .c .F .f

BASEDIR = /home/kuang/web31
SRCDIR  = /nfs/REAS/2010
IODIR   = ${BASEDIR}/ioapi
OBJDIR  = ${BASEDIR}/${BIN}
INSTDIR = /home/kuang/CMAQ/ioapi/Linux2_ia32ifort

# Architecture dependent stuff
# Assumes FC is an f90

include /home/kuang/web31/ioapi/Makeinclude.Linux2_ia32ifort

FFLAGS = -I$(IODIR) -DIOAPICPL $(ARCHFLAGS) $(PARFLAGS) $(FOPTFLAGS) $(ARCHFLAGS)

LDFLAGS = -I$(IODIR) -DIOAPICPL $(DEFINEFLAGS) $(ARCHFLAGS)

#  Incompatibility between netCDF versions before / after v4.1.1:
#  For netCDF v4 and later, you may also need the extra libraries
#  given by netCDF commands
#
#          nc-config --libs
#          nf-config --libs
#
#LIBS = -L${OBJDIR} -lioapi -lnetcdf $(OMPLIBS) $(ARCHLIB) $(ARCHLIBS)

LIBS  = -L${OBJDIR} -lioapi -lnetcdff -lnetcdf $(OMPLIBS) $(ARCHLIB) $(ARCHLIBS)

VPATH = ${OBJDIR}



fSRC = \
txt-nc.f

OBJ = $(fSRC:.f=.o)

EXE = \
txt-nc


#      ----------------------   TOP-LEVEL TARGETS:   ------------------

all: $(EXE)

clean:
	cd ${SRCDIR}; rm $(EXE) $(OBJ)

install: /home/kuang/CMAQ/ioapi/Linux2_ia32ifort
	echo "Installing M3TOOLS in ${INSTDIR}"
	cd ${OBJDIR}; cp $(EXE) $(INSTDIR)

rmexe:
	cd ${OBJDIR}; rm ${EXE}

relink:
	make BIN=${BIN} -i rmexe ; make BIN=${BIN} all

bins:
	make BIN=Linux2_x86_64g95
	make BIN=Linux2_x86_64sun
	make BIN=Linux2_x86_64ifort
	make BIN=Linux2_x86_64g95dbg
	make BIN=Linux2_x86_64sundbg
	make BIN=Linux2_x86_64ifortdbg

binclean:
	make -i BIN=Linux2_x86_64          clean
	make -i BIN=Linux2_x86_64g95       clean
	make -i BIN=Linux2_x86_64sun       clean
	make -i BIN=Linux2_x86_64ifort     clean
	make -i BIN=Linux2_x86_64g95dbg    clean
	make -i BIN=Linux2_x86_64sundbg    clean
	make -i BIN=Linux2_x86_64ifortdbg  clean

binrelink:
	make BIN=Linux2_x86_64         relink
	make BIN=Linux2_x86_64g95      relink
	make BIN=Linux2_x86_64sun      relink
	make BIN=Linux2_x86_64ifort    relink
	make BIN=Linux2_x86_64g95dbg   relink
	make BIN=Linux2_x86_64sundbg   relink
	make BIN=Linux2_x86_64ifortdbg relink

flags:
	echo "BIN=${BIN}"
	echo "FFLAGS=$(FFLAGS)"
	echo "LDFLAGS=$(LDFLAGS)"
	echo "LIBS=$(LIBS)"
	echo "ARCHFLAGS=$(ARCHFLAGS)"
	echo "ARCHLIB=$(ARCHLIB)"
	echo "ARCHLIBS=$(ARCHLIBS)"
	echo "OMPFLAGS=$(OMPFLAGS)"
	echo "OMPLIBS=$(OMPLIBS)"
	echo "FOPTFLAGS=$(FOPTFLAGS)"
	echo "COPTFLAGS=$(COPTFLAGS)"
	echo "PARFLAGS=$(PARFLAGS)"
	echo "PVM_ROOT=$(PVM_ROOT)"
	echo "PVMLIBS=$(PVMLIBS)"


#      -----------------------   RULES:   -------------------------

%.o : %.mod        #  Disable "gmake"s obnoxious implicit Modula-2 rule !!
%.f : %.F          #  Hack for some versions of  "gmake" + "gfortran"

.F.o: m3utilio.mod
	cd ${SRCDIR}; $(FC) $(FPPFLAGS) $(FFLAGS) -c $(SRCDIR)/$<


.f.o:  m3utilio.mod
	cd ${SRCDIR}; $(FC) $(FFLAGS) -c $(SRCDIR)/$<


#  ---------------------------  $(EXE) Program builds:  -----------------



txt-nc:  txt-nc.o
	cd ${SRCDIR}; $(FC) ${LFLAGS} $^ ${LIBS} -o $@

