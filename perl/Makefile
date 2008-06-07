# makefile for perl binding of libstokes
# $Id: Makefile,v 1.7 2008/06/07 03:08:45 kichiki Exp $

SWIG		= /usr/local/bin/swig
PERL		= /usr/local/bin/perl
INC_PERL	= -I/usr/local/lib/perl5/5.8.8/mach/CORE
# for openSUSE
#SWIG		= /usr/bin/swig
#PERL		= /usr/bin/perl
#INC_PERL	= -I/usr/lib/perl5/5.8.8/i586-linux-thread-multi/CORE

# perl -e 'use Config; print $Config{ccflags};'                                 
# at FreeBSD
PERL_OPT        = -DAPPLLIB_EXP="/usr/local/lib/perl5/5.8.8/BSDPAN" -DHAS_FPSETMASK -DHAS_FLOATINGPOINT_H -fno-strict-aliasing -pipe -Wdeclaration-after-statement -I/usr/local/include
# at openSUSE
# PERL_OPT	= -D_REENTRANT -D_GNU_SOURCE -DTHREADS_HAVE_PIDS -DDEBUGGING -fno-strict-aliasing -pipe -Wdeclaration-after-statement -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64

CC		= /usr/bin/cc

RYUON_INC = -I/usr/local/include
RYUON_LIB = -L/usr/local/lib -liter -lstokes
# for openSUSE
#RYUON_INC = -I/home/ichiki/RYUON/include
#RYUON_LIB = -L/home/ichiki/RYUON/lib -liter -lstokes

GSL_INC = `pkg-config --cflags gsl`
GSL_LIB = `pkg-config --libs gsl`

NITSOL_INC = -I/usr/local/include
NITSOL_LIB = -L/usr/local/lib -lnitsol
# for SUSE-Linux
#NITSOL_INC = -I/home/ichiki/RYUON/include
#NITSOL_LIB = -L/home/ichiki/RYUON/lib -lnitsol

GUILE_INC = `guile-config compile`
GUILE_LIB = `guile-config link`

#ATLAS_LIBS	= \
#	-lg2c -lblas -llapack\
#	-lcblas -lf77blas -latlas
ATLAS_LIBS	= \
	-L/usr/local/lib/gcc-4.2.4 -lgfortran \
	-L/usr/local/lib -lalapack -lf77blas -lcblas -latlas\
	-larpack
# for SUSE-Linux
# for openSUSE
#ATLAS_LIBS	= \
#	-lgfortran -lblas -llapack\
#	-larpack
#ATLAS_LIBS	= \
#	-llapack -lptcblas -lptf77blas -latlas\
#	-larpack

LIB_STOKES	= \
	-lm $(RYUON_LIB)	\
	$(ATLAS_LIBS)		\
	-lnetcdf		\
	$(GUILE_LIB)		\
	$(GSL_LIB)		\
	$(NITSOL_LIB)

CFLAGS = -Wall		\
	$(RYUON_INC)

LIBS		= $(LIB_STOKES)


all:		stokes.so


stokes_wrap.c:	../src/stokes.i
	$(SWIG) -perl		\
	$(CFLAGS)		\
	-o stokes_wrap.c	\
	../src/stokes.i

stokes_wrap.o:	stokes_wrap.c
	$(CC) $(CFLAGS)	\
	$(INC_PERL)	\
        $(PERL_OPT)     \
	-c stokes_wrap.c

stokes.so:	stokes_wrap.o
	$(CC) -Wall -shared	\
	-o stokes.so		\
	stokes_wrap.o		\
	$(LIBS)

clean:
	rm -f		\
	stokes.so	\
	*.o		\
	stokes.pm	\
	stokes_wrap.*	\
	test-stokes.res-3f.nc

test:	stokes.so
	rm -f test-stokes.res-3f.nc
	$(PERL) test-stokes.pl
