# makefile for guile binding of libstokes
# $Id: Makefile,v 1.6 2008/06/05 03:52:29 kichiki Exp $

SWIG		= /usr/local/bin/swig
GUILE		= /usr/local/bin/guile
# for openSUSE
#SWIG		= /usr/bin/swig
#GUILE		= /usr/bin/guile

CC		= /usr/bin/cc

RYUON_INC = -I/usr/local/include
RYUON_LIB = -L/usr/local/lib -liter -lstokes
# for openSUSE
#RYUON_INC = 
#RYUON_LIB = -liter -lstokes

GSL_LIB = `pkg-config --libs gsl`
NITSOL_LIB = -L/usr/local/lib -lnitsol
GUILE_LIB = `guile-config link`

#ATLAS_LIBS	= \
#	-lg2c -lblas -llapack\
#	-lcblas -lf77blas -latlas
ATLAS_LIBS	= \
	-L/usr/local/lib/gcc-4.2.4 -lgfortran \
	-L/usr/local/lib -lalapack -lf77blas -lcblas -latlas	\
	-larpack
# for SUSE-Linux
# for openSUSE
#ATLAS_LIBS	= \
#	-lgfortran -lblas -llapack\
#	-larpack

LIB_STOKES	= \
	-lm $(RYUON_LIB)	\
	$(ATLAS_LIBS)		\
	-lnetcdf		\
	$(GUILE_LIB)		\
	$(GSL_LIB)		\
	$(NITSOL_LIB)

CFLAGS		= -Wall -I/usr/local/include

LIBS		= $(LIB_STOKES)


all:		stokes.so


stokes_wrap.c:	../src/stokes.i
	$(SWIG) -guile -emit-setters	\
	$(CFLAGS)			\
	-o stokes_wrap.c		\
	../src/stokes.i

stokes_wrap.o:	stokes_wrap.c
	$(CC) $(CFLAGS)	\
	-c stokes_wrap.c

stokes.so:	stokes_wrap.o
	$(CC) -Wall -shared	\
	-o stokes.so		\
	stokes_wrap.o		\
	$(LIBS)

clean:
	rm -f		\
	stokes.so	\
	stokes_wrap.*	\
	test-stokes.res-3f.nc

test:	stokes.so
	rm -f test-stokes.res-3f.nc
	$(GUILE) -s test-stokes.scm
