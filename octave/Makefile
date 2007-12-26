# makefile for octave binding of libstokes
# $Id: Makefile,v 1.5 2007/12/26 06:55:08 kichiki Exp $

MKOCTFILE	= /usr/local/bin/mkoctfile
OCTAVE		= /usr/local/bin/octave
# for openSUSE
#MKOCTFILE	= /usr/bin/mkoctfile
#OCTAVE		= /usr/bin/octave

RYUON_INC = -I/usr/local/include
RYUON_LIB = -L/usr/local/lib -liter -lstokes
# for openSUSE
#RYUON_INC = 
#RYUON_LIB = -liter -lstokes

GSL_LIB = `pkg-config --libs gsl`
GUILE_LIB = `guile-config link`

#ATLAS_LIBS	= \
#	-lg2c -lblas -llapack\
#	-lcblas -lf77blas -latlas
ATLAS_LIBS	= \
	-L/usr/local/lib/gcc-4.2.3 -lgfortran \
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
	$(GSL_LIB)

LIBS		= $(LIB_STOKES)


all:		stokes_res_3f.oct

stokes_res_3f.oct:	stokes_res_3f.cc
	$(MKOCTFILE)	\
	stokes_res_3f.cc	\
	$(LIBS)

clean:
	rm -f			\
	stokes_res_3f.oct	\
	*.o			\
	*~

test:	stokes_res_3f.oct
	$(OCTAVE) test-stokes.m
