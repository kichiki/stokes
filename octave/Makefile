# makefile for octave binding of libstokes
# $Id: Makefile,v 1.4 2007/12/05 04:11:40 kichiki Exp $
MKOCTFILE	= /usr/local/bin/mkoctfile
OCTAVE		= /usr/local/bin/octave

GSL_LIB = `pkg-config --libs gsl`

#ATLAS_LIBS	= \
#	-lg2c -lblas -llapack\
#	-lcblas -lf77blas -latlas
ATLAS_LIBS	= \
	-L/usr/local/lib/gcc-4.2.3 -lgfortran \
	-L/usr/local/lib -lalapack -lf77blas -lcblas -latlas	\
	-larpack

LIB_STOKES	= \
		-lm -L/usr/local/lib -liter -lstokes	\
		$(ATLAS_LIBS)	\
		-lnetcdf	\
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
	*.o

test:	stokes_res_3f.oct
	$(OCTAVE) test-stokes.m
