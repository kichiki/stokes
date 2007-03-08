# makefile for octave binding of libstokes
# $Id: Makefile,v 1.2 2007/03/08 00:53:08 kichiki Exp $
MKOCTFILE	= /usr/local/bin/mkoctfile
OCTAVE		= /usr/local/bin/octave

#ATLAS_LIBS	= \
#	-lg2c -lblas -llapack\
#	-lcblas -lf77blas -latlas
ATLAS_LIBS	= \
	-L/usr/local/lib/gcc-4.2.0 -lgfortran \
	-L/usr/local/lib -lalapack -lf77blas -lcblas -latlas

LIB_STOKES	= \
		-lm -L/usr/local/lib -liter -lstokes	\
		$(ATLAS_LIBS)\
		-lnetcdf
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
