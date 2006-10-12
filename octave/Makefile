# makefile for octave binding of libstokes
# $Id: Makefile,v 1.1 2006/10/12 16:07:28 ichiki Exp $
MKOCTFILE	= /usr/local/bin/mkoctfile
OCTAVE		= /usr/local/bin/octave
LIB_STOKES	= \
		-lm -L/usr/local/lib -liter -lstokes	\
		-lg2c -lblas -llapack			\
		-lcblas -lf77blas -latlas		\
		-lnetcdf
LIBS		= $(LIB_STOKES)


all:		stokes_res_ewald_3f.oct

stokes_res_ewald_3f.oct:	stokes_res_ewald_3f.cc
	$(MKOCTFILE)	\
	stokes_res_ewald_3f.cc	\
	$(LIBS)

clean:
	rm -f			\
	stokes_res_ewald_3f.oct	\
	*.o

test:	stokes_res_ewald_3f.oct
	$(OCTAVE) test-stokes.m
