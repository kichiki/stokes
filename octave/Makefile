# makefile for octave binding of libstokes
# $Id: Makefile,v 1.6 2008/06/07 03:06:54 kichiki Exp $

MKOCTFILE	= /usr/local/bin/mkoctfile
OCTAVE		= /usr/local/bin/octave
# for openSUSE
#MKOCTFILE	= /usr/bin/mkoctfile
#OCTAVE		= /usr/bin/octave

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
#GUILE_LIB = `guile-config link`
# mkoctfile does not accept "-pthread", so cut it as ...
GUILE_LIB = \
	-I/usr/local/include	\
	-D_THREAD_SAFE		\
	-L/usr/local/lib	\
	-lguile			\
	-lltdl			\
	-lgmp			\
	-lcrypt			\
	-lm

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
