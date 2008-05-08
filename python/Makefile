# makefile for python binding of libstokes
# $Id: Makefile,v 1.7 2008/05/08 03:12:32 kichiki Exp $

SWIG		= /usr/local/bin/swig
PYTHON		= /usr/local/bin/python
INC_PYTHON	= -I/usr/local/include/python2.5
# for openSUSE
#SWIG		= /usr/bin/swig
#PYTHON		= /usr/bin/python
#INC_PYTHON	= -I/usr/include/python2.4

CC		= /usr/bin/cc

RYUON_INC = -I/usr/local/include
RYUON_LIB = -L/usr/local/lib -liter -lstokes
# for openSUSE
#RYUON_INC = -I/home/ichiki/RYUON/include
#RYUON_LIB = -L/home/ichiki/RYUON/lib -liter -lstokes

GSL_LIB = `pkg-config --libs gsl`
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
#ATLAS_LIBS	= \
#	-llapack -lptcblas -lptf77blas -latlas\
#	-larpack

LIB_STOKES	= \
	-lm $(RYUON_LIB)	\
	$(ATLAS_LIBS)		\
	-lnetcdf		\
	$(GUILE_LIB)		\
	$(GSL_LIB)

CFLAGS = -Wall		\
	$(RYUON_INC)

LIBS		= $(LIB_STOKES)


all:		_stokes.so


stokes_wrap.c:	../src/stokes.i
	$(SWIG) -python		\
	$(CFLAGS)		\
	-o stokes_wrap.c	\
	../src/stokes.i

stokes_wrap.o:	stokes_wrap.c
	$(CC) -fPIC	\
	$(CFLAGS)	\
	$(INC_PYTHON)	\
	-c stokes_wrap.c

_stokes.so:	stokes_wrap.o
	$(CC) -Wall -shared	\
	-o _stokes.so		\
	stokes_wrap.o		\
	$(LIBS)

clean:
	rm -f		\
	_stokes.so	\
	*~		\
	*.o		\
	*.pyc		\
	stokes.py*	\
	stokes_wrap.*	\
	test-stokes.res-3f.nc

test:	_stokes.so
	rm -f test-stokes.res-3f.nc
	$(PYTHON) test-stokes.py
