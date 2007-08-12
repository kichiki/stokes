# makefile for python binding of libstokes
# $Id: Makefile,v 1.4 2007/08/12 19:12:38 kichiki Exp $
SWIG		= /usr/local/bin/swig
PYTHON		= /usr/local/bin/python
INC_PYTHON	= -I/usr/local/include/python2.5
CC		= /usr/bin/cc

GUILE_LIB = `guile-config link`

#ATLAS_LIBS	= \
#	-lg2c -lblas -llapack\
#	-lcblas -lf77blas -latlas
ATLAS_LIBS	= \
	-L/usr/local/lib/gcc-4.2.1 -lgfortran \
	-L/usr/local/lib -lalapack -lf77blas -lcblas -latlas

LIB_STOKES	= \
		-lm -L/usr/local/lib -liter -lstokes	\
		$(ATLAS_LIBS)\
		-lnetcdf\
		$(GUILE_LIB)

CFLAGS		= -Wall -I/usr/local/include

LIBS		= $(LIB_STOKES)


all:		_stokes.so


stokes_wrap.c:	../src/stokes.i
	$(SWIG) -python	\
	$(CFLAGS)	\
	-o stokes_wrap.c\
	../src/stokes.i

stokes_wrap.o:	stokes_wrap.c
	$(CC) $(CFLAGS)	\
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
	*.o		\
	stokes.py*	\
	stokes_wrap.*	\
	test-stokes.res-3f.nc

test:	_stokes.so
	rm -f test-stokes.res-3f.nc
	$(PYTHON) test-stokes.py
