# makefile for guile binding of libstokes
# $Id: Makefile,v 1.1 2006/10/12 15:52:07 ichiki Exp $
SWIG		= /usr/local/bin/swig
GUILE		= /usr/local/bin/guile
CC		= /usr/bin/cc
LIB_STOKES	= \
		-lm -L/usr/local/lib -liter -lstokes	\
		-lg2c -lblas -llapack			\
		-lcblas -lf77blas -latlas		\
		-lnetcdf

#CFLAGS		= -Wall -I/usr/local/include
CFLAGS		= -I/usr/local/include

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
