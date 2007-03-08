# makefile for ruby binding of libstokes
# $Id: Makefile,v 1.2 2007/03/08 00:50:33 kichiki Exp $
INC_RUBY	= -I/usr/local/lib/ruby/1.8/i386-freebsd6

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

CFLAGS		= -Wall -I/usr/local/include

LIBS      = $(LIB_STOKES)


all:		stokes.so


stokes_wrap.c:	../src/stokes.i
	swig -ruby		\
	$(CFLAGS)		\
	-o stokes_wrap.c	\
	../src/stokes.i

stokes_wrap.o:	stokes_wrap.c
	cc $(CFLAGS)	\
	$(INC_RUBY)	\
	-c stokes_wrap.c

stokes.so:	stokes_wrap.o
	cc -Wall -shared	\
	-o stokes.so		\
	stokes_wrap.o		\
	$(LIBS)

clean:
	rm -f		\
	stokes.so	\
	stokes_wrap.*	\
	test-stokes.res-3f.nc

test:	stokes.so
	ruby test-stokes.rb
