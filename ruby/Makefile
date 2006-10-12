# makefile for ruby binding of libstokes
# $Id: Makefile,v 1.1 2006/10/12 15:46:42 ichiki Exp $
INC_RUBY	= -I/usr/local/lib/ruby/1.8/i386-freebsd6

LIB_STOKES	= \
		-lm -L/usr/local/lib -liter -lstokes	\
		-lg2c -lblas -llapack			\
		-lcblas -lf77blas -latlas		\
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
