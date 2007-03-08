# makefile for perl binding of libstokes
# $Id: Makefile,v 1.2 2007/03/08 00:45:33 kichiki Exp $
SWIG		= /usr/local/bin/swig
PERL		= /usr/local/bin/perl
INC_PERL	= -I/usr/local/lib/perl5/5.8.8/mach/CORE
CC		= /usr/bin/cc

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

#CFLAGS		= -Wall -I/usr/local/include
CFLAGS		= -I/usr/local/include

LIBS		= $(LIB_STOKES)


all:		stokes.so


stokes_wrap.c:	../src/stokes.i
	$(SWIG) -perl		\
	$(CFLAGS)		\
	-o stokes_wrap.c	\
	../src/stokes.i

stokes_wrap.o:	stokes_wrap.c
	$(CC) $(CFLAGS)	\
	$(INC_PERL)	\
	-c stokes_wrap.c

stokes.so:	stokes_wrap.o
	$(CC) -Wall -shared	\
	-o stokes.so		\
	stokes_wrap.o		\
	$(LIBS)

clean:
	rm -f		\
	stokes.so	\
	*.o		\
	stokes.pm	\
	stokes_wrap.*	\
	test-stokes.res-3f.nc

test:	stokes.so
	rm -f test-stokes.res-3f.nc
	$(PERL) test-stokes.pl
