# $Id: Makefile,v 1.2 2006/10/19 18:49:36 ichiki Exp $

CC           = cc
LD           = $(CC)

GUILE_INC = `guile-config compile`
GUILE_LIB = `guile-config link`

#CFLAGS	     = -Wall -I/usr/local/include -O3 -mtune=pentium
CFLAGS	     = -Wall\
		-pipe -O3 -march=prescott -fomit-frame-pointer\
		-mfpmath=sse -ffast-math\
		$(GUILE_INC)
LDFLAGS      =
LIBS         = -lm -L/usr/local/lib -liter -lstokes\
		-lg2c -lblas -llapack\
		-lcblas -lf77blas -latlas\
		-lnetcdf\
		$(GUILE_LIB)

all:		test-stokes

OBJ =	test-stokes.o

test-stokes:	$(OBJ)
	$(LD) $(LDFLAGS) -o test-stokes $(OBJ) $(LIBS)

clean:	
	rm -f *~ *.o test-stokes
