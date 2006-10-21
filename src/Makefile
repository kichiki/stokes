# $Id: Makefile,v 1.3 2006/10/21 19:07:21 kichiki Exp $

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

all:		test-stokes	\
		xi3

# test-stokes
test-stokes-OBJ =	test-stokes.o

test-stokes:	$(test-stokes-OBJ)
	$(LD) $(LDFLAGS) -o test-stokes	\
	$(test-stokes-OBJ)		\
	$(LIBS)

# xi3
xi3-OBJ =	xi3.o

xi3:	$(xi3-OBJ)
	$(LD) $(LDFLAGS) -o xi3	\
	$(xi3-OBJ)		\
	$(LIBS)

clean:	
	rm -f *~ *.o	\
	test-stokes	\
	xi3
