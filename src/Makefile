# $Id: Makefile,v 1.4 2006/10/23 17:09:59 kichiki Exp $

CC           = cc
LD           = $(CC)

GUILE_INC = `guile-config compile`
GUILE_LIB = `guile-config link`

# for libstokes/src/bench.h
RYUON_INC = -I../../libstokes/src

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

all:	test-stokes	\
	xi3		\
	bench3		\
	stokes3


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


## bench3
bench3.o:	bench3.c
	$(CC) $(CFLAGS) $(RYUON_INC) -c bench3.c

bench3-OBJ =	\
		bench3.o\
		configs.o

bench3:	$(bench3-OBJ)
	$(LD) $(LDFLAGS) -o bench3	\
	$(bench3-OBJ)			\
	$(LIBS)

## stokes3
stokes3-OBJ =	stokes3.o

stokes3:	$(stokes3-OBJ)
	$(LD) $(LDFLAGS) -o stokes3	\
	$(stokes3-OBJ)			\
	$(LIBS)

clean:	
	rm -f *~ *.o	\
	test-stokes	\
	xi3		\
	bench3		\
	stokes3
