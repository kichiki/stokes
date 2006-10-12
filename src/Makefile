CC           = cc
LD           = $(CC)

#CFLAGS	     = -Wall -I/usr/local/include -O3 -mtune=pentium
CFLAGS	     = -Wall -I/usr/local/include\
               -pipe -O3 -march=prescott -fomit-frame-pointer\
               -mfpmath=sse -ffast-math
LDFLAGS      =
LIBS         = -lm -L/usr/local/lib -liter -lstokes\
		-lg2c -lblas -llapack\
		-lcblas -lf77blas -latlas\
		-lnetcdf

all:		test-stokes

OBJ =	test-stokes.o

test-stokes:	$(OBJ)
	$(LD) $(LDFLAGS) -o test-stokes $(OBJ) $(LIBS)

clean:	
	rm -f *~ *.o test-stokes
