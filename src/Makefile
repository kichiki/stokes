# $Id: Makefile,v 1.17 2008/11/01 05:58:26 kichiki Exp $

#------------------------------------------------------------------------------
# for FreeBSD
CC = cc
RYUON_PREFIX = /usr/local
# for SUSE-Linux
#RYUON_PREFIX = /home/ichiki/RYUON

#OPTS = -pipe -O3 -fomit-frame-pointer -ffast-math
OPTS = -pipe -O3 -fomit-frame-pointer -ffast-math\
	-march=prescott -mfpmath=sse

# for gfortran
F2C_LIBS = -L/usr/local/lib/gcc-4.2.5 -lgfortran
# plain g2c library
#F2C_LIBS = -lg2c

#------------------------------------------------------------------------------
# for glacier
#CC = /global/software/intel-9/cc/bin/icc
#RYUON_PREFIX = /global/home/ichiki/RYUON
#OPTS = -O3
#F2C_LIBS = -L/global/software/intel-9/fc/lib -lifcore

#------------------------------------------------------------------------------
# for SUSE-Linux
#CC = cc
#RYUON_PREFIX = /usr/local
#OPTS = \
#	-mtune=opteron -m64 \
#	-pipe -O3 -fomit-frame-pointer -ffast-math
#F2C_LIBS = -lg2c

#------------------------------------------------------------------------------
# for i1200 (openSUSE)
#CC = cc
#RYUON_PREFIX = /usr/local
#OPTS = \
#	-march=pentium3m \
#	-pipe -O3 -fomit-frame-pointer -ffast-math
#F2C_LIBS = -lg2c

#------------------------------------------------------------------------------
# for robson gcc3/g77 (32bit)
#CC = /usr/bin/gcc
#RYUON_PREFIX = /home/ichiki/RYUON3-32
#OPTS = -O3 -mcpu=970 -mtune=970 -maltivec
#F2C_LIBS = -lg2c

# for robson gcc3/g77 (64bit)
#CC = /usr/bin/gcc -m64
#RYUON_PREFIX = /home/ichiki/RYUON3
#OPTS = -O3 -mcpu=970 -mtune=970 -maltivec
#F2C_LIBS = -lg2c

# for robson gcc4/gfortran (32bit)
#CC = /usr/local/bin/gcc4
#RYUON_PREFIX = /home/ichiki/RYUON4-32
#OPTS = -O3 -mcpu=970 -mtune=970 -maltivec
#F2C_LIBS = -L/usr/local/gcc-4.2.2/lib/gcc/powerpc64-suse-linux/4.2.2 -lgfortran

# for robson gcc4/gfortran (64bit)
#CC = /usr/local/bin/gcc4 -m64
#RYUON_PREFIX = /home/ichiki/RYUON4
#OPTS = -O3 -mcpu=970 -mtune=970 -maltivec
#F2C_LIBS = -L/usr/local/gcc-4.2.2/lib/gcc/powerpc64-suse-linux/4.2.2 -lgfortran


RYUON_BIN = $(RYUON_PREFIX)/bin
RYUON_INC = -I$(RYUON_PREFIX)/include
RYUON_LIB = -L$(RYUON_PREFIX)/lib -liter -lstokes

LD = $(CC)
# for SUSE-Linux
#LD = $(CC) -Wl,-R$(RYUON_LIB)


# for header file besides libstokes.h
RYUON_SRC = -I../../libstokes/src


#------------------------------------------------------------------------------
# for FreeBSD with atlas
MATH_LIBS = \
	$(F2C_LIBS)\
	-L/usr/local/lib -lalapack -lf77blas -lcblas -latlas\
	-larpack
# for glacier
#MATH_LIBS = \
#	$(F2C_LIBS)\
#	-L/usr/lib -llapack -lblas \
#	$(RYUON_LIB) -larpack

# for SUSE-Linux
# for openSUSE
#MATH_LIBS = \
#	$(F2C_LIBS)\
#	-lblas -llapack\
#	-larpack
# for robson with local libraries
#MATH_LIBS = \
#	$(F2C_LIBS)\
#	$(RYUON_LIB) -llapack -lblas\
#	-larpack
# for robson (gcc/g77, 64bit)
#MATH_LIBS = \
#	$(F2C_LIBS)\
#	-L/usr/local/lib64 -lcblas -latlas\
#	-larpack
# for robson (gcc4/gfortran, 64bit)
#MATH_LIBS = \
#	$(F2C_LIBS)\
#	-L/usr/local/lib64 -larpack -lessl -lxlf90_r

# for glacier
#ATLAS_LIBS = \
#	-L/usr/lib -llapack -lblas \
#	$(RYUON_LIB) -larpack


#------------------------------------------------------------------------------
# other libraries
GSL_INC = `pkg-config --cflags gsl`
GSL_LIB = `pkg-config --libs gsl`
# for robson
#GSL_INC = `$(RYUON_BIN)/gsl-config --cflags`
#GSL_LIB = `$(RYUON_BIN)/gsl-config --libs`
# for glacier
#GSL_INC = `/usr/bin/gsl-config --cflags`
#GSL_LIB = `/usr/bin/gsl-config --libs`

NITSOL_INC = -I/usr/local/include
NITSOL_LIB = -L/usr/local/lib -lnitsol
# for SUSE-Linux
#NITSOL_INC = -I/home/ichiki/RYUON/include
#NITSOL_LIB = -L/home/ichiki/RYUON/lib -lnitsol
# for robson
# for glacier
#NITSOL_INC = $(RYUON_INC)
#NITSOL_LIB = $(RYUON_LIB) -lnitsol

GUILE_INC = `guile-config compile`
GUILE_LIB = `guile-config link`
# for robson
#GUILE_INC = `$(RYUON_BIN)/guile-config compile`
#GUILE_LIB = `$(RYUON_BIN)/guile-config link`
# for glacier
#GUILE_INC = `/usr/bin/guile-config compile`
#GUILE_LIB = `/usr/bin/guile-config link`

NETCDF_INC = 
NETCDF_LIB = -lnetcdf
# for robson
# for glacier
#NETCDF_INC = $(RYUON_INC)
#NETCDF_LIB = $(RYUON_LIB) -lnetcdf


## for OpenGL, GLUT
GL_CFLAGS = -I/usr/X11R6/include
GL_LIBS   = -L/usr/X11R6/lib -lglut
# for SUSE-Linux
# for openSUSE
#GL_CFLAGS =
#GL_LIBS   = -lglut
# for robson
# for glacier
GL_CFLAGS = -I/usr/include
GL_LIBS   = -L/usr/lib -lglut


CFLAGS = -Wall		\
	$(OPTS)		\
	$(RYUON_INC)	\
	$(NETCDF_INC)	\
	$(GUILE_INC)	\
	$(GSL_INC)	\
	$(NITSOL_INC)

LDFLAGS =

LIB_STOKES = \
	-lm		\
	$(RYUON_LIB)	\
	$(NETCDF_LIB)	\
	$(GUILE_LIB)	\
	$(GSL_LIB)	\
	$(NITSOL_LIB)	\
	$(MATH_LIBS)

LIBS = $(LIB_STOKES)


all:	test-stokes	\
	xi3		\
	bench3		\
	stokes3		\
	stvis-glut	\
	stnc


# test-stokes
test-stokes-OBJ = test-stokes.o

test-stokes:	$(test-stokes-OBJ)
	$(LD) $(LDFLAGS) -o test-stokes	\
	$(test-stokes-OBJ)		\
	$(LIBS)


# xi3
xi3-OBJ = xi3.o

xi3.o: xi3.c
	$(CC) $(CFLAGS) $(RYUON_SRC) -c xi3.c

xi3:	$(xi3-OBJ)
	$(LD) $(LDFLAGS) -o xi3	\
	$(xi3-OBJ)		\
	$(LIBS)


## bench3
bench3-OBJ = \
	bench3.o\
	configs.o

bench3:	$(bench3-OBJ)
	$(LD) $(LDFLAGS)	\
	-o bench3		\
	$(bench3-OBJ)		\
	$(LIBS)

## stokes3
stokes3-OBJ =	\
	stokes3.o	\
	file.o

stokes3:	$(stokes3-OBJ)
	$(LD) $(LDFLAGS)	\
	-o stokes3		\
	$(stokes3-OBJ)		\
	$(LIBS)

## stvis-glut
stvis-glut-OBJ =	stvis-glut.o

stvis-glut.o:	stvis-glut.c
	$(CC) $(CFLAGS) $(GL_CFLAGS) -c stvis-glut.c

stvis-glut:	$(stvis-glut-OBJ)
	$(LD) $(LDFLAGS)	\
	-o stvis-glut		\
	$(stvis-glut-OBJ)	\
	$(LIBS)			\
	-L/usr/X11R6/lib -lglut

## stnc
stnc-OBJ =	\
	stnc.o	\
	file.o

stnc:	$(stnc-OBJ)
	$(LD) $(LDFLAGS)\
	-o stnc		\
	$(stnc-OBJ)	\
	$(LIBS)

clean:	
	rm -f *~ *.o	\
	test-stokes	\
	xi3		\
	bench3		\
	stokes3		\
	stvis-glut	\
	stnc
