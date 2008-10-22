# $Id: Makefile,v 1.16 2008/10/22 06:52:11 kichiki Exp $

CC = cc
RYUON_PREFIX = /usr/local
# for SUSE-Linux
#RYUON_PREFIX = /home/ichiki/RYUON

# for robson gcc3/g77 (32bit)
#CC = /usr/bin/gcc
#RYUON_PREFIX = /home/ichiki/RYUON3-32

# for robson gcc3/g77 (64bit)
#CC = /usr/bin/gcc -m64
#RYUON_PREFIX = /home/ichiki/RYUON3

# for robson gcc4/gfortran (32bit)
#CC = /usr/local/bin/gcc4
#RYUON_PREFIX = /home/ichiki/RYUON4-32

# for robson gcc4/gfortran (64bit)
#CC = /usr/local/bin/gcc4 -m64
#RYUON_PREFIX = /home/ichiki/RYUON4

RYUON_BIN = $(RYUON_PREFIX)/bin
RYUON_INC = -I$(RYUON_PREFIX)/include
RYUON_LIB = -L$(RYUON_PREFIX)/lib -liter -lstokes

LD = $(CC)
# for SUSE-Linux
#LD = $(CC) -Wl,-R$(RYUON_LIB)

#OPTS = -pipe -O3 -fomit-frame-pointer -ffast-math
# for FreeBSD
OPTS = -pipe -O3 -fomit-frame-pointer -ffast-math\
	-march=prescott -mfpmath=sse
# for SUSE-Linux
#OPTS = \
#	-mtune=opteron -m64 \
#	-pipe -O3 -fomit-frame-pointer -ffast-math
# for i1200 (openSUSE)
#OPTS = \
#	-march=pentium3m \
#	-pipe -O3 -fomit-frame-pointer -ffast-math
# for robson
#OPTS = -O3 -mcpu=970 -mtune=970 -maltivec

# for header file besides libstokes.h
RYUON_SRC = -I../../libstokes/src

GSL_INC = `pkg-config --cflags gsl`
GSL_LIB = `pkg-config --libs gsl`
# for robson
#GSL_INC = `$(RYUON_BIN)/gsl-config --cflags`
#GSL_LIB = `$(RYUON_BIN)/gsl-config --libs`

NITSOL_INC = -I/usr/local/include
NITSOL_LIB = -L/usr/local/lib -lnitsol
# for SUSE-Linux
#NITSOL_INC = -I/home/ichiki/RYUON/include
#NITSOL_LIB = -L/home/ichiki/RYUON/lib -lnitsol
# for robson
#NITSOL_INC = $(RYUON_INC)
#NITSOL_LIB = $(RYUON_LIB) -lnitsol

GUILE_INC = `guile-config compile`
GUILE_LIB = `guile-config link`
# for robson
#GUILE_INC = `$(RYUON_BIN)/guile-config compile`
#GUILE_LIB = `$(RYUON_BIN)/guile-config link`

NETCDF_INC = 
NETCDF_LIB = -lnetcdf
# for robson
#NETCDF_INC = $(RYUON_INC)
#NETCDF_LIB = $(RYUON_LIB) -lnetcdf

#ATLAS_LIBS	= \
#	-lg2c -lblas -llapack\
#	-lcblas -lf77blas -latlas\
#	-larpack
# for FreeBSD
ATLAS_LIBS = \
	-L/usr/local/lib/gcc-4.2.5 -lgfortran \
	-L/usr/local/lib -lalapack -lf77blas -lcblas -latlas\
	-larpack
# for SUSE-Linux
# for openSUSE
#ATLAS_LIBS = \
#	-lgfortran -lblas -llapack\
#	-larpack
#ATLAS_LIBS = \
#	-llapack -lptcblas -lptf77blas -latlas\
#	-larpack
# for robson (gcc/g77, 64bit)
#ATLAS_LIBS = -lg2c -L/usr/local/lib64 -lcblas -latlas -L/usr/local/lib64 -larpack
# for robson (gcc/g77, 64bit or 32bit)
#ATLAS_LIBS = \
#	-lg2c \
#	$(RYUON_LIB) -larpack -llapack -lblas
# for robson (gcc4/gfortran, 64bit)
#ATLAS_LIBS = -L/usr/local/lib64 -larpack -lessl -lxlf90_r
# for robson (gcc4/gfortran, 64bit or 32bit)
#ATLAS_LIBS = \
#	-L/usr/local/gcc-4.2.2/lib/gcc/powerpc64-suse-linux/4.2.2 -lgfortran \
#	$(RYUON_LIB) -larpack -llapack -lblas
#ATLAS_LIBS = \
#	-L/usr/local/gcc-4.2.2/lib/gcc/powerpc64-suse-linux/4.2.2 -lgfortran \
#	$(RYUON_LIB) -larpack -llapack -lblas


## for OpenGL, GLUT
GL_CFLAGS = -I/usr/X11R6/include
GL_LIBS   = -L/usr/X11R6/lib -lglut
# for SUSE-Linux
# for openSUSE
#GL_CFLAGS =
#GL_LIBS   = -lglut

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
	$(ATLAS_LIBS)

LIBS = $(LIB_STOKES)


all:	test-stokes	\
	xi3		\
	bench3		\
	stokes3		\
	stvis-glut	\
	stnc


# test-stokes
test-stokes-OBJ =	test-stokes.o

test-stokes:	$(test-stokes-OBJ)
	$(LD) $(LDFLAGS) -o test-stokes	\
	$(test-stokes-OBJ)		\
	$(LIBS)


# xi3
xi3-OBJ =	xi3.o

xi3.o: xi3.c
	$(CC) $(CFLAGS) $(RYUON_SRC) -c xi3.c

xi3:	$(xi3-OBJ)
	$(LD) $(LDFLAGS) -o xi3	\
	$(xi3-OBJ)		\
	$(LIBS)


## bench3
bench3-OBJ =	\
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
