# $Id: Makefile,v 1.15 2008/10/08 03:30:40 kichiki Exp $

CC           = cc
LD           = $(CC)
# for SUSE-Linux
#LD           = $(CC) -Wl,-R/home/ichiki/RYUON/lib

OPTS = \
	-march=prescott -mfpmath=sse \
	-pipe -O3 -fomit-frame-pointer -ffast-math
# for bue9 (SUSE-Linux)
#OPTS = \
#	-mtune=opteron -m64 \
#	-pipe -O3 -fomit-frame-pointer -ffast-math
# for i1200 (openSUSE)
#OPTS = \
#	-march=pentium3m \
#	-pipe -O3 -fomit-frame-pointer -ffast-math


# for header file besides libstokes.h
RYUON_SRC = -I../../libstokes/src

RYUON_INC = -I/usr/local/include
RYUON_LIB = -L/usr/local/lib -liter -lstokes
# for openSUSE
#RYUON_INC = -I/home/ichiki/RYUON/include
#RYUON_LIB = -L/home/ichiki/RYUON/lib -liter -lstokes

GSL_INC = `pkg-config --cflags gsl`
GSL_LIB = `pkg-config --libs gsl`

NITSOL_INC = -I/usr/local/include
NITSOL_LIB = -L/usr/local/lib -lnitsol
# for SUSE-Linux
#NITSOL_INC = -I/home/ichiki/RYUON/include
#NITSOL_LIB = -L/home/ichiki/RYUON/lib -lnitsol

GUILE_INC = `guile-config compile`
GUILE_LIB = `guile-config link`

#ATLAS_LIBS	= \
#	-lg2c -lblas -llapack\
#	-lcblas -lf77blas -latlas
ATLAS_LIBS	= \
	-L/usr/local/lib/gcc-4.2.5 -lgfortran \
	-L/usr/local/lib -lalapack -lf77blas -lcblas -latlas\
	-larpack
# for SUSE-Linux
# for openSUSE
#ATLAS_LIBS	= \
#	-lgfortran -lblas -llapack\
#	-larpack
#ATLAS_LIBS	= \
#	-llapack -lptcblas -lptf77blas -latlas\
#	-larpack

## for OpenGL, GLUT
GL_CFLAGS	= -I/usr/X11R6/include
GL_LIBS		= -L/usr/X11R6/lib -lglut
# for SUSE-Linux
# for openSUSE
#GL_CFLAGS	=
#GL_LIBS		= -lglut

CFLAGS = -Wall		\
	$(OPTS)		\
	$(RYUON_INC)	\
	$(GUILE_INC)	\
	$(GSL_INC)	\
	$(NITSOL_INC)

LDFLAGS      =

LIB_STOKES	= \
	-lm $(RYUON_LIB)	\
	$(ATLAS_LIBS)		\
	-lnetcdf		\
	$(GUILE_LIB)		\
	$(GSL_LIB)		\
	$(NITSOL_LIB)

LIBS		= $(LIB_STOKES)


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
