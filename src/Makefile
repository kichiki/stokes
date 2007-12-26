# $Id: Makefile,v 1.11 2007/12/26 06:56:18 kichiki Exp $

CC           = cc
LD           = $(CC)

GSL_INC = `pkg-config --cflags gsl`
GSL_LIB = `pkg-config --libs gsl`

GUILE_INC = `guile-config compile`
GUILE_LIB = `guile-config link`

RYUON_INC = -I/usr/local/include
RYUON_LIB = -L/usr/local/lib -liter -lstokes
# for openSUSE
#RYUON_INC = 
#RYUON_LIB = -liter -lstokes

# for header file besides libstokes.h
RYUON_SRC = -I../../libstokes/src

## for OpenGL, GLUT
GL_CFLAGS	= -I/usr/X11R6/include
GL_LIBS		= -L/usr/X11R6/lib -lglut
# for openSUSE
#GL_CFLAGS	=
#GL_LIBS		= -lglut

#CFLAGS	     = -Wall -I/usr/local/include -O3 -mtune=pentium
CFLAGS	     = -Wall -I/usr/local/include \
	-pipe -O3 -march=prescott -fomit-frame-pointer \
	-mfpmath=sse -ffast-math \
	$(GUILE_INC)	\
	$(GSL_INC)
# for i1200
#CFLAGS	     = -Wall\
#	-march=pentium3m \
#	-pipe -O3 -fomit-frame-pointer -ffast-math \
#	$(GUILE_INC)	\
#	$(GSL_INC)

#ATLAS_LIBS	= \
#	-lg2c -lblas -llapack\
#	-lcblas -lf77blas -latlas
ATLAS_LIBS	= \
	-L/usr/local/lib/gcc-4.2.3 -lgfortran \
	-L/usr/local/lib -lalapack -lf77blas -lcblas -latlas	\
	-larpack
# for SUSE-Linux
# for openSUSE
#ATLAS_LIBS	= \
#	-lgfortran -lblas -llapack\
#	-larpack

LDFLAGS      =

LIB_STOKES	= \
	-lm $(RYUON_LIB)	\
	$(ATLAS_LIBS)		\
	-lnetcdf		\
	$(GUILE_LIB)		\
	$(GSL_LIB)

LIBS		= $(LIB_STOKES)


all:	test-stokes	\
	xi3		\
	bench3		\
	stokes3		\
	stvis-glut


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

clean:	
	rm -f *~ *.o	\
	test-stokes	\
	xi3		\
	bench3		\
	stokes3		\
	stvis-glut
