/* 3D Display Program by OpenGL
 * Copyright (C) 1997-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: stvis-glut.c,v 1.2 2007/05/15 07:31:47 kichiki Exp $
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
#include <math.h>
#include <stdio.h>	/* printf() fprintf() */
#include <string.h>
#include <netcdf.h> // nc_get_var1_double()
#include <libstokes.h> // struct stokes_nc

#include <GL/glut.h>
#include <stdlib.h>


/** global variables **/
GLenum clearMask = GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT;
GLenum smooth = GL_FALSE;
GLenum lighting = GL_TRUE;
GLenum depth = GL_TRUE;
GLenum animate = GL_FALSE;
GLenum perspective = GL_TRUE; // perspective

struct stokes_nc *nc;
double lattice[3];
int nstep; // frequency to display
int step;  // current step
double *x;
double *xf0;


/*  Initialize material property and light source.
 */
void
Init (void)
{
  // OpenGL is in left-handed coordinate system
  // +x : right
  // +y : up
  // -z : looking at (+z is camera's position)
  gluLookAt(0.5f, 0.5f, 2.0f,  // camera position
	    0.5f, 0.5f, 0.0f,  // position the camera is looking at
	    0.0f, 1.0f, 0.0f); // normalized up vector

  GLfloat light_ambient[]  = {0.0, 0.0, 0.0, 1.0};
  GLfloat light_diffuse[]  = {1.0, 1.0, 1.0, 1.0};
  GLfloat light_specular[] = {1.0, 1.0, 1.0, 1.0};
  /* light_position is NOT default value */
  GLfloat light_position[] = {1.0, 1.0, 1.0, 1.0f};
  //GLfloat light_position[] = {1.0, 1.0, 1.0, 0.0f};
  /* horizontal(x),vertical(z),tangential-horizontal(y),distance */

  glLightfv (GL_LIGHT0, GL_AMBIENT,  light_ambient);
  glLightfv (GL_LIGHT0, GL_DIFFUSE,  light_diffuse);
  glLightfv (GL_LIGHT0, GL_SPECULAR, light_specular);
  glLightfv (GL_LIGHT0, GL_POSITION, light_position);
    
  glFrontFace (GL_CW);
  glEnable (GL_LIGHTING);
  glEnable (GL_LIGHT0);
  glDepthFunc (GL_LESS);
  glEnable (GL_DEPTH_TEST);
  /*glShadeModel (GL_FLAT);*/
  glShadeModel (GL_SMOOTH);

  glEnable (GL_AUTO_NORMAL);
  glEnable (GL_NORMALIZE);
}

void
Reshape (int width, int height)
{
  glViewport (0, 0, (GLint)width, (GLint)height);

  glMatrixMode (GL_PROJECTION);
  glLoadIdentity ();

  if (perspective) // GL_TRUE
    {
      gluPerspective (45.0f,        // fovy (angle [degree])
		      (GLfloat)width / (GLfloat)height, // aspect ratio
		      1.0f,         // near clipping plane
		      100.0f);      // far clipping plane
    }
  else if (width <= height)
    {
      glOrtho (-1.0, // left
	       1.0,  // right
	       -1.0*(GLfloat)height/(GLfloat)width, // bottom
	       1.0*(GLfloat)height/(GLfloat)width,  // top
	       -10.0, // near ? (-z coordinate?)
	       10.0); // far  ? (-z coordinate?)
    }
  else 
    {
      glOrtho (-1.0*(GLfloat)width/(GLfloat)height, // left
	       1.0*(GLfloat)width/(GLfloat)height,  // right
	       -1.0,  // bottom
	       1.0,   // top
	       -10.0, // near ? (-z coordinate?)
	       10.0); // far  ? (-z coordinate?)
    }

  glMatrixMode (GL_MODELVIEW);
}

void
Animate (void)
{
  extern struct stokes_nc *nc;
  extern int nstep;
  extern int step;
  extern double *x;
  extern double *xf0;

  int i;

  double l;
  double lx, ly, lz;
  double rad;


  if (step >= nc->ntime) // beyond the end
    {
      step = nc->ntime -1;
      // halt the animation
      glutIdleFunc(NULL);
    }
  if (step < 0) // before the beginning
    {
      step = 0;
      // halt the animation
      glutIdleFunc(NULL);
    }

  lx = lattice[0];
  ly = lattice[1];
  lz = lattice[2];
  if (lx > ly) l = lx;
  else         l = ly;
  rad = 1.0 / l;

  stokes_nc_get_data (nc, "x", step, x);

  GLfloat amb_ruby[]  = {0.1745, 0.01175, 0.01175,1.0};
  GLfloat dif_ruby[]  = {0.61424, 0.04136, 0.04136,1.0};
  GLfloat spec_ruby[] = {0.727811, 0.626959, 0.626959,1.0};
  GLfloat shine_ruby  = 0.6;

  /*
  GLfloat amb_gold[]  = {0.24725, 0.1995, 0.0745,1.0};
  GLfloat dif_gold[]  = {0.75164, 0.60648, 0.22648,1.0};
  GLfloat spec_gold[] = {0.628281, 0.555802, 0.366065,1.0};
  GLfloat shine_gold  = 0.4;
  */

  GLfloat amb_silver[]  = {0.19225, 0.19225, 0.19225,1.0};
  GLfloat dif_silver[]  = {0.50754, 0.50754, 0.50754,1.0};
  GLfloat spec_silver[] = {0.508273, 0.508273, 0.508273,1.0};
  GLfloat shine_silver  = 0.4;

  GLfloat amb_emerald[]  = {0.0215, 0.1745, 0.0215,1.0};
  GLfloat dif_emerald[]  = {0.07568, 0.61424, 0.07568,1.0};
  GLfloat spec_emerald[] = {0.633, 0.727811, 0.633,1.0};
  GLfloat shine_emerald  = 0.6;

  /*
  GLfloat amb_cyan[]  = {0.0, 0.1, 0.06,1.0};
  GLfloat dif_cyan[]  = {0.0, 0.50980392, 0.50980392,1.0};
  GLfloat spec_cyan[] = {0.50196078, 0.50196078, 0.50196078,1.0};
  GLfloat shine_cyan  = 0.25;
  */

  glClear(clearMask);

  /* frame plane (lx-lz)*/
  glPushMatrix ();
  glBegin(GL_QUADS);
  glMaterialfv (GL_FRONT, GL_AMBIENT, amb_silver);
  glMaterialfv (GL_FRONT, GL_DIFFUSE, dif_silver);
  glMaterialfv (GL_FRONT, GL_SPECULAR, spec_silver);
  glMaterialf (GL_FRONT, GL_SHININESS, shine_silver*128.0);
  glNormal3f (0.0, 0.0,  1.0);
  glVertex3f (0.0, lz/l, -ly/l);
  glVertex3f (1.0, lz/l, -ly/l);
  glVertex3f (1.0, 0.0,  -ly/l);
  glVertex3f (0.0, 0.0,  -ly/l);
  glEnd();
  glPopMatrix ();

  /* frame plane (lx-ly)*/
  glPushMatrix ();
  glBegin(GL_QUADS);
  glMaterialfv (GL_FRONT, GL_AMBIENT, amb_silver);
  glMaterialfv (GL_FRONT, GL_DIFFUSE, dif_silver);
  glMaterialfv (GL_FRONT, GL_SPECULAR, spec_silver);
  glMaterialf (GL_FRONT, GL_SHININESS, shine_silver*128.0);
  glNormal3f (0.0, 0.0,  1.0);
  glVertex3f (0.0, 0.0,  0.0);
  glVertex3f (1.0, 0.0,  0.0);
  glVertex3f (1.0, 0.0, -ly/l);
  glVertex3f (0.0, 0.0, -ly/l);
  glEnd();
  glPopMatrix ();

  /* frame
     glPushMatrix();
     glDisable (GL_LIGHTING);
     glColor3f (0.0, 1.0, 1.0);
     glTranslatef (0.5, 0.5*lz/l, 0.5*ly/l);
     auxWireBox(1.0, lz/l, ly/l);
     glEnable (GL_LIGHTING);
     glPopMatrix ();
  */

  // Note that simulation is done in the right-handed coordinate system
  // +x : right                == +x in OpenGL
  // +y : looking at direction == -z in OpenGL
  // +z : up                   == +y in OpenGL

  /* mobile particles */
  glMaterialfv (GL_FRONT, GL_AMBIENT, amb_emerald);
  glMaterialfv (GL_FRONT, GL_DIFFUSE, dif_emerald);
  glMaterialfv (GL_FRONT, GL_SPECULAR, spec_emerald);
  glMaterialf (GL_FRONT, GL_SHININESS, shine_emerald*128.0);
  // nc->np is only "mobile" particles
  for (i = 0; i < nc->np; i++)
    {
      glPushMatrix ();
      glTranslatef (+x[i*3+0]/l, // +x component
		    +x[i*3+2]/l, // +z component
		    -x[i*3+1]/l);// -y component
      glutSolidSphere (rad, 20, 20);
      glPopMatrix ();
    }

  /* fixed particles */
  if (nc->flag_xf0 != 0)
    {
      glMaterialfv (GL_FRONT, GL_AMBIENT, amb_ruby);
      glMaterialfv (GL_FRONT, GL_DIFFUSE, dif_ruby);
      glMaterialfv (GL_FRONT, GL_SPECULAR, spec_ruby);
      glMaterialf (GL_FRONT, GL_SHININESS, shine_ruby*128.0);
      // nc->npf is only "fixed" particles
      for (i = 0; i < nc->npf; i++)
	{
	  glPushMatrix ();
	  glTranslatef (+xf0[i*3+0]/l, // +x component
			+xf0[i*3+2]/l, // +z component
			-xf0[i*3+1]/l);// -y component
	  glutSolidSphere (rad, 20, 20);
	  glPopMatrix ();
	}
    }

  //glFlush (); // this is done automatically by glutSwapBuffers()
  glutSwapBuffers ();

  // for the next step, skip (nstep-1) steps
  step += nstep;
}

void idle (void)
{
  glutPostRedisplay ();
}

void
keyboard (unsigned char key, int x, int y)
{
  switch (key)
    {
    case '\033':  /* '\033' は ESC の ASCII コード */
      exit(0);
      break;
    case ' ': // start/stop
      animate = !animate;
      if (animate) // TRUE
	glutIdleFunc(idle);
      else // FALSE
	glutIdleFunc(NULL);
      break;
    case 'b': // backward play
      nstep *= -1;
      break;
    case 'e': // end
      glutIdleFunc(NULL);
      break;
    case 'r': // rewind
      step = 0;
      break;
    case 's': /* s -> Toggle smooth/flat */
      smooth = !smooth;
      if (smooth) {
	glShadeModel(GL_SMOOTH);
      } else {
	glShadeModel(GL_FLAT);
      }
      if (!animate) // halting
	{
	  step -= nstep; // to prevent advancing
	  Animate();
	}
      break;
    case 'l': /* l -> lighting */
      lighting = !lighting;
      if (lighting) {
	Init();
	/*
	glEnable (GL_LIGHTING);
	glEnable (GL_LIGHT0);
	glEnable (GL_COLOR_MATERIAL);
	*/
      } else {
	glDisable (GL_LIGHTING);
	glDisable (GL_LIGHT0);
	glDisable (GL_COLOR_MATERIAL);
      }
      if (!animate) // halting
	{
	  step -= nstep; // to prevent advancing
	  Animate();
	}
      break;
    case 'd': /* d -> depth */
      depth = !depth;
      if (depth) {
	glEnable(GL_DEPTH_TEST);
	clearMask |= GL_DEPTH_BUFFER_BIT;
      } else {
	glDisable(GL_DEPTH_TEST);
	clearMask &= ~GL_DEPTH_BUFFER_BIT;
	}
      if (!animate) // halting
	{
	  step -= nstep; // to prevent advancing
	  Animate();
	}
      break;
    case 'i': /* i -> zoom in */
      glTranslatef(0.5, 0.5, 0.0);
      glScalef(1.1, 1.1, 1.1);
      glTranslatef(-0.5, -0.5, 0.0);
      if (!animate) // halting
	{
	  step -= nstep; // to prevent advancing
	  Animate();
	}
      break;
    case 'o': /* o -> zoom out */
      glTranslatef(0.5, 0.5, 0.0);
      glScalef(0.9, 0.9, 0.9);
      glTranslatef(-0.5, -0.5, 0.0);
      if (!animate) // halting
	{
	  step -= nstep; // to prevent advancing
	  Animate();
	}
      break;
    case 'p': /* perspective, orthographic */
      perspective = !perspective;
      Reshape (glutGet (GLUT_WINDOW_WIDTH),
	       glutGet (GLUT_WINDOW_HEIGHT));
      if (!animate) // halting
	{
	  step -= nstep; // to prevent advancing
	  Animate();
	}
      break;
    default:
      return;
    }
}

void
specialkey (int key, int x, int y)
{
  switch (key)
    {
    case GLUT_KEY_LEFT:
      glTranslatef (0.5, 0.5, 0.0);
      glRotatef (-3.0, 0, 1, 0);
      glTranslatef (-0.5, -0.5, 0.0);
      if (!animate) // halting
	{
	  step -= nstep; // to prevent advancing
	  Animate();
	}
      break;
    case GLUT_KEY_RIGHT:
      glTranslatef (0.5, 0.5, 0.0);
      glRotatef (3.0, 0, 1, 0);
      glTranslatef (-0.5, -0.5, 0.0);
      if (!animate) // halting
	{
	  step -= nstep; // to prevent advancing
	  Animate();
	}
      break;
    case GLUT_KEY_DOWN:
      glTranslatef (0.5, 0.5, 0.0);
      glRotatef (-3.0, 1, 0, 0);
      glTranslatef (-0.5, -0.5, 0.0);
      if (!animate) // halting
	{
	  step -= nstep; // to prevent advancing
	  Animate();
	}
      break;
    case GLUT_KEY_UP:
      glTranslatef (0.5, 0.5, 0.0);
      glRotatef (3.0, 1, 0, 0);
      glTranslatef (-0.5, -0.5, 0.0);
      if (!animate) // halting
	{
	  step -= nstep; // to prevent advancing
	  Animate();
	}
      break;
    }
}


void
usage (const char * argv0)
{
  fprintf (stderr, "VERSION:\n");
  fprintf (stderr, "$Id: stvis-glut.c,v 1.2 2007/05/15 07:31:47 kichiki Exp $\n\n");
  fprintf (stderr, "USAGE:\n");
  fprintf (stderr, "%s [OPTIONS]\n", argv0);
  fprintf (stderr, "OPTIONS:\n");
  fprintf (stderr, "\t-f  : stokes-nc file to show\n");
  fprintf (stderr, "\t-ns : display every ns steps"
	   " (default: 1, meaning all steps)\n\n");
  fprintf (stderr, "KEY BINDINGS\n");
  fprintf (stderr, "\tESC       : quit\n");
  fprintf (stderr, "\tSPC       : start/stop\n");
  fprintf (stderr, "\tb         : forward/backward\n");
  fprintf (stderr, "\tr         : rewind to 0 step\n");
  fprintf (stderr, "\ti         : zoom-in\n");
  fprintf (stderr, "\to         : zoom-out\n");
  fprintf (stderr, "\tLEFT/RIGHT: rotate around y (vertical) axis\n");
  fprintf (stderr, "\tUP/DOWN   : rotate around x (horizontal) axis\n");
  fprintf (stderr, "\tl         : lighting (on/off)\n");
  fprintf (stderr, "\td         : depth (on/off)\n");
}

/*  Main Loop
 *  Open window with initial window size, title bar, 
 *  RGBA display mode, and handle input events.
 */
int
main (int argc, char** argv)
{
  extern struct stokes_nc *nc;
  extern double lattice[3];
  extern int nstep;

  int i;

  nc = NULL;
  nstep = 1;
  //nm = 0;
  for (i = 1; i < argc; i++)
    {
      if (strcmp (argv[i], "-f") == 0)
	{
	  if (i + 1 < argc)
	    {
	      nc = stokes_nc_open (argv[++i]);
	    }
	}
      else if (strcmp (argv[i], "-n") == 0)
	{
	  if (i + 1 < argc)
	    {
	      nstep = atoi (argv[++i]);
	    }
	}
      else
	{
	  usage (argv[0]);
	  exit (1);
	}
    }

  //if (nm == 0)
  if (nc == NULL)
    {
      fprintf (stderr, "You Must Believe In Spring!\n");
      exit (1);
    }

  if (nstep <= 0)
    {
      fprintf (stderr, "Invalid n\n");
      exit (1);
    }

  //stokes_nc_print_actives (nc, stdout);
  x = (double *) malloc (sizeof (double) * nc->np * nc->nvec);
  if (x == NULL)
    {
      fprintf (stderr, "allocation error\n");
      exit (1);
    }
  if (nc->flag_xf0 != 0)
    {
      xf0 = (double *) malloc (sizeof (double) * nc->npf * nc->nvec);
      if (xf0 == NULL)
	{
	  fprintf (stderr, "allocation error\n");
	  exit (1);
	}
      stokes_nc_get_data0 (nc, "xf0", xf0);
    }
  stokes_nc_get_array1d (nc, "l", lattice);

  glutInitWindowPosition (0, 0);
  glutInitWindowSize (500, 500);

  glutInit (&argc, argv);
  glutInitDisplayMode (GLUT_DEPTH | GLUT_RGBA | GLUT_DOUBLE);

  glutCreateWindow ("TOY - OpenGL");

  Init ();

  glutDisplayFunc (Animate);
  glutReshapeFunc (Reshape);
  //glutMouseFunc (mouse);
  glutKeyboardFunc(keyboard);
  glutSpecialFunc (specialkey);

  glutMainLoop ();

  free (x);
  if (xf0 != NULL) free (xf0);
  stokes_nc_free (nc);

  return 0;
}
