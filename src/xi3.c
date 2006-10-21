/* tuning program of xi for stokes simulator in 3D for F/FT/FTS versions
 * Copyright (C) 1997-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: xi3.c,v 1.1 2006/10/21 19:03:16 kichiki Exp $
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
#include <stdio.h> /* printf() fprintf() */
#include <stdlib.h> /* exit() */
#include <string.h> /* strcmp() */

#include <libstokes.h> /* struct stokes */
#include <libguile.h> // scm_init_guile()


void
do_xi (struct stokes * sys,
       double tratio, double cutlim,
       int flag_notbl,
       int flag_mat,
       int version)
{
  int np;
  int n;

  double xi;

  double * fts;
  double * uoe;

  double avu;
  double avo;
  double ave;

  int i, j;


  np = sys->np;
  /* allocate memories */
  if      (version == 0) n = np * 3;  // F version
  else if (version == 1) n = np * 6;  // FT version
  else                   n = np * 11; // FT version
  fts = (double *) malloc (sizeof (double) * n);
  uoe = (double *) malloc (sizeof (double) * n);
  if (fts == NULL
      || uoe == NULL)
    {
      fprintf (stderr, "allocation error in main().\n");
      exit (1);
    }

  /* set fts[] */
  for (i = 0; i < n; i ++)
    {
      fts [i] = 1.0;
    }

  xi = xi_by_tratio (sys, tratio);
  stokes_set_xi (sys, xi, cutlim);
  sys->version = version;

  /* call atimes routine */
  if (flag_notbl == 0)
    {
      if (flag_mat != 0)
	{
	  atimes_ewald_3all_matrix (n, fts, uoe, (void *) sys);
	}
      else
	{
	  atimes_ewald_3all (n, fts, uoe, (void *) sys);
	}
    }
  else
    {
      if (flag_mat != 0)
	{
	  atimes_ewald_3all_matrix_notbl (n, fts, uoe, (void *) sys);
	}
      else
	{
	  atimes_ewald_3all_notbl (n, fts, uoe, (void *) sys);
	}
    }

  /* analysis part */
  if (version == 0) // F version
    {
      avu = 0.0;
      for (j = 0; j < np*3; ++j)
	{
	  avu += uoe [j];
	}
      avu /= (double) (np * 3);

      if (flag_notbl == 0)
	{
	  fprintf (stdout, "%f %f %.3f %.3f %.3f %.17e %d %d %d %d\n",
		   tratio, xi,
		   sys->cpu1,
		   sys->cpu2,
		   sys->cpu3,
		   avu,
		   (2*sys->rmaxx+1)*(2*sys->rmaxy+1)*(2*sys->rmaxz+1),
		   (2*sys->kmaxx+1)*(2*sys->kmaxy+1)*(2*sys->kmaxz+1),
		   sys->nr, sys->nk);
	}
      else
	{
	  fprintf (stdout, "%f %f %.3f %.3f %.3f %.17e %d %d\n",
		   tratio, xi,
		   sys->cpu1,
		   sys->cpu2,
		   sys->cpu3,
		   avu,
		   (2*sys->rmaxx+1)*(2*sys->rmaxy+1)*(2*sys->rmaxz+1),
		   (2*sys->kmaxx+1)*(2*sys->kmaxy+1)*(2*sys->kmaxz+1));
	}
    }
  else if (version == 1) // FT version
    {
      avu = 0.0;
      avo = 0.0;
      for (j = 0; j < np; ++j)
	{
	  avu += uoe [j*6 + 0];
	  avu += uoe [j*6 + 1];
	  avu += uoe [j*6 + 2];
	  avo += uoe [j*6 + 3];
	  avo += uoe [j*6 + 4];
	  avo += uoe [j*6 + 5];
	}
      avu /= (double) (np * 3);
      avo /= (double) (np * 3);

      if (flag_notbl == 0)
	{
	  fprintf (stdout, "%f %f %.3f %.3f %.3f %.17e %.17e %d %d %d %d\n",
		   tratio, xi,
		   sys->cpu1,
		   sys->cpu2,
		   sys->cpu3,
		   avu, avo,
		   (2*sys->rmaxx+1)*(2*sys->rmaxy+1)*(2*sys->rmaxz+1),
		   (2*sys->kmaxx+1)*(2*sys->kmaxy+1)*(2*sys->kmaxz+1),
		   sys->nr, sys->nk);
	}
      else
	{
	  fprintf (stdout, "%f %f %.3f %.3f %.3f %.17e %.17e %d %d\n",
		   tratio, xi,
		   sys->cpu1,
		   sys->cpu2,
		   sys->cpu3,
		   avu, avo,
		   (2*sys->rmaxx+1)*(2*sys->rmaxy+1)*(2*sys->rmaxz+1),
		   (2*sys->kmaxx+1)*(2*sys->kmaxy+1)*(2*sys->kmaxz+1));
	}
    }
  else // FTS version
    {
      avu = 0.0;
      avo = 0.0;
      ave = 0.0;
      for (j = 0; j < np; ++j)
	{
	  avu += uoe [j*11 + 0];
	  avu += uoe [j*11 + 1];
	  avu += uoe [j*11 + 2];
	  avo += uoe [j*11 + 3];
	  avo += uoe [j*11 + 4];
	  avo += uoe [j*11 + 5];

	  ave += uoe [j*11 + 6];
	  ave += uoe [j*11 + 7];
	  ave += uoe [j*11 + 8];
	  ave += uoe [j*11 + 9];
	  ave += uoe [j*11 + 10];
	}
      avu /= (double) (np * 3);
      avo /= (double) (np * 3);
      ave /= (double) (np * 5);

      if (flag_notbl == 0)
	{
	  fprintf (stdout, "%f %f %.3f %.3f %.3f %.17e %.17e %.17e %d %d %d %d\n",
		   tratio, xi,
		   sys->cpu1,
		   sys->cpu2,
		   sys->cpu3,
		   avu, avo, ave,
		   (2*sys->rmaxx+1)*(2*sys->rmaxy+1)*(2*sys->rmaxz+1),
		   (2*sys->kmaxx+1)*(2*sys->kmaxy+1)*(2*sys->kmaxz+1),
		   sys->nr, sys->nk);
	}
      else
	{
	  fprintf (stdout, "%f %f %.3f %.3f %.3f %.17e %.17e %.17e %d %d\n",
		   tratio, xi,
		   sys->cpu1,
		   sys->cpu2,
		   sys->cpu3,
		   avu, avo, ave,
		   (2*sys->rmaxx+1)*(2*sys->rmaxy+1)*(2*sys->rmaxz+1),
		   (2*sys->kmaxx+1)*(2*sys->kmaxy+1)*(2*sys->kmaxz+1));
	}
    }

  free (fts);
  free (uoe);
}

void
usage (const char *argv0)
{
  fprintf (stderr, "USAGE\n");
  fprintf (stderr, "%s init-file\n", argv0);
  fprintf (stderr, "\twhere init-file is a SCM file"
	   " (default: xi3.scm)\n\n");
  fprintf (stderr, "Parameters in the init-file:\n");
  fprintf (stderr, "\tversion    : \"F\", \"FT\", or \"FTS\"\n");
  fprintf (stderr,
	   "\tflag-mat   : t for matrix-scheme\n"
	   "\t           : nil for atimes-scheme\n");
  fprintf (stderr,
	   "\tflag-notbl : calculation scheme for the ewald-summation\n"
	   "\t           : t for no-table\n"
	   "\t           : nil for with table\n");
  fprintf (stderr, "\tlattice    : dimensions of the periodic box"
	   " (list or vector of length 3)\n");
  fprintf (stderr, "\tewald-eps  : tolerance value"
	   " for ewald-summation cut-off\n");
  fprintf (stderr, "\tnp         : number of particles\n");
  fprintf (stderr, "\tx          : particle configuration"
	   " (list or vector with length 3*np)\n");
  fprintf (stderr,
	   "\tewald-trs  : (optional) list of tratio parameters\n "
	   "\t           : (list or vector with any length)\n");
}

/* main program */
int
main (int argc, char** argv)
{
  struct stokes * sys = NULL;
  double lat[3];

  int version;
  int flag_mat;
  int flag_notbl;

  int i;
  int np;
  int len;

  double cutlim;
  double tratio;

  char init_file [256];


  /* option analysis */
  strcpy (init_file, "xi3.scm"); // default
  for (i = 1; i < argc; i++)
    {
      if (strcmp (argv [i], "-h") == 0 ||
	  strcmp (argv [i], "--help") == 0)
	{
	  usage (argv[0]);
	  exit (1);
	}
      else
	{
	  strcpy (init_file, argv [i]);
	}
    }
  
  /* parameter set */
  scm_init_guile(); // start the Guile interpreter
  scm_c_primitive_load (init_file); // load initialize script

  // version
  version = 0;
  char * str_version;
  str_version = guile_get_string ("version");
  if (strcmp (str_version, "F") == 0)
    {
      version = 0; // F version
    }
  else if (strcmp (str_version, "FT") == 0)
    {
      version = 1; // FT version
    }
  else if (strcmp (str_version, "FTS") == 0)
    {
      version = 2; // FTS version
    }
  else
    {
      fprintf (stderr, "invalid version %s", str_version);
      exit (1);
    }
  free (str_version);

  // flag-mat
  flag_mat = 0;
  if (guile_get_bool ("flag-mat") != 0) // TRUE
    {
      flag_mat = 1;
    }

  // flag-notbl
  flag_notbl = 0;
  if (guile_get_bool ("flag-notbl") != 0) // TRUE
    {
      flag_notbl = 1;
    }

  if (guile_get_doubles ("lattice", 3, lat) != 1) // FALSE
    {
      fprintf (stderr, "lattice is not defined\n");
      exit (1);
    }

  // np
  np         = guile_get_int    ("np",         0);
  // ewald-eps
  cutlim     = guile_get_double ("ewald-eps",  0.0);


  /* initialization */
  sys = stokes_init ();
  stokes_set_np (sys, np, np);
  stokes_set_l (sys, lat[0], lat[1], lat[2]);

  // x -- configuration of particles
  if (guile_get_doubles ("x", np*3, sys->pos) != 1) // FALSE
    {
      fprintf (stderr, "x is not defined\n");
      exit (1);
    }


  /* write header */
  if (version == 0)      fprintf (stdout, "# F version");
  else if (version == 1) fprintf (stdout, "# FT version");
  else                   fprintf (stdout, "# FTS version");

  if (flag_notbl == 0)   fprintf (stdout, " table");
  else                   fprintf (stdout, " no-table");

  if (flag_mat == 0)     fprintf (stdout, " atimes\n");
  else                   fprintf (stdout, " matrix\n");


  /* main loop */
  // ewald-trs
  double *trs = NULL;
  trs = guile_get_doubles_ ("ewald-trs");
  if (trs == NULL)
    {
      /*tratio = 1.0;*/
      tratio = 0.1;
      for (i = 1; i < 100; i++)
	{
	  tratio *= 1.1;

	  do_xi (sys, tratio, cutlim,
		 flag_notbl, flag_mat, version);
	}
    }
  else
    {
      len = guile_get_length ("ewald-trs");
      for (i = 0; i < len; i ++)
	{
	  do_xi (sys, trs[i], cutlim,
		 flag_notbl, flag_mat, version);
	}
      free (trs);
    }

  stokes_free (sys);

  return 0;
}
