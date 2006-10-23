/* stokesian dynamics simulator under the periodic boundary condition
 * Copyright (C) 1997-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: stokes3.c,v 1.2 2006/10/23 17:12:37 kichiki Exp $
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

#include <libiter.h> /* iter_init() */
#include <libstokes.h> /* struct stokes */
#include <netcdf.h> // nc_sync()
#include <libguile.h> // scm_init_guile()


void
usage (const char *argv0)
{
  fprintf (stderr, "USAGE\n");
  fprintf (stderr, "%s init-file\n", argv0);
  fprintf (stderr, "\twhere init-file is a SCM file"
	   " (default: stokes3.scm)\n\n");
  fprintf (stderr, "Parameters in the init-file:\n");
  fprintf (stderr, "\toutfile    : filename for NetCDF output\n");
  fprintf (stderr, "\tversion    : \"F\", \"FT\", or \"FTS\"\n");
  fprintf (stderr,
	   "\tflag-mat   : #t for matrix-scheme\n"
	   "\t           : #f for atimes-scheme\n");
  fprintf (stderr,
	   "\tflag-lub   : #t for with-lubrication\n"
	   "\t           : #f for no-lubrication\n");
  fprintf (stderr, "\tnp         : number of ALL particles\n");
  fprintf (stderr, "\tnm         : number of mobile particles\n");
  fprintf (stderr, "\tnloop      : number of loops for (10*dt)\n");
  fprintf (stderr, "\tdt         : time interval\n");
  fprintf (stderr, "\tstokes     : effective stokes number\n");
  fprintf (stderr, "\tewald-tr   : time ratio Tr/Tk for ewald summation\n");
  fprintf (stderr, "\tewald-eps  : tolerance value"
	   " for ewald-summation cut-off\n");
  fprintf (stderr, "\tUi         : imposed translational velocity"
	   " (list or vector of length 3)\n");
  fprintf (stderr, "\tOi         : imposed angular velocity"
	   " (list or vector of length 3)\n");
  fprintf (stderr, "\tEi         : imposed strain"
	   " (list or vector of length 5)\n");
  fprintf (stderr, "\tF0         : applied force"
	   " (list or vector of length 3)\n");
  fprintf (stderr, "\tT0         : applied torque"
	   " (list or vector of length 3)\n");
  fprintf (stderr, "\tlattice    : dimensions of the periodic box"
	   " (list or vector of length 3)\n");
  fprintf (stderr, "\tx          : particle configuration"
	   " (list or vector with length 3*np)\n");
}

static void
solve_mix_3all (struct stokes * sys,
		int flag_lub,
		int flag_mat,
		const double *fm, const double *tm, const double *em,
		const double *uf, const double *of, const double *ef,
		double *um, double *om, double *sm,
		double *ff, double *tf, double *sf)
{
  if (sys->version == 0) // F version
    {
      if (flag_lub == 0)
	{
	  if (flag_mat == 0)
	    {
	      solve_mix_ewald_3f
		(sys,
		 fm, uf,
		 um, ff);
	    }
	  else
	    {
	      solve_mix_ewald_3f_matrix
		(sys,
		 fm, uf,
		 um, ff);
	    }
	}
      else
	{
	  if (flag_mat == 0)
	    {
	      solve_mix_lub_ewald_3f
		(sys,
		 fm, uf,
		 um, ff);
	    }
	  else
	    {
	      solve_mix_lub_ewald_3f_matrix
		(sys,
		 fm, uf,
		 um, ff);
	    }
	}
    }
  else if (sys->version == 1) // FT version
    {
      if (flag_lub == 0)
	{
	  solve_mix_ewald_3ft
	    (sys,
	     fm, tm, uf, of,
	     um, om, ff, tf);
	}
      else
	{
	  solve_mix_lub_ewald_3ft
	    (sys,
	     fm, tm, uf, of,
	     um, om, ff, tf);
	}
    }
  else // FTS version
    {
      if (flag_lub == 0)
	{
	  solve_mix_ewald_3fts
	    (sys,
	     fm, tm, em, uf, of, ef,
	     um, om, sm, ff, tf, sf);
	}
      else
	{
	  solve_mix_lub_ewald_3fts
	    (sys,
	     fm, tm, em, uf, of, ef,
	     um, om, sm, ff, tf, sf);
	}
    }
}

/* main program */
int
main (int argc, char** argv)
{
  struct stokes * sys = NULL;
  double Ui[3];
  double Oi[3];
  double Ei[5];
  double F0[3];
  double T0[3];
  double lat[3];
  double xi;
  double ewald_eps;
  double ewald_tr;

  int i;
  int i3, i5;
  int np, np3;
  int nm, nm3, nm5;
  int nf, nf3, nf5;
  int nloop;

  int j, l;

  double dt;
  double t;
  double stokes;

  struct stokes_nc * nc = NULL;

  double * fm, * tm, * sm;
  double * um, * om, * em;
  double *ff = NULL;
  double *tf = NULL;
  double *sf = NULL;
  double *uf = NULL;
  double *of = NULL;
  double *ef = NULL;

  double * x; /* position to write into the output file */
  double * angle; /* angle to write into the output file */

  /* variables for non-zero 'stokes' */
  double * v = NULL; /* velocity for non-zero 'stokes' */
  double * av = NULL; /* angular velocity for non-zero 'stokes' */
  double dts = 0.0;
  //double dtse2;
  //double dte2;
  //double dtse;
  //double dte;
  double edts = 0.0;
  double edts2 = 0.0;

  char init_file [256];
  int version;
  int flag_lub;
  int flag_mat;


  /* option analysis */
  strcpy (init_file, "stokes3.scm"); // default
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

  // flag-lub
  flag_lub = 0;
  if (guile_get_bool ("flag-lub") != 0) // TRUE
    {
      flag_lub = 1;
    }


  // outfile
  char * out_file;
  out_file = guile_get_string ("outfile");

  // other parameters
  np         = guile_get_int    ("np",         0);
  nm         = guile_get_int    ("nm",         0);
  nloop      = guile_get_int    ("nloop",      0);
  dt         = guile_get_double ("dt",         0.0);
  stokes     = guile_get_double ("stokes",     0.0);
  ewald_tr   = guile_get_double ("ewald-tr",   0.0);
  ewald_eps  = guile_get_double ("ewald-eps",  0.0);

  /* imposed flow */
  // Ui
  if (guile_get_doubles ("Ui", 3, Ui) != 1) // FALSE
    {
      fprintf (stderr, "Ui is not defined\n");
      exit (1);
    }
  // Oi
  if (guile_get_doubles ("Oi", 3, Oi) != 1) // FALSE
    {
      fprintf (stderr, "Oi is not defined\n");
      exit (1);
    }
  // Ei
  if (guile_get_doubles ("Ei", 5, Ei) != 1) // FALSE
    {
      fprintf (stderr, "Ei is not defined\n");
      exit (1);
    }

  /* applied force */
  // F0
  if (guile_get_doubles ("F0", 3, F0) != 1) // FALSE
    {
      fprintf (stderr, "F0 is not defined\n");
      exit (1);
    }
  // T0
  if (guile_get_doubles ("T0", 3, T0) != 1) // FALSE
    {
      fprintf (stderr, "T0 is not defined\n");
      exit (1);
    }
  // lattice
  if (guile_get_doubles ("lattice", 3, lat) != 1) // FALSE
    {
      fprintf (stderr, "lattice is not defined\n");
      exit (1);
    }

  /* initialization */
  sys = stokes_init ();
  sys->version = version;
  stokes_set_np (sys, np, nm);
  stokes_set_l (sys, lat[0], lat[1], lat[2]);
  xi = xi_by_tratio (sys, ewald_tr);
  stokes_set_xi (sys, xi, ewald_eps);

  sys->lubcut = 2.0000000001;
  stokes_set_iter (sys, "gmres", 2000, 20, 1.0e-6, 1, stderr);

  np3 = np * 3;
  nm3 = nm * 3;
  nm5 = nm * 5;
  nf = np - nm;
  nf3 = nf * 3;
  nf5 = nf * 5;

  stokes_set_Ui (sys, Ui[0], Ui[1], Ui[2]);
  stokes_set_Oi (sys, Oi[0], Oi[1], Oi[2]);
  stokes_set_Ei (sys, Ei[0], Ei[1], Ei[2], Ei[3], Ei[4]);

  x     = (double *) malloc (sizeof (double) * np3);
  angle = (double *) malloc (sizeof (double) * np3);

  fm = (double *) malloc (sizeof (double) * nm3);
  tm = (double *) malloc (sizeof (double) * nm3);
  sm = (double *) malloc (sizeof (double) * nm5);
  um = (double *) malloc (sizeof (double) * nm3);
  om = (double *) malloc (sizeof (double) * nm3);
  em = (double *) malloc (sizeof (double) * nm5);
  if (nf > 0)
    {
      ff = (double *) malloc (sizeof (double) * nf3);
      tf = (double *) malloc (sizeof (double) * nf3);
      sf = (double *) malloc (sizeof (double) * nf5);
      uf = (double *) malloc (sizeof (double) * nf3);
      of = (double *) malloc (sizeof (double) * nf3);
      ef = (double *) malloc (sizeof (double) * nf5);
    }

  // x -- particle configuration
  if (guile_get_doubles ("x", np3, x) != 1) // FALSE
    {
      fprintf (stderr, "x is not defined\n");
      exit (1);
    }

  if (stokes > 0.0)
    {
      dts   = dt / stokes;
      //dtse2 = dts * exp (dts / 2.0);
      //dte2  = dt * exp (- dts / 2.0);
      //dtse  = dts * exp (dts);
      //dte   = dt * exp (- dts);
      edts  = exp (- dts);
      edts2 = - stokes * (edts - 1.0);

      v  = (double *) malloc (sizeof (double) * np3);
      av = (double *) malloc (sizeof (double) * np3);

      // v -- initial velocity
      if (guile_get_doubles ("v", np3, x) != 1) // FALSE
	{
	  for (i = 0; i < np3; i ++)
	    {
	      v [i] = 0.0;
	    }
	}
    }

  /* define given parameter */
  for (i = 0; i < nm; ++i)
    {
      i3 = i * 3;
      i5 = i * 5;

      fm [i3 + 0] = F0[0];
      fm [i3 + 1] = F0[1];
      fm [i3 + 2] = F0[2];

      tm [i3 + 0] = T0[0];
      tm [i3 + 1] = T0[1];
      tm [i3 + 2] = T0[2];

      em [i5 + 0] = 0.0;
      em [i5 + 1] = 0.0;
      em [i5 + 2] = 0.0;
      em [i5 + 3] = 0.0;
      em [i5 + 4] = 0.0;
    }
  for (i = 0; i < nf; ++i)
    {
      i3 = i * 3;
      i5 = i * 5;

      uf [i3 + 0] = 0.0;
      uf [i3 + 1] = 0.0;
      uf [i3 + 2] = 0.0;

      of [i3 + 0] = 0.0;
      of [i3 + 1] = 0.0;
      of [i3 + 2] = 0.0;

      ef [i5 + 0] = 0.0;
      ef [i5 + 1] = 0.0;
      ef [i5 + 2] = 0.0;
      ef [i5 + 3] = 0.0;
      ef [i5 + 4] = 0.0;
    }


  /* initialize NetCDF for output */
  if (nf == 0)
    {
      if (version == 0) // F version
	{
	  nc = stokes_nc_mob_f_init (out_file, np);
	}
      else if (version == 1) // FT version
	{
	  nc = stokes_nc_mob_ft_init (out_file, np);
	}
      else // FTS version
	{
	  nc = stokes_nc_mob_fts_init (out_file, np);
	}
    }
  else
    {
      if (version == 0) // F version
	{
	  nc = stokes_nc_mob_fix_f_i0_init (out_file, nm, nf);
	}
      else if (version == 1) // FT version
	{
	  nc = stokes_nc_mob_fix_ft_i0_init (out_file, nm, nf);
	}
      else // FTS version
	{
	  nc = stokes_nc_mob_fix_fts_i0_init (out_file, nm, nf);
	}
      stokes_nc_set_xf0 (nc, x + nm*3);
    }
  free (out_file);

  // set parameters for NetCDF
  stokes_nc_set_l (nc, lat);
  if (version == 0) // F version
    {
      stokes_nc_set_f0 (nc, fm);
      if (nf > 0)
	{
	  stokes_nc_set_uf0 (nc, uf);
	}
    }
  else if (version == 0) // FT version
    {
      stokes_nc_set_f0 (nc, fm);
      stokes_nc_set_t0 (nc, tm);
      if (nf > 0)
	{
	  stokes_nc_set_uf0 (nc, uf);
	  stokes_nc_set_of0 (nc, of);
	}
    }
  else // FTS version
    {
      stokes_nc_set_f0 (nc, fm);
      stokes_nc_set_t0 (nc, tm);
      stokes_nc_set_e0 (nc, em);
      if (nf > 0)
	{
	  stokes_nc_set_uf0 (nc, uf);
	  stokes_nc_set_of0 (nc, of);
	  stokes_nc_set_ef0 (nc, ef);
	}
    }


  /* mail loop */
  t = 0.0;
  for (l = 0; l < nloop; l++)
    {
      fprintf (stdout, "%d steps\n", l);

      for (j = 0; j < 10; j++)
	{
	  stokes_set_pos (sys, x);

	  solve_mix_3all (sys, flag_lub, flag_mat,
			  fm, tm, em, uf, of, ef,
			  um, om, sm, ff, tf, sf);

	  /* Euler method */
	  if (stokes == 0.0)
	    {
	      for (i = 0; i < nm3; ++i)
		{
		  x [i] += um [i] * dt;
		  angle [i] += om [i] * dt;
		}
	    }
	  else
	    {
	      for (i = 0; i < nm3; ++i)
		{
		  x [i] += um [i] * dt + (v [i] - um [i]) * edts2;
		  angle [i] += om [i] * dt + (av [i] - om [i]) * edts2;

		  v [i] = um [i] + (v [i] - um [i]) * edts;
		  av [i] = om [i] + (av [i] - om [i]) * edts;
		}

	      collide_particles (sys, x, v, 1.0);
	      collide_wall_z (sys, x, v, 1.0, 0.0, 0.0);
	    }

	  /* 2D trick
	  for (i = 0; i < nm; i ++)
	    {
	      iy = i * 3 + 1;
	      x [iy] = 1.0;
	    } */

	  check_periodic (sys, x);
	  check_angle (sys, angle);

	  t += dt;
	}

      /* output the results */
      stokes_nc_set_time (nc, l, t);
      stokes_nc_set_x (nc, l, x);
      if (version == 0) // F version
	{
	  stokes_nc_set_u (nc, l, um);
	  if (nf > 0)
	    {
	      stokes_nc_set_ff (nc, l, ff);
	    }
	}
      else if (version == 0) // FT version
	{
	  stokes_nc_set_u (nc, l, um);
	  stokes_nc_set_o (nc, l, om);
	  if (nf > 0)
	    {
	      stokes_nc_set_ff (nc, l, ff);
	      stokes_nc_set_tf (nc, l, tf);
	    }
	}
      else // FTS version
	{
	  stokes_nc_set_u (nc, l, um);
	  stokes_nc_set_o (nc, l, om);
	  stokes_nc_set_s (nc, l, sm);
	  if (nf > 0)
	    {
	      stokes_nc_set_ff (nc, l, ff);
	      stokes_nc_set_tf (nc, l, tf);
	      stokes_nc_set_sf (nc, l, sf);
	    }
	}
      /* flush the data */
      nc_sync(nc->id);
    }

  printf ("Normaly Terminated !\n");

  stokes_nc_free (nc);


  free (x);
  free (angle);

  free (v);
  free (av);

  free (fm);
  free (tm);
  free (sm);
  free (um);
  free (om);
  free (em);
  free (ff);
  free (tf);
  free (sf);
  free (uf);
  free (of);
  free (ef);

  stokes_free (sys);

  return 0;
}
