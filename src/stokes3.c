/* stokesian dynamics simulator for both periodic and non-periodic systems
 * Copyright (C) 1997-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: stokes3.c,v 1.7 2007/05/11 02:10:14 kichiki Exp $
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
#include <stdio.h> /* printf() fprintf() */
#include <stdlib.h> /* exit() */
#include <string.h> /* strcmp() */
#include "memory-check.h"

#include <libiter.h> /* iter_init() */
#include <libstokes.h> /* struct stokes */
#include <netcdf.h> // nc_sync()
#include <libguile.h> // scm_init_guile()

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>


void
usage (const char *argv0)
{
  fprintf (stderr, "Stokesian dynamics simulator\n");
  fprintf (stderr, "$Id: stokes3.c,v 1.7 2007/05/11 02:10:14 kichiki Exp $\n\n");
  fprintf (stderr, "USAGE\n");
  fprintf (stderr, "%s init-file\n", argv0);
  fprintf (stderr, "\twhere init-file is a SCM file"
	   " (default: stokes3.scm)\n\n");
  fprintf (stderr, "Parameters in the init-file:\n");
  fprintf (stderr, "* output parameters\n");
  fprintf (stderr, "\toutfile    : filename for NetCDF output\n");
  fprintf (stderr, "\tdt         : time interval for outputs\n");
  fprintf (stderr, "\tnloop      : number of loops for dt\n");
  fprintf (stderr, "* core libstokes parameters\n");
  fprintf (stderr, "\tversion    : \"F\", \"FT\", or \"FTS\"\n");
  fprintf (stderr,
	   "\tflag-mat   : #t for matrix-scheme,\n"
	   "\t           : #f for atimes-scheme.\n");
  fprintf (stderr,
	   "\tflag-lub   : #t for with-lubrication,\n"
	   "\t           : #f for no-lubrication.\n");
  fprintf (stderr,
	   "\tperiodic   : #f for non-periodic systems,\n"
	   "\t           : #t for periodic systems.\n"
	   "\t             set the next three parameters"
	   " for the periodic case.\n");
  fprintf (stderr, "\tewald-tr   : time ratio Tr/Tk for ewald summation"
	   " (see xi3 in details.)\n");
  fprintf (stderr, "\tewald-eps  : tolerance value"
	   " for ewald-summation cut-off\n");
  fprintf (stderr, "\tlattice    : size of the periodic box"
	   " (list or vector of length 3)\n");
  fprintf (stderr, "* ODE parameters\n");
  fprintf (stderr, "\tode-solver : GSL ODE solver\n"
	   "\t\t\"rk2\"    Embedded Runge-Kutta (2, 3) method.\n"
	   "\t\t\"rk4\"    4th order (classical) Runge-Kutta.\n"
	   "\t\t\"rkf45\"  Embedded Runge-Kutta-Fehlberg (4, 5) method.\n"
	   "\t\t\"rkck\"   Embedded Runge-Kutta Cash-Karp (4, 5) method.\n"
	   "\t\t\"rk8pd\"  Embedded Runge-Kutta Prince-Dormand (8,9) method.\n"
	   "\t\t\"rk2imp\" Implicit 2nd order Runge-Kutta at Gaussian points.\n"
	   "\t\t\"rk4imp\" Implicit 4th order Runge-Kutta at Gaussian points.\n"
	   "\t\t\"gear1\"  M=1 implicit Gear method.\n"
	   "\t\t\"gear2\"  M=2 implicit Gear method.\n");
  fprintf (stderr, "\tode-eps    : GSL ODE control parameter eps\n");
  fprintf (stderr, "* system parameters\n");
  fprintf (stderr, "\tnp         : number of ALL (mobile and fixed) particles\n");
  fprintf (stderr, "\tnm         : number of mobile particles\n");
  fprintf (stderr, "\tx          : particle configuration"
	   " (list or vector with length 3*np)\n");
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
  fprintf (stderr, "\tstokes     : effective stokes number\n");
  fprintf (stderr, "\tncol       : frequency of collision check in dt"
	   " for stokes != 0\n");
  fprintf (stderr, "* bond parameters (for chains)\n");
  fprintf (stderr, "\tbonds      : bonds among particles,"
	   " list in the following form\n"
	   "\t\t(define bonds\n"
	   "\t\t  '(\n"
	   "\t\t    (; bond 1\n"
	   "\t\t     1.0 ; spring const\n"
	   "\t\t     2.1 ; natural distance\n"
	   "\t\t     ((0 1) ; list of pairs\n"
	   "\t\t      (1 2)\n"
	   "\t\t      (2 3)))\n"
	   "\t\t    (; bond 2\n"
	   "\t\t     1.0 ; spring const\n"
	   "\t\t     2.5 ; natural distance\n"
	   "\t\t     ((4 5) ; list of pairs\n"
	   "\t\t      (5 6)\n"
	   "\t\t      (6 7)))\n"
	   "\t\t    )\n"
	   "\t\t  )\n"
	   );
  fprintf (stderr,
	   "\tflag_relax : #f stokesian dynamics,\n"
	   "\t           : #t relaxation dynamics for bonds.\n");
  fprintf (stderr, "\tgamma      : friction coefficient for the relaxation dynamics\n");
}




/* main program */
int
main (int argc, char** argv)
{
  /* option analysis */
  char init_file [256];
  strcpy (init_file, "stokes3.scm"); // default

  int i;
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

  // outfile
  char * out_file;
  out_file = guile_get_string ("outfile");


  // flag-relax
  int flag_relax;
  flag_relax = 0;
  if (guile_get_bool ("flag-relax") != 0) // TRUE
    {
      flag_relax = 1;
    }

  // version
  int version;
  version = 0;
  char * str_version = NULL;
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
  int flag_mat;
  flag_mat = 0;
  if (guile_get_bool ("flag-mat") != 0) // TRUE
    {
      flag_mat = 1;
    }

  // flag-lub
  int flag_lub;
  flag_lub = 0;
  if (guile_get_bool ("flag-lub") != 0) // TRUE
    {
      flag_lub = 1;
    }

  // ode-solver
  char * str_ode_solver = NULL;
  str_ode_solver = guile_get_string ("ode-solver");
  // ode-eps
  double ode_eps = guile_get_double ("ode-eps", 1.0e-6);

  // other parameters
  int nloop = guile_get_int    ("nloop",  0);
  double dt = guile_get_double ("dt",     0.0);
  double st = guile_get_double ("stokes", 0.0);
  int ncol  = guile_get_int    ("ncol",   10);

  /* imposed flow */
  // Ui
  double Ui[3];
  if (guile_get_doubles ("Ui", 3, Ui) != 1) // FALSE
    {
      fprintf (stderr, "Ui is not defined\n");
      exit (1);
    }
  // Oi
  double Oi[3];
  if (guile_get_doubles ("Oi", 3, Oi) != 1) // FALSE
    {
      fprintf (stderr, "Oi is not defined\n");
      exit (1);
    }
  // Ei
  double Ei[5];
  if (guile_get_doubles ("Ei", 5, Ei) != 1) // FALSE
    {
      fprintf (stderr, "Ei is not defined\n");
      exit (1);
    }

  /* applied force */
  // F0
  double F0[3];
  if (guile_get_doubles ("F0", 3, F0) != 1) // FALSE
    {
      fprintf (stderr, "F0 is not defined\n");
      exit (1);
    }
  // T0
  double T0[3];
  if (guile_get_doubles ("T0", 3, T0) != 1) // FALSE
    {
      fprintf (stderr, "T0 is not defined\n");
      exit (1);
    }


  int np = guile_get_int ("np",         0);
  int nm = guile_get_int ("nm",         0);

  int np3 = np * 3;
  int nm3 = nm * 3;

  // position of ALL particles (BOTH mobile and fixed)
  // this is used to write into the output file, too.
  double *x = (double *)malloc (sizeof (double) * np3);
  CHECK_MALLOC (x, "main");

  // x -- particle configuration
  if (guile_get_doubles ("x", np3, x) != 1) // FALSE
    {
      fprintf (stderr, "x is not defined\n");
      exit (1);
    }

  // bonds
  struct bonds *bonds = guile_get_bonds ("bonds");
  if (bonds == NULL) // FALSE
    {
      fprintf (stderr, "main: fail to parse bonds\n");
      exit (1);
    }

  // relaxation dynamics only with bond interaction
  double gamma = guile_get_double ("gamma", 1.0);


  // initialize struct stokes *sys
  struct stokes *sys = stokes_init ();
  sys->version = version;
  stokes_set_np (sys, np, nm);

  double lubmin = guile_get_double ("lub-min", 2.0000000001);
  sys->lubmin2 = lubmin * lubmin;
  sys->lubmax = guile_get_double ("lub-max", 4.0);
  stokes_set_iter (sys, "gmres", 2000, 20, 1.0e-6, 0, stderr);

  stokes_set_Ui (sys, Ui[0], Ui[1], Ui[2]);
  stokes_set_Oi (sys, Oi[0], Oi[1], Oi[2]);
  stokes_set_Ei (sys, Ei[0], Ei[1], Ei[2], Ei[3], Ei[4]);


  // radius of ALL particles (BOTH mobile and fixed)
  int flag_poly = 0; // for stokes_nc_init()
  double *a = (double *)malloc (sizeof (double) * np);
  CHECK_MALLOC (a, "main");
  if (guile_get_doubles ("a", np, a) != 1) // FALSE
    {
      // "a" is not given, so that system is monodisperse
      // sys->a is NULL by default, so do nothing here
    }
  else
    {
      // the system is polydisperse
      flag_poly = 1;
      stokes_set_radius (sys, a);
    }
  free (a);


  // periodic
  double lat[3];
  if (guile_get_bool ("periodic") != 0) // TRUE
    {
      // periodic
      sys->periodic = 1;

      if (guile_get_doubles ("lattice", 3, lat) != 1) // FALSE
	{
	  fprintf (stderr, "lattice is not defined\n");
	  exit (1);
	}
      stokes_set_l (sys, lat[0], lat[1], lat[2]);

      double ewald_tr  = guile_get_double ("ewald-tr",  0.0);
      double ewald_eps = guile_get_double ("ewald-eps", 0.0);
      double xi = xi_by_tratio (sys, ewald_tr);
      stokes_set_xi (sys, xi, ewald_eps);
    }
  else
    {
      // non-periodic
      sys->periodic = 0;
    }



  // initialize the dependent variable for ODE y[], where
  // for st==0, y[nm3]   =  x[nm3], positions of mobile particles,
  // for st!=0, y[nm3*2] = (x[nm3],U[nm3]), for st != 0
  int n; // dimension for ODE integrator
  if (st <= 0.0)
    {
      // positions of mobile particles
      n = nm3;
    }
  else
    {
      // position and velocity of mobile particles
      n = nm3 * 2;
    }
  double *y = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (y, "main");

  // set the initial configuration
  for (i = 0; i < nm3; i ++)
    {
      y[i] = x[i];
    }
  if (st > 0.0)
    {
      // initial velocity
      for (i = 0; i < nm3; i ++)
	{
	  y[nm3 + i] = 0.0;
	}
    }


  // set constant parameters for ode_params
  int nm5 = nm * 5;
  double *F = (double *)malloc (sizeof (double) * nm3);
  double *T = (double *)malloc (sizeof (double) * nm3);
  double *E = (double *)malloc (sizeof (double) * nm5);
  CHECK_MALLOC (F, "main");
  CHECK_MALLOC (T, "main");
  CHECK_MALLOC (E, "main");
  int i3, i5;
  for (i = 0; i < nm; i ++)
    {
      i3 = i * 3;
      i5 = i * 5;

      F [i3  ] = F0 [0];
      F [i3+1] = F0 [1];
      F [i3+2] = F0 [2];

      T [i3  ] = T0 [0];
      T [i3+1] = T0 [1];
      T [i3+2] = T0 [2];

      E [i5  ] = 0.0;
      E [i5+1] = 0.0;
      E [i5+2] = 0.0;
      E [i5+3] = 0.0;
      E [i5+4] = 0.0;
    }
  // prepare variables for fixed particles
  //int nf = sys->np - nm;
  int nf = np - nm;
  int nf3 = nf * 3;
  int nf5 = nf * 5;
  double *uf = NULL;
  double *of = NULL;
  double *ef = NULL;
  if (nf > 0)
    {
      uf = (double *)malloc (sizeof (double) * nf3);
      of = (double *)malloc (sizeof (double) * nf3);
      ef = (double *)malloc (sizeof (double) * nf5);
      CHECK_MALLOC (uf, "main");
      CHECK_MALLOC (of, "main");
      CHECK_MALLOC (ef, "main");
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

  // set ode_params
  double *pos_fixed = NULL;
  if (np > nm)
    {
      pos_fixed = x + nm*3;
    }

  struct ode_params *ode_params
    = ode_params_init (sys, pos_fixed,
		       F, T, E,
		       uf, of, ef,
		       flag_lub, flag_mat,
		       st,
		       bonds,
		       gamma);
  CHECK_MALLOC (ode_params, "main");


  int (*f_dydt)(double, const double *, double *, void *) = NULL;
  // note that dydt_hydro() and dydt_hydro_st() can handle bonds.
  if (st > 0.0)
    {
      f_dydt = dydt_hydro_st;
    }
  else
    {
      // st == 0
      f_dydt = dydt_hydro;
    }

  gsl_odeiv_system GSL_ODE_SYSTEM
    = {f_dydt,                // function dy/dt
       NULL,                // jacobian (optional for simpler solver)
       n,                   // size_t dimension
       (void *)ode_params}; // void * params

  // GSL ODE integrator
  const gsl_odeiv_step_type *GSL_ODE_TYPE = NULL;
  if (strcmp (str_ode_solver, "rk2") == 0)
    {
      GSL_ODE_TYPE = gsl_odeiv_step_rk2;
    }
  else if (strcmp (str_ode_solver, "rk4") == 0)
    {
      GSL_ODE_TYPE = gsl_odeiv_step_rk4;
    }
  else if (strcmp (str_ode_solver, "rkf45") == 0)
    {
      GSL_ODE_TYPE = gsl_odeiv_step_rkf45;
    }
  else if (strcmp (str_ode_solver, "rkck") == 0)
    {
      GSL_ODE_TYPE = gsl_odeiv_step_rkck;
    }
  else if (strcmp (str_ode_solver, "rk8pd") == 0)
    {
      GSL_ODE_TYPE = gsl_odeiv_step_rk8pd;
    }
  else if (strcmp (str_ode_solver, "rk2imp") == 0)
    {
      GSL_ODE_TYPE = gsl_odeiv_step_rk2imp;
    }
  else if (strcmp (str_ode_solver, "rk4imp") == 0)
    {
      GSL_ODE_TYPE = gsl_odeiv_step_rk4imp;
    }
  else if (strcmp (str_ode_solver, "gear1") == 0)
    {
      GSL_ODE_TYPE = gsl_odeiv_step_gear1;
    }
  else if (strcmp (str_ode_solver, "gear2") == 0)
    {
      GSL_ODE_TYPE = gsl_odeiv_step_gear2;
    }
  else
    {
      fprintf (stderr, "invalid ode-solver %s", str_ode_solver);
      exit (1);
    }
  free (str_ode_solver);

  gsl_odeiv_step *GSL_ODE_STEP;
  GSL_ODE_STEP = gsl_odeiv_step_alloc (GSL_ODE_TYPE, n);

  gsl_odeiv_control *GSL_ODE_CONTROL;
  GSL_ODE_CONTROL = gsl_odeiv_control_y_new (ode_eps, 0.0);

  gsl_odeiv_evolve *GSL_ODE_EVOLVE;
  GSL_ODE_EVOLVE = gsl_odeiv_evolve_alloc (n);


  // initialize NetCDF and set the constant parameters
  struct stokes_nc *nc
    = stokes_nc_init (out_file, np, nf,
		      version, flag_poly, 0);
  if (nf == 0)
    {
      // mobility problem (no fixed particles)
      if (flag_poly != 0)
	{
	  stokes_nc_set_a (nc, sys->a);
	}

      if (version == 0)
	{
	  // F version
	  stokes_nc_set_f0 (nc, F);
	}
      else if (version == 1)
	{
	  // FT version
	  stokes_nc_set_f0 (nc, F);
	  stokes_nc_set_t0 (nc, T);
	}
      else
	{
	  // FTS version
	  stokes_nc_set_f0 (nc, F);
	  stokes_nc_set_t0 (nc, T);
	  stokes_nc_set_e0 (nc, E);
	}
    }
  else
    {
      // mix problem (with fixed particles)
      if (flag_poly != 0)
	{
	  stokes_nc_set_a (nc, sys->a);
	  stokes_nc_set_af (nc, sys->a + sys->nm);
	}

      if (version == 0)
	{
	  // F version
	  stokes_nc_set_f0 (nc, F);
	  stokes_nc_set_uf0 (nc, uf);
	}
      else if (version == 1)
	{
	  // FT version
	  stokes_nc_set_f0 (nc, F);
	  stokes_nc_set_t0 (nc, T);
	  stokes_nc_set_uf0 (nc, uf);
	  stokes_nc_set_of0 (nc, of);
	}
      else
	{
	  // FTS version
	  stokes_nc_set_f0 (nc, F);
	  stokes_nc_set_t0 (nc, T);
	  stokes_nc_set_e0 (nc, E);
	  stokes_nc_set_uf0 (nc, uf);
	  stokes_nc_set_of0 (nc, of);
	  stokes_nc_set_ef0 (nc, ef);
	}
      stokes_nc_set_xf0 (nc, x + nm*3);
    }
  free (out_file);

  // non-periodic system
  if (sys->periodic == 1)
    {
      stokes_nc_set_l (nc, lat);
    }



  /* mail loop */
  double t = 0.0;
  double h = dt / (double)ncol; // initial time step for ODE integrator
  double ddt = dt / (double)ncol; // time step for collision check for st != 0
  int l;
  double ptime0 = ptime_ms_d();
  for (l = 0; l < nloop; l++)
    {
      fprintf (stdout, "%d steps\n", l);

      // integrate from t to t_out
      double t_out = t + dt;
      while (t < t_out)
	{
	  if (st == 0.0)
	    {
	      // (st==0) no collision checks
	      gsl_odeiv_evolve_apply (GSL_ODE_EVOLVE,
				      GSL_ODE_CONTROL,
				      GSL_ODE_STEP,
				      &GSL_ODE_SYSTEM,
				      &t, t_out,
				      &h, y);

	      // check periodicity
	      if (sys->periodic == 1)
		{
		  check_periodic (sys, y);
		}
	    }
	  else
	    {
	      // (st!=0) ncol times collision checks in dt step
	      for (i = 0; i < ncol; i ++)
		{
		  double tt_out = t + ddt;
		  gsl_odeiv_evolve_apply (GSL_ODE_EVOLVE,
					  GSL_ODE_CONTROL,
					  GSL_ODE_STEP,
					  &GSL_ODE_SYSTEM,
					  &t, tt_out,
					  &h, y);

		  // set pos for mobile particles
		  // note that x[np3] has both mobile and fixed particles
		  // and collide_* and check_periodic need that
		  int j;
		  for (j = 0; j < nm3; j ++)
		    {
		      x [j] = y [j];
		    }

		  // check the collision
		  collide_particles (sys,
				     x, y + nm3,
				     1.0);
		  collide_wall_z (sys,
				  x, y + nm3,
				  1.0, 0.0, 0.0);

		  // check periodicity
		  if (sys->periodic == 1)
		    {
		      check_periodic (sys, y);
		    }
		}
	    }
	}

      // output the results
      stokes_nc_set_time (nc, l, t);
      stokes_nc_set_x (nc, l, y);

      // flush the data
      nc_sync(nc->id);
    }
  double ptime1 = ptime_ms_d();

  printf ("Normaly Terminated !\n");
  printf ("CPU time : %.3f [m sec]\n", ptime1 - ptime0);


  free (ode_params);
  stokes_nc_free (nc);
  free (x);
  free (y);

  stokes_free (sys);
  bonds_free (bonds);

  free (F);
  free (T);
  free (E);
  if (nf > 0)
    {
      free (uf);
      free (of);
      free (ef);
    }

  return 0;
}
