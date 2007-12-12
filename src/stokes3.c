/* stokesian dynamics simulator for both periodic and non-periodic systems
 * Copyright (C) 1997-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: stokes3.c,v 1.20 2007/12/12 06:36:25 kichiki Exp $
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

#include <libstokes.h>
#include <netcdf.h> // nc_sync()

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#include "file.h" // check_file()

void
usage (const char *argv0)
{
  fprintf (stderr, "Stokesian dynamics simulator\n");
  fprintf (stderr, "$Id: stokes3.c,v 1.20 2007/12/12 06:36:25 kichiki Exp $\n\n");
  fprintf (stderr, "USAGE\n");
  fprintf (stderr, "%s [OPTIONS] init-file\n", argv0);
  fprintf (stderr, "\t-h or --help     : this message.\n");
  fprintf (stderr, "\t-o or --output   : overwrite the file name for output\n");
  fprintf (stderr, "\t-c or --continue : give nloop for continuation\n");
  fprintf (stderr, "\tinit-file : SCM file (default: stokes3.scm)\n\n");
  fprintf (stderr, "Parameters in the init-file:\n");
  fprintf (stderr, "* output parameters\n");
  fprintf (stderr, "\toutfile    : filename for NetCDF output\n");
  fprintf (stderr, "\tdt         : time interval for outputs\n");
  fprintf (stderr, "\tnloop      : number of loops for dt\n");
  fprintf (stderr,
	   "\tflag-Q     : #t output quaternion,\n"
	   "\t           : #f no quaternion in the output.\n");
  fprintf (stderr, "* core libstokes parameters\n");
  fprintf (stderr, "\tversion    : \"F\", \"FT\", or \"FTS\"\n");
  fprintf (stderr,
	   "\tflag-mat   : #t for matrix-scheme,\n"
	   "\t           : #f for atimes-scheme.\n");
  fprintf (stderr,
	   "\tflag-lub   : #t for with-lubrication,\n"
	   "\t           : #f for no-lubrication.\n");
  fprintf (stderr,
	   "\trmin       : parameter for mininum distance in"
	   " (ai+aj) * rmin.\n");
  fprintf (stderr,
	   "\tlub-min    : mininum cut-off distance for lubrication\n");
  fprintf (stderr,
	   "\tlub-max    : maximum cut-off distance for lubrication\n");
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
  fprintf (stderr, "* libiter parameters (used if \"flag-mat\" is #f)\n");
  fprintf (stderr,
	   "\tIT-solver  : solver for libiter\n"
	   "\t\t\"cg\"       conjugate gradient method\n"
	   "\t\t\"cgs\"      conjugate gradient squared (Weiss' Algorithm 11)\n"
	   "\t\t\"bicgstab\" Bi-CG stabilized (Weiss' Algorithm 12)\n"
	   "\t\t\"sta\"      Bi-CG stab method (another implementation)\n"
	   "\t\t\"sta2\"     Bi-CG dtab2 method\n"
	   "\t\t\"gpb\"      gpbi-cg method\n"
	   "\t\t\"otmk\"     orthomin method\n"
	   "\t\t\"gmres\"    generalized minimum residual method\n");
  fprintf (stderr, "\tIT-max     : max number of iteraction\n");
  fprintf (stderr, "\tIT-n       : restart number\n");
  fprintf (stderr, "\tIT-eps     : accuracy of the solution\n");
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
  fprintf (stderr, "\tnp     : number of ALL (mobile and fixed) particles\n");
  fprintf (stderr, "\tnm     : number of mobile particles\n");
  fprintf (stderr, "\tx      : particle configuration"
	   " (list or vector with length 3*np)\n");
  fprintf (stderr, "\ta      : radius of particles"
	   " (list or vector with length np)\n"
	   "\t\tby default (if not given), monodisperse system\n");
  fprintf (stderr, "\tslip   : slip length of particles"
	   " (list or vector with length np)\n"
	   "\t\tby default (if not given), no-slip particles\n");
  fprintf (stderr, "\tUi     : imposed translational velocity"
	   " (list or vector of length 3)\n");
  fprintf (stderr, "\tOi     : imposed angular velocity"
	   " (list or vector of length 3)\n");
  fprintf (stderr, "\tEi     : imposed strain"
	   " (list or vector of length 5)\n");
  fprintf (stderr, "\tF0     : applied force"
	   " (list or vector of length 3)\n");
  fprintf (stderr, "\tT0     : applied torque"
	   " (list or vector of length 3)\n");
  fprintf (stderr, "\tstokes : effective stokes number\n");
  fprintf (stderr, "\tncol   : frequency of collision check in dt"
	   " for stokes != 0\n");
  fprintf (stderr, "* bond parameters (for chains)\n");
  fprintf (stderr, "\tbonds      : bonds among particles,"
	   " list in the following form\n"
           "\t(define bonds '(\n"
           "\t  (; bond 1\n"
           "\t   0       ; 1) spring type\n"
           "\t   (       ; 2) spring parameters (list with 3 elements)\n"
           "\t    0      ;    fene = 0 means (p1, p2) = (A^{sp}, L_{s})\n"
           "\t    1.0    ;"
	   "    p1   = A^{sp}, scaled spring constant  (for fene == 0)\n"
           "\t    2.1)   ;"
	   "    p2   = L_{s} / a, scaled max extension (for fene == 0)\n"
           "\t   ((0 1)  ; 3) list of pairs\n"
           "\t    (1 2)\n"
           "\t    (2 3))\n"
           "\t    -1)    ; 4) number of exclusion for lubrication\n"
           "\t           ;    negative means all particles in the chain is excluded.\n"
           "\t  (; bond 2\n"
           "\t   2       ; 1) spring type\n"
           "\t   (       ; 2) spring parameters (list with 3 elements)\n"
           "\t    1      ;    fene = 1 means (p1, p2) = (N_{K,s}, b_{K})\n"
           "\t    19.8   ;"
	   "    p1 = N_{K,s}, the Kuhn steps for a spring (for fene = 1)\n"
           "\t    106.0) ;"
	   "    p2 = b_{K} [micro m], the Kuhn length     (for fene = 1)\n"
           "\t   ((4 5)  ; 3) list of pairs\n"
           "\t    (5 6)\n"
           "\t    (6 7))\n"
           "\t     1)    ; 4) number of exclusion for lubrication\n"
           "\t ))\n"
           "\twhere spring types are\n"
	   "\t  0 : Hookean spring (Asp * (r - Ls)\n"
	   "\t  1 : wormlike chain (WLC)\n"
	   "\t  2 : inverse Langevin chain (ILC)\n"
	   "\t  3 : Cohen's Pade approximation\n"
	   "\t  4 : Warner spring\n"
	   "\t  5 : Hookean spring (Asp * r / Ls)\n"
	   );
  fprintf (stderr,
	   "\tflag_relax : #f stokesian dynamics,\n"
	   "\t           : #t relaxation dynamics for bonds.\n");
  fprintf (stderr,
	   "\tgamma      : friction coefficient for the relaxation dynamics\n");
  fprintf (stderr, "* Excluded-Volume parameters\n");
  fprintf (stderr,
	   "\tev-v       : v [(micro m)^3] for each chain type\n"
	   "\t             (list or vector with length of bonds' length)\n");
  fprintf (stderr,
	   "\tev-lim     : maximum distance for EV interaction [micro m]\n");
  fprintf (stderr, "* Brownian dynamics' parameters\n");
  fprintf (stderr,
	   "\tpeclet     : peclet number (negative means no Brownian force)\n");
  fprintf (stderr,
	   "\tlength     : unit of the length scale [micro m]\n");
  fprintf (stderr,
	   "\tBD-seed    : seed for Brownian force random number generator\n");
  fprintf (stderr,
	   "\tn-cheb-minv: number of Chebyshev coefficients for M^{-1}\n");
  fprintf (stderr,
	   "\tn-cheb-lub : number of Chebyshev coefficients for L\n"
	   "\t             (if zero is given, use Cholesky decomposition)\n");
  fprintf (stderr, "\tBD-scheme : Brownian dynamics time integration scheme\n"
	   "\t\t\"mid-point\"        The mid-point algorithm.\n"
	   "\t\t\"BanchioBrady03\"   Banchio and Brady (2003).\n"
	   "\t\t\"BallMelrose97\"    Ball and Melrose (1997).\n"
	   "\t\t\"JendrejackEtal00\" Jendrejack et al (2000).\n");
  fprintf (stderr,
	   "\tBB-n        : step parameter for Banchio-Brady03 algorithm.\n");
  fprintf (stderr,
	   "\tdt-lim      : lower bound to shrink dt to prevent overlaps.\n"
	   "\t\tset \"dt\" if you don't want to adjust \"dt\" "
	   "but just reject it.\n");
}


/* main program */
int
main (int argc, char** argv)
{
  /**
   * option analysis
   */
  char init_file [256];
  strcpy (init_file, "stokes3.scm"); // default

  int i;
  int nloop_arg = 0;
  char *out_file = NULL;
  for (i = 1; i < argc; i++)
    {
      if (strcmp (argv [i], "-h") == 0 ||
	  strcmp (argv [i], "--help") == 0)
	{
	  usage (argv[0]);
	  exit (1);
	}
      else if (strcmp (argv [i], "-c") == 0 ||
	  strcmp (argv [i], "--continue") == 0)
	{
	  nloop_arg = atoi (argv [++i]);
	}
      else if (strcmp (argv [i], "-o") == 0 ||
	  strcmp (argv [i], "--output") == 0)
	{
	  i++;
	  int len = strlen (argv [i]);
	  out_file = (char *)malloc (sizeof (char) * (len + 1));
	  strcpy (out_file, argv [i]);
	}
      else
	{
	  strcpy (init_file, argv [i]);
	}
    }
  
  /**
   * parameter set
   */
  guile_load (init_file);

  // outfile
  if (out_file == NULL)
    {
      out_file = guile_get_string ("outfile");
    }
  else
    {
      fprintf (stderr, "outfile is overwritten by %s\n", out_file);
    }

  // flag-relax
  int flag_relax = 0;
  if (guile_get_bool ("flag-relax") != 0) // TRUE
    {
      flag_relax = 1;
    }

  // version
  int version = 0;
  char *str_version = guile_get_string ("version");
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
  int flag_mat = 0;
  if (guile_get_bool ("flag-mat") != 0) // TRUE
    {
      flag_mat = 1;
    }

  // flag-lub
  int flag_lub = 0;
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

  /**
   * imposed flow
   */
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

  /**
   * applied force
   */
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


  int np = guile_get_int ("np", 0);
  int nm = guile_get_int ("nm", 0);

  int np3 = np * 3;
  int nm3 = nm * 3;

  /**
   * position of ALL particles (BOTH mobile and fixed)
   * this is used to write into the output file, too.
   */
  double *x = (double *)malloc (sizeof (double) * np3);
  CHECK_MALLOC (x, "main");

  // x -- particle configuration
  if (guile_get_doubles ("x", np3, x) != 1) // FALSE
    {
      fprintf (stderr, "x is not defined\n");
      exit (1);
    }

  /**
   * Brownian dynamics parameters
   */
  // peclet number
  double peclet = guile_get_double ("peclet", -1.0);
  double length = guile_get_double ("length", 1.0);
  int BD_seed = guile_get_int ("BD-seed", 0);

  // chebyshev polynomials
  int n_minv = guile_get_int ("n-cheb-minv", 0);
  int n_lub  = guile_get_int ("n-cheb-lub", 0);
  // time-integration scheme for Brownian dynamics
  char *str_BD_scheme = guile_get_string ("BD-scheme");
  int BD_scheme;
  if (strcmp (str_BD_scheme, "mid-point") == 0)
    {
      BD_scheme = 0;
    }
  else if (strcmp (str_BD_scheme, "BanchioBrady03") == 0)
    {
      BD_scheme = 1;
    }
  else if (strcmp (str_BD_scheme, "BallMelrose97") == 0)
    {
      BD_scheme = 2;
    }
  else if (strcmp (str_BD_scheme, "JendrejackEtal00") == 0)
    {
      BD_scheme = 3;
    }
  else
    {
      fprintf (stderr, "invalid BD-scheme %s", str_BD_scheme);
      exit (1);
    }
  free (str_BD_scheme);
  // step parameter for BB03 algorithm
  double BB_n = guile_get_double ("BB-n", 100.0);
  // lower bound to shrink dt to prevent overlaps
  double dt_lim = guile_get_double ("dt-lim", 1.0e-12);


  /**
   * bonds
   */
  struct bonds *bonds = guile_get_bonds ("bonds");
  if (bonds == NULL) // FALSE
    {
      fprintf (stderr, "main: fail to parse bonds\n");
      exit (1);
    }
  bonds_set_FENE (bonds, length, peclet);
  /* check
  for (i = 0; i < bonds->n; i ++)
    {
      fprintf (stdout, "# bonds [%d] %d %d %e %e\n",
	       i, bonds->type[i], bonds->fene[i],
	       bonds->p1[i],
	       bonds->p2[i]);
    }
  */

  // relaxation dynamics only with bond interaction
  double gamma = guile_get_double ("gamma", 1.0);

  /**
   * excluded volume
   */
  struct EV *ev = NULL;
  double ev_lim = guile_get_double ("ev-lim", 1.0);
  ev_lim /= length; // scale by length
  double ev_r2 = ev_lim * ev_lim;
  if (bonds->n > 0)
    {
      double *ev_v = (double *)malloc (sizeof (double) * bonds->n);
      CHECK_MALLOC (ev_v, "main");
      if (guile_get_doubles ("ev-v", bonds->n, ev_v) != 1) // FALSE
	{
	  fprintf (stderr, "excluded-volume is not defined\n");
	  //exit (1);
	}
      else
	{
	  ev = EV_init (bonds, length, peclet,
			ev_r2, ev_v, np);
	}
      free (ev_v);
    }

  // initialize struct stokes *sys
  struct stokes *sys = stokes_init ();
  sys->version = version;
  stokes_set_np (sys, np, nm);

  sys->rmin = guile_get_double ("rmin", 0.0);
  double lubmin = guile_get_double ("lub-min", 2.0000000001);
  sys->lubmin2 = lubmin * lubmin;
  sys->lubmax = guile_get_double ("lub-max", 4.0);

  /* set exclusion list for lub by bonds */
  list_ex_set_by_bonds (sys->ex_lub, bonds);

  // iterative solver
  char *str_it_solver = NULL;
  str_it_solver = guile_get_string ("IT-solver");
  int it_max = guile_get_int ("IT-max", 2000);
  int it_n   = guile_get_int ("IT-n", 20);
  double it_eps = guile_get_double ("IT-eps", 1.0e-6);
  int it_debug = guile_get_int ("IT-debug", 0);
  stokes_set_iter (sys, str_it_solver, it_max, it_n, it_eps, it_debug, stderr);
  free (str_it_solver);

  stokes_set_Ui (sys, Ui[0], Ui[1], Ui[2]);
  stokes_set_Oi (sys, Oi[0], Oi[1], Oi[2]);
  stokes_set_Ei (sys, Ei[0], Ei[1], Ei[2], Ei[3], Ei[4]);


  // radius of ALL particles (BOTH mobile and fixed)
  //int flag_poly = 0; // for stokes_nc_init()
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
      //flag_poly = 1;
      stokes_set_radius (sys, a);
    }
  free (a);


  // slip length of ALL particles (BOTH mobile and fixed)
  double *slip = (double *)malloc (sizeof (double) * np);
  CHECK_MALLOC (slip, "main");
  if (guile_get_doubles ("slip", np, slip) != 1) // FALSE
    {
      // "slip" is not given, so that system is no-slip
      // sys->slip is NULL by default, so do nothing here
    }
  else
    {
      // set slip parameters
      stokes_set_slip (sys, slip);
    }
  free (slip);


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

  // quaternion
  int flag_Q = 0;
  if (guile_get_bool ("flag-Q") != 0) // TRUE
    {
      if (version > 1) flag_Q = 1;
      // set flag_Q = 1 for F version
    }

  /* initialize the dependent variable for ODE y[], where
   * for flag_Q == 0,
   *   for st==0, y[nm3  ] =  x[nm3], positions of mobile particles,
   *   for st!=0, y[nm3*2] = (x[nm3],U[nm3]), for st != 0
   * for flag_Q != 0,
   *   for st==0, y[nm3  +nm4], quaternion[nm4] is added
   *   for st!=0, y[nm3*2+nm4], quaternion[nm4] is added
   */
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
  int nq = 0; // defined for flag_Q!=0 by nm3 (for st==0), nm3*2 (for st!=0)
  if (flag_Q != 0)
    {
      // follow the angles by quaternion
      nq = n;
      n += nm * 4;
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
  if (flag_Q != 0)
    {
      for (i = 0; i < nm; i ++)
	{
	  y[nq + i*4+0] = 0.0;
	  y[nq + i*4+1] = 0.0;
	  y[nq + i*4+2] = 0.0;
	  y[nq + i*4+3] = 1.0;
	}
    }


  /**
   * set constant parameters for ode_params
   */
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

      if (sys->a == NULL)
	{
	  F [i3  ] = F0 [0];
	  F [i3+1] = F0 [1];
	  F [i3+2] = F0 [2];

	  T [i3  ] = T0 [0];
	  T [i3+1] = T0 [1];
	  T [i3+2] = T0 [2];
	}
      else
	{
	  double rad = sys->a[i];
	  double a3 = rad*rad*rad;
	  F [i3  ] = F0 [0] *a3;
	  F [i3+1] = F0 [1] *a3;
	  F [i3+2] = F0 [2] *a3;

	  double a4 = a3 * rad;
	  T [i3  ] = T0 [0] *a4;
	  T [i3+1] = T0 [1] *a4;
	  T [i3+2] = T0 [2] *a4;
	}

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

  /* set position for the fixed particles in sys */
  if (np > nm)
    {
      stokes_set_pos_fixed (sys, x + nm*3);
    }

  /**
   * set ode_params
   */
  struct ode_params *ode_params = NULL;
  int (*f_dydt)(double, const double *, double *, void *) = NULL;
  const gsl_odeiv_step_type *GSL_ODE_TYPE = NULL;
  gsl_odeiv_system GSL_ODE_SYSTEM;
  gsl_odeiv_step    *GSL_ODE_STEP    = NULL;
  gsl_odeiv_control *GSL_ODE_CONTROL = NULL;
  gsl_odeiv_evolve  *GSL_ODE_EVOLVE  = NULL;

  struct BD_params *BD_params = NULL;
  struct BD_imp *BDimp = NULL;
  if (peclet < 0.0) // no Brownian force
    {
      ode_params = ode_params_init (sys,
				    F, T, E,
				    uf, of, ef,
				    flag_lub, flag_mat,
				    st,
				    bonds,
				    gamma);
      CHECK_MALLOC (ode_params, "main");


      // note that dydt_hydro() and dydt_hydro_st() can handle bonds.
      if (flag_Q == 0)
	{
	  if (st > 0.0) f_dydt = dydt_hydro_st;
	  else          f_dydt = dydt_hydro;
	}
      else
	{
	  if (st > 0.0) f_dydt = dydt_Q_hydro_st;
	  else          f_dydt = dydt_Q_hydro;
	}

      // asign params for gsl_odeiv_system GSL_ODE_SYSTEM
      GSL_ODE_SYSTEM.function  = f_dydt;
      GSL_ODE_SYSTEM.jacobian  = NULL;
      GSL_ODE_SYSTEM.dimension = n;
      GSL_ODE_SYSTEM.params    = ode_params;

      // GSL ODE integrator
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

      GSL_ODE_STEP = gsl_odeiv_step_alloc (GSL_ODE_TYPE, n);
      GSL_ODE_CONTROL = gsl_odeiv_control_y_new (ode_eps, 0.0);
      GSL_ODE_EVOLVE = gsl_odeiv_evolve_alloc (n);
    }
  else // Brownian dynamics
    {
      BD_params = BD_params_init (sys,
				  BD_seed,
				  F, T, E,
				  uf, of, ef,
				  flag_lub, flag_mat,
				  st,
				  bonds,
				  gamma,
				  ev,
				  flag_Q,
				  peclet,
				  ode_eps,
				  n_minv, // n of chebyshev for minv
				  n_lub,  // n of chebyshev for lub
				  BD_scheme,
				  BB_n,
				  dt_lim);
      CHECK_MALLOC (BD_params, "main");

      // implicit scheme
      if (BD_scheme == 3)
	{
	  BDimp = BD_imp_init (BD_params,
			       1000,  // itmax
			       ode_eps
			       );
	  CHECK_MALLOC (BDimp, "main");
	}
    }


  /**
   * initialize NetCDF and set the constant parameters
   */
  double t;
  int l0;
  struct stokes_nc *nc = NULL;
  if (nloop_arg > 0)
    {
      // continuation
      if (check_file (out_file) == 0)
	{
	  fprintf (stderr, "result file %s does not exist.\n", out_file);
	  exit (1);
	}
      nc = stokes_nc_reopen (out_file);
      CHECK_MALLOC (nc, "main");

      if (stokes_nc_check_params (nc, sys,
				  flag_Q,
				  Ui, Oi, Ei, F, T, E,
				  uf, of, ef, x + nm3, // xf
				  lat,
				  1.0e-16)
	  != 0)
	{
	  fprintf (stderr, "result file %s does not match to"
		   " the init script %s\n",
		   out_file, init_file);
	  exit (1);
	}

      // set the loop parameters
      l0 = nc->ntime;
      /* l0 is the starting index because index starts from 0
       * so that we need to access "l0-1" to get the last step infos
       */
      nloop = l0 + nloop_arg;
      t  = stokes_nc_get_time_step (nc, l0-1);

      // set the configuration at the current time
      stokes_nc_get_data (nc, "x", l0-1, y);
      if (flag_Q != 0)
	{
	  stokes_nc_get_data (nc, "q", l0-1, y + nq);
	}
    }
  else
    {
      // create new stokes_nc
      nc = stokes_nc_set_by_params (out_file,
				    sys,
				    flag_Q,
				    Ui, Oi, Ei, F, T, E,
				    uf, of, ef, x + nm3, // xf
				    lat);
      CHECK_MALLOC (nc, "main");

      // set the loop parameters
      l0 = 0;
      // nloop is given by the script
      t = 0.0;
    }
  free (out_file);


  /**
   * mail loop
   */
  double h = dt / (double)ncol; // initial time step for ODE integrator
  double ddt = dt / (double)ncol; // time step for collision check for st != 0
  int l;
  double ptime0 = ptime_ms_d();
  for (l = l0; l < nloop; l++)
    {
      fprintf (stdout, "%d steps\n", l + 1);

      // integrate from t to t_out
      double t_out = t + dt;
      while (t < t_out)
	{
	  if (st == 0.0)
	    {
	      if (peclet < 0.0) // no Brownian force
		{
		  // (st==0) no collision checks
		  gsl_odeiv_evolve_apply (GSL_ODE_EVOLVE,
					  GSL_ODE_CONTROL,
					  GSL_ODE_STEP,
					  &GSL_ODE_SYSTEM,
					  &t, t_out,
					  &h, y);
		}
	      else // Brownian dynamics
		{
		  if (BD_params->scheme == 3)
		    {
		      BD_imp_ode_evolve (BDimp,
					 &t, t_out,
					 &h, y);
		    }
		  else
		    {
		      BD_ode_evolve (BD_params,
				     &t, t_out,
				     &h, y);
		    }
		}

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
      if (flag_Q != 0)
	{
	  stokes_nc_set_q (nc, l, y + nq);
	}

      // flush the data
      nc_sync(nc->id);
    }
  double ptime1 = ptime_ms_d();

  printf ("Normaly Terminated !\n");
  printf ("CPU time : %.3f [m sec]\n", ptime1 - ptime0);


  ode_params_free (ode_params);
  BD_params_free (BD_params);
  BD_imp_free (BDimp);
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
