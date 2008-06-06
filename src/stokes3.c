/* stokesian dynamics simulator for both periodic and non-periodic systems
 * Copyright (C) 1997-2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: stokes3.c,v 1.33 2008/06/06 04:23:47 kichiki Exp $
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
  fprintf (stdout, "Stokesian dynamics simulator\n");
  fprintf (stdout, "$Id: stokes3.c,v 1.33 2008/06/06 04:23:47 kichiki Exp $\n\n");
  fprintf (stdout, "USAGE\n");
  fprintf (stdout, "%s [OPTIONS] init-file\n", argv0);
  fprintf (stdout, "\t-h or --help     : this message.\n");
  fprintf (stdout, "\t-o or --output   : overwrite the file name for output\n");
  fprintf (stdout, "\t-c or --continue : give nloop for continuation\n");
  fprintf (stdout, "\tinit-file : SCM file (default: stokes3.scm)\n\n");
  fprintf (stdout, "Parameters in the init-file:\n");
  fprintf (stdout, "* output parameters\n");
  fprintf (stdout, "\toutfile    : filename for NetCDF output\n");
  //fprintf (stdout, "\tdt         : time interval for outputs\n");
  fprintf (stdout, "\tdt         : time interval for the integrator\n");
  fprintf (stdout, "\tnout       : frequency of data outputs\n");
  fprintf (stdout,
	   "\t   NOTE for GSL integrator (non-Brownian), this is just an initial one,\n"
	   "\t\tfor non-zero Stokes, \"dt\" is the upper limit,\n"
	   "\t\tfor zero Stokes, \"dt * nout\" is the upper limit,\n"
	   "\t\tfor Brownian, \"dt\" is fixed by this value.\n");
  fprintf (stdout, "\tnloop      : number of output steps (time duration is \"dt*nout*nloop\")\n");
  fprintf (stdout,
	   "\tflag-Q     : #t output quaternion,\n"
	   "\t           : #f no quaternion in the output.\n");
  fprintf (stdout, "* core libstokes parameters\n");
  fprintf (stdout, "\tversion    : \"F\", \"FT\", or \"FTS\"\n");
  fprintf (stdout,
	   "\tflag-noHI  : #t no hydrodynamic interaction (only self term),\n"
	   "\t           : #f with hydrodynamic interactions.\n");
  fprintf (stdout,
	   "\tflag-mat   : #t for matrix-scheme,\n"
	   "\t           : #f for atimes-scheme.\n");
  fprintf (stdout,
	   "\tflag-lub   : #t for with-lubrication,\n"
	   "\t           : #f for no-lubrication.\n");
  fprintf (stdout,
	   "\trmin       : parameter for mininum distance in"
	   " (ai+aj) * rmin.\n");
  fprintf (stdout,
	   "\tlub-min    : mininum cut-off distance for lubrication\n");
  fprintf (stdout,
	   "\tlub-max    : maximum cut-off distance for lubrication\n");
  fprintf (stdout,
	   "\tperiodic   : #f for non-periodic systems,\n"
	   "\t           : #t for periodic systems.\n"
	   "\t             set the next three parameters"
	   " for the periodic case.\n");
  fprintf (stdout, "\tewald-tr   : time ratio Tr/Tk for ewald summation"
	   " (see xi3 in details.)\n");
  fprintf (stdout, "\tewald-eps  : tolerance value"
	   " for ewald-summation cut-off\n");
  fprintf (stdout, "\tlattice    : size of the periodic box"
	   " (list or vector of length 3)\n");
  fprintf (stdout, "* libiter parameters (used if \"flag-mat\" is #f)\n");
  fprintf (stdout,
	   "\tIT-solver  : solver for libiter\n"
	   "\t\t\"cg\"       conjugate gradient method\n"
	   "\t\t\"cgs\"      conjugate gradient squared (Weiss' Algorithm 11)\n"
	   "\t\t\"bicgstab\" Bi-CG stabilized (Weiss' Algorithm 12)\n"
	   "\t\t\"sta\"      Bi-CG stab method (another implementation)\n"
	   "\t\t\"sta2\"     Bi-CG dtab2 method\n"
	   "\t\t\"gpb\"      gpbi-cg method\n"
	   "\t\t\"otmk\"     orthomin method\n"
	   "\t\t\"gmres\"    generalized minimum residual method\n");
  fprintf (stdout, "\tIT-max     : max number of iteraction\n");
  fprintf (stdout, "\tIT-n       : restart number\n");
  fprintf (stdout, "\tIT-eps     : accuracy of the solution\n");
  fprintf (stdout, "* ODE parameters\n");
  fprintf (stdout, "\tode-solver : GSL ODE solver\n"
	   "\t\t\"rk2\"    Embedded Runge-Kutta (2, 3) method.\n"
	   "\t\t\"rk4\"    4th order (classical) Runge-Kutta.\n"
	   "\t\t\"rkf45\"  Embedded Runge-Kutta-Fehlberg (4, 5) method.\n"
	   "\t\t\"rkck\"   Embedded Runge-Kutta Cash-Karp (4, 5) method.\n"
	   "\t\t\"rk8pd\"  Embedded Runge-Kutta Prince-Dormand (8,9) method.\n"
	   "\t\t\"rk2imp\" Implicit 2nd order Runge-Kutta at Gaussian points.\n"
	   "\t\t\"rk4imp\" Implicit 4th order Runge-Kutta at Gaussian points.\n"
	   "\t\t\"gear1\"  M=1 implicit Gear method.\n"
	   "\t\t\"gear2\"  M=2 implicit Gear method.\n");
  fprintf (stdout, "\tode-eps    : GSL ODE control parameter eps\n");
  fprintf (stdout, "* system parameters\n");
  fprintf (stdout, "\tnp     : number of ALL (mobile and fixed) particles\n");
  fprintf (stdout, "\tnm     : number of mobile particles\n");
  fprintf (stdout, "\tx      : particle configuration"
	   " (list or vector with length 3*np)\n");
  fprintf (stdout, "\ta      : radius of particles"
	   " (list or vector with length np)\n"
	   "\t\tby default (if not given), monodisperse system\n");
  fprintf (stdout, "\tslip   : slip length of particles"
	   " (list or vector with length np)\n"
	   "\t\tby default (if not given), no-slip particles\n");
  fprintf (stdout, "\tUi     : imposed translational velocity"
	   " (list or vector of length 3)\n");
  fprintf (stdout, "\tOi     : imposed angular velocity"
	   " (list or vector of length 3)\n");
  fprintf (stdout, "\tEi     : imposed strain"
	   " (list or vector of length 5)\n");
  fprintf (stdout,
	   "\tshear-mode : 0 imposed flow is given by Ui, Oi, Ei (default)\n"
	   "\t             1 for simple shear (x = flow dir, y = grad dir)\n"
	   "\t             2 for simple shear (x = flow dir, z = grad dir)\n"
	   "\t             NOTE: (Ui,Oi,Ei) is overwritten for shear-mode != 0\n");
  fprintf (stdout, "\tshear-rate : for the shear-mode != 0\n");
  fprintf (stdout, "\tshear-shift: the initial cell-shift for the shear-mode != 0\n");
  fprintf (stdout, "\tF0     : applied force"
	   " (list or vector of length 3)\n");
  fprintf (stdout, "\tT0     : applied torque"
	   " (list or vector of length 3)\n");
  fprintf (stdout, "\tstokes : effective stokes number\n");
  fprintf (stdout, "* bond parameters (for chains)\n");
  fprintf (stdout, "\tbonds      : bonds among particles,"
	   " list in the following form\n"
           "\t(define bonds '(\n"
           "\t  (; bond 1\n"
           "\t   0       ; 1) spring type\n"
           "\t   (       ; 2) spring parameters (list with 3 elements)\n"
           "\t    0      ;    fene = 0 means (p1, p2) = (A^{sp}, L_{s})\n"
           "\t    1.0    ;    p1   = A^{sp}, scaled spring constant\n"
           "\t    2.1)   ;    p2   = L_{s} / a, scaled max extension\n"
           "\t   ((0 1)  ; 3) list of pairs\n"
           "\t    (1 2)\n"
           "\t    (2 3))\n"
           "\t    -1)    ; 4) number of exclusion for lubrication\n"
           "\t           ;    negative means all particles in the chain is excluded.\n"
           "\t  (; bond 2\n"
           "\t   2       ; 1) spring type\n"
           "\t   (       ; 2) spring parameters (list with 3 elements)\n"
           "\t    1      ;    fene = 1 means (p1, p2) = (N_{K,s}, b_{K})\n"
           "\t    19.8   ;    p1 = N_{K,s}, the Kuhn steps for a spring\n"
           "      106.0) ;    p2 = b_{K} [nm], the Kuhn length\n"
           "\t           ;    note that, for dWLC (type == 6),\n"
           "             ;    (p1, p2) = (k, r0 [nm]), where the potential is\n"
           "\t           ;    (k/2) * (kT / r0^2) * (r-r0)^2\n"
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
           "\t  6 : Hookean spring for dWLC\n"
	   );
  fprintf (stdout,
	   "\tflag_relax : #f stokesian dynamics,\n"
	   "\t           : #t relaxation dynamics for bonds.\n");
  fprintf (stdout,
	   "\tgamma      : friction coefficient for the relaxation dynamics\n");
  fprintf (stdout, "* Excluded-Volume parameters\n");
  fprintf (stdout,
	   "\tev-v       : v [nm^3] for each chain type\n"
	   "\t             (list or vector with length of bonds' length)\n");
  fprintf (stdout,
	   "\tev-lim     : maximum distance for EV interaction [nm]\n");
  fprintf (stdout, "* angle parameters (for chains)\n");
  fprintf (stdout, "\tangles      : angles among particles,"
	   " list in the following form\n"
           "\t(define angles '(\n"
           "\t  (; angle type 1\n"
           "\t   10.0    ; 1) constant (k^{angle})\n"
           "\t   0.0     ; 2) angle in degree (theta_0)\n"
           "\t   0       ; 3) scale flag (0 == scaled)\n"
           "\t           ;    in this case, the above value for k is just used.\n"
           "\t   ((0 1 2); 4) list of triplets\n"
           "\t    (1 2 3)\n"
           "\t    (2 3 4)\n"
           "\t   )\n"
           "\t  )\n"
           "\t  (; angle type 2\n"
           "\t   20.0    ; 1) constant (k^{angle})\n"
           "\t   90.0    ; 2) angle in degree (theta_0)\n"
           "\t   1       ; 3) scale flag (1 == not scaled yet)\n"
           "\t           ;    in this case, the potential is given by\n"
           "\t           ;    (k/2) * kT * (theta - theta_0)^2\n"
           "\t   ((3 4 5); 4) list of triplets\n"
           "\t    (4 5 6)\n"
           "\t   )\n"
           "\t  )\n"
           "\t))\n"
	   );
  fprintf (stdout, "* Excluded-Volume Debye-Huckel parameters\n");
  fprintf (stdout, "\tev-dh       : list in the following form\n"
           "\t(define ev-dh '(\n"
           "\t  ; system parameters\n"
	   "\t  1.0e-6   ; 1) epsilon for the cut-off distance of EV_DH interaction\n"
           "\t  298.0    ; 2) temperature [K]\n"
           "\t  80.0     ; 3) dielectric constant of the solution\n"
           "\t  3.07     ; 4) Debye length [nm]\n"
           "\t  (        ; 5) list of chain types\n"
           "\t   (; chain type 1\n"
           "\t    2.43    ; 1) nu [e/nm]\n"
           "\t    5.00    ; 2) l0 [nm]\n"
           "\t    (0 1 2) ; 3) list of particles\n"
           "\t   )\n"
           "\t   (; chain type 2\n"
           "\t    2.00    ; 1) nu [e/nm]\n"
           "\t    4.00    ; 2) l0 [nm]\n"
           "\t    (3 4)   ; 3) list of particles\n"
           "\t   )\n"
           "\t  )\n"
           "\t))\n"
	   );
  fprintf (stdout, "* Excluded-Volume Lennard-Jones parameters\n");
  fprintf (stdout, "\tev-LJ       : list in the following form\n"
           "\t(define ev-LJ '(\n"
           "\t (; LJ type 1\n"
           "\t  10.0 ; 1) LJ parameter epsilon in kT (so this is dimensionless value)\n"
           "\t  1.0  ; 2) LJ parameter r0 in \"length\" (so this is dimensionless value)\n"
           "\t  (    ; 3) list of particles\n"
           "\t   0 1 2\n"
           "\t  )\n"
           "\t )\n"
           "\t (; LJ type 2\n"
           "\t  8.0  ; 1) LJ parameter epsilon in kT (so this is dimensionless value)\n"
           "\t  2.0  ; 2) LJ parameter r0 in \"length\" (so this is dimensionless value)\n"
           "\t  (    ; 3) list of particles\n"
           "\t   3 4\n"
           "\t  )\n"
           "\t )\n"
           "\t))\n"
	   );
  fprintf (stdout, "* Confinement force parameters\n");
  fprintf (stdout, "\tconfinement : list in the following form\n"
           "\tfor spherical confinement,\n"
           "\t (define confinement '(\n"
           "\t   10.0 ;; LJ parameter epsilon in kT (so this is dimensionless value)\n"
           "\t   1.0  ;; LJ parameter r0 in \"length\" (so this is dimensionless value)\n"
           "\t   \"sphere\"\n"
           "\t   10.0 ;; radius of the cavity at (0, 0, 0)\n"
           "\t ))\n"
           "\tfor spherical confinement with a hole,\n"
           "\t (define confinement '(\n"
           "\t   10.0 ;; LJ parameter epsilon in kT (so this is dimensionless value)\n"
           "\t   1.0  ;; LJ parameter r0 in \"length\" (so this is dimensionless value)\n"
           "\t   \"sphere+hole\"\n"
           "\t   10.0 ;; radius of the cavity at (0, 0, 0)\n"
           "\t   1.0  ;; radius of the hole at (0, 0, 1) direction\n"
           "\t ))\n"
           "\tfor cylindrical confinement,\n"
           "\t (define confinement '(\n"
           "\t   10.0 ;; LJ parameter epsilon in kT (so this is dimensionless value)\n"
           "\t   1.0  ;; LJ parameter r0 in \"length\" (so this is dimensionless value)\n"
           "\t   \"cylinder\" ;; the cylinder center goes through (0,0,0) and (x,y,z)\n"
           "\t   10.0       ;; radius of the cylinder\n"
           "\t   1.0  0.0  0.0 ;; direction vector (x, y, z) of the cylinder\n"
           "\t ))\n"
           "\tfor dumbbell confinement,\n"
           "\t (define confinement '(\n"
           "\t   10.0 ;; LJ parameter epsilon in kT (so this is dimensionless value)\n"
           "\t   1.0  ;; LJ parameter r0 in \"length\" (so this is dimensionless value)\n"
           "\t   \"dumbbell\" ;; the origin is at the center of the cylinder\n"
           "\t   10.0       ;; left cavity radius centered at (center1, 0, 0)\n"
           "\t   10.0       ;; right cavity radius centered at (center2, 0, 0)\n"
           "\t   2.0        ;; length of the cylinder\n"
           "\t   1.0        ;; cylinder radius\n"
           "\t ))\n"
           "\tfor 2D hexagonal confinement with cylinder pipe,\n"
           "\t (define confinement '(\n"
           "\t   10.0 ;; LJ parameter epsilon in kT (so this is dimensionless value)\n"
           "\t   1.0  ;; LJ parameter r0 in \"length\" (so this is dimensionless value)\n"
           "\t   \"hex2d\"\n"
           "\t   10.0    ;; cavity radius\n"
           "\t   1.0     ;; cylinder radius\n"
           "\t   12.0    ;; lattice spacing\n"
           "\t ))\n"
           "\tfor porous media (outside of the 3D hexagonal particle array)\n"
           "\t (define confinement '(\n"
           "\t   10.0 ;; LJ parameter epsilon in kT (so this is dimensionless value)\n"
           "\t   1.0  ;; LJ parameter r0 in \"length\" (so this is dimensionless value)\n"
           "\t   \"porous\"\n"
           "\t   10.0    ;; particle radius\n"
           "\t   20.0    ;; lattice spacing in x (2R for touching case)\n"
           "\t ))\n"
	   );
  fprintf (stdout, "* Brownian dynamics' parameters\n");
  fprintf (stdout,
	   "\tpeclet       : peclet number (negative means no Brownian force)\n");
  fprintf (stdout,
	   "\tlength       : unit of the length scale [nm]\n");
  fprintf (stdout,
	   "\tBD-seed      : seed for Brownian force random number generator\n");
  fprintf (stdout,
	   "\tn-cheb-minv  : number of Chebyshev coefficients for M^{-1}\n");
  fprintf (stdout,
	   "\tn-cheb-lub   : number of Chebyshev coefficients for L\n"
	   "\t\t(if zero is given, use Cholesky decomposition)\n");
  fprintf (stdout,
	   "\tBD-scheme    : Brownian dynamics time integration scheme\n"
	   "\t\t\"mid-point\"        The mid-point algorithm.\n"
	   "\t\t\"BanchioBrady03\"   Banchio and Brady (2003).\n"
	   "\t\t\"BallMelrose97\"    Ball and Melrose (1997).\n"
	   "\t\t\"JendrejackEtal00\" Jendrejack et al (2000).\n"
	   "\t\t\"semi-implicit-PC\" semi-implicit predictor-corrector.\n");
  fprintf (stdout,
	   "\tBB-n         : step parameter for Banchio-Brady03 algorithm.\n");
  fprintf (stdout,
	   "\tBD-nl-solver : nonlinear solver for implicit schemes\n"
	   "\t\t\"GSL\"    GSL multiroot solver\n"
	   "\t\t\"NITSOL\" Newton-GMRES solver by Pernice and Walker (1998)\n");
  fprintf (stdout, "* dt-ajustment parameters for Brownian dynamics\n"
	   "\tNOTE: if rmin is defined by non-zero, the process is skipped.\n");
  fprintf (stdout,
	   "\tBD-rmin : overlap-param for dt-adjustment process in BD.\n"
	   "\t          the condition is (r2 < rmin * a2).\n");
  fprintf (stdout,
	   "\tdt-lim  : lower bound to shrink dt to prevent overlaps.\n"
	   "\t          set equal to \"dt\" if you don't want to adjust \"dt\"\n"
	   "\t          but just reject it.\n");
  fprintf (stdout, "NOTE on the unit of length:\n"
	   "\tin the above, we take [nm] for the unit of length.\n"
	   "\thowever, we can use [micro m] if we replace everything, that is,\n"
	   "\t\t\"length\"\n"
	   "\t\tp2 for \"bonds\" with fene=1\n"
	   "\t\t\"ev-v\"\n"
	   "\t\t\"ev-lim\"\n"
	   "\t\tmax. distance, Debye length, nu, and l0 in \"ev-dh\"\n"
	   );
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

  // flag-noHI
  int flag_noHI = 0;
  if (guile_get_bool ("flag-noHI") != 0) // TRUE
    {
      flag_noHI = 1;
    }

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
  int nout  = guile_get_int    ("nout",   10);

  /**
   * imposed flow Ui, Oi, Ei
   */
  double Ui[3];
  double Oi[3];
  double Ei[5];

  int shear_mode = guile_get_int ("shear-mode", 0);
  double shear_rate = guile_get_double ("shear-rate", 0.0);
  double shear_shift = guile_get_double ("shear-shift", 0.0);
  if (shear_mode == 0)
    {
      if (guile_get_doubles ("Ui", 3, Ui) != 1) // FALSE
	{
	  fprintf (stderr, "Ui is not defined\n");
	  exit (1);
	}
      if (guile_get_doubles ("Oi", 3, Oi) != 1) // FALSE
	{
	  fprintf (stderr, "Oi is not defined\n");
	  exit (1);
	}
      if (guile_get_doubles ("Ei", 5, Ei) != 1) // FALSE
	{
	  fprintf (stderr, "Ei is not defined\n");
	  exit (1);
	}
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

  int flag_BD = 0;
  if (peclet >= 0.0) flag_BD = 1;


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
  else if (strcmp (str_BD_scheme, "semi-implicit-PC") == 0)
    {
      BD_scheme = 4;
    }
  else
    {
      fprintf (stderr, "invalid BD-scheme %s", str_BD_scheme);
      exit (1);
    }
  free (str_BD_scheme);
  // step parameter for BB03 algorithm
  double BB_n = guile_get_double ("BB-n", 100.0);

  // nonlinear solver for implicit schemes
  char *str_BD_nl_solver = guile_get_string ("BD-nl-solver");
  int BD_nl_solver;
  if (strcmp (str_BD_scheme, "GSL") == 0)
    {
      BD_nl_solver = 0;
    }
  else if (strcmp (str_BD_scheme, "NITSOL") == 0)
    {
      BD_nl_solver = 1;
    }
  else
    {
      fprintf (stderr, "invalid BD-nl-solver %s", str_BD_nl_solver);
      exit (1);
    }
  free (str_BD_nl_solver);

  // factor for overlap check for dt-adjustment process in BD
  double BD_rmin = guile_get_double ("BD-rmin", 1.0);
  // lower bound to shrink dt to prevent overlaps
  double dt_lim = guile_get_double ("dt-lim", 1.0e-12);

  /**
   * bonds
   */
  struct bonds *bonds = bonds_guile_get ("bonds");
  if (bonds != NULL)
    {
      bonds_set_FENE (bonds, length, peclet);
    }
  else
    {
      bonds = bonds_init ();
    }
  // bonds == NULL means no bond interaction

  // relaxation dynamics only with bond interaction
  double gamma = guile_get_double ("gamma", 1.0);

  /**
   * excluded volume
   */
  double ev_lim = guile_get_double ("ev-lim", 1.0);
  ev_lim /= length; // scale by length
  double ev_r2 = ev_lim * ev_lim;
  struct EV *ev = EV_guile_get ("ev-v", bonds, length, peclet, ev_r2, np);
  // ev == NULL means no EV interaction

  /**
   * angles
   */
  struct angles *ang = angles_guile_get ("angles");
  // ang == NULL means no angle interaction
  if (ang != NULL)
    {
      angles_scale_k (ang, length, peclet);
    }

  /**
   * EV_DH
   */
  struct EV_DH *ev_dh = EV_DH_guile_get ("ev-dh",
					 length, peclet, np);
  // ev_dh == NULL means no EV-DH interaction

  /**
   * EV_LJ
   */
  struct EV_LJ *ev_LJ = EV_LJ_guile_get ("ev-LJ", np);
  if (ev_LJ != NULL)
    {
      EV_LJ_scale (ev_LJ, peclet);
    }
  // ev_LJ == NULL means no EV-LJ interaction

  /**
   * confinement
   */
  struct confinement *cf = CF_guile_get ("confinement");
  if (cf != NULL)
    {
      CF_set (cf, peclet);
    }
  // cf == NULL means no confinement force


  // initialize struct stokes *sys
  struct stokes *sys = stokes_init ();
  CHECK_MALLOC (sys, "main");
  sys->version = version;
  stokes_set_np (sys, np, nm);

  sys->rmin = guile_get_double ("rmin", 0.0);
  double lubmin = guile_get_double ("lub-min", 2.0000000001);
  sys->lubmin2 = lubmin * lubmin;
  sys->lubmax = guile_get_double ("lub-max", 4.0);

  /* set exclusion list for lub by bonds */
  if (flag_noHI == 0 && flag_lub != 0)
    {
      /* currently, it is only used in the lubrication calculation */
      list_ex_set_by_bonds (sys->ex_lub, bonds);
    }

  // iterative solver
  char *str_it_solver = NULL;
  str_it_solver = guile_get_string ("IT-solver");
  int it_max = guile_get_int ("IT-max", 2000);
  int it_n   = guile_get_int ("IT-n", 20);
  double it_eps = guile_get_double ("IT-eps", 1.0e-6);
  int it_debug = guile_get_int ("IT-debug", 0);
  stokes_set_iter (sys, str_it_solver, it_max, it_n, it_eps, it_debug, stderr);
  free (str_it_solver);


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

  // imposed flow
  if (shear_mode == 0)
    {
      stokes_set_Ui (sys, Ui[0], Ui[1], Ui[2]);
      stokes_set_Oi (sys, Oi[0], Oi[1], Oi[2]);
      stokes_set_Ei (sys, Ei[0], Ei[1], Ei[2], Ei[3], Ei[4]);
    }
  else
    {
      stokes_set_shear (sys, shear_mode, shear_rate);
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

      /*
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
      */
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
  if (flag_BD == 0) // no Brownian force
    {
      ode_params = ode_params_init (sys,
				    F, T, E,
				    uf, of, ef,
				    flag_noHI,
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
				  flag_noHI,
				  flag_lub, flag_mat,
				  st,
				  bonds,
				  gamma,
				  ev,
				  ang,
				  ev_dh,
				  ev_LJ,
				  cf,
				  flag_Q,
				  peclet,
				  ode_eps,
				  n_minv, // n of chebyshev for minv
				  n_lub,  // n of chebyshev for lub
				  BD_scheme,
				  BB_n,
				  BD_rmin,
				  dt_lim);
      CHECK_MALLOC (BD_params, "main");

      // implicit scheme
      if (BD_scheme > 2) // either JendrejackEtal00 or semi-implicit-PC
	{
	  BDimp = BD_imp_init (BD_params,
			       BD_nl_solver,
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
      l0 = nc->ntime - 1;
      /* l0 is the starting index because index starts from 0
       * so that we need to access "l0-1" to get the last step infos
       */
      //nloop = l0 + nloop_arg + 1;
      nloop = l0 + nloop_arg; // I guess we don't need "+1" here.
      t  = stokes_nc_get_time_step (nc, l0);

      // set the configuration at the current time
      stokes_nc_get_data (nc, "x", l0, y);
      if (flag_Q != 0)
	{
	  stokes_nc_get_data (nc, "q", l0, y + nq);
	}
      if (sys->shear_mode != 0)
	{
	  int status = nc_get_var1_double
	    (nc->id, nc->shear_shift_id, NULL, &shear_shift);
	  if (status != NC_NOERR)
	    {
	      fprintf (stderr,
		       "at nc_get_var1_double() for shear_shift in main\n");
	    }
	}
      if (flag_BD != 0) // Brownian
	{
	  stokes_nc_get_rng (nc, l0, BD_params->rng);
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
				    lat,
				    shear_mode, shear_rate,
				    flag_BD);
      CHECK_MALLOC (nc, "main");

      // set the loop parameters
      l0 = 0;
      // nloop is given by the script
      t = 0.0;

      // output the configuration at t = 0.0 (l0 = 0)
      // NOTE: for continuation, the last configuration at l0 step
      stokes_nc_set_time (nc, l0, t);
      stokes_nc_set_x (nc, l0, y);
      if (flag_Q != 0)
	{
	  stokes_nc_set_q (nc, l0, y + nq);
	}
      if (sys->shear_mode != 0)
	{
	  stokes_nc_set_shear_shift (nc, l0, shear_shift);
	}
      if (flag_BD != 0) // Brownian
	{
	  stokes_nc_set_rng (nc, l0, BD_params->rng);
	}

      // flush the data
      nc_sync(nc->id);

    }
  free (out_file);


  /**
   * mail loop
   */
  double h = dt; // initial time step for ODE integrator
  double dt_out = dt * (double)nout;
  int l;
  double ptime0 = ptime_ms_d();
  for (l = l0; l < nloop; l++)
    {
      fprintf (stdout, "%d steps\n", l + 1);

      double t_out = t + dt_out;

      // set reference for the shear_shift at t = t (not t_out)
      if (flag_BD == 0) // non Brownian
	{
	  ode_set_shear_shift_ref (ode_params, t, shear_shift);
	}
      else // Brownian
	{
	  BD_set_shear_shift_ref (BD_params, t, shear_shift);
	}

      // loop for dt_out
      if (flag_BD == 0) // non Brownian
	{
	  // GSL integrator is adaptive scheme
	  while (t < t_out)
	    {
	      if (st == 0.0)
		{
		  // no collision checks
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
	      else // st != 0.0
		{
		  // collision checks for nout times in dt_out (every dt step)
		  for (i = 0; i < nout; i ++)
		    {
		      double t_end = t + dt;
		      gsl_odeiv_evolve_apply (GSL_ODE_EVOLVE,
					      GSL_ODE_CONTROL,
					      GSL_ODE_STEP,
					      &GSL_ODE_SYSTEM,
					      &t, t_end,
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
		    } // end of nout loop
		} // end if st == 0
	    }
	}
      else // Brownian dynamics (peclet >= 0)
	{
	  // discretize dt_out by nout
	  for (i = 0; i < nout; i ++)
	    {
	      double t_end = t + dt;

	      if (BD_params->scheme > 2)
		{
		  // either JendrejackEtal00 or semi-implicit-PC
		  BD_imp_ode_evolve (BDimp,
				     &t, t_end,
				     &h, y);
		}
	      else
		{
		  // explicit schemes
		  BD_ode_evolve (BD_params,
				 &t, t_end,
				 &h, y);
		}
	      // check periodicity
	      if (sys->periodic == 1)
		{
		  check_periodic (sys, y);
		}
	    }
	}
      // end of the loop for dt_out

      // output the results at t (should be t_out) with (l+1) step
      stokes_nc_set_time (nc, l+1, t);
      stokes_nc_set_x (nc, l+1, y);
      if (flag_Q != 0)
	{
	  stokes_nc_set_q (nc, l+1, y + nq);
	}
      if (sys->shear_mode != 0)
	{
	  if (flag_BD == 0) // non Brownian
	    {
	      shear_shift
		= stokes_get_shear_shift (sys, t,
					  ode_params->t0, ode_params->s0);
	    }
	  else // Brownian
	    {
	      shear_shift
		= stokes_get_shear_shift (sys, t,
					  BD_params->t0, BD_params->s0);
	    }
	  stokes_nc_set_shear_shift (nc, l+1, shear_shift);
	}
      if (flag_BD != 0) // Brownian
	{
	  stokes_nc_set_rng (nc, l+1, BD_params->rng);
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
  angles_free (ang);

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
