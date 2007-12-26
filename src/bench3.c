/* benchmark code for stokes simulator in 3D for F/FT/FTS versions
 * Copyright (C) 1997-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bench3.c,v 1.7 2007/12/26 06:47:19 kichiki Exp $
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

#include <libiter.h> /* struct iter */
#include <libstokes.h> /* struct stokes */
#include "configs.h"


/* make label
 */
static void
make_label (int flag_prob,
	    int flag_lub,
	    int flag_mat,
	    char *iter_solver,
	    int version,
	    int flag_periodic,
	    char *label)
{
  // version
  char l_ver[4];
  switch (version)
    {
    case 0:
      strcpy (l_ver, "F");
      break;
    case 1:
      strcpy (l_ver, "FT");
      break;
    case 2:
      strcpy (l_ver, "FTS");
      break;
    default:
      strcpy (l_ver, "unknown");
      break;
    }
  // periodic
  char l_per[3];
  if (flag_periodic == 0)
    {
      strcpy (l_per, "NP");
    }
  else
    {
      strcpy (l_per, "P");
    }
  // problem
  char l_prb[4];
  if (flag_prob == 0)
    {
      strcpy (l_prb, "mob");
    }
  else if (flag_prob == 1)
    {
      strcpy (l_prb, "mix");
    }
  else
    {
      strcpy (l_prb, "res");
    }
  // lub
  char l_lub[6];
  if (flag_lub == 0)
    {
      strcpy (l_lub, "nolub");
    }
  else
    {
      strcpy (l_lub, "lub");
    }
  // mat
  if (flag_mat == 0)
    {
      sprintf (label, "%s %s %s %s %s",
	       l_ver, l_per, l_prb, l_lub, iter_solver);
    }
  else
    {
      sprintf (label, "%s %s %s %s %s",
	       l_ver, l_per, l_prb, l_lub, "matrix");
    }
}

/* solve the problem
 * OUTPUT
 *  returned value : CPU time
 */
static double
solve_problem (struct stokes *sys,
	       int flag_prob, int flag_lub, int flag_mat, int version,
	       double *f,
	       double *t,
	       double *e,
	       double *uf,
	       double *of,
	       double *ef,
	       double *u,
	       double *o,
	       double *s,
	       double *ff,
	       double *tf,
	       double *sf)
{
  double t0 = 0.0; // for compiler warning
  double t1 = 0.0;

  /* call atimes routine */
  if (flag_prob != 2) // mobility or mixed problem
    {
      if (flag_lub == 0) // no lubrication
	{
	  if (flag_mat == 0)
	    {
	      switch (version)
		{
		case 0: // F version
		  t0 = ptime_ms_d();
		  solve_mix_3f
		    (sys, f, uf,
		     u, ff);
		  t1 = ptime_ms_d();
		  break;
		case 1: // FT version
		  t0 = ptime_ms_d();
		  solve_mix_3ft
		    (sys, f, t, uf, of,
		     u, o, ff, tf);
		  t1 = ptime_ms_d();
		  break;
		case 2: // FTS version
		  t0 = ptime_ms_d();
		  solve_mix_3fts
		    (sys, f, t, e, uf, of, ef,
		     u, o, s, ff, tf, sf);
		  t1 = ptime_ms_d();
		  break;
		default:
		  fprintf (stderr, "invalid version %d\n", version);
		  break;
		}
	    }
	  else
	    {
	      switch (version)
		{
		case 0: // F version
		  t0 = ptime_ms_d();
		  solve_mix_3f_matrix
		    (sys, f, uf,
		     u, ff);
		  t1 = ptime_ms_d();
		  break;
		case 1: // FT version
		  t0 = ptime_ms_d();
		  solve_mix_3ft_matrix
		    (sys, f, t, uf, of,
		     u, o, ff, tf);
		  t1 = ptime_ms_d();
		  break;
		case 2: // FTS version
		  t0 = ptime_ms_d();
		  solve_mix_3fts_matrix
		    (sys, f, t, e, uf, of, ef,
		     u, o, s, ff, tf, sf);
		  t1 = ptime_ms_d();
		  break;
		default:
		  fprintf (stderr, "invalid version %d\n", version);
		  break;
		}
	    }
	}
      else // with lubrication
	{
	  if (flag_mat == 0)
	    {
	      switch (version)
		{
		case 0: // F version
		  t0 = ptime_ms_d();
		  solve_mix_lub_3f
		    (sys, f, uf,
		     u, ff);
		  t1 = ptime_ms_d();
		  break;
		case 1: // FT version
		  t0 = ptime_ms_d();
		  solve_mix_lub_3ft
		    (sys, f, t, uf, of,
		     u, o, ff, tf);
		  t1 = ptime_ms_d();
		  break;
		case 2: // FTS version
		  t0 = ptime_ms_d();
		  solve_mix_lub_3fts
		    (sys, f, t, e, uf, of, ef,
		     u, o, s, ff, tf, sf);
		  t1 = ptime_ms_d();
		  break;
		default:
		  fprintf (stderr, "invalid version %d\n", version);
		  break;
		}
	    }
	  else
	    {
	      switch (version)
		{
		case 0: // F version
		  t0 = ptime_ms_d();
		  solve_mix_lub_3f_matrix
		    (sys, f, uf,
		     u, ff);
		  t1 = ptime_ms_d();
		  break;
		case 1: // FT version
		  t0 = ptime_ms_d();
		  solve_mix_lub_3ft_matrix
		    (sys, f, t, uf, of,
		     u, o, ff, tf);
		  t1 = ptime_ms_d();
		  break;
		case 2: // FTS version
		  t0 = ptime_ms_d();
		  solve_mix_lub_3fts_matrix
		    (sys, f, t, e, uf, of, ef,
		     u, o, s, ff, tf, sf);
		  t1 = ptime_ms_d();
		  break;
		default:
		  fprintf (stderr, "invalid version %d\n", version);
		  break;
		}
	    }
	}
    }
  else // resistance problem
    {
      if (flag_lub == 0) // no lubrication
	{
	  if (flag_mat == 0)
	    {
	      switch (version)
		{
		case 0: // F version
		  t0 = ptime_ms_d();
		  solve_res_3f
		    (sys, u,
		     f);
		  t1 = ptime_ms_d();
		  break;
		case 1: // FT version
		  t0 = ptime_ms_d();
		  solve_res_3ft
		    (sys, u, o,
		     f, t);
		  t1 = ptime_ms_d();
		  break;
		case 2: // FTS version
		  t0 = ptime_ms_d();
		  solve_res_3fts
		    (sys, u, o, e,
		     f, t, s);
		  t1 = ptime_ms_d();
		  break;
		default:
		  fprintf (stderr, "invalid version %d\n", version);
		  break;
		}
	    }
	  else
	    {
	      switch (version)
		{
		case 0: // F version
		  t0 = ptime_ms_d();
		  solve_res_3f_matrix
		    (sys, u,
		     f);
		  t1 = ptime_ms_d();
		  break;
		case 1: // FT version
		  t0 = ptime_ms_d();
		  solve_res_3ft_matrix
		    (sys, u, o,
		     f, t);
		  t1 = ptime_ms_d();
		  break;
		case 2: // FTS version
		  t0 = ptime_ms_d();
		  solve_res_3fts_matrix
		    (sys, u, o, e,
		     f, t, s);
		  t1 = ptime_ms_d();
		  break;
		default:
		  fprintf (stderr, "invalid version %d\n", version);
		  break;
		}
	    }
	}
      else // with lubrication
	{
	  if (flag_mat == 0)
	    {
	      switch (version)
		{
		case 0: // F version
		  t0 = ptime_ms_d();
		  solve_res_lub_3f
		    (sys, u,
		     f);
		  t1 = ptime_ms_d();
		  break;
		case 1: // FT version
		  t0 = ptime_ms_d();
		  solve_res_lub_3ft
		    (sys, u, o,
		     f, t);
		  t1 = ptime_ms_d();
		  break;
		case 2: // FTS version
		  t0 = ptime_ms_d();
		  solve_res_lub_3fts
		    (sys, u, o, e,
		     f, t, s);
		  t1 = ptime_ms_d();
		  break;
		default:
		  fprintf (stderr, "invalid version %d\n", version);
		  break;
		}
	    }
	  else
	    {
	      switch (version)
		{
		case 0: // F version
		  t0 = ptime_ms_d();
		  solve_res_lub_3f_matrix
		    (sys, u,
		     f);
		  t1 = ptime_ms_d();
		  break;
		case 1: // FT version
		  t0 = ptime_ms_d();
		  solve_res_lub_3ft_matrix
		    (sys, u, o,
		     f, t);
		  t1 = ptime_ms_d();
		  break;
		case 2: // FTS version
		  t0 = ptime_ms_d();
		  solve_res_lub_3fts_matrix
		    (sys, u, o, e,
		     f, t, s);
		  t1 = ptime_ms_d();
		  break;
		default:
		  fprintf (stderr, "invalid version %d\n", version);
		  break;
		}
	    }
	}
    }

  return (t1 - t0);
}

/*
 * INPUT
 *  n_sc : number of particles in one-direction for SC lattice.
 *         so that the total number of particles is (n_sc)^3.
 *  phi  : volume fraction for the configuration
 *  flag_prob : 0 == mobility problem
 *              1 == mixed problem
 *              2 == resistance problem
 */
void
bench_SC (int n_sc, double phi,
	  double ewald_tr,
	  double ewald_eps,
	  int flag_prob,
	  int flag_lub,
	  int flag_mat,
	  char *iter_solver,
	  int version,
	  int flag_periodic)
{
  struct stokes *sys = stokes_init ();
  sys->version = version;

  int np, nm, nf;
  np = n_sc * n_sc * n_sc;
  if (flag_prob == 1)
    {
      // mixed problem
      nm = np - np/2;
      nf = np - nm;
    }
  else
    {
      // no fixed particles
      nm = np;
      nf = 0;
    }
  int np3, nm3, nf3;
  int np5, nm5, nf5;
  np3 = np * 3;
  nm3 = nm * 3;
  nf3 = nf * 3;
  np5 = np * 5;
  nm5 = nm * 5;
  nf5 = nf * 5;

  stokes_set_np (sys, np, nm);

  double lx, ly, lz;
  init_config_SC (phi, n_sc, n_sc, n_sc,
		  sys->pos, &lx, &ly, &lz);
  if (flag_periodic == 0)
    {
      // non-periodic boundary condition
      sys->periodic = 0;
    }
  else
    {
      // periodic boundary condition
      sys->periodic = 1;
      stokes_set_l (sys, lx, ly, lz);

      double xi = xi_by_tratio (sys, ewald_tr);
      stokes_set_xi (sys, xi, ewald_eps);
    }

  sys->lubmin2 = 4.0000000001;
  sys->lubmax = 4.0;
  stokes_set_iter (sys, iter_solver, 2000, 20, 1.0e-6, 1, stderr);

  double *f = NULL;
  double *t = NULL;
  double *e = NULL;
  double *uf = NULL;
  double *of = NULL;
  double *ef = NULL;
  double *u = NULL;
  double *o = NULL;
  double *s = NULL;
  double *ff = NULL;
  double *tf = NULL;
  double *sf = NULL;
  int i;

  if (flag_prob != 2) // mobility or mixed problem
    {
      switch (version)
	{
	case 0: // F version
	  f  = (double *)malloc (sizeof (double) * nm3);
	  u  = (double *)malloc (sizeof (double) * nm3);
	  uf = (double *)malloc (sizeof (double) * nf3);
	  ff = (double *)malloc (sizeof (double) * nf3);
	  for (i = 0; i < nm3; i ++)
	    {
	      f[i] = 1.0;
	    }
	  for (i = 0; i < nf3; i ++)
	    {
	      uf[i] = 0.0;
	    }
	  break;
	case 1: // FT version
	  f  = (double *)malloc (sizeof (double) * nm3);
	  t  = (double *)malloc (sizeof (double) * nm3);
	  u  = (double *)malloc (sizeof (double) * nm3);
	  o  = (double *)malloc (sizeof (double) * nm3);
	  uf = (double *)malloc (sizeof (double) * nf3);
	  of = (double *)malloc (sizeof (double) * nf3);
	  ff = (double *)malloc (sizeof (double) * nf3);
	  tf = (double *)malloc (sizeof (double) * nf3);
	  for (i = 0; i < nm3; i ++)
	    {
	      f[i] = 1.0;
	      t[i] = 1.0;
	    }
	  for (i = 0; i < nf3; i ++)
	    {
	      uf[i] = 0.0;
	      of[i] = 0.0;
	    }
	  break;
	case 2: // FTS version
	  f  = (double *)malloc (sizeof (double) * nm3);
	  t  = (double *)malloc (sizeof (double) * nm3);
	  s  = (double *)malloc (sizeof (double) * nm5);
	  u  = (double *)malloc (sizeof (double) * nm3);
	  o  = (double *)malloc (sizeof (double) * nm3);
	  e  = (double *)malloc (sizeof (double) * nm5);
	  uf = (double *)malloc (sizeof (double) * nf3);
	  of = (double *)malloc (sizeof (double) * nf3);
	  ef = (double *)malloc (sizeof (double) * nf5);
	  ff = (double *)malloc (sizeof (double) * nf3);
	  tf = (double *)malloc (sizeof (double) * nf3);
	  sf = (double *)malloc (sizeof (double) * nf5);
	  for (i = 0; i < nm3; i ++)
	    {
	      f[i] = 1.0;
	      t[i] = 1.0;
	    }
	  for (i = 0; i < nm5; i ++)
	    {
	      e[i] = 1.0;
	    }
	  for (i = 0; i < nf3; i ++)
	    {
	      uf[i] = 0.0;
	      of[i] = 0.0;
	    }
	  for (i = 0; i < nf5; i ++)
	    {
	      ef[i] = 1.0;
	    }
	  break;
	default:
	  fprintf (stderr, "invalid version %d\n", version);
	  break;
	}
    }
  else // resistance problem
    {
      switch (version)
	{
	case 0: // F version
	  f  = (double *)malloc (sizeof (double) * np3);
	  u  = (double *)malloc (sizeof (double) * np3);
	  for (i = 0; i < np3; i ++)
	    {
	      u[i] = 1.0;
	    }
	  break;
	case 1: // FT version
	  f  = (double *)malloc (sizeof (double) * np3);
	  t  = (double *)malloc (sizeof (double) * np3);
	  u  = (double *)malloc (sizeof (double) * np3);
	  o  = (double *)malloc (sizeof (double) * np3);
	  for (i = 0; i < np3; i ++)
	    {
	      u[i] = 1.0;
	      o[i] = 1.0;
	    }
	  break;
	case 2: // FTS version
	  f  = (double *)malloc (sizeof (double) * np3);
	  t  = (double *)malloc (sizeof (double) * np3);
	  s  = (double *)malloc (sizeof (double) * np5);
	  u  = (double *)malloc (sizeof (double) * np3);
	  o  = (double *)malloc (sizeof (double) * np3);
	  e  = (double *)malloc (sizeof (double) * np5);
	  for (i = 0; i < np3; i ++)
	    {
	      u[i] = 1.0;
	      o[i] = 1.0;
	    }
	  for (i = 0; i < np5; i ++)
	    {
	      e[i] = 1.0;
	    }
	  break;
	default:
	  fprintf (stderr, "invalid version %d\n", version);
	  break;
	}
    }

  double dt
    = solve_problem (sys,
		     flag_prob, flag_lub, flag_mat, version,
		     f, t, e, uf, of, ef,
		     u, o, s, ff, tf, sf);

  /* make label */
  char label[80];
  make_label (flag_prob,
	      flag_lub,
	      flag_mat,
	      iter_solver,
	      version,
	      flag_periodic,
	      label);


  /* analysis */
  double av_u;
  double av_ff;
  double av_o;
  double av_tf;
  double av_e;
  double av_ef;

  double av_f;
  double av_t;
  double av_s;

  if (flag_prob == 0) // mobility problem
    {
      switch (version)
	{
	case 0: // F version
	  av_u  = 0.0;
	  for (i = 0; i < nm3; i ++)
	    {
	      av_u  += u[i];
	    }
	  av_u  /= (double) nm3;

	  av_ff = 0.0;
	  for (i = 0; i < nf3; i ++)
	    {
	      av_ff += ff[i];
	    }
	  if (nf > 0) av_ff /= (double) nf3;

	  fprintf (stdout, "%s %d %d %.3f %e %e\n",
		   label,
		   np, nm, dt, av_u, av_ff);

	  free (f);
	  free (u);
	  free (uf);
	  free (ff);
	  break;
	case 1: // FT version
	  av_u  = 0.0;
	  av_ff = 0.0;
	  av_o  = 0.0;
	  av_tf = 0.0;
	  for (i = 0; i < nm3; i ++)
	    {
	      av_u  += u[i];
	      av_o  += o[i];
	    }
	  av_u  /= (double) nm3;
	  av_o  /= (double) nm3;

	  for (i = 0; i < nf3; i ++)
	    {
	      av_ff += ff[i];
	      av_tf += tf[i];
	    }
	  if (nf > 0)
	    {
	      av_ff /= (double) nf3;
	      av_tf /= (double) nf3;
	    }

	  fprintf (stdout, "%s %d %d %.3f %e %e %e %e\n",
		   label,
		   np, nm, dt, av_u, av_ff, av_o, av_tf);

	  free (f);
	  free (t);
	  free (u);
	  free (o);
	  free (uf);
	  free (of);
	  free (ff);
	  free (tf);
	  break;
	case 2: // FTS version
	  av_u  = 0.0;
	  av_o  = 0.0;
	  for (i = 0; i < nm3; i ++)
	    {
	      av_u  += u[i];
	      av_o  += o[i];
	    }
	  av_u  /= (double) nm3;
	  av_o  /= (double) nm3;

	  av_ff = 0.0;
	  av_tf = 0.0;
	  for (i = 0; i < nf3; i ++)
	    {
	      av_ff += ff[i];
	      av_tf += tf[i];
	    }
	  if (nf > 0)
	    {
	      av_ff /= (double) nf3;
	      av_tf /= (double) nf3;
	    }

	  av_e  = 0.0;
	  for (i = 0; i < nm5; i ++)
	    {
	      av_e  += e[i];
	    }
	  av_e /= (double) nm5;

	  av_ef = 0.0;
	  for (i = 0; i < nf5; i ++)
	    {
	      av_ef += ef[i];
	    }
	  if (nf > 0) av_ef /= (double) nf5;

	  fprintf (stdout, "%s %d %d %.3f %e %e %e %e %e %e\n",
		   label,
		   np, nm, dt, av_u, av_ff, av_o, av_tf, av_e, av_ef);

	  free (f);
	  free (t);
	  free (s);
	  free (u);
	  free (o);
	  free (e);
	  free (uf);
	  free (of);
	  free (ef);
	  free (ff);
	  free (tf);
	  free (sf);
	  break;
	default:
	  fprintf (stderr, "invalid version %d\n", version);
	  break;
	}
    }
  else // resistance problem
    {
      switch (version)
	{
	case 0: // F version
	  av_f  = 0.0;
	  for (i = 0; i < np3; i ++)
	    {
	      av_f  += f[i];
	    }
	  av_f  /= (double) np3;

	  fprintf (stdout, "%s %d %d %.3f %e\n",
		   label,
		   np, nm, dt, av_f);

	  free (f);
	  free (u);
	  break;
	case 1: // FT version
	  av_f  = 0.0;
	  av_t  = 0.0;
	  for (i = 0; i < np3; i ++)
	    {
	      av_f  += f[i];
	      av_t  += t[i];
	    }
	  av_f  /= (double) np3;
	  av_t  /= (double) np3;

	  fprintf (stdout, "%s %d %d %.3f %e %e\n",
		   label,
		   np, nm, dt, av_f, av_t);

	  free (f);
	  free (t);
	  free (u);
	  free (o);
	  break;
	case 2: // FTS version
	  av_f  = 0.0;
	  av_t  = 0.0;
	  av_s  = 0.0;
	  for (i = 0; i < np3; i ++)
	    {
	      av_f  += f[i];
	      av_t  += t[i];
	    }
	  av_f  /= (double) np3;
	  av_t  /= (double) np3;

	  for (i = 0; i < np5; i ++)
	    {
	      av_s  += s[i];
	    }
	  av_s  /= (double) np5;

	  fprintf (stdout, "%s %d %d %.3f %e %e %e\n",
		   label,
		   np, nm, dt, av_f, av_t, av_s);

	  free (f);
	  free (t);
	  free (s);
	  free (u);
	  free (o);
	  free (e);
	  break;
	default:
	  fprintf (stderr, "invalid version %d\n", version);
	  break;
	}
    }


  stokes_free (sys);
}


void
usage (char * argv0)
{
  fprintf (stdout, "Benchmark Test for libstokes\n");
  fprintf (stdout, "$Id: bench3.c,v 1.7 2007/12/26 06:47:19 kichiki Exp $\n\n");
  fprintf (stdout, "USAGE\n");
  fprintf (stdout, "%s [options]\n", argv0);
  fprintf (stdout, "OPTIONS\n");
  fprintf (stdout, "\t-f, -ft, -fts : specify version (default F)\n");
  fprintf (stdout, "\t-prob         : 0 for mobility problem (default)\n"
	   "\t              : 1 for mixed problem\n"
	   "\t              : 2 for resistance problem\n");
  fprintf (stdout, "\t-lub          : calc with lub"
	   " (default: no lub)\n");
  fprintf (stdout, "\t-mat          : use matrix scheme"
	   " (default: atimes)\n");
  fprintf (stdout, "\t-open         : calc under non-periodic B.C."
	   " (default: periodic B.C.)\n");
  fprintf (stdout, "\t-tr           : ewald_tr (default: 10.0)\n");
  fprintf (stdout, "\t-eps          : cut-off limit for ewald-sum"
	   " (default: 1.0e-12)\n");
  fprintf (stdout, "\t-iter         : iterative scheme. ignored for -mat.\n");
  fprintf (stdout, "\t\t0 = cg\n");
  fprintf (stdout, "\t\t1 = cgs\n");
  fprintf (stdout, "\t\t2 = bicgstab\n");
  fprintf (stdout, "\t\t3 = sta (another bicgstab)\n");
  fprintf (stdout, "\t\t4 = sta2\n");
  fprintf (stdout, "\t\t5 = gpb\n");
  fprintf (stdout, "\t\t6 = otmk\n");
  fprintf (stdout, "\t\t7 = gmres (default)\n");
  fprintf (stdout, "\t-all : calc all possible options. "
	   "(above options are just ignored.)\n");
}


/* main program */
int
main (int argc, char** argv)
{
  /* option analysis */
  int version   = 0;
  int flag_mat  = 0;
  int flag_lub  = 0;
  int flag_prob = 0;
  int flag_iter = 7; // gmres
  int flag_periodic = 1; // periodic B.C.
  int flag_all = 0;
  double ewald_eps = 1.0e-12;
  double ewald_tr  = 10.0;
  int i;
  for (i = 1; i < argc; i++)
    {
      if (strcmp (argv [i], "-mat") == 0)
	{
	  flag_mat = 1;
	}
      else if (strcmp (argv [i], "-lub") == 0)
	{
	  flag_lub = 1;
	}
      else if (strcmp (argv [i], "-prob") == 0)
	{
	  flag_prob = atoi (argv [++i]);
	}
      else if (strcmp (argv [i], "-f") == 0)
	{
	  version = 0;
	}
      else if (strcmp (argv [i], "-ft") == 0)
	{
	  version = 1;
	}
      else if (strcmp (argv [i], "-fts") == 0)
	{
	  version = 2;
	}
      else if (strcmp (argv [i], "-tr") == 0)
	{
	  ewald_tr = atof (argv [++i]);
	}
      else if (strcmp (argv [i], "-eps") == 0)
	{
	  ewald_eps = atof (argv [++i]);
	}
      else if (strcmp (argv [i], "-iter") == 0)
	{
	  flag_iter = atoi (argv [++i]);
	}
      else if (strcmp (argv [i], "-open") == 0)
	{
	  flag_periodic = 0;
	}
      else if (strcmp (argv [i], "-all") == 0)
	{
	  flag_all = 1;
	}
      else
	{
	  usage (argv[0]);
	  exit (1);
	}
    }

  char iter_solver [80];
  switch (flag_iter)
    {
    case 0:
      strcpy (iter_solver, "cg");
      break;
    case 1:
      strcpy (iter_solver, "cgs");
      break;
    case 2:
      strcpy (iter_solver, "bicgstab");
      break;
    case 3:
      strcpy (iter_solver, "sta");
      break;
    case 4:
      strcpy (iter_solver, "sta2");
      break;
    case 5:
      strcpy (iter_solver, "gpb");
      break;
    case 6:
      strcpy (iter_solver, "otmk");
      break;
    case 7:
      strcpy (iter_solver, "gmres");
      break;
    default:
      fprintf (stderr, "invalid -iter option\n");
      usage (argv[0]);
      exit (1);
    }

  /* write header
  if (version == 0)
    {
      fprintf (stdout, "# F version");
    }
  else if (version == 1)
    {
      fprintf (stdout, "# FT version");
    }
  else
    {
      fprintf (stdout, "# FTS version");
    }

  if (flag_prob == 0)
    {
      fprintf (stdout, " mob");
      if (flag_fix == 1)
	{
	  fprintf (stdout, " with fix");
	}
    }
  else
    {
      fprintf (stdout, " res");
    }

  if (flag_lub == 0)
    {
      fprintf (stdout, " no-lub");
    }
  else
    {
      fprintf (stdout, " with lub");
    }

  if (flag_mat == 0)
    {
      fprintf (stdout, " iterative(%s)\n", iter_solver);
    }
  else
    {
      fprintf (stdout, " matrix\n");
    }
  */

  /* main loop */
  double phi = 0.3;
  for (i = 2; i < 100; i ++)
    {
      if (flag_all == 0)
	{
	  bench_SC (i, phi,
		    ewald_tr,
		    ewald_eps,
		    flag_prob,
		    flag_lub,
		    flag_mat,
		    iter_solver,
		    version,
		    flag_periodic);
	}
      else
	{
	  for (flag_periodic = 0; flag_periodic < 2; flag_periodic ++)
	    {
	      for (flag_prob = 0; flag_prob < 3; flag_prob ++)
		{
		  for (flag_lub = 0; flag_lub < 2; flag_lub ++)
		    {
		      for (version = 0; version < 3; version ++)
			{
			  char label[80];
			  make_label (flag_prob,
				      flag_lub,
				      1,  //flag_mat,
				      "", //iter_solver,
				      version,
				      flag_periodic,
				      label);
			  fprintf (stderr, "# %s\n", label);

			  // matrix
			  flag_mat = 1;
			  bench_SC (i, phi,
				    ewald_tr,
				    ewald_eps,
				    flag_prob,
				    flag_lub,
				    flag_mat,
				    iter_solver,
				    version,
				    flag_periodic);
			  // atimes
			  flag_mat = 0;
			  // cg
			  bench_SC (i, phi,
				    ewald_tr,
				    ewald_eps,
				    flag_prob,
				    flag_lub,
				    flag_mat,
				    "cg",
				    version,
				    flag_periodic);
			  // cgs
			  bench_SC (i, phi,
				    ewald_tr,
				    ewald_eps,
				    flag_prob,
				    flag_lub,
				    flag_mat,
				    "cgs",
				    version,
				    flag_periodic);
			  // bicgstab
			  bench_SC (i, phi,
				    ewald_tr,
				    ewald_eps,
				    flag_prob,
				    flag_lub,
				    flag_mat,
				    "bicgstab",
				    version,
				    flag_periodic);
			  // sta
			  bench_SC (i, phi,
				    ewald_tr,
				    ewald_eps,
				    flag_prob,
				    flag_lub,
				    flag_mat,
				    "sta",
				    version,
				    flag_periodic);
			  // sta2
			  bench_SC (i, phi,
				    ewald_tr,
				    ewald_eps,
				    flag_prob,
				    flag_lub,
				    flag_mat,
				    "sta2",
				    version,
				    flag_periodic);
			  // gpb
			  bench_SC (i, phi,
				    ewald_tr,
				    ewald_eps,
				    flag_prob,
				    flag_lub,
				    flag_mat,
				    "gpb",
				    version,
				    flag_periodic);
			  // otmk
			  bench_SC (i, phi,
				    ewald_tr,
				    ewald_eps,
				    flag_prob,
				    flag_lub,
				    flag_mat,
				    "otmk",
				    version,
				    flag_periodic);
			  // gmres
			  bench_SC (i, phi,
				    ewald_tr,
				    ewald_eps,
				    flag_prob,
				    flag_lub,
				    flag_mat,
				    "gmres",
				    version,
				    flag_periodic);
			}
		    }
		}
	    }
	}
    }

  return 0;
}
