/* benchmark code for stokes simulator in 3D for F/FT/FTS versions
 * Copyright (C) 1997-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bench3.c,v 1.5 2007/05/04 02:33:38 kichiki Exp $
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


/*
 * INPUT
 *  n_sc : number of particles in one-direction for SC lattice.
 *         so that the total number of particles is (n_sc)^3.
 *  phi  : volume fraction for the configuration
 */
void
bench_SC (int n_sc, double phi,
	  double ewald_tr,
	  double ewald_eps,
	  int flag_prob,
	  int flag_fix,
	  int flag_lub,
	  int flag_mat,
	  char * iter_solver,
	  int version)
{
  struct stokes *sys = stokes_init ();
  sys->version = version;

  int np, nm, nf;
  np = n_sc * n_sc * n_sc;
  if (flag_fix != 0)
    {
      nm = np - np/2;
      nf = np - nm;
    }
  else
    {
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

  sys->periodic = 1; // periodic boundary condition
  double lx, ly, lz;
  init_config_SC (phi, n_sc, n_sc, n_sc,
		  sys->pos, &lx, &ly, &lz);
  stokes_set_l (sys, lx, ly, lz);

  double xi = xi_by_tratio (sys, ewald_tr);
  stokes_set_xi (sys, xi, ewald_eps);

  //sys->lubmin = 2.0000000001;
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

  if (flag_prob == 0) // mobility problem
    {
      switch (version)
	{
	case 0: // F version
	  f  = (double *)calloc (nm3, sizeof (double));
	  u  = (double *)calloc (nm3, sizeof (double));
	  uf = (double *)calloc (nf3, sizeof (double));
	  ff = (double *)calloc (nf3, sizeof (double));
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
	  f  = (double *)calloc (nm3, sizeof (double));
	  t  = (double *)calloc (nm3, sizeof (double));
	  u  = (double *)calloc (nm3, sizeof (double));
	  o  = (double *)calloc (nm3, sizeof (double));
	  uf = (double *)calloc (nf3, sizeof (double));
	  of = (double *)calloc (nf3, sizeof (double));
	  ff = (double *)calloc (nf3, sizeof (double));
	  tf = (double *)calloc (nf3, sizeof (double));
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
	  f  = (double *)calloc (nm3, sizeof (double));
	  t  = (double *)calloc (nm3, sizeof (double));
	  s  = (double *)calloc (nm5, sizeof (double));
	  u  = (double *)calloc (nm3, sizeof (double));
	  o  = (double *)calloc (nm3, sizeof (double));
	  e  = (double *)calloc (nm5, sizeof (double));
	  uf = (double *)calloc (nf3, sizeof (double));
	  of = (double *)calloc (nf3, sizeof (double));
	  ef = (double *)calloc (nf5, sizeof (double));
	  ff = (double *)calloc (nf3, sizeof (double));
	  tf = (double *)calloc (nf3, sizeof (double));
	  sf = (double *)calloc (nf5, sizeof (double));
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
	  f  = (double *)calloc (np3, sizeof (double));
	  u  = (double *)calloc (np3, sizeof (double));
	  for (i = 0; i < np3; i ++)
	    {
	      u[i] = 1.0;
	    }
	  break;
	case 1: // FT version
	  f  = (double *)calloc (np3, sizeof (double));
	  t  = (double *)calloc (np3, sizeof (double));
	  u  = (double *)calloc (np3, sizeof (double));
	  o  = (double *)calloc (np3, sizeof (double));
	  for (i = 0; i < np3; i ++)
	    {
	      u[i] = 1.0;
	      o[i] = 1.0;
	    }
	  break;
	case 2: // FTS version
	  f  = (double *)calloc (np3, sizeof (double));
	  t  = (double *)calloc (np3, sizeof (double));
	  s  = (double *)calloc (np5, sizeof (double));
	  u  = (double *)calloc (np3, sizeof (double));
	  o  = (double *)calloc (np3, sizeof (double));
	  e  = (double *)calloc (np5, sizeof (double));
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

  double t0 = 0.0;
  double t1 = 0.0;
  /* call atimes routine */
  if (flag_prob == 0) // mobility problem
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

	  fprintf (stdout, "%d %d %f %e %e\n",
		   np, nm, t1 - t0, av_u, av_ff);

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

	  fprintf (stdout, "%d %f %e %e %e %e\n",
		   np, t1 - t0, av_u, av_ff, av_o, av_tf);

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

	  fprintf (stdout, "%d %f %e %e %e %e %e %e\n",
		   np, t1 - t0, av_u, av_ff, av_o, av_tf, av_e, av_ef);

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

	  fprintf (stdout, "%d %f %e\n",
		   np, t1 - t0, av_f);

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

	  fprintf (stdout, "%d %f %e %e\n",
		   np, t1 - t0, av_f, av_t);

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

	  fprintf (stdout, "%d %f %e %e %e\n",
		   np, t1 - t0, av_f, av_t, av_s);

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
  fprintf (stderr, "Benchmark Test for libstokes\n");
  fprintf (stderr, "$Id: bench3.c,v 1.5 2007/05/04 02:33:38 kichiki Exp $\n\n");
  fprintf (stderr, "USAGE\n");
  fprintf (stderr, "%s [options]\n", argv0);
  fprintf (stderr, "OPTIONS\n");
  fprintf (stderr, "\t-f, -ft, -fts : specify version (default F)\n");
  fprintf (stderr, "\t-res          : calc res problem"
	   " (default: mob problem)\n");
  fprintf (stderr, "\t-fix          : make half of particles fixed"
	   " (default: no fixed)\n"
	   "\t\tNOTE this is ignored for res problem.\n");
  fprintf (stderr, "\t-lub          : calc with lub"
	   " (default: no lub)\n");
  fprintf (stderr, "\t-mat          : use matrix scheme"
	   " (default: atimes)\n");
  fprintf (stderr, "\t-tr           : ewald_tr (default: 10.0)\n");
  fprintf (stderr, "\t-eps          : cut-off limit for ewald-sum"
	   " (default: 1.0e-12)\n");
  fprintf (stderr, "\t-iter         : iterative scheme. ignored for -mat.\n");
  fprintf (stderr, "\t\t0 = cg\n");
  fprintf (stderr, "\t\t1 = cgs\n");
  fprintf (stderr, "\t\t2 = bicgstab\n");
  fprintf (stderr, "\t\t3 = sta (another bicgstab)\n");
  fprintf (stderr, "\t\t4 = sta2\n");
  fprintf (stderr, "\t\t5 = gpb\n");
  fprintf (stderr, "\t\t6 = otmk\n");
  fprintf (stderr, "\t\t7 = gmres (default)\n");
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
  int flag_fix  = 0;
  int flag_iter = 7; // gmres
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
      else if (strcmp (argv [i], "-res") == 0)
	{
	  flag_prob = 1;
	}
      else if (strcmp (argv [i], "-fix") == 0)
	{
	  flag_fix = 1;
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

  /* write header */
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

  /* main loop */
  double phi = 0.3;
  for (i = 2; i < 100; i ++)
    {
      bench_SC (i, phi,
		ewald_tr,
		ewald_eps,
		flag_prob,
		flag_fix,
		flag_lub,
		flag_mat,
		iter_solver,
		version);
    }

  return 0;
}
