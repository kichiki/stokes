/* stokes-netcdf utility tool
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: stnc.c,v 1.2 2008/10/22 06:51:48 kichiki Exp $
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

#include "file.h" // check_file()

void
usage (const char *argv0)
{
  fprintf (stdout, "stokes-netcdf utility tool\n");
  fprintf (stdout, "$Id: stnc.c,v 1.2 2008/10/22 06:51:48 kichiki Exp $\n\n");
  fprintf (stdout, "USAGE\n");
  fprintf (stdout, "%s [OPTIONS]\n", argv0);
  fprintf (stdout, "\t-h or --help     : this message.\n");
  fprintf (stdout, "\t-i or --input    : input stokes-netcdf file\n");
  fprintf (stdout, "\t-o or --output   : output stokes-netcdf file\n");
  fprintf (stdout, "\t-i0 or --begin   : first step (default: 0)\n");
  fprintf (stdout, "\t-i1 or --end     : end step (default: the last step)\n");
}


/* main program */
int
main (int argc, char** argv)
{
  /**
   * option analysis
   */
  int i;
  char *in_file = NULL;
  char *out_file = NULL;
  int i0 = 0;
  int i1 = -1;
  for (i = 1; i < argc; i++)
    {
      if (strcmp (argv [i], "-i") == 0 ||
	  strcmp (argv [i], "--input") == 0)
	{
	  i++;
	  int len = strlen (argv [i]);
	  in_file = (char *)malloc (sizeof (char) * (len + 1));
	  strcpy (in_file, argv [i]);
	}
      else if (strcmp (argv [i], "-o") == 0 ||
	  strcmp (argv [i], "--output") == 0)
	{
	  i++;
	  int len = strlen (argv [i]);
	  out_file = (char *)malloc (sizeof (char) * (len + 1));
	  strcpy (out_file, argv [i]);
	}
      else if (strcmp (argv [i], "-i0") == 0 ||
	  strcmp (argv [i], "--begin") == 0)
	{
	  i0 = atoi (argv [++i]);
	}
      else if (strcmp (argv [i], "-i1") == 0 ||
	  strcmp (argv [i], "--end") == 0)
	{
	  i1 = atoi (argv [++i]);
	}
      else
	{
	  usage (argv[0]);
	  exit (1);
	}
    }
  
  if (in_file == NULL ||
      out_file == NULL)
    {
      fprintf (stderr, "you must give input and output files.\n");
      exit (1);
    }

  if (check_file (in_file) != 1) // not found
    {
      fprintf (stderr, "input file %s does not exist.\n", in_file);
      exit (1);
    }
  if (check_file (out_file) == 1) // found
    {
      fprintf (stderr, "output file %s already exist.\n", out_file);
      exit (1);
    }


  /**
   * initialize NetCDF and set the constant parameters
   */
  double t;
  //struct stokes_nc *nc_in = stokes_nc_reopen (in_file);
  struct stokes_nc *nc_in = stokes_nc_open (in_file);
  CHECK_MALLOC (nc_in, "main");
  free (in_file);


  // check
  fprintf (stdout, "active entries:\n");
  stokes_nc_print_actives (nc_in, stdout);
  fprintf (stdout, "\nit contains %d steps (from 0 to %d).\n\n",
	   nc_in->ntime, nc_in->ntime - 1);


  struct stokes *sys = stokes_init ();
  CHECK_MALLOC (sys, "main");
  int flag_q;
  double Ui[3];
  double Oi[3];
  double Ei[5];
  double lat[3];
  int shear_mode;
  double shear_rate;
  int flag_rng;

  int nm = nc_in->np;
  double *x = (double *)malloc (sizeof (double) * nm * 3);
  CHECK_MALLOC (x, "main");
  double *q = (double *)malloc (sizeof (double) * nm * 4);
  CHECK_MALLOC (q, "main");
  double *F = (double *)malloc (sizeof (double) * nm * 3);
  double *T = (double *)malloc (sizeof (double) * nm * 3);
  double *E = (double *)malloc (sizeof (double) * nm * 5);
  CHECK_MALLOC (F, "main");
  CHECK_MALLOC (T, "main");
  CHECK_MALLOC (E, "main");


  double *uf = NULL;
  double *of = NULL;
  double *ef = NULL;
  double *xf = NULL;
  int nf = nc_in->npf;
  if (nf > 0)
    {
      uf = (double *)malloc (sizeof (double) * nf * 3);
      of = (double *)malloc (sizeof (double) * nf * 3);
      ef = (double *)malloc (sizeof (double) * nf * 5);
      xf = (double *)malloc (sizeof (double) * nf * 3);
      CHECK_MALLOC (uf, "main");
      CHECK_MALLOC (of, "main");
      CHECK_MALLOC (ef, "main");
      CHECK_MALLOC (xf, "main");
    }

  stokes_nc_get_params (nc_in,
			sys,
			&flag_q,
			Ui, Oi, Ei,
			F, T, E,
			uf, of, ef,
			xf,
			lat,
			&shear_mode,
			&shear_rate,
			&flag_rng);

  struct stokes_nc *nc_out
    = stokes_nc_set_by_params (out_file,
			       sys,
			       flag_q,
			       Ui, Oi, Ei, F, T, E,
			       uf, of, ef, xf,
			       lat,
			       shear_mode, shear_rate,
			       flag_rng);
  CHECK_MALLOC (nc_out, "main");
  free (out_file);


  struct KIrand *rng = NULL;
  if (flag_rng != 0)
    {
      rng = KIrand_init ();
      CHECK_MALLOC (rng, "main");
    }


  /**
   * mail loop
   */
  if (i1 < 0)
    {
      i1 = nc_in->ntime - 1;
    }
  fprintf (stdout, "output steps from %d to %d.\n\n", i0, i1);


  int l;
  double shear_shift;
  size_t index = 0;
  for (l = i0; l <= i1; l++)
    {
      //fprintf (stdout, "%d steps\n", l);


      // get information at "l" step
      t  = stokes_nc_get_time_step (nc_in, l);

      // set the configuration at the current time
      stokes_nc_get_data (nc_in, "x", l, x);
      if (flag_q != 0)
	{
	  stokes_nc_get_data (nc_in, "q", l, q);
	}
      if (shear_mode != 0)
	{
	  int status = nc_get_var1_double
	    (nc_in->id, nc_in->shear_shift_id, &index, &shear_shift);
	  if (status != NC_NOERR)
	    {
	      fprintf (stderr,
		       "at nc_get_var1_double() for shear_shift in main\n");
	    }
	}
      if (flag_rng != 0) // Brownian
	{
	  stokes_nc_get_rng (nc_in, l, rng);
	}



      // output the results at t (should be t_out) with (l) step
      stokes_nc_set_time (nc_out, l, t);
      stokes_nc_set_x (nc_out, l, x);
      if (flag_q != 0)
	{
	  stokes_nc_set_q (nc_out, l, q);
	}
      if (shear_mode != 0)
	{
	  stokes_nc_set_shear_shift (nc_out, l, shear_shift);
	}
      if (flag_rng != 0) // Brownian
	{
	  stokes_nc_set_rng (nc_out, l, rng);
	}

      // flush the data
      nc_sync(nc_out->id);
    }


  // house keeping
  stokes_nc_free (nc_in);
  stokes_free (sys);
  free (x);
  free (q);
  free (F);
  free (T);
  free (E);
  if (uf != NULL) free (uf);
  if (of != NULL) free (of);
  if (ef != NULL) free (ef);
  if (xf != NULL) free (xf);
  stokes_nc_free (nc_out);
  if (rng != NULL) KIrand_free (rng);


  return 0;
}
