/* test code for libstokes
 * Copyright (C) 2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: test-stokes.c,v 1.3 2006/10/19 18:48:53 ichiki Exp $
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
#include <stdio.h>
#include <stdlib.h>

#include <libstokes.h>
#include <libiter.h>


/* main program */
int
main (int argc, char** argv)
{
  struct stokes * sys = NULL;

  sys = stokes_init ();

  int np, nm;
  np = 8;
  nm = 8;
  stokes_set_np (sys, np, nm);

  double lx, ly, lz;
  lx = 10.0;
  ly = 10.0;
  lz = 10.0;
  stokes_set_l (sys, lx, ly, lz);

  double tratio, cutlim, xi;
  tratio = 60.25;
  xi = xi_by_tratio (sys, tratio);

  cutlim = 1.0e-12;
  stokes_set_xi (sys, xi, cutlim);

  fprintf (stdout, "xi = %f\n", xi);

  sys->lubcut = 2.0000000001;
  stokes_set_iter (sys, "gmres", 2000, 20, 1.0e-6, 1, stdout);

  int i;
  double * pos;
  double * u;
  double * f;
  pos = (double *) calloc (np * 3, sizeof (double));
  u   = (double *) calloc (np * 3, sizeof (double));
  f   = (double *) calloc (np * 3, sizeof (double));

  pos[ 0] = 0.0; // x component
  pos[ 1] = 0.0; // y component
  pos[ 2] = 0.0; // z component
  pos[ 3] = 5.0; pos[ 4] = 0.0; pos[ 5] = 0.0;
  pos[ 6] = 0.0; pos[ 7] = 5.0; pos[ 8] = 0.0;
  pos[ 9] = 0.0; pos[10] = 0.0; pos[11] = 5.0;
  pos[12] = 5.0; pos[13] = 5.0; pos[14] = 0.0;
  pos[15] = 0.0; pos[16] = 5.0; pos[17] = 5.0;
  pos[18] = 5.0; pos[19] = 0.0; pos[20] = 5.0;
  pos[21] = 5.0; pos[22] = 5.0; pos[23] = 5.0;

  for (i = 0; i < np * 3; i ++)
    {
      u[i] = 1.0;
    }

  fprintf (stdout, "pos:\n");
  for (i = 0; i < np; i ++)
    {
      fprintf (stdout, "%d %f %f %f\n",
	       i,
	       pos[i*3 + 0],
	       pos[i*3 + 1],
	       pos[i*3 + 2]);
    }

  fprintf (stdout, "u:\n");
  for (i = 0; i < np; i ++)
    {
      fprintf (stdout, "%d %f %f %f\n",
	       i,
	       u[i*3 + 0],
	       u[i*3 + 1],
	       u[i*3 + 2]);
    }

  stokes_set_pos (sys, pos);

  solve_res_ewald_3f (sys, u, f);

  fprintf (stdout, "f:\n");
  for (i = 0; i < np; i ++)
    {
      fprintf (stdout, "%d %f %f %f\n",
	       i,
	       f[i*3 + 0],
	       f[i*3 + 1],
	       f[i*3 + 2]);
    }

  return 0;
}
