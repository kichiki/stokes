/* octave wrapper of calc_res_3f()
 * Copyright (C) 2006-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: stokes_res_3f.cc,v 1.2 2007/04/20 02:18:26 kichiki Exp $
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

#include <octave/oct.h>

extern "C" {
#include <stdlib.h>
#include <libiter.h>
#include <libstokes.h>
}


DEFUN_DLD (stokes_res_3f, args, ,
  "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {} stokes_res_3f (@var{pos}, @var{u}, @var{l}, @var{tz})\n\
@cindex Stokes resistance solver of F version in 3D with periodic boundary\n\
This function is a solver of Stokesian dynamics for resistance problem\n\
in F version in libstokes (http://ryuon.sourceforge.net/libstokes.html).\n\
Position @var{pos} and velocity @var{u} of particles are given and force\n\
@var{f} is returned. Lattice vector of the periodic boundary is given by\n\
@var{l}. An optional parameter @var{tz} is for the splitting of real and\n\
reciprocal spaces in Ewald summation.\n\
\n\
@example\n\
l   = [10, 10, 10];\n\
pos = [0, 0, 0, 5, 5, 5];\n\
u   = [1, 1, 1, 1, 1, 1];\n\
f = stokes_res_3f (pos, u, l)\n\
@end example\n\
@end deftypefn")
{
  octave_value retval;
  int nargin = args.length();

  if (nargin < 3 || nargin > 4)
    {
      print_usage ("stokes");
      return retval;
    }

  // pos
  ColumnVector pos (args(0).vector_value ());
  // u
  ColumnVector u   (args(1).vector_value ());
  // l
  ColumnVector l   (args(2).vector_value ());

  // z
  double ewald_tr;
  if (nargin == 4)
    {
      octave_value tr_arg = args(3);
      if (tr_arg.is_scalar_type())
	{
	  ewald_tr = tr_arg.double_value ();
	}
      else
	{
	  print_usage ("stokes_res_3f");
	  return retval;
	}
    }
  else
    {
      // default value
      ewald_tr = 10.0;
    }

  int n = pos.length();
  if (n != u.length())
    {
      error ("stokes_res_3f: wrong lengths for pos and u");
      return retval;
    }
  if (n % 3 != 0)
    {
      error ("stokes_res_3f: wrong length of pos");
      return retval;
    }
  int np = n / 3;

  if (l.length() != 3)
    {
      error ("stokes_res_3f: wrong lengths of l");
      return retval;
    }

  struct stokes *sys;
  sys = stokes_init ();
  stokes_set_np (sys, np, np);

  sys->periodic = 1; // periodic boundary condition
  double lx = l(0);
  double ly = l(1);
  double lz = l(2);
  stokes_set_l (sys, lx, ly, lz);

  double xi = xi_by_tratio (sys, ewald_tr);
  double ewald_eps = 1.0e-12;
  stokes_set_xi (sys, xi, ewald_eps);

  sys->lubmin = 2.0000000001;
  sys->lubmax = 4.0;
  stokes_set_iter (sys, "gmres", 2000, 20, 1.0e-6, 1, stdout);

  int i;
  double *d_pos = NULL;
  double *d_u   = NULL;
  double *d_f   = NULL;
  d_pos = (double *)calloc (n, sizeof (double));
  d_u   = (double *)calloc (n, sizeof (double));
  d_f   = (double *)calloc (n, sizeof (double));
  for (i = 0; i < n; i ++)
    {
      d_pos[i] = pos(i);
      d_u[i]   = u(i);
    }

  stokes_set_pos (sys, d_pos);
  solve_res_3f (sys, d_u, d_f);

  ColumnVector f (n);
  for (i = 0; i < n; i ++)
    {
      f(i) = d_f[i];
    }

  free (d_pos);
  sys->pos = NULL;

  free (d_u);
  free (d_f);

  stokes_free (sys);

  return octave_value (f);
}
