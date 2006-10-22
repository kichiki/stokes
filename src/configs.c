/* utility routines for Ewald-summation code in 3D
 * Copyright (C) 2001-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: configs.c,v 1.1 2006/10/22 03:53:41 kichiki Exp $
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
#include <stdlib.h> // srand48(), drand48()

#include "configs.h"


/* check overlap
 * INPUT
 * OUTPUT (return value)
 *  0 : no-overlap
 *  1 : overlap
 */
static int
check_overlap (int np, double lx, double ly, double lz, double * pos)
{
  int i;
  int ix, iy, iz;
  int j;
  int jx, jy, jz;
  int kx, ky, kz;

  double rr;
  double x, y, z;


  for (i = 0; i < np; ++i)
    {
      ix = i * 3;
      iy = ix + 1;
      iz = ix + 2;
      for (j = i + 1; j < np; ++j)
	{
	  jx = j * 3;
	  jy = jx + 1;
	  jz = jx + 2;
	  for (kx = -1; kx <= 1; ++kx)
	    {
	      for (ky = -1; ky <= 1; ++ky)
		{
		  for (kz = -1; kz <= 1; ++kz)
		    {
		      x = pos [ix] - pos [jx] + (double) kx * lx;
		      y = pos [iy] - pos [jy] + (double) ky * ly;
		      z = pos [iz] - pos [jz] + (double) kz * lz;
		      rr = x * x
			+ y * y
			+ z * z;
		      if (rr <= 4.0)
			{
			  return 1;
			}
		    }
		}
	    }
	}
    }
  return 0;
}

/* initialize configuration and the primary cell of simple cubic lattic
 * INPUT
 *  phi        : volume fraction of particles
 *  nx, ny, nz : # particles in each direction
 * OUTPUT
 *  pos [(nx * ny * nz) * 3] : positions of particles
 *  lx, ly, lz : geometry of the primary cell
 */
void
init_config_SC (double phi, int nx, int ny, int nz,
		double *pos, double *lx, double *ly, double *lz)
{
  int j;
  int ix, iy, iz;


  (*lx) = (*ly) = (*lz) = pow (4.0 * M_PI / phi / 3.0, 1.0 / 3.0);

  /* extend the primary cell with 1 particle to that with "n" particles*/
  j = 0;
  for (ix = 0; ix < nx; ix ++)
    {
      for (iy = 0; iy < ny; iy ++)
	{
	  for (iz = 0; iz < nz; iz ++)
	    {
	      pos [j * 3 + 0] = (double) ix * (*lx);
	      pos [j * 3 + 1] = (double) iy * (*ly);
	      pos [j * 3 + 2] = (double) iz * (*lz);
	      j ++;
	    }
	}
    }
  (*lx) *= (double) nx;
  (*ly) *= (double) ny;
  (*lz) *= (double) nz;
}

/* initialize configuration and the primary cell of BCC
 * INPUT
 *  phi        : volume fraction of particles
 *  nx, ny, nz : # single cell in each direction
 *             : so that total # particles is 2*nx*ny*nz
 * OUTPUT
 *  pos [(nx * ny * nz) 2 * * 3] : positions of particles
 *  lx, ly, lz : geometry of the primary cell
 */
void
init_config_BCC (double phi, int nx, int ny, int nz,
		 double *pos, double *lx, double *ly, double *lz)
{
  int j;
  int ix, iy, iz;


  (*lx) = (*ly) = (*lz) = pow (8.0 * M_PI / phi / 3.0, 1.0 / 3.0);

  /* extend the primary cell with 1 particle to that with "n" particles*/
  j = 0;
  for (ix = 0; ix < nx; ix ++)
    {
      for (iy = 0; iy < ny; iy ++)
	{
	  for (iz = 0; iz < nz; iz ++)
	    {
	      pos [j * 3 + 0] = (double) ix * (*lx);
	      pos [j * 3 + 1] = (double) iy * (*ly);
	      pos [j * 3 + 2] = (double) iz * (*lz);
	      j ++;
	      pos [j * 3 + 0] = ((double) ix + 0.5) * (*lx);
	      pos [j * 3 + 1] = ((double) iy + 0.5) * (*ly);
	      pos [j * 3 + 2] = ((double) iz + 0.5) * (*lz);
	      j ++;
	    }
	}
    }
  (*lx) *= (double) nx;
  (*ly) *= (double) ny;
  (*lz) *= (double) nz;
}

/* initialize configuration and the primary cell of FCC
 * INPUT
 *  phi        : volume fraction of particles
 *  nx, ny, nz : # single cell in each direction
 *             : so that total # particles is 4*nx*ny*nz
 * OUTPUT
 *  pos [(nx * ny * nz) 4 * * 3] : positions of particles
 *  lx, ly, lz : geometry of the primary cell
 */
void
init_config_FCC (double phi, int nx, int ny, int nz,
		 double *pos, double *lx, double *ly, double *lz)
{
  int j;
  int ix, iy, iz;


  (*lx) = (*ly) = (*lz) = pow (16.0 * M_PI / phi / 3.0, 1.0 / 3.0);

  /* extend the primary cell with 1 particle to that with "n" particles*/
  j = 0;
  for (ix = 0; ix < nx; ix ++)
    {
      for (iy = 0; iy < ny; iy ++)
	{
	  for (iz = 0; iz < nz; iz ++)
	    {
	      pos [j * 3 + 0] = (double) ix * (*lx);
	      pos [j * 3 + 1] = (double) iy * (*ly);
	      pos [j * 3 + 2] = (double) iz * (*lz);
	      j ++;
	      pos [j * 3 + 0] = ((double) ix + 0.0) * (*lx);
	      pos [j * 3 + 1] = ((double) iy + 0.5) * (*ly);
	      pos [j * 3 + 2] = ((double) iz + 0.5) * (*lz);
	      j ++;
	      pos [j * 3 + 0] = ((double) ix + 0.5) * (*lx);
	      pos [j * 3 + 1] = ((double) iy + 0.0) * (*ly);
	      pos [j * 3 + 2] = ((double) iz + 0.5) * (*lz);
	      j ++;
	      pos [j * 3 + 0] = ((double) ix + 0.5) * (*lx);
	      pos [j * 3 + 1] = ((double) iy + 0.5) * (*ly);
	      pos [j * 3 + 2] = ((double) iz + 0.0) * (*lz);
	      j ++;
	    }
	}
    }
  (*lx) *= (double) nx;
  (*ly) *= (double) ny;
  (*lz) *= (double) nz;
}

/* initialize random configuration and the primary cell
 * INPUT
 *  seed : seed for random
 *  phi : volume fraction of particles
 *  np  : # particles
 * OUTPUT
 *  pos [np * 3] : positions of particles
 *  lx, ly, lz   : geometry of the primary cell, where lz = 2.0 is fixed.
 */
void
init_config_random (long seed,
		    double phi, int np,
		    double *pos, double *lx, double *ly, double *lz)
{
  int i;
  int ix, iy, iz;
  int i_tor = 1000;
  int j;


  srand48 (seed);

  (*lx) = (*ly) = (*lz)
    = pow (4.0 * M_PI * (double) np / phi / 3.0,
	   1.0 / 3.0);

  j = 0;
  for (i = 0; i < np; i ++)
    {
    retry_init_config_random:
      ix = i * 3;
      iy = ix + 1;
      iz = ix + 2;

      pos [ix] = (*lx) * drand48 ();
      pos [iy] = (*ly) * drand48 ();
      pos [iz] = (*lz) * drand48 ();

      if (check_overlap (i + 1, (*lx), (*ly), (*lz), pos) != 0)
	{
	  ++j;
	  if (j > i_tor)
	    {
	      j = 0;
	      if (i > 0)
		--i;
	    }
	  goto retry_init_config_random;
	}
    }
}

