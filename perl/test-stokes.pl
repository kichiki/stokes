# test code for libstokes
# Copyright (C) 2006-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
# $Id: test-stokes.pl,v 1.7 2007/08/12 20:03:25 kichiki Exp $
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

use stokes;

$sys = stokes::stokes_init();

$np = 8;
$nm = 8;
stokes::stokes_set_np($sys, $np, $nm);
# you must call stokes_set_np() because
# this also allocate the memory for pos.

$sys->{periodic} = 1; # periodic boundary condition
$lx = 10.0;
$ly = 10.0;
$lz = 10.0;
stokes::stokes_set_l($sys, $lx, $ly, $lz);

$ewald_tr = 60.25;
$xi = stokes::xi_by_tratio($sys, $ewald_tr);
$ewald_eps = 1.0e-12;
stokes::stokes_set_xi($sys, $xi, $ewald_eps);

print "xi = ", $xi, "\n";

$sys->{lubmin2} = 4.0000000001;
$sys->{lubmax}  = 4.0;
stokes::stokes_set_iter($sys, "gmres", 2000, 20, 1.0e-6,
			1, stokes::get_stdout());

$pos = new stokes::darray($np*3);
$u   = new stokes::darray($np*3);
$f   = new stokes::darray($np*3);

$pos->setitem( 0, 0.0); # x component
$pos->setitem( 1, 0.0); # y component
$pos->setitem( 2, 0.0); # z component
$pos->setitem( 3, 5.0); $pos->setitem( 4, 0.0); $pos->setitem( 5, 0.0);
$pos->setitem( 6, 0.0); $pos->setitem( 7, 5.0); $pos->setitem( 8, 0.0);
$pos->setitem( 9, 0.0); $pos->setitem(10, 0.0); $pos->setitem(11, 5.0);
$pos->setitem(12, 5.0); $pos->setitem(13, 5.0); $pos->setitem(14, 0.0);
$pos->setitem(15, 0.0); $pos->setitem(16, 5.0); $pos->setitem(17, 5.0);
$pos->setitem(18, 5.0); $pos->setitem(19, 0.0); $pos->setitem(20, 5.0);
$pos->setitem(21, 5.0); $pos->setitem(22, 5.0); $pos->setitem(23, 5.0);

for ($i=0; $i < $np*3; $i++){
    $u->setitem($i, 1.0);
}

print "pos:\n";
for ($i=0; $i < $np; $i++){
    printf ("%d %f %f %f\n",
	    $i,
	    $pos->getitem($i*3),
	    $pos->getitem($i*3+1),
	    $pos->getitem($i*3+2));
}
print "u:\n";
for ($i=0; $i < $np; $i++){
    printf ("%d %f %f %f\n",
	    $i,
	    $u->getitem($i*3),
	    $u->getitem($i*3+1),
	    $u->getitem($i*3+2));
}

#$sys->{pos} = $pos;
stokes::stokes_set_pos($sys, $pos);
stokes::solve_res_3f($sys, $u, $f);

$nc_f = stokes::stokes_nc_init("test-stokes.res-3f.nc",
			       $np,
			       0, # nf
			       0, # version
			       0, # flag_poly
			       0, # flag_Q
			       0);# flag_it (time-dependent imposed flow)
# f0, x, u are active
stokes::stokes_nc_set_f0($nc_f, $f);
stokes::stokes_nc_set_time($nc_f, 0, 0.0);
stokes::stokes_nc_set_x($nc_f, 0, $pos);
stokes::stokes_nc_set_u($nc_f, 0, $u);

stokes::stokes_nc_free($nc_f);


print "f:\n";
for ($i=0; $i < $np; $i++){
    printf ("%d %f %f %f\n",
	    $i,
	    $f->getitem($i*3),
	    $f->getitem($i*3+1),
	    $f->getitem($i*3+2));
}
