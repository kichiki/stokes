/* SWIG interface for libstokes
 * Copyright (C) 2006-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: stokes.i,v 1.3 2007/11/30 06:36:09 kichiki Exp $
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
%module stokes
%{
#include <stdio.h>

#include <libiter.h>
#include <libstokes.h>
%}

%inline%{
  FILE * get_stdin (void) {
    return (stdin);
  }
  FILE * get_stdout (void) {
    return (stdout);
  }
  FILE * get_stderr (void) {
    return (stderr);
  }
%}

%include "carrays.i"
%array_class(double, darray);
%array_class(int, iarray);

%include <libiter.h>
%include <libstokes.h>
