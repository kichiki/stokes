/* some I/O utility routines
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: file.h,v 1.1 2007/05/15 07:19:28 kichiki Exp $
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
#ifndef	_FILE_H_
#define	_FILE_H_


/* return the last place of '/' in name
 * OUTPUT
 *  returned value : -1 when no '/' is presented
 */
int
get_dir (const char * name);

/* check wheterh the file exists or not (from man opendir())
 * INPUT
 * OUTPUT
 *  returned value : 0 NOT_FOUND (false)
 *                   1 FOUND     (true)
 */
int
check_file (const char * name);


#endif /* !_FILE_H_ */
