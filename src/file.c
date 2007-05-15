/* some I/O utility routines
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: file.c,v 1.1 2007/05/15 07:19:06 kichiki Exp $
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
#include <string.h> // strlen()
#include <sys/types.h>
#include <dirent.h>


/* return the last place of '/' in name
 * OUTPUT
 *  returned value : -1 when no '/' is presented
 */
int
get_dir (const char * name)
{
  int len = strlen (name);
  int i;
  for (i = 0; i < len; i ++)
    {
      if (name [len-1-i] == '/')
	{
	  /* name [0 ~ len-1-i] is the path including '/' */
	  /* name [len-i ~ len-1] is file */
	  return (len-1-i);
	}
    }
  return (-1);
}

/* check wheterh the file exists or not (from man opendir())
 * INPUT
 * OUTPUT
 *  returned value : 0 NOT_FOUND (false)
 *                   1 FOUND     (true)
 */
int
check_file (const char * name)
{
  int l = strlen (name);
  int i = get_dir (name);

  char *path;
  char *file;
  if (i < 0)
    {
      path = (char *)malloc (sizeof (char) * 2);
      path [0] = '.';
      path [1] = '\0';

      file = (char *)malloc (sizeof (char) * (l+1));
      strcpy (file, name);
    }
  else
    {
      /* -- like 'xxx/yyy/zzz' */
      path = (char *)malloc (sizeof (char) * (i+1));
      strncpy (path, name, i);
      path [i] = '\0';
      file = (char *)malloc (sizeof (char) * (l-i));
      strncpy (file, name + i + 1, l-(i+1));
      file [l-(i+1)] = '\0';
    }
  /*
  fprintf (stderr, "whole = %s\n", name);
  fprintf (stderr, "path = %s\n", path);
  fprintf (stderr, "file = %s\n", file);
  exit (0);*/
  int len = strlen (file);

  DIR *dirp = opendir (path);
  if (dirp == NULL)
    {
      fprintf (stderr, "cannot open directory %s\n", path);
      return 0;
    }
  free (path);

  struct dirent *dp = readdir (dirp);
  while (dp != NULL)
    {
      if (
#ifdef _D_EXACT_NAMLEN
	  /* for LINUX */
	  _D_EXACT_NAMLEN(dp) == len
#else
	  dp->d_namlen == len
#endif
	  && !strcmp (dp->d_name, file))
	{
	  closedir (dirp);
	  free (file);
	  return 1;/* FOUND */
	}
    }
  closedir (dirp);
  free (file);
  return 0; /* NOT_FOUND */
}

