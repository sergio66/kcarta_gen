/*** see http://www.koders.com/c/

/* Copyright (C) 1997 DJ Delorie, see COPYING.DJ for details */
/* Copyright (C) 1996 DJ Delorie, see COPYING.DJ for details */
/* Copyright (C) 1995 DJ Delorie, see COPYING.DJ for details */

#include <stdio.h>
#include <errno.h>
#include <fcntl.h>
#include "/home/sergio/KCARTA/UTILITY/FSEEK/Src/file.h"
#include "/home/sergio/KCARTA/UTILITY/FSEEK/Src/io.h"

int fseek_local(FILE *f, long offset, int ptrname)
{
  long p = -1;			/* can't happen? */
  if ( f == NULL ) {
  __set_errno (EINVAL);
       return -1;
  }
  
  f->_flag &= ~_IOEOF;
  if (!OPEN4WRITING(f))
  {
    if (f->_base && !(f->_flag & _IONBF))
    {
      p = ftell(f);
      if (ptrname == SEEK_CUR)
      {
      offset += p;
      ptrname = SEEK_SET;
      }
      /* check if the target position is in the buffer and
        optimize seek by moving inside the buffer */
      if (ptrname == SEEK_SET && (f->_flag & (_IOUNGETC|_IOREAD|_IOWRT )) == 0
      && p-offset <= f->_ptr-f->_base && offset-p <= f->_cnt)
      {
        f->_ptr+=offset-p;
        f->_cnt+=p-offset;
        return 0;
      }
    }

 
    p = lseek(fileno(f), offset, ptrname);
    f->_cnt = 0;
    f->_ptr = f->_base;
    f->_flag &= ~_IOUNGETC;
  }
  else 
  {
    p = fflush(f);
    return lseek(fileno(f), offset, ptrname) == -1 || p == EOF ?
      -1 : 0;
  }
  return p==-1 ? -1 : 0;

}

/* An ftell() function that works around platform bugs.
   Copyright (C) 2007 Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along
   with this program; if not, write to the Free Software Foundation,
   Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.  */

#include "/home/sergio/KCARTA/UTILITY/FSEEK/Src/config.h"

/* Specification.  */
#include <stdio.h>

#include <errno.h>
/* Get off_t.  */
#include <unistd.h>

long
ftell (FILE *fp)
{
  /* Use the replacement ftello function with all its workarounds.  */
  off_t offset = ftello (fp);
  if (offset == (long)offset)
    return (long)offset;
  else
    {
      errno = EOVERFLOW;
      return -1;
    }
}

