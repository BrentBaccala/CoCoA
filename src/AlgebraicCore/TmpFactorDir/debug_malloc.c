//   Copyright (c)  1997-2006  John Abbott

//   This file is part of the source of CoCoALib, the CoCoA Library.

//   CoCoALib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.

//   CoCoALib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.

//   You should have received a copy of the GNU General Public License
//   along with CoCoALib.  If not, see <http://www.gnu.org/licenses/>.

#include <stdlib.h>
#include <stdio.h>

/* debugging allocation functions for GMP */

void *debug_malloc(size_t sz)
{
  void *ans;

  ans = malloc(sz);
  printf("ALLOC %lx (%ld bytes)\n", (long)ans, (long)sz);
  return ans;
}

void debug_free(void *ptr, size_t sz)
{
  printf("FREEING %lx (%ld bytes)\n", (long)ptr, (long)sz);
  free(ptr);
}

void *debug_realloc(void *ptr, size_t old_sz, size_t new_sz)
{
  void *ans;
  ans = realloc(ptr, new_sz);
  if (ans == ptr) return ans;
  printf("FREEING %lx (%ld bytes) realloc\n", (long)ptr, (long)old_sz);
  printf("ALLOC %lx (%ld bytes) realloc\n", (long)ans, (long)new_sz);
  return ans;
}
