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

#ifndef DMPZinteger_content_h
#define DMPZinteger_content_h

#include "DMPZ.h"

/* Compute integer content of a DMPZ polynomial and put result in content. */
/* The algorithm is not terribly clever.                                   */
/* If f = 0 result is zero, otherwise it is positive.                      */

void DMPZinteger_content(mpz_t content, DMPZ f);

#endif
