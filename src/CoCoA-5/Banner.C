//   Copyright (c) 2012 Anna M. Bigatti
//
//   This file is part of the source of CoCoALib, the CoCoA Library.
//
//   CoCoALib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   CoCoALib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with CoCoALib.  If not, see <http://www.gnu.org/licenses/>.

#include "Banner.H"

using std::string;

namespace CoCoA
{

  std::string CoCoA5Banner()
  {
    return 
      "   ______      ______      ___         ______\n"
      "  / ____/___  / ____/___  /   |       / ____/\n"
      " / /   / __ `/ /   / __ `/ /| |______/___ `  \n"
      "/ /___/ /_/ / /___/ /_/ / ___ /_____/___/ /  \n"
      "`____/`____/`____/`____/_/  |_|    /_____/   ";
  }

  std::string CoCoA5BannerNonFixedWidthFonts()
  {
    return
      "\n"
      "/******      C o C o A - 5      ******/";
  }

} // end of namespace CoCoA

