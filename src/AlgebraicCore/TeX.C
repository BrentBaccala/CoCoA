//   Copyright (c)  2015  Brent Baccala

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


#include "CoCoA/TeX.H"

#include <iostream>
using std::ostream;

namespace CoCoA
{

  //---------------------------------------------------------------------------
  // Standard library's 'iword' feature is used to add an ios output flag

  long& TeX_flag_ref(std::ios_base& s) {
    static int TeX_index = std::ios_base::xalloc();
    return s.iword(TeX_index);
  }

  bool TeX_mode(std::ios_base& s) {
    return TeX_flag_ref(s);
  }

  void set_TeX_flag(std::ios_base& s, long n) {
    TeX_flag_ref(s) = n;
  }

  std::ostream& TeX(std::ostream& os) {
    set_TeX_flag(os, true);
    return os;
  }

  std::ostream& TeXoff(std::ostream& os) {
    set_TeX_flag(os, false);
    return os;
  }

  ios_manipulator TeX(bool mode) {
    if (mode) return TeX;
    else return TeXoff;
  }


} // end of namespace CoCoA
