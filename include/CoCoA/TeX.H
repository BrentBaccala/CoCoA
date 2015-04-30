#ifndef CoCoA_TeX_H
#define CoCoA_TeX_H

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


#include <iosfwd>


namespace CoCoA
{

  //---------------------------------------------------------------------------
  // Standard library's 'iword' feature is used to add an ios output flag

  bool TeX_mode(std::ios_base& s);

  typedef std::ostream& (* ios_manipulator) (std::ostream&);

  ios_manipulator TeX(bool mode = true);

} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/OpenMath.H,v 1.10 2012/05/04 16:59:43 abbott Exp $
// $Log: OpenMath.H,v $
// Revision 1.10  2012/05/04 16:59:43  abbott
// For future compatibility with auto conversion from bool to 3-way bool,
// removed sending of MachineInt; instead added operators for sending
// (unsigned) int and (unsigned) long [but no other integral types for now].
//
// Revision 1.9  2012/02/02 15:55:48  abbott
// Replaced include of MachineInt by fwd decl.
//
// Revision 1.8  2011/11/09 13:47:56  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.7  2011/08/14 15:52:18  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.6  2011/08/12 15:48:31  abbott
// Added virtual mem fns for BigInt type.
//
// Revision 1.5  2011/03/11 14:49:08  abbott
// Changed size_t into long.
//
// Revision 1.4  2010/03/18 13:53:48  abbott
// Minor rationalization to OpenMath implementation: moved op<<(OMOut,ZZ) to ZZ files.
//
// Revision 1.3  2008/12/16 21:10:32  abbott
// Replaced the various output fns for different sort of machine integers by a
// single one for MachineInt.
//
// Revision 1.2  2007/10/30 17:14:12  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.6  2006/12/06 17:14:12  cocoa
// -- removed #include "config.H"
//
// Revision 1.5  2006/11/27 13:04:35  cocoa
// Added explicit OpenMath output operators for all integer types.
//
// Revision 1.4  2006/11/23 17:10:52  cocoa
// -- changed: OpenMathOutput and OpenMathInput are now a class (instead of typedef)
//
// Revision 1.3  2006/11/02 13:25:44  cocoa
// Simplification of header files: the OpenMath classes have been renamed.
// Many minor consequential changes.
//
// Revision 1.2  2006/10/06 14:04:15  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//


#endif
