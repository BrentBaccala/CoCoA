//   Copyright (c)  2005,2007,2012  John Abbott

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


#include "CoCoA/BuildInfo.H"
#include "CoCoA/config.H"

#include <limits>
#include <iostream>
using std::ostream;
using std::endl;
#include <string>
using std::string;

namespace CoCoA
{

  namespace BuildInfo
  {
    // Build info is actually passed in via three preprocessor variables:
    // COCOA_VERSION, COCOA_CXX, COCOA_CXXFLAGS -- all three are quoted strings.

    const std::string& version()
    {
      static const string info(COCOA_VERSION);
      return info;
    }

    const std::string& compiler()
    {
      static const string info(COCOA_CXX);
      return info;
    }

    const std::string& CompilationFlags()
    {
      static const string info(COCOA_CXXFLAGS);
      return info;
    }

    const std::string& CompilationDefines()
    {
      static const string info(COCOA_CXXFLAGS_DEFINES);
      return info;
    }


    void PrintAll(std::ostream& out)
    {
    // Below we use string literal juxtaposition, so that the build
    // info can readily be extracted from libcocoa.a by the command
    // strings libcocoa.a | egrep "CoCoA::BuildInfo"
      out << endl
          << "CoCoA::BuildInfo Summary of build information for CoCoALib:" << endl
          << "CoCoA::BuildInfo CoCoALib Version: " COCOA_VERSION << endl
          << "CoCoA::BuildInfo Compiler: " COCOA_CXX << endl
          << "CoCoA::BuildInfo Compilation Flags: " COCOA_CXXFLAGS << endl
          << "CoCoA::BuildInfo Compilation Flags Defines: " COCOA_CXXFLAGS_DEFINES << endl
          << "CoCoA::BuildInfo NumBits int = " << std::numeric_limits<unsigned int>::digits << endl
          << "CoCoA::BuildInfo NumBits long = " << std::numeric_limits<unsigned long>::digits << endl
          << "CoCoA::BuildInfo NumBits SmallExponent_t = " << std::numeric_limits<SmallExponent_t>::digits << endl
          << endl;
    }

  } // end of namespace CoCoA::BuildInfo

} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/BuildInfo.C,v 1.8 2012/10/24 12:10:45 abbott Exp $
// $Log: BuildInfo.C,v $
// Revision 1.8  2012/10/24 12:10:45  abbott
// Revised PrintAll: now prints also NumBits of SmallExponent_t
//
// Revision 1.7  2012/10/02 11:45:09  bigatti
// -- added CompilationDefines
//
// Revision 1.6  2012/10/02 10:36:12  abbott
// Revised interface to BuildInfo information strings.
// Several consequential changes.
//
// Revision 1.5  2012/06/20 12:26:26  bigatti
// -- now exporting also compiler and compilation flags
//
// Revision 1.4  2011/03/11 17:38:21  bigatti
// -- added missing include
//
// Revision 1.3  2011/03/11 17:18:55  bigatti
// -- added info on int and long size
//
// Revision 1.2  2007/10/30 17:14:08  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.4  2007/02/10 18:44:03  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.3  2007/02/08 22:34:22  cocoa
// Changed BuildInfo: only BuildInfo version is publicly visible.
// Only BuildInfo needs complicated compile-time flags, so several changes
// to Makefiles etc.  Added a new example: ex-BuildInfo.C
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.2  2005/10/12 15:52:09  cocoa
// Completed test-RingFp1 and corrected/cleaned the SmallFp*
// and RingFp* files.
//
// Some minor tidying elsewhere.
//
// Revision 1.1  2005/10/06 16:36:42  cocoa
// Added the capability find out build information at run-time.
// The Makefiles should be a little tidier too.
//
