#ifndef CoCoA_OnlineHelp_H
#define CoCoA_OnlineHelp_H

//   Copyright (c)  2010  Anna Bigatti

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

#include <string>
//using std::string;
#include <vector>
//using std::vector;


namespace CoCoA
{
namespace OnlineHelp
{

  void PrintMan(std::string keyword, std::ostream &out);
  void ReloadMan(std::ostream &out);
  void PrintAllExamples(std::ostream &out);
  void PrintAllExamplesWithoutOutput(std::ostream &out);
  void PrintWordlist(std::ostream &out);
  const char* XMLFileName();

} // namespace OnlineHelp
} // namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/CoCoA-5/OnlineHelp.H,v 1.11 2014/03/26 11:50:43 abbott Exp $
// $Log: OnlineHelp.H,v $
// Revision 1.11  2014/03/26 11:50:43  abbott
// Summary: Changed signature of PrintMan (1st arg was const string&, now string)
// Author: JAA
//
// Revision 1.10  2012/06/04 09:34:05  bigatti
// -- added PrintWordlist
//
// Revision 1.9  2012/04/04 13:56:35  bigatti
// -- added PrintAllExamplesWithoutOutput
//
// Revision 1.8  2012/02/24 13:10:12  bigatti
// -- added ReloadMan
//
// Revision 1.7  2011/02/16 16:14:08  bigatti
// -- class deefinition moved into .C file
// -- added storing of types, and function for <commands_and_functions_for>
// -- cleaning up
//
// Revision 1.6  2011/02/14 10:10:07  bigatti
// -- fixed &amp; &apos; <backslash/>  ...
// -- added function PrintAllExamples
//
// Revision 1.5  2010/09/01 13:24:48  bigatti
// -- moved all manual functions into CoCoA::OnlineHelp namespace
//
// Revision 1.4  2010/09/01 12:27:52  lagorio
// *** empty log message ***
//
// Revision 1.3  2010/09/01 07:46:19  lagorio
// *** empty log message ***
//
// Revision 1.2  2010/08/31 14:55:58  bigatti
// -- added cvs log
//

#endif
