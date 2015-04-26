//   Copyright (c) 2009 Giovanni Lagorio (lagorio@disi.unige.it)
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

#include <iostream>
#include <boost/tr1/memory.hpp>
#include <cstring>
#include <cctype>
#include <boost/iostreams/filter/newline.hpp>
#include "LineProviders.H"

namespace CoCoA {
namespace LexerNS {

using namespace std;
using namespace boost;
using namespace boost::iostreams;
using namespace CoCoA;

ReferenceCountedObject::~ReferenceCountedObject() {}

FileLineProvider::FileLineProvider(const std::string &file)
{
  filename = file;
	input.push(newline_filter(newline::posix));
	const file_source fs(filename.c_str());
	if (!fs.is_open())
		throw IOException("Cannot open file \""+filename+"\" for reading.");
	input.push(fs);
}

  FileRegionLineProvider::FileRegionLineProvider(const std::string &file, long FromLine, long FromChar, long ToLine, long ToChar) :
      myFromLine(FromLine),
      myFromChar(FromChar),
      myToLine(ToLine),
      myToChar(ToChar),
      myCurrLine(0)
  {
    filename = file;
    if (FromLine <= 0 || FromChar <= 0 || ToLine <= 0 || ToChar <= 0)
      throw IOException("Line number and char position must be positive");
    if (FromLine > ToLine || (FromLine == ToLine && FromChar > ToChar))
      throw IOException("Impossible region");
    input.push(newline_filter(newline::posix));
    const file_source fs(filename.c_str());
    if (!fs.is_open())
      throw IOException("Cannot open file \""+filename+"\" for reading.");
    input.push(fs);
}

string InteractiveLineProvider::promptSuffix("#");

string InteractiveLineProvider::prompt(const LexerStatus &ls, const ParserNS::ParserStatus &ps) {
	if (ps.isInRecoveryMode())
		throw ParserNS::AskingForNewInteractiveInputDuringRecoveryException();
	string s = ps.promptHeader();
	if (ps.isInTheMiddleOfAStatement())
		s.append("... ; ");
	else if (ls.isInMultiLineComment())
		s.append("/* ... ");
	else if (ls.isInStringLiteral())
		s.append("\" ... ");
	s.append(promptSuffix);
	s.append(" ");
	return s;
}

bool GetlineLineProvider::doReadNextLine(const LexerStatus &ls, const ParserNS::ParserStatus &ps, string &chars) {
	if (cin.eof())
		return false;
// Line below does automatic prompt suppression if cin.syn_with_stdio(false) has been called (commented out first line of main() in Main.C)
//        if (cin.rdbuf()->in_avail() == 0) cout << prompt(ls, ps);
	cout << prompt(ls, ps); // always print prompt
	getline(cin, chars);
	return true;
}

bool FileLineProvider::doReadNextLine(const LexerStatus &, const ParserNS::ParserStatus &, string &chars) {
	if (input.eof())
		return false;
	getline(input, chars);
	if (input.bad())
		throw IOException("Cannot read from \""+filename+"\".");
#ifdef C5IDE
	this->wholeFile += chars;
	this->wholeFile += '\n';
#endif // #ifdef C5IDE
	return true;
}

bool FileRegionLineProvider::doReadNextLine(const LexerStatus &, const ParserNS::ParserStatus &, string &chars)
{
  ++myCurrLine;  // Become index of the line we are about to read...
  const bool EndOfRegion = (myCurrLine > myToLine) ||
                           (myCurrLine == myToLine && myToChar == 1);
  if (input.eof() || EndOfRegion) return false;
  getline(input, chars);
  if (input.bad())
    throw IOException("Cannot read from \""+filename+"\".");
  if (myCurrLine < myFromLine) chars.clear();
  if (myCurrLine == myFromLine && static_cast<size_t>(myFromChar-1) > chars.size())
    throw IOException("SourceRegion FromChar beyond EOL");
  if (myCurrLine == myToLine && static_cast<size_t>(myToChar-1) > chars.size())
    throw IOException("SourceRegion ToChar beyond EOL");
  if (myCurrLine == myToLine) chars = chars.erase(myToChar-1);       // This line first;
  if (myCurrLine == myFromLine) chars = chars.substr(myFromChar-1);  // this line second!
#ifdef C5IDE
  this->wholeFile += chars;
  this->wholeFile += '\n';
#endif // #ifdef C5IDE
  return true;
}

} // namespace LexerNS
} // namespace CoCoA
