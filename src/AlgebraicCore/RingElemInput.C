//   Copyright (c)  2014  Anna Bigatti

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


// Source code for input of polynomials and polynomial expressions

#include "CoCoA/RingElemInput.H"

#include "CoCoA/BigInt.H"
#include "CoCoA/error.H"
#include "CoCoA/ring.H"
#include "CoCoA/symbol.H"

#include <istream>
using std::istream;
//#include <iostream> // just for debugging
#include <sstream>
using std::istringstream;
#include <string>
using std::string;


namespace CoCoA
{

  char WhatsNext(std::istream& in)
  { // recognized characters
    in >> std::ws; // skip white spaces
    //    std::cout << " WhatsNext: testing " << char(in.peek()) << std::endl;
    if (in.eof()) return '\0';
    if (isdigit(in.peek())) return 'd'; // d for digit
    if (isalpha(in.peek())) return 'a'; // a for alpha
    if (in.peek()=='+') return '+';
    if (in.peek()=='-') return '-';
    if (in.peek()=='*') return '*';
    if (in.peek()=='/') return '/';
    if (in.peek()=='^') return '^';
    if (in.peek()=='(') return '(';
    if (in.peek()==')') return ')';
    if (in.peek()==';') return ';';
    CoCoA_ERROR("Illegal char \'"+string(1,char(in.peek()))+"\'","WhatsNext");
    return 'o';
  }


  void CoCoA_ERROR_NextChar(std::istream& in, const string& FuncName)
  {
    in.clear();
    CoCoA_ERROR("Unexpected \'"+ string(1,char(in.peek())) +"\'", FuncName);
  }


  RingElem ReadFactor(const ring& P, std::istream& in)
  {
    symbol s("temp");
    BigInt n;
    RingElem tmp(P);

    switch ( WhatsNext(in) )
    {
    case 'a':
      in >> s;  if (!in) CoCoA_ERROR_NextChar(in, "ReadFactor: symbol");
      tmp = RingElem(P, s);
      break;
    case 'd':
      in >> n;  if (!in) CoCoA_ERROR_NextChar(in, "ReadFactor: BigInt");
      tmp = RingElem(P, n);
      break;
    case '(':
      in.ignore(); // '('
      tmp = ReadExpr(P, in);
      if ( WhatsNext(in) != ')' )  CoCoA_ERROR_NextChar(in, "ReadFactor: ()");
      in.ignore(); // ')'
      break;
    default: CoCoA_ERROR_NextChar(in, "ReadFactor");
    }
    if ( WhatsNext(in) != '^' )  return tmp;
    in.ignore(); // '^'
    in >> n;  if (!in) CoCoA_ERROR_NextChar(in, "ReadFactor: exponent");
    return power(tmp, n);
  }


  RingElem ReadProduct(const ring& P, std::istream& in)
  {
    char w = '*';
    RingElem resProd(one(P));

    while (true)
    {
      if ( w == '*' )  resProd *= ReadFactor(P, in);
      else             resProd /= ReadFactor(P, in); // ( w == '/' )
      if ( WhatsNext(in) != '*' &&  WhatsNext(in) != '/' )  break;
      // uncomment the following to prevent a/b*c, a/b/c
      // if ( w == '/' ) 
      //   CoCoA_ERROR_NextChar(in, "ReadProduct: ambiguous after \'/\'");
      in >> w; // '*' or '/'
    }
    return resProd;
  }


  RingElem ReadExpr(const ring& P, std::istream& in)
  {
    char w = WhatsNext(in);
    RingElem f(P);

    while (true)
    {
      if ( w == '-' || w == '+' )  in.ignore(); // '+/-'
      if ( w == '-' )  f -= ReadProduct(P,in);  else  f += ReadProduct(P,in);
      w = WhatsNext(in);
      if ( (w == ';') || (w == ')') || (w == '\0') )  break;
      if ( (w != '-') && (w != '+') )  CoCoA_ERROR_NextChar(in, "ReadExpr");
    }
    return f;
  }


  RingElem ReadExprSemicolon(const ring& P, std::istream& in)
  {
    RingElem f = ReadExpr(P, in);
    if ( WhatsNext(in) != ';' )  CoCoA_ERROR_NextChar(in, "ReadExprSemicolon");
    in.ignore(); // ';'
    return f;
  }


  //-------- input from string --------------------
  RingElem ReadExprSemicolon(const ring& P, const std::string& s)
  {
    istringstream is(s);
    return ReadExprSemicolon(P, is);
  }


  RingElem ReadExpr(const ring& P, const std::string& s)
  {
    istringstream is(s);
    return ReadExpr(P, is);
  }


} // end of namespace CoCoA

// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/RingElemInput.C,v 1.3 2014/03/21 15:57:19 bigatti Exp $
// $Log: RingElemInput.C,v $
// Revision 1.3  2014/03/21 15:57:19  bigatti
// -- Ring is now first argument of ReadExpr
//
// Revision 1.2  2014/01/30 16:21:53  bigatti
// -- removed restriction about 1/3*x
//
// Revision 1.1  2014/01/30 15:15:33  bigatti
// -- first import
//
