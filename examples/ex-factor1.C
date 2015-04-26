// Copyright (c)  2009  John Abbott
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"
using namespace CoCoA;
using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This example shows how to interpret the result of a factorization. \n";

const string LongDescription =
  "This example shows how to interpret the result of a factorization. \n"
  "It creates a ring element (belonging to a polynomial ring), and    \n"
  "factorizes it.  The result is a \"factorization object\".  We show \n"
  "how to access the various parts of this object.                    \n";
//----------------------------------------------------------------------

void program()
{
  GlobalManager CoCoAFoundations;

  SparsePolyRing P = NewPolyRing(RingQQ(), 2); // QQ[x,y];
  const RingElem& x = indet(P,0);
  const RingElem& y = indet(P,1);
  RingElem f = power(x,96) - power(y,96); // f = x^96-y^96

  const factorization<RingElem> FacInfo = factor(f);

  // These are convenient aliases for the 3 fields in the factorization:
  const RingElem&         content   = FacInfo.myRemainingFactor();
  const vector<RingElem>& IrredFacs = FacInfo.myFactors();
  const vector<long>&     mult      = FacInfo.myMultiplicities();

  cout << "The factors of " << f << " are:" << endl;
  if (!IsOne(content))
    cout << "content: " << content << endl;
  const int NumIrredFacs = len(IrredFacs);
  for (int i = 0; i != NumIrredFacs; ++i)
  {
    cout << IrredFacs[i] << "  with multiplicity  " << mult[i] << endl;
  }
}

//----------------------------------------------------------------------
// Use main() to handle any uncaught exceptions and warn the user about them.
int main()
{
  try
  {
    program();
    return 0;
  }
  catch (const CoCoA::ErrorInfo& err)
  {
    cerr << "***ERROR***  UNCAUGHT CoCoA error";
    ANNOUNCE(cerr, err);
  }
  catch (const std::exception& exc)
  {
    cerr << "***ERROR***  UNCAUGHT std::exception: " << exc.what() << endl;
  }
  catch(...)
  {
    cerr << "***ERROR***  UNCAUGHT UNKNOWN EXCEPTION" << endl;
  }
  return 1;
}

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-factor1.C,v 1.5 2014/03/24 12:09:20 abbott Exp $
// $Log: ex-factor1.C,v $
// Revision 1.5  2014/03/24 12:09:20  abbott
// Summary: Major revision to public interface of factorization template class
// Author: JAA
//
// Revision 1.4  2012/10/05 09:29:43  abbott
// Changed myExponents into myMultiplicities.
//
// Revision 1.3  2012/02/08 17:52:17  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.2  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.1  2009/09/23 14:09:36  abbott
// First example for (polynomial) factorization.
//
//
