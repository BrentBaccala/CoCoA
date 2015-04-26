// Copyright (c) 2011  John Abbott
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"
using namespace CoCoA;
using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This program makes simple use of the decimal conversion function, \n"
  "MantissaAndExponent10";

const string LongDescription =
  "Example of use of MantissaAndExponent10 and MantExp10 structure.          \n"
  "Also use of FloatStr to convert a integer/rational into a decimal string. \n";
//----------------------------------------------------------------------


void program()
{
  GlobalManager CoCoAFoundations;

  cout << ShortDescription << endl;

  const BigRat v = -(1-BigRat(5,100000)); // -0.99995
  cout << "True rational value is " << v << endl << endl;
  for (int sigfig = 2; sigfig <= 6; ++sigfig)
  {
    cout << "Asking for " << sigfig << " significant figures:" << endl;
    const MantExp10 ME = MantissaAndExponent10(v, sigfig);
    cout << "MantExp = " << ME << endl;
    cout << endl;
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
