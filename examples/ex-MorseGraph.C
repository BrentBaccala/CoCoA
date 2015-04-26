// Copyright (c) 2013  Mario Albert
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace CoCoA;
using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "In this file we explain how to compute a free Resolution via Pommaret Basis in CoCoA.  \n";

const string LongDescription =
  " JBResolution, JBBettiDiagramm, JBMinimalResolution. \n";
//----------------------------------------------------------------------



void program()
{
  GlobalManager CoCoAFoundations;

  cout << ShortDescription << endl;
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing_DMPI(Q,3);
  const vector<RingElem>& x = indets(polyRing);
  vector<RingElem> input;
  input.push_back(   power(x[0],2) + 4*power(x[1],2) - 4);
  input.push_back(   - x[0] + 2*power(x[1],2));


  //////////////////////////////////////////////////////////////////
  std::cout << "Resolution of ideal" << std::endl;

  JBMill mill1(ExtendedJanetBasis(input));
  std::vector<matrix> res1(JBResolution(mill1));
  const long LenRes = len(res1);
  for (long i=0; i < LenRes; ++i)
  {
    std::cout << "F[" << i << "] == " << res1[i] << std::endl;
  }
  

  //////////////////////////////////////////////////////////////////
  std::cout << "Betti Diagram of homogenized Ideal" << std::endl;

  input = JBReturnGB(mill1);
  for (std::vector<RingElem>::iterator i = input.begin(); i != input.end(); ++i)
  {
    *i = homog(*i, x[2]);
  }
  JBMill mill2 = ExtendedJanetBasis(input);
  std::cout << JBBettiDiagram(mill2) << std::endl;


  //////////////////////////////////////////////////////////////////
  std::cout << "Minimal Free Resolution of homogenized Ideal" << std::endl;

  std::vector<matrix> res2(JBMinimalResolution(mill2));
  for (long i = 0; i < len(res2); ++i)
  {
    std::cout << "F[" << i << "] == " << res2[i] << std::endl;
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
