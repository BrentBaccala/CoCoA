// Copyright (c) 2013  John Abbott -- Author: Mario Albert
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"
#include <bitset>
using namespace CoCoA;
using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "In this file we explain how to compute the Janet Basis in CoCoA.  \n";

const string LongDescription =
  " We explain the function JanetBasis and show whicht options we can choose. \n";
//----------------------------------------------------------------------

void output(vector<RingElem> vec)
{
  for(vector<RingElem>::iterator iter = vec.begin(); iter != vec.end(); ++iter)
    {
      cout << *iter << endl;
    }
}


void program()
{
  GlobalManager CoCoAFoundations;

  cout << ShortDescription << endl;
  ring Q = RingQQ();

  cout <<"/////////////////////////////////////Cyclic 5//////////////////////////////////////////" << endl;
  SparsePolyRing polyRingCyc5 = NewPolyRing(Q, SymbolRange("x",0,4), StdDegLex);
  const vector<RingElem> x = indets(polyRingCyc5);
  vector<RingElem> cyclic5;
  cyclic5.push_back(x[0]+x[1]+x[2]+x[3]+x[4]);
  cyclic5.push_back(x[0]*x[1] + x[2]*x[1] + x[2]*x[3] + x[3]*x[4] + x[0]*x[4]);
  cyclic5.push_back(x[1]*x[2]*x[0] + x[2]*x[3]*x[1] + x[4]*x[2]*x[3] + x[0]*x[4]*x[1] + x[0]*x[4]*x[3]);
  cyclic5.push_back(x[0]*x[1]*x[2]*x[3]+x[0]*x[1]*x[2]*x[4]+x[0]*x[1]*x[3]*x[4]+x[0]*x[2]*x[3]*x[4]+x[2]*x[3]*x[1]*x[4]);
  cyclic5.push_back(x[0]*x[1]*x[2]*x[3]*x[4] -1 );
  ideal I = ideal(cyclic5);
  double timing = CpuTime();
  vector<RingElem> GroebnerBasis = TidyGens(I);
  cout << "Time taken (TinyGens) is " << CpuTime() - timing << endl;
  cout << "size of GroebnerBasis = " << GroebnerBasis.size() << endl << endl;


  cout << "---------computing the Janet Base with a list of polynomials(TQBlockLow):------------------------" << endl;
  vector<RingElem> result;
  cout << "-----------------------------first two criteria + GB---------------------------------" << endl;
  result = JanetBasis(cyclic5, GB, TQBlockLow);
  std::cout << "result.size() = " << result.size() << std::endl;

  cout << "-----------------------------first two criteria + JB---------------------------------" << endl;
  result = JanetBasis(cyclic5, TQBlockLow);
  std::cout << "result.size() = " << result.size() << std::endl;

  cout << "-----------------------------first + third criteria + GB-----------------------------" << endl;
  result = JanetBasis(cyclic5, std::bitset<3>(4), GB, TQBlockLow);
  std::cout << "result.size() = " << result.size() << std::endl;

  cout << "-----------------------------first + third criteria + JB-----------------------------" << endl;
  result = JanetBasis(cyclic5, bitset<3>(4), TQBlockLow);
  std::cout << "result.size() = " << result.size() << std::endl << endl;

  cout << "---------computing the Janet Base with a list of polynomials (TQBlockHigh):------------------------" << endl;
  cout << "-----------------------------first two criteria + GB---------------------------------" << endl;
  result = JanetBasis(cyclic5, GB, TQBlockHigh);
  std::cout << "result.size() = " << result.size() << std::endl;

  cout << "-----------------------------first two criteria + JB---------------------------------" << endl;
  result = JanetBasis(cyclic5, TQBlockHigh);
  std::cout << "result.size() = " << result.size() << std::endl;

  cout << "-----------------------------first + third criteria + GB-----------------------------" << endl;
  result = JanetBasis(cyclic5, std::bitset<3>(4), GB, TQBlockHigh);
  std::cout << "result.size() = " << result.size() << std::endl;

  cout << "-----------------------------first + third criteria + JB-----------------------------" << endl;
  result = JanetBasis(cyclic5, bitset<3>(4), TQBlockHigh);
  std::cout << "result.size() = " << result.size() << std::endl << endl;

  cout << "---------computing the Janet Base with a list of polynomials (TQDegree):------------------------" << endl;
  cout << "-----------------------------first two criteria + GB---------------------------------" << endl;
  result = JanetBasis(cyclic5, GB, TQDegree);
  std::cout << "result.size() = " << result.size() << std::endl;

  cout << "-----------------------------first two criteria + JB---------------------------------" << endl;
  result = JanetBasis(cyclic5, TQDegree);
  std::cout << "result.size() = " << result.size() << std::endl;

  cout << "-----------------------------first + third criteria + GB-----------------------------" << endl;
  result = JanetBasis(cyclic5, std::bitset<3>(4), GB, TQDegree);
  std::cout << "result.size() = " << result.size() << std::endl;

  cout << "-----------------------------first + third criteria + JB-----------------------------" << endl;
  result = JanetBasis(cyclic5, bitset<3>(4), TQDegree);
  std::cout << "result.size() = " << result.size() << std::endl << endl;


  cout << "---------now the same again, but with an ideal instead of a vector of polynomials----" << endl;
  ideal cyc5Ideal(cyclic5);


  cout << "-----------------------------first two criteria + GB---------------------------------" << endl;
  result = JanetBasis(cyc5Ideal, GB); 
  std::cout << "result.size() = " << result.size() << std::endl;

  cout << "-----------------------------first two criteria + JB---------------------------------" << endl;
  result = JanetBasis(cyc5Ideal);
  std::cout << "result.size() = " << result.size() << std::endl;

  cout << "-----------------------------first + third  criteria + GB----------------------------" << endl;
  result = JanetBasis(cyc5Ideal, std::bitset<3>(3), GB);
  std::cout << "result.size() = " << result.size() << std::endl;

  cout << "-----------------------------first + third  criteria + GB----------------------------" << endl;
  result = JanetBasis(cyc5Ideal, bitset<3>(3));
  std::cout << "result.size() = " << result.size() << std::endl << endl << endl ;

  cout << "//////////////////////Example over Fp7//////////////////////////////////////////////" << endl;
  ring Fp7 = NewZZmod(7);
  SparsePolyRing polyRing = NewPolyRing(Q, SymbolRange("x",0,2), NewStdDegLexOrdering(3));
  const vector<RingElem> y = indets(polyRing);
  vector<RingElem> gEx;
  gEx.push_back(2*y[0]*y[1]*power(y[2],3) - power(y[0],4));
  gEx.push_back(7*power(y[1],5)*y[2] - 6*power(y[1],4));
  gEx.push_back(4*power(y[1],6));
  gEx.push_back(2*power(y[0],4)*y[1]*y[1]);

  vector<RingElem> result2;
  cout << "-----------------------------first two criteria + GB---------------------------------" << endl;
  result2 = JanetBasis(gEx, GB); 
  std::cout << "result2.size() = " << result2.size() << std::endl << endl;
  cout << "-----------------------output JanetBasis----------------------------------------------" << endl;
  output(result2);
  cout << endl;


  cout << "---------------------------first two criteria + JB-----------------------------------" << endl;
  result2 = JanetBasis(gEx);
  std::cout << "result2.size() = " << result2.size() << std::endl << endl;


  cout << "----------------------------third and first criteria + GB----------------------------" << endl;
  result2 = JanetBasis(gEx, std::bitset<3>(4), GB);
  std::cout << "result2.size() = " << result2.size() << std::endl << endl;


  cout << "-------------------------first and third criteria + JB-------------------------------" << endl;
  result2 = JanetBasis(gEx, bitset<3>(4));
  std::cout << "result2.size() = " << result2.size() << std::endl << endl;


  cout << "------------------------------TidyGens-----------------------------------------------" << endl;
  ideal gIdeal(gEx);
  result2 = TidyGens(gIdeal);
  std::cout << "result2.size() = " << result2.size() << std::endl << endl;

  cout << "-----------------------------output Tidy Gens----------------------------------------" << endl;
  output(result2);
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
