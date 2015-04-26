// Copyright (c) 2013  John Abbott
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"
using namespace CoCoA;
using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This program shows how to reconstruct a rational from several     \n"
  "residue-modulus pairs (allowing for some of them to be wrong).    \n";

const string LongDescription =
  "This program shows how to reconstruct a rational from several        \n"
  "residue-modulus pairs (allowing for some of them to be wrong).       \n"
  "CoCoALib offers two methods: one using continued fractions, the other\n"
  "using lattice reduction.  The actual reconstruction is performed     \n"
  "inside ContFracTrial or LatticeTrial.                                \n";

//----------------------------------------------------------------------

long PossiblyFaultyResidue(const BigRat& TrueAns, long p, double ErrorRate)
{
  if (den(TrueAns)%p == 0 || RandomBiasedBool(ErrorRate))
    return RandomLong(0,p-1); // BAD residue

  return ((num(TrueAns)%p)*InvMod(den(TrueAns)%p, p))%p;  // GOOD residue
}

struct ResModPair
{
  ResModPair(long r, long m): res(r), mod(m) {}
  long res;
  long mod;
};

vector<ResModPair> MakeResidues(const BigRat& TrueAns, double ErrorRate)
{
  const long NumPairs = 500;
  vector<ResModPair> ans; ans.reserve(NumPairs);
  long p = 101;
  for (long i=0; i < NumPairs; ++i)
  {
    p = NextPrime(p);
    const long res = PossiblyFaultyResidue(TrueAns, p, ErrorRate);
    ans.push_back(ResModPair(res, p));
  }
  return ans;
}


void ContFracTrial(const vector<ResModPair>& ModularImages, const BigRat& TrueAns)
{
  cout << "ContFrac method: ";
  RatReconstructByContFrac reconstructor(0); // 0 --> use default threshold

  for (long NumIters=0; NumIters < len(ModularImages); ++NumIters)
  {
    const long res = ModularImages[NumIters].res;
    const long mod = ModularImages[NumIters].mod;
    reconstructor.myAddInfo(res,mod);
    if (IsConvincing(reconstructor))
    {
      if (ReconstructedRat(reconstructor) == TrueAns)
      {
        cout << " success after " << 1+NumIters << " iters!" << endl;
        return;
      }
      cout << "false-positive after " << 1+NumIters << " iters;  ";
    }
  }
  cout << "FAILED after " << len(ModularImages) << " iters!" << endl;
}


void LatticeTrial(const vector<ResModPair>& ModularImages, const BigRat& TrueAns)
{
  cout << "Lattice method: ";
  RatReconstructByLattice reconstructor(0); // 0 --> use default SafetyFactor

  for (long NumIters=0; NumIters < len(ModularImages); ++NumIters)
  {
    const long res = ModularImages[NumIters].res;
    const long mod = ModularImages[NumIters].mod;
    reconstructor.myAddInfo(res,mod);
    if (IsConvincing(reconstructor))
    {
      if (ReconstructedRat(reconstructor) == TrueAns)
      {
        cout << " success after " << 1+NumIters << " iters!" << endl;
        return;
      }
      cout << "false-positive after " << 1+NumIters << " iters;  ";
    }
  }
  cout << "FAILED after " << len(ModularImages) << " iters!" << endl;
}


void program()
{
  GlobalManager CoCoAFoundations;

  BigInt N,D;
  vector<ResModPair> ModularImages;

  cout << "ContFrac and lattice method are similar if result is BALANCED:" << endl;

  N = RandomBigInt(-power(2,350), power(2,350));
  D = RandomBigInt(1, power(2, 350));
  BigRat BalancedRat(N,D);

  cout << "Trial: with error rate = 0.1" << endl;
  ModularImages = MakeResidues(BalancedRat, 0.1);
  ContFracTrial(ModularImages,BalancedRat);
  LatticeTrial(ModularImages, BalancedRat);

  cout << "Trial: with error rate = 0.3" << endl;
  ModularImages = MakeResidues(BalancedRat, 0.3);
  ContFracTrial(ModularImages,BalancedRat);
  LatticeTrial(ModularImages, BalancedRat);


  cout << endl
       << "ContFrac works better if result is UNBALANCED:" << endl;

  N = RandomBigInt(-power(2,500), power(2,500));
  D = RandomBigInt(1, power(2, 200));
  BigRat UnbalancedRat(N,D);

  cout << "Trial: with error rate = 0.1" << endl;
  ModularImages = MakeResidues(UnbalancedRat, 0.1);
  ContFracTrial(ModularImages,UnbalancedRat);
  LatticeTrial(ModularImages, UnbalancedRat);

  cout << "Trial: with error rate = 0.3" << endl;
  ModularImages = MakeResidues(UnbalancedRat, 0.3);
  ContFracTrial(ModularImages,UnbalancedRat);
  LatticeTrial(ModularImages, UnbalancedRat);

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

