#include "CoCoA/TmpResolutionMinimization.H"
#include "CoCoA/matrix.H"
#include "CoCoA/utils.H"

using std::make_pair;

namespace CoCoA
{

  std::pair<long, long> ResolutionMinimization::myFindPivot(ConstMatrixView m) const
  {
    const long nrows = NumRows(m);
    const long ncols = NumCols(m);
    for (long j = 0; j < ncols; ++j)
    {
      for (long i = 0; i < nrows; ++i)
      {
        const RingElemAlias mij = m->myEntry(i,j);
//        if (!myRing->myIsZero(raw(mij)) && myRing->myIsConstant(raw(mij)))  // SLIGHTLY FASTER
        if (myRing->myIsInvertible(raw(mij)))
//        if (IsInvertible(m(i, j)))
        {
          return make_pair(i, j);
        }
      }
    }
    return make_pair(-1, -1);
  }


  void ResolutionMinimization::myManipulateMatrix(matrix& m, long r, const std::vector<RingElem>& PivotColumn) const
  {
    CoCoA_ASSERT(0 <= r && r < NumRows(m));
    CoCoA_ASSERT(len(PivotColumn) == NumRows(m));
    CoCoA_ASSERT(IsInvertible(PivotColumn[r]));
    // m behaves as a pointer -> no need to give something back
    const RingElem ScaleFactor = (-1)/PivotColumn[r];
    for (long i = 0; i < NumRows(m); ++i)
    {
      if (i != r)
      {
        m->myAddRowMul(i, r, ScaleFactor*PivotColumn[i]);
      }
    }
  }


  void ResolutionMinimization::myMinimization()
  {
    std::vector<matrix>::iterator begin(myResolution.begin());
    ++begin;
    for (std::vector<matrix>::iterator it = begin; it != myResolution.end(); ++it)
    {
      std::pair<long, long> PivotElement(myFindPivot(*it));
      while (PivotElement.first != -1)
      {
        std::vector<RingElem> PivotColumn(myGetColumn(*it, PivotElement.second));
        DeleteCol(*it, PivotElement.second);
        myManipulateMatrix(*it, PivotElement.first, PivotColumn);
        DeleteRow(*it, PivotElement.first);
        if (it != myResolution.begin())
        {
          std::vector<matrix>::iterator TmpIter(it);
          --TmpIter;
          DeleteCol(*TmpIter, PivotElement.first);
        }
        if (it != (--myResolution.end()))
        {
          std::vector<matrix>::iterator TmpIter(it);
          ++TmpIter;
          DeleteRow(*TmpIter, PivotElement.second);
        }
        PivotElement = myFindPivot(*it);
      }
    }
  }


  std::vector<RingElem> ResolutionMinimization::myGetColumn(ConstMatrixView m, long col) const
  {
    CoCoA_ASSERT(0 <= col && col < NumCols(m));
    std::vector<RingElem> res;
    const long nrows = NumRows(m);
    res.reserve(nrows);
    for (long i = 0; i < nrows; ++i)
    {
      res.push_back(m(i, col));
    }
    return res;
  }

} // end of namespace CoCoA
