#ifdef CoCoA_WITH_NORMALIZ
//   Copyright (c)  2010-2014 Anna M. Bigatti, Christof Soeger

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


// Source code for Normaliz integration

#include "CoCoA/ExternalLibs-Normaliz.H"
#include "CoCoA/BigInt.H"
#include "CoCoA/BigRat.H"
#include "CoCoA/IntOperations.H"
#include "CoCoA/TmpPPVector.H"
#include "CoCoA/convert.H"
#include "CoCoA/error.H"
#include "CoCoA/matrix.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/symbol.H"
#include "CoCoA/QuasiPoly.H"

#include "CoCoA/VectorOperations.H"  // just for debugging

#include "gmpxx.h"
#include "libnormaliz/cone.h"
#include "libnormaliz/map_operations.h"
using libnormaliz::Cone;
using libnormaliz::ConeProperties;



// #include <list>
// using std::list;
#include <map>
using std::map;
#include <vector>
using std::vector;
#include <ostream>
using std::ostream;
using std::endl;
#include <sstream>  // for ErrorMessage


namespace CoCoA
{


  class TemporarilySet
  {
  public:
    TemporarilySet(bool& var, bool value): myVar(var), myOrigValue(var) { myVar = value; }
    ~TemporarilySet() { myVar = myOrigValue; }
  private: // data members
    bool& myVar;
    bool myOrigValue;
  };

  namespace Normaliz
  {

    libnormaliz::InputType ToInputType(const std::string& TypeString)
    {
      return libnormaliz::to_type(TypeString);
    }

    class ConeImpl: protected IntrusiveReferenceCount
    {
      friend class SmartPtrIRC<const ConeImpl>; // Morally "friend Cone", so it can alter reference count.

      public:
        ConeImpl(const std::vector<std::vector<BigInt> >& v, libnormaliz::Type::InputType InputType);
        ConeImpl(const std::map< libnormaliz::InputType, std::vector<std::vector<BigInt> > >& m);
        ConeImpl(const PPVector& v, libnormaliz::Type::InputType InputType);
        friend std::vector<std::vector<BigInt> > HilbertBasis(const cone& c);
        friend std::vector<std::vector<BigInt> > ModuleGenerators(const cone& c);
        friend std::vector<std::vector<BigInt> > Generators(const cone& c);
        friend std::vector<std::vector<BigInt> > ExtremeRays(const cone& c);
        friend std::vector<std::vector<BigInt> > VerticesOfPolyhedron(const cone& c);
        friend std::vector<std::vector<BigInt> > Deg1Elements(const cone& c);
        friend std::vector<std::vector<BigInt> > GeneratorsOfToricRing(const cone& c);
        friend std::vector<std::vector<BigInt> > SupportHyperplanes(const cone& c);
        friend std::vector<std::vector<BigInt> > Equations(const cone& c);
        friend std::vector<std::vector<BigInt> > Congruences(const cone& c);
        friend std::vector<std::vector<BigInt> > ExcludedFaces(const cone& c);
        //friend std::vector<BigInt> HVector(const cone& c);
        friend HPSeries HilbertSeries(const cone& c);
        friend RingElem HilbertPoly(const cone& c);
        friend QuasiPoly HilbertQuasiPoly(const cone& c);
        friend BigRat multiplicity(const cone& c);
        friend std::vector<BigRat> grading(const cone& c);
        friend bool IsPointed(const cone& c);
        friend bool IsInhomogeneous(const cone& c);
        friend bool IsIntegrallyClosed(const cone& c);
        friend bool IsDeg1HilbertBasis(const cone& c);
        friend long EmbeddingDim(const cone& c);
        friend long rank(const cone& c);
        friend long RecessionRank(const cone& c);
        friend long AffineDim(const cone& c);
        friend long ModuleRank(const cone& c);
        friend std::vector<BigInt> dehomogenization(const cone& c);
        friend BigInt shift(const cone& c);

      public:
        void myComputation(const libnormaliz::ConeProperties& CPs) const;
        void myComputation(libnormaliz::ConeProperty::Enum CP) const;    
        void myComputation(libnormaliz::ConeProperty::Enum CP1, libnormaliz::ConeProperty::Enum CP2) const;    
        void myComputation() const;   //default: compute everything possible 
        bool isComputed(libnormaliz::ConeProperty::Enum CP) const;    
        bool isComputed(libnormaliz::ConeProperty::Enum CP1, libnormaliz::ConeProperty::Enum CP2) const;    

      private:
        mutable bool mySmallConeIsGood;
        mutable libnormaliz::Cone<long> myConeLong; // be careful!
        mutable libnormaliz::Cone<mpz_class> myConeMPZ; // be careful!
    };


    // implementation of cone

    cone::cone(const std::vector<std::vector<BigInt> >& v, libnormaliz::Type::InputType InputType) : mySmartPtr(new ConeImpl(v, InputType)) {}

    cone::cone(const std::map< libnormaliz::InputType, std::vector<std::vector<BigInt> > >& m) : mySmartPtr(new ConeImpl(m)) {}

    cone::cone(ConstMatrixView M, libnormaliz::Type::InputType InputType) : mySmartPtr(new ConeImpl(MatrixToVecVecBigInt(M), InputType)) {}

    cone::cone(const ConeImpl* ptr): mySmartPtr(ptr) {}
    cone::cone(const cone& c) : mySmartPtr(c.mySmartPtr) {}
    cone::~cone() {}


    // printing
    ostream& operator<< (ostream& out, const cone& C)
    {
      using namespace libnormaliz;
      out << "cone(EmbeddingDim = " << EmbeddingDim(C);
      if (C.isComputed(ConeProperty::ExtremeRays)) {
        out << ", rank = " << rank(C);
        out << ", pointed = " << IsPointed(C);
      }
      if (C.isComputed(ConeProperty::Generators))
        out << ", NumGenerators = " << Generators(C).size();
      if (C.isComputed(ConeProperty::SupportHyperplanes))
        out << ", NumSupportHyperplanes = " << SupportHyperplanes(C).size();
      out << ")";
      return out;
    }
    
    // computations
    void cone::myComputation(const libnormaliz::ConeProperties& CPs) const { mySmartPtr->myComputation(CPs); }
    void cone::myComputation(libnormaliz::ConeProperty::Enum CP) const { mySmartPtr->myComputation(CP); }
    void cone::myComputation(libnormaliz::ConeProperty::Enum CP1, libnormaliz::ConeProperty::Enum CP2) const { mySmartPtr->myComputation(CP1, CP2); }
    void cone::myComputation() const { mySmartPtr->myComputation(); }
    bool cone::isComputed(libnormaliz::ConeProperty::Enum CP) const { return mySmartPtr->isComputed(CP); }
    bool cone::isComputed(libnormaliz::ConeProperty::Enum CP1, libnormaliz::ConeProperty::Enum CP2) const { return mySmartPtr->isComputed(CP1, CP2); }
    const ConeImpl* cone::operator->() const { return mySmartPtr.operator->(); }

    namespace  // conversion functions, implementation at the end of the file
    {
      std::vector<BigInt> LongVToBigIntV(const std::vector<long>& VIn);
      std::vector<BigInt> MPZ_classVToBigIntV(const std::vector<mpz_class>& VIn);
      std::vector<BigRat> LongVToBigRatV(const std::vector<long>& VIn, const BigInt& denom);
      std::vector<BigRat> MPZ_classVToBigRatV(const std::vector<mpz_class>& VIn, const BigInt& denom);
      std::vector<long> BigIntVToLongV(const std::vector<BigInt>& VIn);
      std::vector<mpz_class> BigIntVToMPZ_classV(const std::vector<BigInt>& VIn);
      void ConvertFromNormaliz(std::vector<std::vector<BigInt> >& v, const std::vector<std::vector<long> >& l);
      void ConvertFromNormaliz(std::vector<std::vector<BigInt> >& v, const std::vector<std::vector<mpz_class> >& l);
      std::vector<std::vector<long> >  ReturnVecVecLong(const std::vector<std::vector<BigInt> >& v);
      map< libnormaliz::InputType, vector<vector<long> > >  ReturnMapVecVecLong(const map< libnormaliz::InputType, vector<vector<BigInt> > > m);
      std::vector<std::vector<mpz_class> >  ReturnVecVecMPZ_class(const std::vector<std::vector<BigInt> >& v);
      map< libnormaliz::InputType, vector<vector<mpz_class> > >  ReturnMapVecVecMPZ_class(const map< libnormaliz::InputType, vector<vector<BigInt> > > m);
      std::vector<std::vector<BigInt> > PPVectorToVecVecBigInt(const PPVector& ppv);
      void VecVecBigIntToPPVector(PPVector& ppv, const std::vector<std::vector<BigInt> >& M);
    } // end of anonymous namespace


    // implementation of ConeImpl
    ConeImpl::ConeImpl(const std::vector<std::vector<BigInt> >& v, libnormaliz::Type::InputType InputType):
          myConeLong(ReturnVecVecLong(v), InputType),
          myConeMPZ(ReturnVecVecMPZ_class(v), InputType)
    {
      mySmallConeIsGood = (myConeLong.getBasisChange().get_rank() > 0);
    }

    ConeImpl::ConeImpl(const std::map< libnormaliz::InputType, std::vector<std::vector<BigInt> > >& m):
      myConeLong(ReturnMapVecVecLong(m)),
      myConeMPZ(ReturnMapVecVecMPZ_class(m))
    {
      mySmallConeIsGood = (myConeLong.getBasisChange().get_rank() > 0);
    }

    void ConeImpl::myComputation(const libnormaliz::ConeProperties& CPs) const
    {
      try
      {
        if (mySmallConeIsGood)
        {
          TemporarilySet save(libnormaliz::test_arithmetic_overflow, true);
          libnormaliz::ConeProperties missing = myConeLong.compute(CPs);
          if (missing.any())
          {
            std::ostringstream os;
            os << missing;
            CoCoA_ERROR("Normaliz cannot compute "+ os.str(), "myComputation long");
            return;
          }
        }
      }
      catch (libnormaliz::ArithmeticException&) { mySmallConeIsGood = false; } 
      if ( !mySmallConeIsGood )
      {
        TemporarilySet save(libnormaliz::test_arithmetic_overflow, false);
        libnormaliz::ConeProperties missing = myConeMPZ.compute(CPs);
        if (missing.any())
        {
          std::ostringstream os;
          os << missing;
          CoCoA_ERROR("Normaliz cannot compute "+ os.str(), "myComputation MPZ");
        }
      }
    }

    void ConeImpl::myComputation(libnormaliz::ConeProperty::Enum CP) const
    {
      myComputation(ConeProperties(CP));
    }

    void ConeImpl::myComputation(libnormaliz::ConeProperty::Enum CP1,
                                 libnormaliz::ConeProperty::Enum CP2) const
    {
      myComputation(ConeProperties(CP1, CP2));
    }

    void ConeImpl::myComputation() const
    {
      myComputation(ConeProperties(libnormaliz::ConeProperty::DefaultMode));
    }

    bool ConeImpl::isComputed(libnormaliz::ConeProperty::Enum CP) const
    {
      if (mySmallConeIsGood)
      {
        return myConeLong.isComputed(CP);
      }
      else
      {
        return myConeMPZ.isComputed(CP);
      }
    }

    bool ConeImpl::isComputed(libnormaliz::ConeProperty::Enum CP1,
                              libnormaliz::ConeProperty::Enum CP2) const
    {
        return isComputed(CP1) && isComputed(CP2);
    }

    // friend functions which are the interface to the user
    std::vector<std::vector<BigInt> > HilbertBasis(const cone& c)
    {
      std::vector<std::vector<BigInt> > v;
      c->myComputation(libnormaliz::ConeProperty::HilbertBasis);
      if (c->mySmallConeIsGood)
        ConvertFromNormaliz(v, c->myConeLong.getHilbertBasis());
      else
        ConvertFromNormaliz(v, c->myConeMPZ.getHilbertBasis());
      return v;
    }

    std::vector<std::vector<BigInt> > ModuleGenerators(const cone& c)
    {
      std::vector<std::vector<BigInt> > v;
      c->myComputation(libnormaliz::ConeProperty::ModuleGenerators);
      if (c->mySmallConeIsGood)
        ConvertFromNormaliz(v, c->myConeLong.getModuleGenerators());
      else
        ConvertFromNormaliz(v, c->myConeMPZ.getModuleGenerators());
      return v;
    }

    std::vector<std::vector<BigInt> > Deg1Elements(const cone& c)
    { 
      std::vector<std::vector<BigInt> > v;
      c->myComputation(libnormaliz::ConeProperty::Deg1Elements);
      if (c->mySmallConeIsGood)
        ConvertFromNormaliz(v, c->myConeLong.getDeg1Elements());
      else
        ConvertFromNormaliz(v, c->myConeMPZ.getDeg1Elements());
      return v;
    }

    std::vector<std::vector<BigInt> > Generators(const cone& c)
    {
      std::vector<std::vector<BigInt> > v;
      c->myComputation(libnormaliz::ConeProperty::Generators);
      if (c->mySmallConeIsGood)
        ConvertFromNormaliz(v, c->myConeLong.getGenerators());
      else
        ConvertFromNormaliz(v, c->myConeMPZ.getGenerators());
      return v;
    }

    std::vector<std::vector<BigInt> > ExtremeRays(const cone& c)
    {
      std::vector<std::vector<BigInt> > v;
      c->myComputation(libnormaliz::ConeProperty::ExtremeRays);
      if (c->mySmallConeIsGood)
        ConvertFromNormaliz(v, c->myConeLong.getExtremeRays());
      else
        ConvertFromNormaliz(v, c->myConeMPZ.getExtremeRays());
      return v;
    }

    std::vector<std::vector<BigInt> > VerticesOfPolyhedron(const cone& c)
    {
      std::vector<std::vector<BigInt> > v;
      c->myComputation(libnormaliz::ConeProperty::VerticesOfPolyhedron);
      if (c->mySmallConeIsGood)
        ConvertFromNormaliz(v, c->myConeLong.getVerticesOfPolyhedron());
      else
        ConvertFromNormaliz(v, c->myConeMPZ.getVerticesOfPolyhedron());
      return v;
    }

    std::vector<std::vector<BigInt> > GeneratorsOfToricRing(const cone& c)
    {
      std::vector<std::vector<BigInt> > v;
      c->myComputation(libnormaliz::ConeProperty::GeneratorsOfToricRing);
      if (c->mySmallConeIsGood)
        ConvertFromNormaliz(v, c->myConeLong.getGeneratorsOfToricRing());
      else
        ConvertFromNormaliz(v, c->myConeMPZ.getGeneratorsOfToricRing());
      return v;
    }


    // The following constraints depend all on ConeProperty::SupportHyperplanes
    std::vector<std::vector<BigInt> > SupportHyperplanes(const cone& c)
    { 
      std::vector<std::vector<BigInt> > v;
      c->myComputation(libnormaliz::ConeProperty::SupportHyperplanes);
      if (c->mySmallConeIsGood)
        ConvertFromNormaliz(v, c->myConeLong.getSupportHyperplanes());
      else
        ConvertFromNormaliz(v, c->myConeMPZ.getSupportHyperplanes());
      return v;
    }

    std::vector<std::vector<BigInt> > Equations(const cone& c)
    { 
      std::vector<std::vector<BigInt> > v;
      c->myComputation(libnormaliz::ConeProperty::SupportHyperplanes);
      if (c->mySmallConeIsGood)
        ConvertFromNormaliz(v, c->myConeLong.getEquations());
      else
        ConvertFromNormaliz(v, c->myConeMPZ.getEquations());
      return v;
    }

    std::vector<std::vector<BigInt> > Congruences(const cone& c)
    { 
      std::vector<std::vector<BigInt> > v;
      c->myComputation(libnormaliz::ConeProperty::SupportHyperplanes);
      if (c->mySmallConeIsGood)
        ConvertFromNormaliz(v, c->myConeLong.getCongruences());
      else
        ConvertFromNormaliz(v, c->myConeMPZ.getCongruences());
      return v;
    }

    // excluded faces are only available if they were explicitly given
    std::vector<std::vector<BigInt> > ExcludedFaces(const cone& c)
    {
      std::vector<std::vector<BigInt> > v;
      c->myComputation(libnormaliz::ConeProperty::ExcludedFaces);
      if (c->mySmallConeIsGood)
        ConvertFromNormaliz(v, c->myConeLong.getExcludedFaces());
      else
        ConvertFromNormaliz(v, c->myConeMPZ.getExcludedFaces());
      return v;
    }



/* REPLACED by the better function HilbertSeries
    std::vector<BigInt> HVector(const cone& c)
    { 
      c->myComputation(libnormaliz::ConeProperty::HilbertSeries);
      if (c->mySmallConeIsGood)
        return MPZ_classVToBigIntV(c->myConeLong.getHilbertSeries().getNum());
      else
        return MPZ_classVToBigIntV(c->myConeMPZ.getHilbertSeries().getNum());
    }
*/
    HPSeries HilbertSeries(const cone& c)
    { 
      c->myComputation(libnormaliz::ConeProperty::HilbertSeries);
      const libnormaliz::HilbertSeries& HS = c->mySmallConeIsGood ? c->myConeLong.getHilbertSeries() : c->myConeMPZ.getHilbertSeries();
      return HPSeries(MPZ_classVToBigIntV(HS.getNum()), 
                      libnormaliz::to_vector(HS.getDenom()));
    }

    RingElem HilbertPoly(const cone& c)
    { 
      c->myComputation(libnormaliz::ConeProperty::HilbertSeries);
      const libnormaliz::HilbertSeries& HS = c->mySmallConeIsGood ? c->myConeLong.getHilbertSeries() : c->myConeMPZ.getHilbertSeries();
      if (HS.getPeriod() != 1)
      {
        CoCoA_ERROR("Hilbert function is a quasi-polynomial of periode > 1. This function works for regular polynomials only.", "HilbertPoly");
      }
      vector<BigRat> coeffs = MPZ_classVToBigRatV( HS.getHilbertQuasiPolynomial()[0],
                                                   BigInt(HS.getHilbertQuasiPolynomialDenom().get_mpz_t()));
      PolyRing QQt = RingQQt(1);
      const RingElem t = indet(QQt,0);
      RingElem tpower = one(QQt);
      RingElem hpoly(QQt);
      for (long i=0; i<len(coeffs); ++i)
      {
          hpoly += coeffs[i] * tpower;   // NB tpower = t^i
          tpower *= t;
      }

      return hpoly;
    }

    QuasiPoly HilbertQuasiPoly(const cone& c)
    { 
      c->myComputation(libnormaliz::ConeProperty::HilbertSeries);
      const libnormaliz::HilbertSeries& HS = c->mySmallConeIsGood ? c->myConeLong.getHilbertSeries() : c->myConeMPZ.getHilbertSeries();
      const long period = HS.getPeriod();
      if (period < 1)
      {
        CoCoA_ERROR("Hilbert function not computed.", "HilbertQuasiPoly");
      }

      const PolyRing QQt = RingQQt(1);
      const RingElem t = indet(QQt,0);
      RingElem tpower;
      vector<BigRat> coeffs;
      vector<RingElem> qp = vector<RingElem>(period, RingElem(QQt));
      for (long j=0; j<period; ++j)
      {
        coeffs = MPZ_classVToBigRatV( HS.getHilbertQuasiPolynomial()[j],
                             BigInt(HS.getHilbertQuasiPolynomialDenom().get_mpz_t()));
        tpower = one(QQt);
        for (long i=0; i<len(coeffs); ++i)
        {
            qp[j]  += coeffs[i] * tpower;   // NB tpower = t^i
            tpower *= t;
        }
      }

      return QuasiPoly(qp);
    }

    BigRat multiplicity(const cone& c)
    {
      c->myComputation(libnormaliz::ConeProperty::Multiplicity);
      if (c->mySmallConeIsGood)
        return BigRat(c->myConeLong.getMultiplicity().get_mpq_t());
      else
        return BigRat(c->myConeMPZ.getMultiplicity().get_mpq_t());
    }



    bool IsPointed(const cone& c)
    {
      c->myComputation(libnormaliz::ConeProperty::ExtremeRays);
      if (c->mySmallConeIsGood)
        return c->myConeLong.isPointed();
      else
        return c->myConeMPZ.isPointed();
    }

    bool IsInhomogeneous(const cone& c)
    {
      // is always known
      if (c->mySmallConeIsGood)
        return c->myConeLong.isInhomogeneous();
      else
        return c->myConeMPZ.isInhomogeneous();
    }

    bool IsIntegrallyClosed(const cone& c)
    {
      c->myComputation(libnormaliz::ConeProperty::HilbertBasis);
      if (c->mySmallConeIsGood)
        return c->myConeLong.isIntegrallyClosed();
      else
        return c->myConeMPZ.isIntegrallyClosed();
    }

    bool IsDeg1HilbertBasis(const cone& c)
    {
      c->myComputation(libnormaliz::ConeProperty::HilbertBasis, libnormaliz::ConeProperty::Grading);
      if (c->mySmallConeIsGood)
        return c->myConeLong.isDeg1HilbertBasis();
      else
        return c->myConeMPZ.isDeg1HilbertBasis();
    }

   

    // dimension and rank invariants
    long EmbeddingDim(const cone& c)
    {
      // is always known
      if (c->mySmallConeIsGood)
        return c->myConeLong.getEmbeddingDim();
      else
        return c->myConeMPZ.getEmbeddingDim();
    }

    long rank(const cone& c)
    {
      c->myComputation(libnormaliz::ConeProperty::ExtremeRays);
      if (c->mySmallConeIsGood)
        return c->myConeLong.getRank();
      else
        return c->myConeMPZ.getRank();
    }

    std::vector<BigRat> grading(const cone& c)
    {
      c->myComputation(libnormaliz::ConeProperty::Grading);
      if (c->mySmallConeIsGood)
      {
        return LongVToBigRatV(c->myConeLong.getGrading(), BigInt(c->myConeLong.getGradingDenom()));
      }
      else
      {
        return MPZ_classVToBigRatV(c->myConeMPZ.getGrading(), BigInt(c->myConeMPZ.getGradingDenom().get_mpz_t()));
      }
    }

    // only for inhomogeneous case:
    long RecessionRank(const cone& c)
    {
      c->myComputation(libnormaliz::ConeProperty::RecessionRank);
      if (c->mySmallConeIsGood)
        return c->myConeLong.getRecessionRank();
      else
        return c->myConeMPZ.getRecessionRank();
    }
    long AffineDim(const cone& c)
    {
      c->myComputation(libnormaliz::ConeProperty::ExtremeRays);
      if (c->mySmallConeIsGood)
        return c->myConeLong.getRank();
      else
        return c->myConeMPZ.getRank();
    }
    long ModuleRank(const cone& c)
    {
      c->myComputation(libnormaliz::ConeProperty::ModuleRank);
      if (c->mySmallConeIsGood)
        return c->myConeLong.getModuleRank();
      else
        return c->myConeMPZ.getModuleRank();
    }

    BigInt shift(const cone& c)
    {
      c->myComputation(libnormaliz::ConeProperty::Shift);
      if (c->mySmallConeIsGood)
        return BigInt(c->myConeLong.getShift());
      else
        return BigInt(c->myConeMPZ.getShift().get_mpz_t());
    }
 
    std::vector<BigInt> dehomogenization(const cone& c)
    {
      c->myComputation(libnormaliz::ConeProperty::Dehomogenization);
      if (c->mySmallConeIsGood)
        return LongVToBigIntV(c->myConeLong.getDehomogenization());
      else
        return MPZ_classVToBigIntV(c->myConeMPZ.getDehomogenization());
    }

    // not so important in this library
    // size_t getTriangulationSize() const;
    // Integer getTriangulationDetSum() const;


/**************  applications to monomials  **************/

    namespace // anonymous for file local fns
    {
      PPVector HilbertBasis_PPVector(const PPVector& ppv, libnormaliz::InputType t)
      {
        cone c (PPVectorToVecVecBigInt(ppv), t);
        std::vector<std::vector<BigInt> > res = HilbertBasis(c);
        //create an emptyPPVector with the same PPMonoid
        PPVector ppv_res(PPM(ppv), DMR(ppv));
        VecVecBigIntToPPVector(ppv_res,res);
        return ppv_res;
      }

      //only for modes where normaliz adds a component, like polytope or rees_algebra
      PPVector HilbertBasis_PPVector(const PPVector& ppv, long sym_pos, libnormaliz::Type::InputType t)
      { 
        long dim = NumIndets(PPM(ppv));
        if (sym_pos < 0 || sym_pos >= dim)
          CoCoA_ERROR("Symbol position not valid", "HilbertBasis_PPVector");

        vector< vector<BigInt> > vv_input = PPVectorToVecVecBigInt(ppv);
        if (sym_pos != dim-1)
          CoCoA_ERROR("Symbol needs to be the last symbol of the PPMonoid (other implementations are missing)",
                      "HilbertBasis_PPVector (with extra symbol)");
        for (long i=0; i < len(vv_input); ++i) 
        {
          if (vv_input[i][sym_pos] != 0)
            CoCoA_ERROR("Symbol must not be used in the power products", "HilbertBasis_PPVector (with extra symbol)");
          vv_input[i].erase(vv_input[i].begin()+sym_pos);
        }
        if (t != libnormaliz::Type::polytope && t != libnormaliz::Type::rees_algebra)
          CoCoA_ERROR("Invalid InputType for this method", "HilbertBasis_PPVector (with extra symbol)");

        
        cone c (vv_input, t);
        vector<vector<BigInt> > res = HilbertBasis(c);
        if (sym_pos != dim-1)
        {
          for (long i=0; i < dim; ++i) 
          {
            swap(res[i][dim],res[i][sym_pos]);
          }
        }
        //create an emptyPPVector with the same PPMonoid
        PPVector ppv_res(PPM(ppv), DMR(ppv));
        VecVecBigIntToPPVector(ppv_res,res);

        return ppv_res;
      }

      PPVector HilbertBasis_PPVector(const std::vector<std::vector<BigInt> >& Mat, libnormaliz::InputType t, PPMonoid ppm)
      {
        cone c (Mat, t);
        std::vector<std::vector<BigInt> > res = HilbertBasis(c);
        //create an empty PPVector with the same PPMonoid
        PPVector ppv_res(ppm,NewDivMaskNull());
        VecVecBigIntToPPVector(ppv_res,res);
        return ppv_res;
      }

      PPVector HilbertBasis_PPVector(const cone& C, PPMonoid ppm)
      {
        std::vector<std::vector<BigInt> > res = HilbertBasis(C);
        //create an empty PPVector with the same PPMonoid
        PPVector ppv_res(ppm,NewDivMaskNull());
        VecVecBigIntToPPVector(ppv_res,res);
        return ppv_res;
      }
    } // end anonymous namespace

    PPVector NormalToricRing(const PPVector& ppv)
    { 
      return HilbertBasis_PPVector(ppv, libnormaliz::Type::normalization);
    }

    PPVector IntClosureToricRing(const PPVector& ppv)
    {
      return HilbertBasis_PPVector(ppv, libnormaliz::Type::integral_closure);
    }

    PPVector IntClosureMonIdeal(const PPVector& ppv)
    {
      cone c (PPVectorToVecVecBigInt(ppv), libnormaliz::Type::rees_algebra);
      vector<vector<BigInt> > res = HilbertBasis(c);
      vector<vector<BigInt> > res1;
      for (vector<vector<BigInt> >::const_iterator it=res.begin(); it!=res.end(); ++it)
      {
        //BigInt tmp = it->back();
        if (it->back() == 1 ) {
          res1.push_back(*it);
          res1.back().pop_back();
        }
      }

      //create an empty PPVector with the same PPMonoid
      PPVector ppv_res1(PPM(ppv), DMR(ppv));
      VecVecBigIntToPPVector(ppv_res1,res1);
      return ppv_res1;
    }

    /* If you want also the whole normalization of the rees algebra then the
     * ring must have an unused symbol that we can use for homogenization.
     * You have to specify its position.
     */
    PPVector IntClosureMonIdeal(const PPVector& ppv, long sym_pos)
    {
      return HilbertBasis_PPVector(ppv, sym_pos, libnormaliz::Type::rees_algebra);
    }
      
    PPVector EhrhartRing(const PPVector& ppv, long sym_pos)
    {
      return HilbertBasis_PPVector(ppv, sym_pos, libnormaliz::Type::polytope);
    }

/**************  torus invariants and valuation rings  **************/
    PPVector TorusInvariants(const vector< vector<BigInt> >& T, const PPMonoid& ppm)
    {
      if (T.empty())
      {
        CoCoA_ERROR("Matrix should be non-empty", "TorusInvariants");
      }
      if (NumIndets(ppm) != len(T[0]))
      {
        CoCoA_ERROR("Number of columns in matrix does not match number of variables in ring", "TorusInvariants");
      }
      return HilbertBasis_PPVector(T, libnormaliz::Type::equations, ppm);
    }
    
    PPVector FiniteDiagInvariants(const vector< vector<BigInt> >& Cong, const PPMonoid& ppm)
    {
      if (Cong.empty())
      {
        CoCoA_ERROR("Matrix should be non-empty", "FiniteDiagInvariants");
      }
      if (NumIndets(ppm) != len(Cong[0])-1)
      {
        CoCoA_ERROR("Number of columns in matrix -1 does not match number of variables in ring", "FiniteDiagInvariants");
      }
      return HilbertBasis_PPVector(Cong, libnormaliz::Type::congruences, ppm);
    }

    PPVector DiagInvariants(const vector< vector<BigInt> >& T, const vector< vector<BigInt> >& Cong, const PPMonoid& ppm)
    {
      if (T.empty() && Cong.empty())
      {
        CoCoA_ERROR("At least one Matrix should be non-empty", "DiagInvariants");
      }
      if (T.empty())
        return FiniteDiagInvariants(Cong, ppm);
      if (Cong.empty())
        return TorusInvariants(T, ppm);

      if (NumIndets(ppm) != len(T[0]))
      {
        CoCoA_ERROR("Number of columns in matrix does not match number of variables in ring", "DiagInvariants");
      }
      if (NumIndets(ppm) != len(Cong[0])-1)
      {
        CoCoA_ERROR("Number of columns in matrix -1 does not match number of variables in ring", "DiagInvariants");
      }
      std::map< libnormaliz::InputType, std::vector<std::vector<BigInt> > > cone_input;
      cone_input[libnormaliz::Type::equations] = T;
      cone_input[libnormaliz::Type::congruences] = Cong;
      cone C(cone_input);
      return HilbertBasis_PPVector(C, ppm);
    }


    PPVector IntersectionValRings (const vector< vector<BigInt> >& V, const PPMonoid& ppm)
    {
      long dim = len(V[0]);
      std::map< libnormaliz::InputType, std::vector<std::vector<BigInt> > > cone_input;
      cone_input[libnormaliz::Type::inequalities] = V;
      vector< vector<BigInt> > positive_signs = vector< vector<BigInt> >(1,vector<BigInt>(dim,BigInt(1)));
      cone_input[libnormaliz::Type::signs] = positive_signs;
      cone C(cone_input);
      return HilbertBasis_PPVector(C, ppm);
    }


/********************************************************
 ***               conversion functions               ***
 ********************************************************/

    namespace  // conversion functions
    {
      std::vector<BigInt> LongVToBigIntV(const std::vector<long>& VIn)
      {
        std::vector<BigInt> v;
        v.reserve(VIn.size());
        for (vector<long>::const_iterator it=VIn.begin(); it!=VIn.end(); ++it)
          v.push_back(BigInt(*it));
        return v;
      }

      std::vector<BigInt> MPZ_classVToBigIntV(const std::vector<mpz_class>& VIn)
      {
        std::vector<BigInt> v;
        v.reserve(VIn.size());
        for (vector<mpz_class>::const_iterator it=VIn.begin(); it!=VIn.end(); ++it)
          v.push_back(BigInt(it->get_mpz_t()));
        return v;
      }

      std::vector<BigRat> LongVToBigRatV(const std::vector<long>& VIn, const BigInt& denom)
      {
        std::vector<BigRat> v;
        v.reserve(VIn.size());
        for (vector<long>::const_iterator it=VIn.begin(); it!=VIn.end(); ++it)
          v.push_back(BigRat(BigInt(*it),denom));
        return v;
      }

      std::vector<BigRat> MPZ_classVToBigRatV(const std::vector<mpz_class>& VIn, const BigInt& denom)
      {
        std::vector<BigRat> v;
        v.reserve(VIn.size());
        for (vector<mpz_class>::const_iterator it=VIn.begin(); it!=VIn.end(); ++it)
          v.push_back(BigRat(BigInt(it->get_mpz_t()),denom));
        return v;
      }


      std::vector<long> BigIntVToLongV(const std::vector<BigInt>& VIn)
      {
        std::vector<long> v;
        long n;
        v.reserve(VIn.size());
        for (vector<BigInt>::const_iterator it=VIn.begin(); it!=VIn.end(); ++it)
          if (IsConvertible(n, *it))  v.push_back(n);
          else CoCoA_ERROR(ERR::BadConvert, "normaliz: BigIntVToLongV");
        return v;
      }


      std::vector<mpz_class> BigIntVToMPZ_classV(const std::vector<BigInt>& VIn)
      {
        std::vector<mpz_class> v;
        v.reserve(VIn.size());
        for (vector<BigInt>::const_iterator it=VIn.begin(); it!=VIn.end(); ++it)
          v.push_back(mpz_class(mpzref(*it)));
        return v;
      }


      void ConvertFromNormaliz(std::vector<std::vector<BigInt> >& v, const std::vector<std::vector<long> >& l)
      {
        v.clear();      // not exception safe
        v.reserve(l.size());
        //      transform(l.begin(), l.end(), v.begin(), LongVToBigIntV);
        for (vector<vector<long> >::const_iterator it=l.begin(); it!=l.end(); ++it)
          v.push_back(LongVToBigIntV(*it));
      }    


      void ConvertFromNormaliz(std::vector<std::vector<BigInt> >& v, const std::vector<std::vector<mpz_class> >& l)
      {
        v.clear();      // not exception safe
        v.reserve(l.size());
        //      transform(l.begin(), l.end(), v.begin(), LongVToBigIntV);
        for (vector<vector<mpz_class> >::const_iterator it=l.begin(); it!=l.end(); ++it)
          v.push_back(MPZ_classVToBigIntV(*it));
      }    


      std::vector<std::vector<long> >  ReturnVecVecLong(const std::vector<std::vector<BigInt> >& v)
      {
        std::vector<std::vector<long> > l;
        //      transform(v.begin(), v.end(), l.begin(), BigIntVToLongV);
        try
        {
            for (vector<vector<BigInt> >::const_iterator it=v.begin(); it!=v.end(); ++it)
              l.push_back(BigIntVToLongV(*it));
            //      std::clog << l << std::endl;
        } catch (...) {l.clear();}
        return l;
      }


      map< libnormaliz::InputType, vector<vector<long> > >  ReturnMapVecVecLong(const map< libnormaliz::InputType, vector<vector<BigInt> > > m)
      {
        map< libnormaliz::InputType, vector<vector<long> > > ret;
        for (map<libnormaliz::InputType, vector<vector<BigInt> > >::const_iterator it=m.begin(); it!=m.end(); ++it)
          ret.insert(make_pair((*it).first, ReturnVecVecLong((*it).second)));
        return ret;
      }
      
      
      std::vector<std::vector<mpz_class> >  ReturnVecVecMPZ_class(const std::vector<std::vector<BigInt> >& v)
      {
        std::vector<std::vector<mpz_class> > l;
        //      transform(v.begin(), v.end(), l.begin(), BigIntVToLongV);
        for (vector<vector<BigInt> >::const_iterator it=v.begin(); it!=v.end(); ++it)
          l.push_back(BigIntVToMPZ_classV(*it));
        //      std::clog << l << std::endl;
        return l;
      }


      map< libnormaliz::InputType, vector<vector<mpz_class> > >  ReturnMapVecVecMPZ_class(const map< libnormaliz::InputType, vector<vector<BigInt> > > m)
      {
        map< libnormaliz::InputType, vector<vector<mpz_class> > > ret;
        for (map<libnormaliz::InputType, vector<vector<BigInt> > >::const_iterator it=m.begin(); it!=m.end(); ++it)
          ret.insert(make_pair((*it).first, ReturnVecVecMPZ_class((*it).second)));
        return ret;
      }


      std::vector<std::vector<BigInt> > PPVectorToVecVecBigInt(const PPVector& ppv)
      {
        vector<BigInt> tmp;
        const long n =  len(ppv);
        vector<vector<BigInt> > v(n);
        for (long i=0; i<n; ++i)
        {
          BigExponents(tmp,PP(ppv[i]));
          v[i]=tmp;
        }
        return v;
      }

      void VecVecBigIntToPPVector(PPVector& ppv, const std::vector<std::vector<BigInt> >& M)
      {
        PPMonoid ppm = PPM(ppv);
        const long n =  len(M);
        for (long i=0; i < n; ++i)
        {
          ppv.myPushBack(PPMonoidElem(ppm, M[i]));
        }
      }

    } // end of anonymous namespace

    std::vector<std::vector<BigInt> > MatrixToVecVecBigInt(ConstMatrixView M)
    {
      const ErrorInfo ErrMesg("Matrix entries must be integer", "MatrixToVecVecBigInt (cone ctor)");

      vector<vector<BigInt> > v;
      for (long i=0; i<NumRows(M); ++i)
      {
        v.push_back(vector<BigInt>());
        for (long j=0; j<NumCols(M); ++j)
          v[i].push_back(ConvertTo<BigInt>(M(i,j), ErrMesg));
      }
      return v;
    }

    PPVector MonomialsToPPV(const std::vector<RingElem>& v)
    {
      if(!AreMonomials(v)) {
        CoCoA_ERROR("Expected list of monomials","MonomialsToPPV");
      }
      if(v.empty()) {
        CoCoA_ERROR("List of monomials has to be non-empty","MonomialsToPPV");
      }
      //convert it to a PPVector
      PPVector ppv(PPM(owner(v[0])), NewDivMaskNull());
      convert(ppv,v);
      return ppv;
    }


  } // namespace Normaliz
} // namespace CoCoA

// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/ExternalLibs-Normaliz.C,v 1.29 2014/07/31 14:45:17 abbott Exp $
// $Log: ExternalLibs-Normaliz.C,v $
// Revision 1.29  2014/07/31 14:45:17  abbott
// Summary: Merged io.H and UtilsTemplate.H into new header VectorOperations.H
// Author: JAA
//
// Revision 1.28  2014/07/14 10:02:56  abbott
// Summary: Christof has added some fns which use quasi polys
// Author: JAA
//
// Revision 1.27  2014/07/11 10:12:51  abbott
// Summary: Christof added printing for cones
// Author: JAA
//
// Revision 1.26  2014/07/07 12:13:53  abbott
// Summary: Christof's updates (& removed AsSparsePolyRing)
// Author: JAA
//
// Revision 1.25  2014/07/01 15:31:21  bigatti
// -- fix by Christof Soeger
//
// Revision 1.24  2014/06/17 10:06:44  abbott
// Summary: Removed pointless AsPolyRing
// Author: JAA
//
// Revision 1.23  2014/05/12 14:33:51  bigatti
// -- updated from Christof Soeger
//
// Revision 1.22  2014/05/09 14:56:30  bigatti
// -- new fn by Christof Soeger
//
// Revision 1.21  2014/01/30 09:56:49  bigatti
// -- removed DOS newlines
//
// Revision 1.20  2014/01/29 17:38:35  bigatti
// -- added HilbertSeries (by Christof Soeger)
// -- added multiplicity (by Christof Soeger)
//
// Revision 1.19  2014/01/28 09:44:42  abbott
// Tidier impl of MatrixToVecVecBigInt using new ConvertTo syntax.
//
// Revision 1.18  2013/07/12 14:53:32  abbott
// Improved indentation.
//
// Revision 1.17  2013/03/15 17:49:17  abbott
// Minor cleaning; removed an out-of-date TODO comment.
//
// Revision 1.16  2012/10/08 13:52:02  bigatti
// -- more cleaning and updates by Christof Soeger
//
// Revision 1.15  2012/10/05 10:17:04  bigatti
// by Christof Soeger
// * Made the NewCone functions to constructors of cone.
// * Introduced some abbreviation in the method names:
//  IntegralClosure -> IntClosure
//  MonomialIdeal -> MonIdeal
//  Normaliz -> Nmz  (the prefix for CoCoA5 functions)
// * New function IntClosureMonIdeal
//
// Revision 1.14  2012/09/28 13:56:36  abbott
// Cleaned up code with Christof's help -- now cleaner & more symmetrical.
// Also fixed a subtle bug where mySmallConeIsGood was not updated as it should have been.
//
// Revision 1.13  2012/08/03 16:31:22  bigatti
// -- changed: procedural --> functional (by C.Soeger)
//
// Revision 1.12  2012/07/25 12:46:27  bigatti
// -- merged with C.Soeger changes for NormalizComputation (init from map, etc)
//
// Revision 1.11  2012/07/19 17:12:05  abbott
// Added NewCone -- unified pseudo-ctor so user does not have to choose between long and BigInt.
//
// Revision 1.10  2012/06/19 14:44:41  bigatti
// -- changed Ht1 --> Deg1, changed def of HVector (by C.Soeger)
//
// Revision 1.9  2011/11/07 11:09:32  bigatti
// -- new ctors taking ConstRefRingElem
// -- added PPVector NormalToricRing(const PPVector& ppv)
// -- removed void HilbertBasis(std::vector<std::vector<BigInt> >& v, const cone& c)
//
// Revision 1.8  2011/10/04 13:03:02  bigatti
// -- new logo for gui
//
// Revision 1.7  2011/09/30 12:55:44  bigatti
// -- introduced namespace "Normaliz" and removed Normaliz from function names
// -- input of Normaliz functions in CoCoA-5 is now a matrix instead of
//    vector<vector<BigInt>>
//
// Revision 1.6  2011/08/14 15:52:17  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.5  2011/07/20 15:31:21  bigatti
// -- fixed InputType in (undefined) functions
//
// Revision 1.4  2011/07/20 13:49:37  bigatti
// -- added "Normaliz" postfix to Normaliz function calls
//
// Revision 1.3  2011/07/20 12:45:12  bigatti
// -- new normaliz interface (not yet public)
//
// Revision 1.2  2011/02/17 16:50:04  bigatti
// -- getting ready for new official veson on Normaliz: added HVector, removed Triangulation
//
// Revision 1.1  2010/11/05 14:27:18  bigatti
// -- was TmpNormaliz**
//
// Revision 1.6  2010/10/12 11:22:54  bigatti
// -- TmpNormaliz.H simplified:
//    now cone is a smart pointer to ConeBase
//    and the concrete classes are entirely in the .C file
// -- added Ht1Elements, SupportHyperplanes, Triangulation
// -- added some text in ex-Normaliz1 and 2
//
// Revision 1.5  2010/10/08 13:40:37  bigatti
// -- moved inclusion of libnormaliz.cpp and instantiation of types
//    into TmpNormalizTypes.C
//
// Revision 1.4  2010/10/08 10:39:32  bigatti
// -- extended interface for normaliz
//

#endif
