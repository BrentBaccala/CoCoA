//   Copyright (c)  2007,2010,2011  John Abbott

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


#include "CoCoA/GlobalManager.H"
#include "CoCoA/FractionField.H"
#include "CoCoA/MemPool.H"
#include "CoCoA/DenseUPolyRing.H" // for Hilbert-Poincare' series
#include "CoCoA/utils.H"

#include "CoCoA/error.H"
#include "TmpHilbertDir/TmpPoincareCPP.H" // for Hilbert-Poincare' series

#include "gmp.h"

#include <algorithm>
using std::min;
using std::max;
#include <cstdlib>
using std::malloc;
using std::realloc;
using std::free;
#include <iostream>
// using std::cerr & std::endl for serious warning in GlobalManager dtor
#include <cstring>
using std::memcpy;


// These 3 fns must have C linkage to work with GMP's mem mgr setter.
extern "C"
{
  void* CoCoA_GMP_alloc(size_t sz);
  void* CoCoA_GMP_realloc(void* ptr, size_t oldsz, size_t newsz);
  void CoCoA_GMP_free(void* ptr, size_t sz);
}

void* CoCoA_GMP_alloc(size_t sz)
{
  CoCoA_ASSERT(CoCoA::GlobalGMPSliceSize() != 0);
  if (sz <= CoCoA::GlobalGMPSliceSize())
    return CoCoA::GlobalGMPPoolPtr()->alloc();
  return malloc(sz);
}

void* CoCoA_GMP_realloc(void* ptr, size_t oldsz, size_t newsz)
{
  CoCoA_ASSERT(CoCoA::GlobalGMPSliceSize() != 0);
  if (oldsz <= CoCoA::GlobalGMPSliceSize() &&
      newsz <= CoCoA::GlobalGMPSliceSize())
    return ptr;

  if (oldsz > CoCoA::GlobalGMPSliceSize() &&
      newsz > CoCoA::GlobalGMPSliceSize())
    return realloc(ptr, newsz);

  const size_t n = min(oldsz, newsz);
  void* dest = CoCoA_GMP_alloc(newsz);
  memcpy(dest, ptr, n);
  CoCoA_GMP_free(ptr, oldsz);
  return dest;
}

void CoCoA_GMP_free(void* ptr, size_t sz)
{
  CoCoA_ASSERT(CoCoA::GlobalGMPSliceSize() != 0);
  if (sz <= CoCoA::GlobalGMPSliceSize())
    CoCoA::GlobalGMPPoolPtr()->free(ptr);
  else
    free(ptr);
}



namespace CoCoA
{

  // Pseudo-ctors for RingZZ and RingQ.
  ring MakeUniqueInstanceOfRingZZ(); // Defined in RingZZ.C.
  FractionField MakeUniqueInstanceOfRingQQ(const ring&); // Defined in RingQ.C.

  // Checking fns to be called immediately before calling dtors for RingQ and RingZZ
  bool RingZZStillInUse(const ring& ZZ);  // Defined in RingZZ.C
  bool RingQQStillInUse(const FractionField& Q);  // Defined in RingQ.C


  // The static members of GlobalManager -- effectively global variables.
  const GlobalManager* GlobalManager::ourGlobalDataPtr = 0;
  std::size_t GlobalManager::GMPSliceSize = 0; // size in bytes of slices in the MemPool (compile-time constant)
  MemPool* GlobalManager::GMPPoolPtr = 0;
  std::size_t GlobalManager::ourHPMaxPower = 100;  // for Hilbert-Poincare' series

  const GlobalManager* GlobalManager::ptr(const char* FnName)
  {
    if (GlobalManager::ourGlobalDataPtr == 0)
      CoCoA_ERROR(ERR::NoGlobalMgr, FnName);
    return GlobalManager::ourGlobalDataPtr;
  }


  GlobalManager::GMPMemMgr::GMPMemMgr(GlobalSettings::AllocatorSetting choice, std::size_t SliceSize)
  {
    if (choice == GlobalSettings::SystemAllocator) return;

    myPoolPtr.reset(new MemPool(SliceSize, "Global GMP MemPool")); // must do this first to be exception safe
    GlobalManager::GMPPoolPtr = myPoolPtr.get();
    GlobalManager::GMPSliceSize = GlobalManager::GMPPoolPtr->mySliceSize();
    mp_get_memory_functions(&myPrevAlloc, &myPrevRealloc, &myPrevFree);
    mp_set_memory_functions(&CoCoA_GMP_alloc, &CoCoA_GMP_realloc, &CoCoA_GMP_free);
  }


  GlobalManager::GMPMemMgr::~GMPMemMgr()
  {
    if (myPoolPtr.get() == 0) return;

    mp_set_memory_functions(myPrevAlloc, myPrevRealloc, myPrevFree);
    GlobalManager::GMPSliceSize = 0;
    GlobalManager::GMPPoolPtr = 0;
  }


  // ----------------------------------------------------------------------

  const std::size_t GlobalSettings::ourDefaultSliceSize = 2*sizeof(long);
  const GlobalSettings::ResidueSetting GlobalSettings::ourDefaultResidueSetting = GlobalSettings::SymmResidues;
  const GlobalSettings::AllocatorSetting GlobalSettings::ourDefaultAllocatorSetting = GlobalSettings::SystemAllocator;


  GlobalSettings::GlobalSettings():
      myResidueSettingHasBeenSet(false),
      myAllocatorSettingHasBeenSet(false),
      mySliceSizeHasBeenSet(false),
      myResidueSetting(ourDefaultResidueSetting),
      myAllocatorSetting(ourDefaultAllocatorSetting),
      mySliceSize(ourDefaultSliceSize)
  {}

  GlobalSettings& GlobalSettings::mySetResidueSetting(ResidueSetting r)
  {
    CoCoA_ASSERT(!myResidueSettingHasBeenSet);
    myResidueSettingHasBeenSet = true;
    myResidueSetting = r;
    return *this;
  }

  GlobalSettings& GlobalSettings::mySetAllocatorSetting(AllocatorSetting a)
  {
    CoCoA_ASSERT(!myAllocatorSettingHasBeenSet);
    myAllocatorSettingHasBeenSet = true;
    myAllocatorSetting = a;
    return *this;
  }

  GlobalSettings& GlobalSettings::mySetSliceSize(std::size_t SliceSize)
  {
    CoCoA_ASSERT(!mySliceSizeHasBeenSet);
    mySliceSizeHasBeenSet = true;
    mySliceSize = SliceSize;
    return *this;
  }



  GlobalSettings GlobalSettings::operator()(std::size_t SliceSize) const
  {
    CoCoA_ASSERT(!mySliceSizeHasBeenSet && myAllocatorSetting != SystemAllocator);
    GlobalSettings ans(*this);
    return ans.mySetSliceSize(SliceSize);
  }


  GlobalSettings operator+(const GlobalSettings& arg1, const GlobalSettings& arg2)
  {
    GlobalSettings ans;
    if (arg1.myResidueSettingHasBeenSet && arg2.myResidueSettingHasBeenSet)
      CoCoA_ERROR(ERR::BadGlobalSettings, "residue setting");
    if (arg1.myResidueSettingHasBeenSet) ans.mySetResidueSetting(arg1.myResidueSetting);
    if (arg2.myResidueSettingHasBeenSet) ans.mySetResidueSetting(arg2.myResidueSetting);

    if (arg1.myAllocatorSettingHasBeenSet && arg2.myAllocatorSettingHasBeenSet)
      CoCoA_ERROR(ERR::BadGlobalSettings, "allocator setting");
    if (arg1.myAllocatorSettingHasBeenSet) ans.mySetAllocatorSetting(arg1.myAllocatorSetting);
    if (arg2.myAllocatorSettingHasBeenSet) ans.mySetAllocatorSetting(arg2.myAllocatorSetting);

    if (arg1.mySliceSizeHasBeenSet && arg2.mySliceSizeHasBeenSet)
      CoCoA_ERROR(ERR::BadGlobalSettings, "GMPAllocator slice size");
    if (arg1.mySliceSizeHasBeenSet) ans.mySetSliceSize(arg1.mySliceSize);
    if (arg2.mySliceSizeHasBeenSet) ans.mySetSliceSize(arg2.mySliceSize);

    return ans;
  }


  const GlobalSettings UseSymmResidues(GlobalSettings().mySetResidueSetting(GlobalSettings::SymmResidues));
  const GlobalSettings UseNonNegResidues(GlobalSettings().mySetResidueSetting(GlobalSettings::NonNegResidues));
  const GlobalSettings UseSystemAllocatorForGMP(GlobalSettings().mySetAllocatorSetting(GlobalSettings::SystemAllocator));
  const GlobalSettings UseGMPAllocator(GlobalSettings().mySetAllocatorSetting(GlobalSettings::GMPAllocator));


  // ----------------------------------------------------------------------

  GlobalManager::ZZQQMgr::ZZQQMgr():
      myRingZZ(MakeUniqueInstanceOfRingZZ()),
      myRingQQ(MakeUniqueInstanceOfRingQQ(myRingZZ))
  {}

  GlobalManager::ZZQQMgr::~ZZQQMgr()
  {
    if (RingZZStillInUse(myRingZZ) || RingQQStillInUse(myRingQQ))
      std::cerr << std::endl
                << "*************************" << std::endl
                << "*** IMMINENT DISASTER ***" << std::endl
                << "*************************" << std::endl
                << std::endl
                << ">>>  CoCoA::GlobalManager being destroyed while CoCoA objects still live!  <<<" << std::endl
                << std::endl;
  }


  // ----------------------------------------------------------------------

  GlobalManager::GlobalManager(const GlobalSettings& settings):
      myResidueSetting(settings.myResidueSetting),
      myGMPMemMgr(settings.myAllocatorSetting, settings.mySliceSize),
      myZZQQMgr()
  {
// !!!***NOT THREAD SAFE***!!!  Must make next 3 lines atomic.
    // Complain if a GlobalManager object has already been created
    if (ourGlobalDataPtr != 0)
      CoCoA_ERROR(ERR::GlobalManager2, "GlobalManager ctor");
    ourGlobalDataPtr = this;
  }


  GlobalManager::~GlobalManager()
  {
// !!!***NOT THREAD SAFE***!!!  Must make check on ourGlobalDataPtr exclusive!
    if (ourGlobalDataPtr == 0)
      CoCoA_ERROR(ERR::SERIOUS, "GlobalManager dtor");
    // Delete registered globals in reverse order
    for (long i=len(myDtorList)-1; i >= 0; --i)
      myDtorList[i].myDtor(myDtorList[i].myPtr);
    ourGlobalDataPtr = 0; // "deregister" the global data
  }


//   bool DefaultResiduesAreSymm()
//   {
//     return GlobalManager::ptr("DefaultResiduesAreSymm")->myResiduesAreSymm;
//   }

  GlobalSettings::ResidueSetting DefaultResidueSetting()
  {
    return GlobalManager::ptr("DefaultResidueSetting")->myResidueSetting;
  }

  
  GlobalManager::PtrDtor::PtrDtor(void* ptr, void (*dtor)(void*)): myPtr(ptr), myDtor(dtor) {}

  void RegisterDtorForGlobal(void* ptr, void (*dtor)(void*))
  {
    GlobalManager::ptr("RegisterDtorForGlobal")->myDtorList.push_back(GlobalManager::PtrDtor(ptr,dtor));
  }


  //----------------------------------------------------------------------
  // pre-computed power list for univariate Hilbert-Poincare Series
  void MakeGlobalHPPowerList(const DenseUPolyRing& P)
  {  
    MakeHPPowerList(GlobalManager::ptr("MakeGlobalHPPowerList")->myHPPowerList,
                    P,
                    GlobalManager::ourHPMaxPower);
  }


  int HPPowerListMaxDeg()
  {  
    return GlobalManager::ptr("HPPowerListMaxDeg")->ourHPMaxPower;
  }

  
  ConstRefRingElem HPPowerList(int exp)
  {
    if (exp>HPPowerListMaxDeg())
      CoCoA_ERROR(ERR::ArgTooBig, "HPPowerList");
    return GlobalManager::ptr("HPPowerList")->myHPPowerList[exp];
  }
  
    
  void CopyHPPower(RingElem& res, int exp)
  {
    if (exp<=HPPowerListMaxDeg())
      res = HPPowerList(exp);
    else
    {
      res = HPPowerList(HPPowerListMaxDeg());
      const DenseUPolyRing HSRing = owner(res);
      for (long i=HPPowerListMaxDeg(); i<exp; ++i)
        HSRing->myMulBy1MinusXExp(raw(res), 1);
    }    
  }
  
  RandomSource& GlobalRandomSource()
  {
    return GlobalManager::ptr("GlobalRandomSource")->myRandomSource;
  }

} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/GlobalManager.C,v 1.18 2014/07/09 13:01:17 abbott Exp $
// $Log: GlobalManager.C,v $
// Revision 1.18  2014/07/09 13:01:17  abbott
// Summary: Removed AsDenseUPolyRing
// Author: JAA
//
// Revision 1.17  2014/07/01 12:40:53  bigatti
// -- added CopyHPPower and argument check in HPPowerList
//
// Revision 1.16  2014/04/30 16:07:40  abbott
// Summary: Replaced X.size() by len(X)
// Author: JAA
//
// Revision 1.15  2013/06/17 08:54:02  abbott
// Added RegisterDtorForGlobal.
//
// Revision 1.14  2012/10/15 12:35:24  abbott
// Added  std::  prefix.
//
// Revision 1.13  2012/02/08 13:47:16  bigatti
// -- changed Z,Q --> ZZ,QQ
//
// Revision 1.12  2011/05/19 13:54:48  abbott
// Replaced DefaultResiduesAreSymm by DefaultResidueSetting.
//
// Revision 1.11  2011/05/03 10:03:32  abbott
// Added GlobalRandomSource.
// Internally added GlobalManager::ptr to allow neater implementations.
//
// Revision 1.10  2010/11/17 15:52:33  abbott
// Removed out-of-date include of GMPAllocator.H.
//
// Revision 1.9  2010/11/11 17:45:08  abbott
// Moved GMPMemMgr so that it is a nested class inside GlobalManager.
//
// Revision 1.8  2010/10/29 12:06:41  bigatti
// -- added globals for Hilbert-Poincare' series
//
// Revision 1.7  2010/10/27 20:58:45  abbott
// Major reorganization of GlobalManager and GMPAllocator.
//
// Revision 1.6  2010/10/22 14:03:04  abbott
// Major change to GMPAllocator -- it is now set/activated by the GlobalManager.
// This is a Friday afternoon check-in... hope to check in cleaner code in the
// next few days.
//
// Revision 1.5  2010/09/30 14:28:23  abbott
// Replaced auto_ptrs to RingZ and RingQ by direct values; ctor changed accordingly.
//
// Dtor now checks that ref counts in RingZ and RingQ are correct; if not, a rude
// message is printed on cerr (and the program will probably crash after the
// GlobalManager has been destroyed).
//
// Revision 1.4  2009/12/03 17:40:36  abbott
// Added include directives for ZZ.H (as a consequence of removing
// the directive from ring.H).
//
// Revision 1.3  2009/05/14 09:39:29  abbott
// Added possibility to specify "symmetric" or "non-negative" residues
// in quotients of ZZ.  Affects printing of elements in quotients of ZZ
// (also changed printing of elements in general quotient rings).
// Consequent changes in several tests.
//
// Revision 1.2  2007/10/30 17:14:08  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.2  2007/03/05 21:33:13  cocoa
// Improved/cleaned GlobalManager; added doc too.
//
// Revision 1.1  2007/03/03 14:02:11  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.1  2007/03/02 16:46:28  cocoa
// New foundations object which calls ctors and dtors of global objects.
//
