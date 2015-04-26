//   Copyright (c)  2014  John Abbott

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


#include "CoCoA/ProgressReporter.H"
#include "CoCoA/time.H"

#include <cmath>
using std::floor;
#include <iostream>
using std::cout;
using std::endl;

namespace CoCoA
{

  ProgressReporter::ProgressReporter(double interval):
      myIntervalCount(0),
      myCheckCount(1),
      myTotalCount(0),
      myLastCheckTime(CpuTime()),
      myTargetInterval(interval),
      myNextPrintTime(myLastCheckTime+interval)
  {}


  bool ProgressReporter::myIsTimeToPrint()
  {
    myTotalCount += myIntervalCount;
    myIntervalCount = 0;
    const double t = CpuTime();
    const double ratio = (t-myLastCheckTime)/myTargetInterval;
    myLastCheckTime = t;
    if (8*ratio >= 1) { decrease125(myCheckCount); }
    if (32*ratio < 1)
    {
      increase125(myCheckCount);
      // Now adjust myIntervalCount, so next print will be at a multiple of myCheckCount
      const long tmp = myTotalCount%myCheckCount;
      if (tmp != 0)
      {
        myTotalCount += myCheckCount-tmp;
        myIntervalCount = tmp-myCheckCount;
      }
      return false;
    }
    if (t < myNextPrintTime) return false;
    myNextPrintTime = t+myTargetInterval;
    return true;
  }

  void ProgressReporter::myPrintReport()
  {
    cout << "--> Progress count=" << myTotalCount << "   time=" << myLastCheckTime << endl;
  }

  void ProgressReporter::myPrintReport(long arg1)
  {
    cout << "--> Progress at "<< arg1 << "   count=" << myTotalCount << "   time=" << myLastCheckTime << endl;
  }

  void ProgressReporter::myPrintReport(long arg1, long arg2)
  {
    cout << "--> Progress at (" << arg1 << ", " << arg2 << ")   count=" << myTotalCount << "   time=" << myLastCheckTime<< endl;
  }


  // Quick makeshift impl.
  std::ostream& operator<<(std::ostream& out, const ProgressReporter& PR)
  {
    out << "ProgressReporter(intvl=" << PR.myTargetInterval << ")";
    return out;
  }


  void increase125(long& n)
  {
    long pwr10 = 1;
    while (n >= 10) { n /= 10; pwr10 *= 10; }
    if (n == 1) { n = 2*pwr10; return; }
    if (n == 2) { n = 5*pwr10; return; }
    n = 10*pwr10;
  }

  void decrease125(long& n)
  {
    if (n == 1) return;
    long pwr10 = 1;
    while (n >= 10) { n /= 10; pwr10 *= 10; }
    if (n == 1) { n = pwr10/2; return; }
    if (n == 2) { n = pwr10; return; }
    n = 2*pwr10;
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/ProgressReporter.C,v 1.2 2014/04/22 14:53:43 abbott Exp $
// $Log: ProgressReporter.C,v $
// Revision 1.2  2014/04/22 14:53:43  abbott
// Summary: Changed format of progress reports; now respects the chosen interval (rather than aiming for multiples of chosen interval)
// Author: JAA
//
// Revision 1.1  2014/04/22 13:26:01  abbott
// Summary: New class ProgressReporter
// Author: JAA
//
//
