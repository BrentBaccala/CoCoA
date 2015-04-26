//   Copyright (c)  2005  Massimo Caboara

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


#include <iostream>
using std::ostream;
using std::endl;
#include <utility>
using std::pair;
#include <list>
using std::list;
#include "CoCoA/VectorOperations.H"  // for templates which print lists/vectors (in myStampaPairs)

#include "CoCoA/TmpGRStats.H"
#include "CoCoA/IntOperations.H"


namespace CoCoA
{


/*----------------- class DegStats functions ----------------*/

  DegStats::DegStats(const degree& Deg):
    myDeg(Deg)
  {
    myPInserted=0;
    myGMKilled=0;
    myCopKilled=0;
    myBKilled=0;
    myUseful=0;
    myUseless=0;
    myPairsNo=0;
  }//ctor
  
  DegStats::DegStats(const degree& Deg,
           unsigned int PInserted,
	   unsigned int GMKilled,
           unsigned int CopKilled,
	   unsigned int BKilled,
	   unsigned int Useful,
	   unsigned int Useless,
           unsigned int PairsNo):
    myDeg(Deg)
  {
    myPInserted=PInserted;
    myGMKilled=GMKilled;
    myCopKilled=CopKilled;
    myBKilled=BKilled;
    myUseful=Useful;
    myUseless=Useless;
    myPairsNo=PairsNo;
  }//ctor

/*----------------- class Stats functions ----------------*/


  Stats::Stats(unsigned int NumGens, int level)
  {
    mySetLevel(level);// Sets the flags. Has to be there.
    myNumGens=NumGens;
    myPInserted=0;
    myGMKilled=0;
    myCopKilled=0;
    myBKilled=0;
    myBTouched=0;
    myGMTouched=0;
    myUseful=0;
    myUseless=0;
    myPolyDeleted=0;
    myPolyDHed=0;
    myDegDH=0;
    myNReductions=0;
    myReductionTime=0.0;
    myTotalTime=0.0;
   }

// WARN: copy all the data, level may be out of sync with the flags
 Stats::Stats(unsigned int NumGens,
        int level,
        unsigned int PInserted,
        unsigned int GMKilled,
        unsigned int CopKilled,
        unsigned int Useful,
        unsigned int Useless,
        unsigned int BKilled,
        unsigned int BTouched,
        unsigned int GMTouched,
        unsigned int PolyDeleted,
        unsigned int PolyDHed,
        unsigned int DegDH,
	unsigned int NReductions,
        double ReductionTime,
        double TotalTime,
    	bool PrintReduction,
	bool PrintDeg,
	bool PrintGM,
        bool PrintCop,
	bool PrintBC,
	bool PrintNumPair,
	bool PrintFinal,
	bool PrintFinalFull,
	bool PrintFinalSimple,
	bool PrintNewPairs,
	bool PrintPolyDeleted,
	bool PrintPolyDH,
        bool PrintPolyLen,
        bool PrintKill)
{
      myNumGens=NumGens;
      myLevel=level;
      myPInserted=PInserted;
      myGMKilled=GMKilled;
      myCopKilled=CopKilled;
      myBKilled=BKilled;
      myBTouched=BTouched;
      myGMTouched=GMTouched;
      myUseful=Useful;
      myUseless=Useless;
      myPolyDeleted=PolyDeleted;
      myPolyDHed=PolyDHed;
      myDegDH=DegDH;
      myNReductions=NReductions,
      myReductionTime=ReductionTime;
      myTotalTime=TotalTime;
      myPrintReduction=PrintReduction;
      myPrintDeg=PrintDeg;
      myPrintGM=PrintGM;
      myPrintCop=PrintCop;
      myPrintBC=PrintBC;
      myPrintNumPair=PrintNumPair;
      myPrintFinal=PrintFinal;
      myPrintFinalFull=PrintFinalFull;
      myPrintFinalSimple=PrintFinalSimple;
      myPrintNewPairs=PrintNewPairs;
      myPrintPolyDeleted=PrintPolyDeleted;
      myPrintPolyDH=PrintPolyDH;
      myPrintPolyLen=PrintPolyLen;
      myPrintKill=PrintKill;
}

  void Stats::mySetLevel(int level)
  {
    myLevel=level;
    myPrintReduction=false;
    myPrintDeg=false;
    myPrintCop=false;
    myPrintGM=false;
    myPrintBC=false;
    myPrintNumPair=false;
    myPrintFinal=false;
    myPrintFinalFull=false;
    myPrintFinalSimple=false;
    myPrintNewPairs=false;
    myPrintPolyDeleted=false;
    myPrintPolyDH=true;
    myPrintPolyLen=false;
    myPrintKill=false;
    switch (level) {
        case 7:// Reductors and pairs printing - only for real hard nosed debuggers
        // TODO
        case 6:// poly reduction history
          myPrintPolyLen=true;
          myPrintFinalFull=true;
        case 5://
          myPrintGM=true;
          myPrintBC=true;
          myPrintPolyDeleted=true;
          myPrintPolyDH=true;
          myPrintNewPairs=true;
          myPrintCop=true;
        case 4:// pair by pair stats
          myPrintKill=true;
          myPrintReduction=true;
        case 3:// deg by deg and number of pairs
          myPrintNumPair=true;
          myPrintDeg=true;
        case 2:// disable at the moment
          //myPrintFinalFull=true;
        case 1:// some final stats
          myPrintFinal=true;
          break;
	case 0:
          myPrintFinalSimple=true;
          break;
        default:
        case -1:
          break;
    };
  }

  void Stats::myUpgradeDegStats(const degree& new_deg, unsigned int Pairs_no)
  {
    DegStats New(new_deg);
    //  New.Deg=new_deg;
    for (list<DegStats>::const_iterator it=myDegByDeg.begin();
         it!=myDegByDeg.end();
         ++it)
    {
      New.myPInserted+=it->myPInserted;
      New.myGMKilled+=it->myGMKilled;
      New.myCopKilled+=it->myCopKilled;
      New.myUseful+=it->myUseful;
      New.myUseless+=it->myUseless;
      New.myBKilled+=it->myBKilled;
    }
    New.myPInserted=myPInserted-New.myPInserted;
    New.myGMKilled=myGMKilled-New.myGMKilled;
    New.myCopKilled=myCopKilled-New.myCopKilled;
    New.myBKilled=myBKilled-New.myBKilled;
    New.myUseful=myUseful-New.myUseful;
    New.myUseless=myUseless-New.myUseless;
    New.myPairsNo=Pairs_no;

    myDegByDeg.push_back(New);
  }


  void Stats::myStampa(ostream& out)const
  {
    if (myPrintFinalFull)
    {

      if (!myDegByDeg.empty()) {
        out<<"The degree by degree Stats\n";
        for (list<DegStats>::const_iterator it=myDegByDeg.begin();
             it!=myDegByDeg.end();++it) {
          out<<"Degree "<<it->myDeg<<endl;
          out<<"** Pairs Reduced         "<<it->myUseful+it->myUseless
              <<" of which "<<it->myUseful<<" useful and "
              <<it->myUseless<<" useless\n";
          out<<"** Pairs Inserted        "<<it->myPInserted;
          if (it->myPInserted!=0)	
	  out<<" of which "<<it->myPInserted-it->myGMKilled-it->myCopKilled<<" survived";
          out<<endl;
          if (it->myPInserted!=0) {
            out<<"****  of which GMKilled    "<<it->myGMKilled<<endl;
            out<<"****  of which CopKilled   "<<it->myCopKilled<<endl;
            out<<"****  of which BKilled     "<<it->myBKilled<<endl;	
          };
          out<<"Pairs at end of deg      "<<it->myPairsNo<<endl;
          out<<"\n";
        };//for
      };//if
      if (myPrintPolyLen)
      {
        out<< "Poly Lens\n";
        out<<"[[";
        list< pair<unsigned int, unsigned int> >::const_iterator it1;
        for (it1=myPolyLens.begin();it1!=myPolyLens.end();++it1)
          out<<",["<<it1->first<<","<<it1->second<<"]";
        out<<"]"<<endl;
      }
    };// if (myPrintFinalFull)

    if (myPrintFinal)
    {
      out<<"The Stats\n";
      out<<" Poly in Basis  "<<myUseful-myPolyDeleted<<endl;
      out<<" Pairs Reduced        "<<myUseful+myUseless
          <<" = ("<<myUseful+myUseless-myNumGens<<"+"<<myNumGens<<" Gens)"
          <<" of which "<<myUseful<<" useful and "
	  <<myUseless<<" useless\n";
      if (myBKilled==0) out<<" Minimal Pairs"<<myUseful+myUseless+myCopKilled-myNumGens
                           <<" = (Reduced "<<myUseful+myUseless<<" + Coprime "
                           <<myCopKilled<<" - Gens "<<myNumGens<<")"<<endl;
      if (myPolyDeleted!=0)
        out<<" Poly Deleted         "<<myPolyDeleted<<endl;
      if (myPolyDHed!=0)
      {
        out<<" Poly Dehomog'ed      "<<myPolyDHed<<endl;// number of polynomials dehomogenized
        out<<" Degs dropped           "<<myDegDH<<endl;// sum of the degrees dropped
      }
      out<<" Pairs Inserted        "<<myPInserted
          <<" + "<<myNumGens<<" generators"
          <<" of which "<<myUseful+myUseless<<" survived"<<endl;
      out<<"   of which GMKilled     "<<myGMKilled<<endl;
      out<<"   of which CopKilled    "<<myCopKilled<<endl;
      out<<"   of which BKilled      "<<myBKilled<<endl;
      out<<" GM considered pairs   "<<myGMTouched<<endl;

      out<<" Back considered pairs "<<myBTouched<<endl;
     }//if (myPrintFinal)
     if (myPrintFinalSimple)
      {
        out << "[log]  reductions=" << myNReductions;
        out << "\tTotalTime="<<-myTotalTime<< endl;
      }//if (myPrintFinalSimple)
  }//myStampa


  Stats& Stats::operator=(const Stats& rhs)
  {
    if (this == &rhs) return *this;
    myNumGens = rhs.myNumGens;
    myLevel = rhs.myLevel;
    myPInserted = rhs.myPInserted;
    myGMKilled = rhs.myGMKilled;
    myCopKilled = rhs.myCopKilled;
    myBKilled = rhs.myBKilled;
    myBTouched = rhs.myBTouched;
    myGMTouched = rhs.myGMTouched;
    myUseful = rhs.myUseful;
    myUseless = rhs.myUseless;
    myPolyDeleted = rhs.myPolyDeleted;
    myPolyDHed = rhs.myPolyDHed;
    myDegDH = rhs.myDegDH;
    myNReductions = rhs.myNReductions;
    myReductionTime = rhs.myReductionTime;
    myTotalTime = rhs.myTotalTime;
    myPrintReduction = rhs.myPrintReduction;
    myPrintDeg = rhs.myPrintDeg;
    myPrintGM = rhs.myPrintGM;
    myPrintCop = rhs.myPrintCop;
    myPrintBC = rhs.myPrintBC;
    myPrintNumPair = rhs.myPrintNumPair;
    myPrintFinal = rhs.myPrintFinal;
    myPrintFinalFull = rhs.myPrintFinalFull;
    myPrintFinalSimple = rhs.myPrintFinalSimple;
    myPrintNewPairs = rhs.myPrintNewPairs;
    myPrintPolyDeleted = rhs.myPrintPolyDeleted;
    myPrintPolyDH = rhs.myPrintPolyDH;
    myPrintPolyLen = rhs.myPrintPolyLen;
    myPrintKill = rhs.myPrintKill;
    
    return *this;
  }
  


}// end namespace cocoa
		

// RCS header/log on the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/TmpGRStats.C,v 1.8 2014/07/31 14:45:18 abbott Exp $
// $Log: TmpGRStats.C,v $
// Revision 1.8  2014/07/31 14:45:18  abbott
// Summary: Merged io.H and UtilsTemplate.H into new header VectorOperations.H
// Author: JAA
//
// Revision 1.7  2014/04/30 16:15:40  abbott
// Summary: Replaced X.size()==0 by X.empty()
// Author: JAA
//
// Revision 1.6  2013/01/31 11:41:18  bigatti
// -- added operator=
//
// Revision 1.5  2012/05/28 09:18:20  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.4  2010/03/23 14:43:07  bigatti
// -- class GRingInfo estracted from TmpGPoly
//
// Revision 1.3  2008/09/19 14:08:16  bigatti
// -- modified GRStats (M.Caboara)
//
// Revision 1.2  2007/10/30 17:14:07  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1  2007/03/09 18:56:56  bigatti
// -- added Tmp prefix to Groebner related files
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.4  2007/03/08 18:22:29  cocoa
// Just whitespace cleaning.
//
// Revision 1.3  2007/03/07 17:04:31  cocoa
// -- several changes by M.Caboara: more operations on ideals,
//    exception cleaner, coding conventions, WSugar, dynamic
//
// Revision 1.2  2006/12/21 13:48:33  cocoa
// Made all increment/decrement calls prefix (except where the must be postfix).
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.2  2006/01/17 15:44:56  cocoa
// -- chamges by Max for operations with modules
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.1.1.1  2005/05/03 15:47:31  cocoa
// Imported files
//
// Revision 1.2  2005/04/19 14:06:04  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.4  2004/06/29 17:10:22  cocoa
// Partially tidied use of "protected" and "private" in various
// base classes.  Checking in at the end of the day -- it works,
// and I wouldn't want it to be lost next time point's disk
// misbehaves.
//
// Revision 1.3  2004/06/16 16:13:41  cocoa
// Improved I/O facilities with knock-on changes
//
// Revision 1.2  2004/05/27 16:30:15  cocoa
// -- removed ";" at the end of function bodies (g++ gives error on them)
//
// Revision 1.1.1.1  2003/09/24 12:55:43  cocoa
// Imported files
//
// Revision 1.6  2003/06/23 17:09:42  abbott
// Minor cleaning prior to public release.
// Improved the include directives,
//
// Revision 1.5  2003/05/29 16:46:22  bigatti
// - added: myLevel
//
// Revision 1.4  2003/05/14 17:03:34  bigatti
// *** empty log message ***
//

