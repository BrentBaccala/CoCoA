//   Copyright (c)  2007,2009  John Abbott

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


#include "CoCoA/CanonicalHom.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/PolyRing.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/FractionField.H"

// #include <iostream>  // for debugging only

namespace CoCoA
{

  RingHom CanonicalHom(const ring& domain, const ring& codomain)
  {
    if (domain == codomain) return IdentityHom(domain);

    // Two easy cases:
    if (IsZZ(domain)) return ZZEmbeddingHom(codomain);
    if (IsQQ(domain)) return QQEmbeddingHom(codomain); // NB result is only a partial hom!!

    // Check codomain first, as this makes it possible to exploit certain "shortcuts"
    if (IsFractionField(codomain))
    {
      if (domain == BaseRing(codomain))
        return EmbeddingHom(codomain);
      else
	return EmbeddingHom(codomain)(CanonicalHom(domain, BaseRing(codomain)));
    }
    if (IsPolyRing(codomain))
    {
      if (domain == CoeffRing(codomain))
        return CoeffEmbeddingHom(codomain);
      else
	return CoeffEmbeddingHom(codomain)(CanonicalHom(domain, CoeffRing(codomain)));
    }
    if (IsQuotientRing(codomain))
    {
      const QuotientRing QR = codomain;
      if (domain == BaseRing(QR))
        return QuotientingHom(QR);
      else
	return QuotientingHom(QR)(CanonicalHom(domain, BaseRing(QR)));
    }

    CoCoA_ERROR(ERR::CanonicalHomFail, "CanonicalHom(R1,R2)");
    return IdentityHom(codomain); // Never executed; just to keep the compiler quiet.
  }


  RingHom TmpChainCanonicalHom(const ring& domain, const ring& codomain)
  {
    //    std::clog << " TmpChainCanonicalHom ";
    try { return CanonicalHom(domain, codomain); }
    catch (const CoCoA::ErrorInfo& err) {if (err!=ERR::CanonicalHomFail) throw;}
    if (IsFractionField(codomain))
    {
      //      std::clog << " FrF ";
      return EmbeddingHom(codomain) (TmpChainCanonicalHom(domain, BaseRing(codomain)));
    }
    if (IsPolyRing(codomain))
    {
      //      std::clog << " P ";
      return CoeffEmbeddingHom(codomain) (TmpChainCanonicalHom(domain, CoeffRing(codomain)));
    }
    if (IsQuotientRing(codomain))
    {
      //      std::clog << " QR ";
      const QuotientRing QR = codomain;
      return QuotientingHom(QR) (TmpChainCanonicalHom(domain, BaseRing(QR)));
    }
    CoCoA_ERROR(ERR::CanonicalHomFail, "TmpChainCanonicalHom(R1,R2)");
    return IdentityHom(codomain); // Never executed; just to keep the compiler quiet.
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/CanonicalHom.C,v 1.10 2014/07/08 13:14:40 abbott Exp $
// $Log: CanonicalHom.C,v $
// Revision 1.10  2014/07/08 13:14:40  abbott
// Summary: Removed AsQuotientRing; added new defn of BaseRing
// Author: JAA
//
// Revision 1.9  2014/07/08 08:33:18  abbott
// Summary: Removed AsFractionField
// Author: JAA
//
// Revision 1.8  2014/07/07 12:12:17  abbott
// Summary: Removed AsPolyRing
// Author: JAA
//
// Revision 1.7  2012/02/10 10:26:40  bigatti
// -- changed RingZ.H, RingQ.H --> RingZZ.H, RingQQ.H
//
// Revision 1.6  2012/02/08 15:07:08  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.5  2011/02/18 12:56:08  bigatti
// -- added TmpChainCanonicalHom
//
// Revision 1.4  2009/07/24 14:21:03  abbott
// Cleaned up include directives, and added some fwd decls.
//
// Revision 1.3  2009/07/24 12:27:46  abbott
// Added some comments.
//
// Revision 1.2  2007/10/30 17:14:08  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.2  2007/03/07 11:18:01  bigatti
// -- fixed CanonicalHom goto flag
//
// Revision 1.1  2007/03/05 21:25:02  cocoa
// New CanonicalHom pseudo-ctor.
//
//
