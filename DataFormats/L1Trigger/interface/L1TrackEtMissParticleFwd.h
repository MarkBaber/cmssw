#ifndef L1Trigger_L1TrackEtMissParticleFwd_h
#define L1Trigger_L1TrackEtMissParticleFwd_h
// -*- C++ -*-
//
// Package:     L1Trigger
// Class  :     L1TrackEtMissParticleFwd
// 
/**\class L1TrackEtMissParticleRef \file L1TrackEtMissParticleFwd.h DataFormats/L1Trigger/interface/L1TrackEtMissParticleFwd.h \author Werner Sun

 Description: typedefs for L1TrackEtMissParticleRef and associated containers.
*/
//
// Original Author:  Werner Sun
//         Created:  Sat Jul 15 14:28:43 EDT 2006
// $Id: L1TrackEtMissParticleFwd.h,v 1.3 2007/12/18 03:26:49 wsun Exp $
//

// system include files

// user include files

// forward declarations
#include <vector>
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"


namespace l1extra {

   class L1TrackEtMissParticle ;

   //typedef edm::RefProd< L1TrackEtMissParticle > L1TrackEtMissParticleRefProd ;

   typedef std::vector< L1TrackEtMissParticle > L1TrackEtMissParticleCollection ;

   //typedef edm::Ref< L1TrackEtMissParticleCollection > L1TrackEtMissParticleRef ;
   //typedef edm::RefVector< L1TrackEtMissParticleCollection > L1TrackEtMissParticleRefVector ;
   //typedef std::vector< L1TrackEtMissParticleRef > L1TrackEtMissParticleVectorRef ;
}

#endif
