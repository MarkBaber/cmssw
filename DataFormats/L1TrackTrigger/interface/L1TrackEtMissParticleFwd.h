#ifndef L1TrackTrigger_L1TrackEtMissParticleFwd_h
#define L1TrackTrigger_L1TrackEtMissParticleFwd_h
// -*- C++ -*-
//
// Package:     L1Trigger
// Class  :     L1TrackEtMissParticleFwd
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
