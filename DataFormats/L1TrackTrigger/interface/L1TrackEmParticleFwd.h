#ifndef L1TrackTrigger_L1EmParticleFwd_h
#define L1TrackTrigger_L1EmParticleFwd_h

// -*- C++ -*-
//
// Package:     L1Trigger
// Class  :     L1TrackEmParticleFwd
// 

#include <vector>
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"


namespace l1extra {

   class L1TrackEmParticle ;

   typedef std::vector< L1TrackEmParticle > L1TrackEmParticleCollection ;

   typedef edm::Ref< L1TrackEmParticleCollection > L1TrackEmParticleRef ;
   typedef edm::RefVector< L1TrackEmParticleCollection > L1TrackEmParticleRefVector ;
   typedef std::vector< L1TrackEmParticleRef > L1TrackEmParticleVectorRef ;
}

#endif


