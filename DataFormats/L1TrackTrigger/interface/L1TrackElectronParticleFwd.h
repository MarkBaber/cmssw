#ifndef L1TrackTrigger_L1ElectronParticleFwd_h
#define L1TrackTrigger_L1ElectronParticleFwd_h

// -*- C++ -*-
//
// Package:     L1Trigger
// Class  :     L1TrackElectronParticleFwd
// 

#include <vector>
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"


namespace l1extra {

   class L1TrackElectronParticle ;

   typedef std::vector< L1TrackElectronParticle > L1TrackElectronParticleCollection ;

   typedef edm::Ref< L1TrackElectronParticleCollection > L1TrackElectronParticleRef ;
   typedef edm::RefVector< L1TrackElectronParticleCollection > L1TrackElectronParticleRefVector ;
   typedef std::vector< L1TrackElectronParticleRef > L1TrackElectronParticleVectorRef ;
}

#endif


