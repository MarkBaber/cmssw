// -*- C++ -*-
//
// Package:     L1Trigger
// Class  :     L1TrackEmParticle
// 

#include "DataFormats/L1TrackTrigger/interface/L1TrackEmParticle.h"


using namespace l1extra ;


L1TrackEmParticle::L1TrackEmParticle()
{
}

L1TrackEmParticle::L1TrackEmParticle( const LorentzVector& p4,
         const edm::Ref< L1EmParticleCollection >& egRef,
	 float tkisol )
   : LeafCandidate( ( char ) 0, p4 ),
     egRef_ ( egRef ),
     TrkIsol_ ( tkisol ) 
{

}




