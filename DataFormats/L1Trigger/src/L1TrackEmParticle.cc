// -*- C++ -*-
//
// Package:     L1Trigger
// Class  :     L1TrackEmParticle
// 

#include "DataFormats/L1Trigger/interface/L1TrackEmParticle.h"


using namespace l1extra ;


L1TrackEmParticle::L1TrackEmParticle()
{
}

L1TrackEmParticle::L1TrackEmParticle( const LorentzVector& p4,
         const edm::Ref< L1EmParticleCollection >& egRef,
         //const edm::Ref< L1TkTrackCollectionType >& trkRef,
         int bx)
   : LeafCandidate( ( char ) 0, p4 ),
     egRef_ ( egRef ),
     bx_ ( bx )  
{

  if ( egRef_.isNonnull()  ) {
	//type_ = egRef()->type()  ;
	//type_ = egRef_.type() ;
	L1EmParticle eg = 
	L1EmParticle::EmType typ = egRef.type() ;
        if (typ == L1EmParticle::EmType::kIsolated) {
	  type_ = kIsolated;
	}
	else {
	  type_ = kNonIsolated ;
	}
  }

}




