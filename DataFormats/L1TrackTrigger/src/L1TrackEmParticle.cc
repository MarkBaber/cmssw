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
         //const edm::Ptr< L1TkTrackType >& trkPtr,
	 float tkisol ,
	 int bx )
   : LeafCandidate( ( char ) 0, p4 ),
     egRef_ ( egRef ),
     //trkPtr_ ( trkPtr ),
     TrkIsol_ ( tkisol ) ,
     bx_ ( bx )  
{


  if ( egRef_.isNonnull()  ) {
	
        L1EmParticle::EmType typ =  getEGRef() -> type() ;
        if (typ == L1EmParticle::EmType::kIsolated) {
	  type_ = kIsolated;
	}
	else {
	  type_ = kNonIsolated ;
	}
  }

  

}




