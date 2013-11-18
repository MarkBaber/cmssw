// -*- C++ -*-
//
// Package:     L1Trigger
// Class  :     L1TrackEmParticle
// 

#include "DataFormats/L1TrackTrigger/interface/L1TrackElectronParticle.h"


using namespace l1extra ;


L1TrackElectronParticle::L1TrackElectronParticle()
{
}

L1TrackElectronParticle::L1TrackElectronParticle( const LorentzVector& p4,
         const edm::Ref< L1EmParticleCollection >& egRef,
         const edm::Ptr< L1TkTrackType >& trkPtr,
	 float tkisol )
   : L1TrackEmParticle( p4, egRef, tkisol) ,
     trkPtr_ ( trkPtr )

{

 if ( trkPtr_.isNonnull() ) {
	float z = getTrkPtr() -> getVertex().z(); 
	setTrkzVtx( z );
 }
}



