// -*- C++ -*-
//
// Package:     L1Trigger
// Class  :     L1TkEmParticle
// 

#include "DataFormats/L1TrackTrigger/interface/L1TkElectronParticle.h"


using namespace l1extra ;


L1TkElectronParticle::L1TkElectronParticle()
{
}

L1TkElectronParticle::L1TkElectronParticle( const LorentzVector& p4,
         const edm::Ref< L1EmParticleCollection >& egRef,
         const edm::Ptr< L1TkTrackType >& trkPtr,
	 float tkisol )
   : L1TkEmParticle( p4, egRef, tkisol) ,
     trkPtr_ ( trkPtr )

{

 if ( trkPtr_.isNonnull() ) {
	float z = getTrkPtr() -> getVertex().z(); 
	setTrkzVtx( z );

	// temporary stuff:
	float chi2 = getTrkPtr() -> getChi2();
	float ptTra = getTrkPtr() -> getMomentum().perp();
	std::vector< edm::Ptr< L1TkStub_PixelDigi_ > > theStubs = getTrkPtr() -> getStubPtrs();
	int nStubs = (int)theStubs.size();

	chi2_ = chi2;
	ptTra_ = ptTra ;
	nStubs_ = nStubs ;
 }
}



