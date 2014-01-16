#ifndef L1TkTrigger_L1ElectronParticle_h
#define L1TkTrigger_L1ElectronParticle_h

// -*- C++ -*-
//
// Package:     L1Trigger
// Class  :     L1TkEmParticle
// 

#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/Ptr.h"

#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"

#include "DataFormats/L1TrackTrigger/interface/L1TkEmParticle.h"

//#include "SimDataFormats/SLHC/interface/L1TkTrack.h"
#include "SimDataFormats/SLHC/interface/StackedTrackerTypes.h"


namespace l1extra {
         
   class L1TkElectronParticle : public L1TkEmParticle
   {     
         
      public:

  typedef L1TkTrack_PixelDigi_                          L1TkTrackType;
  typedef std::vector< L1TkTrackType >                  L1TkTrackCollectionType;
           
         L1TkElectronParticle();

	 L1TkElectronParticle( const LorentzVector& p4,
			    const edm::Ref< L1EmParticleCollection >& egRef,
			    const edm::Ptr< L1TkTrackType >& trkPtr,
			    float tkisol = -999. );

	virtual ~L1TkElectronParticle() {}

         // ---------- const member functions ---------------------

         const edm::Ptr< L1TkTrackType >& getTrkPtr() const
         { return trkPtr_ ; }

 	 float getTrkzVtx() const { return TrkzVtx_ ; }


         // ---------- member functions ---------------------------

	 void setTrkzVtx(float TrkzVtx)  { TrkzVtx_ = TrkzVtx ; }

	// temporary stuff:
	void setDeltas(float dphi, float dphi_prime, float deta, float dr) { dphi_ = dphi ;
							   dphi_prime_ = dphi_prime;
							   deta_ = deta ;
							   dr_ = dr ; }
      private:

         edm::Ptr< L1TkTrackType > trkPtr_ ;
	 float TrkzVtx_ ;

	// temporary stuff:
	float ptTra_ ;
	int nStubs_ ;
	float dphi_ ;
	float dphi_prime_;
	float deta_ ;
	float dr_ ;
	float chi2_ ;

    };
}

#endif


