#ifndef L1TrackTrigger_L1ElectronParticle_h
#define L1TrackTrigger_L1ElectronParticle_h

// -*- C++ -*-
//
// Package:     L1Trigger
// Class  :     L1TrackEmParticle
// 

#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/Ptr.h"

#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"

#include "DataFormats/L1TrackTrigger/interface/L1TrackEmParticle.h"

//#include "SimDataFormats/SLHC/interface/L1TkTrack.h"
#include "SimDataFormats/SLHC/interface/StackedTrackerTypes.h"


namespace l1extra {
         
   class L1TrackElectronParticle : public L1TrackEmParticle
   {     
         
      public:

  typedef L1TkTrack_PixelDigi_                          L1TkTrackType;
  typedef std::vector< L1TkTrackType >                  L1TkTrackCollectionType;
           
         L1TrackElectronParticle();

	 L1TrackElectronParticle( const LorentzVector& p4,
			    const edm::Ref< L1EmParticleCollection >& egRef,
			    const edm::Ptr< L1TkTrackType >& trkPtr,
			    float tkisol = -999. );

	virtual ~L1TrackElectronParticle() {}

         // ---------- const member functions ---------------------

         const edm::Ptr< L1TkTrackType >& getTrkPtr() const
         { return trkPtr_ ; }

 	 float getTrkzVtx() const { return TrkzVtx_ ; }

         // ---------- member functions ---------------------------

	 void setTrkzVtx(float TrkzVtx)  { TrkzVtx_ = TrkzVtx ; }

      private:

         edm::Ptr< L1TkTrackType > trkPtr_ ;
	 float TrkzVtx_ ;

    };
}

#endif


