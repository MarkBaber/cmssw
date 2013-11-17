#ifndef L1TrackTrigger_L1EmParticle_h
#define L1TrackTrigger_L1EmParticle_h

// -*- C++ -*-
//
// Package:     L1Trigger
// Class  :     L1TrackEmParticle
// 

#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"

//#include "SimDataFormats/SLHC/interface/L1TkTrack.h"

namespace l1extra {
         
   class L1TrackEmParticle : public reco::LeafCandidate
   {     
         
      public:
         enum EmType
         {
            kIsolated,
            kNonIsolated,
            kUndefined, 
            kNumOfEmTypes
         } ;

  //typedef L1TkTrack_PixelDigi_                          L1TkTrackType;
  //typedef std::vector< L1TkTrackType >                  L1TkTrackCollectionType;
           
         L1TrackEmParticle();

	 L1TrackEmParticle( const LorentzVector& p4,
			    const edm::Ref< L1EmParticleCollection >& egRef,
			    //const edm::Ref< L1TkTrackCollectionType >& trkRef,
			    int bx = 0 );

	virtual ~L1TrackEmParticle() {}

         // ---------- const member functions ---------------------

         EmType type() const
         { return type_ ; }

	 const edm::Ref< L1EmParticleCollection >& egRef() const
	 { return egRef_ ; }

	 //const edm::Ref< L1TkTrackCollectionType >& trkRef() const
	 //{ return trkRef_ ; }

	 float getTrkIsol() const { return TrkIsol_ ; } 
 	 float getTrkzVtx() const { return TrkzVtx_ ; }
         int bx() const
         { return bx_ ; } 

         // ---------- member functions ---------------------------

	 void setTrkIsol(float TrkIsol)  { TrkIsol_ = TrkIsol ; }
	 void setTrkzVtx(float TrkzVtx)  { TrkzVtx_ = TrkzVtx ; }

      private:

         EmType type_ ;
	 edm::Ref< L1EmParticleCollection > egRef_ ;
	 //emd::Ref< L1TkTrackCollectionType > trkRef_ ;
	 float TrkIsol_;
	 float TrkzVtx_ ;
	 int bx_ ;
    };
}

#endif


