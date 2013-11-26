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

#include "SimDataFormats/SLHC/interface/StackedTrackerTypes.h"


namespace l1extra {
         
   class L1TrackEmParticle : public reco::LeafCandidate
   {     
         
      public:
           
         L1TrackEmParticle();

	 L1TrackEmParticle( const LorentzVector& p4,
			    const edm::Ref< L1EmParticleCollection >& egRef,
			    float tkisol = -999. );

	virtual ~L1TrackEmParticle() {}

         // ---------- const member functions ---------------------

	 const edm::Ref< L1EmParticleCollection >& getEGRef() const
	 { return egRef_ ; }

	 float getTrkIsol() const { return TrkIsol_ ; } 

         float getTrkIsol_v1() const { return TrkIsol_v1_ ; }
         float getTrkIsol_v2() const { return TrkIsol_v2_ ; }


         // ---------- member functions ---------------------------

	 void setTrkIsol(float TrkIsol)  { TrkIsol_ = TrkIsol ; }

         void setTrkIsol_v1(float TrkIsol)  { TrkIsol_v1_ = TrkIsol ; }
         void setTrkIsol_v2(float TrkIsol)  { TrkIsol_v2_ = TrkIsol ; }

	 int bx() const;

      private:

	 edm::Ref< L1EmParticleCollection > egRef_ ;
	 float TrkIsol_;

	// variants for the isolation :
	float TrkIsol_v1_;
	float TrkIsol_v2_;
    };
}

#endif


