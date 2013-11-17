// -*- C++ -*-
//
//
// dummy producer for a L1TrackEmParticle
// 

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/L1Trigger/interface/L1TrackEmParticle.h"
#include "DataFormats/L1Trigger/interface/L1TrackEmParticleFwd.h"

#include "DataFormats/Math/interface/LorentzVector.h"


// for L1Tracks:
#include "SimDataFormats/SLHC/interface/StackedTrackerTypes.h"


using namespace l1extra ;

//
// class declaration
//

class L1TrackEmParticleProducer : public edm::EDProducer {
   public:

  typedef L1TkTrack_PixelDigi_                          L1TkTrackType;
  typedef std::vector< L1TkTrackType >                               L1TkTrackCollectionType;

      explicit L1TrackEmParticleProducer(const edm::ParameterSet&);
      ~L1TrackEmParticleProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      //virtual void beginRun(edm::Run&, edm::EventSetup const&);
      //virtual void endRun(edm::Run&, edm::EventSetup const&);
      //virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      //virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------
	edm::InputTag L1EGammaLabel;
	edm::InputTag L1TrackLabel;

} ;


//
// constructors and destructor
//
L1TrackEmParticleProducer::L1TrackEmParticleProducer(const edm::ParameterSet& iConfig)
{

   L1EGammaLabel = iConfig.getParameter<edm::InputTag>("L1EGammaLabel") ;
   L1TrackLabel = iConfig.getParameter<edm::InputTag>("L1TrackLabel");

   produces<L1TrackEmParticleCollection>();
}

L1TrackEmParticleProducer::~L1TrackEmParticleProducer() {
}

// ------------ method called to produce the data  ------------
void
L1TrackEmParticleProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

 std::auto_ptr<L1TrackEmParticleCollection> result(new L1TrackEmParticleCollection);

 edm::Handle<L1EmParticleCollection> EGammaHandle;
 iEvent.getByLabel(L1EGammaLabel,EGammaHandle);
 std::vector<L1EmParticle>::const_iterator egIter ;

 edm::Handle<L1TkTrackCollectionType> L1TkTrackHandle;
 //iEvent.getByLabel("L1Tracks","Level1TkTracks",L1TkTrackHandle);
 iEvent.getByLabel(L1TrackLabel, L1TkTrackHandle);
 L1TkTrackCollectionType::const_iterator trackIter;


 int ieg = 0;
 for (egIter = EGammaHandle->begin();  egIter != EGammaHandle->end(); ++egIter) {

    edm::Ref< L1EmParticleCollection > EGammaRef( EGammaHandle, ieg );
    ieg ++;

	// match the L1EG object with a L1Track
	// here dummy : I take the hightest PT track

	float ptmax = 0;
	int itr = 0;
	int itrack = -1;
	for (trackIter = L1TkTrackHandle->begin(); trackIter != L1TkTrackHandle->end(); ++trackIter) {
		float pt = trackIter->getMomentum().perp();	
	        if ( pt > ptmax) {
		  ptmax = pt;
		  itrack = itr;
		}
		itr ++ ;
	}
	edm::Ref< L1TkTrackCollectionType > L1TrackRef( L1TkTrackHandle, itrack) ;
	float pt = L1TrackRef -> getMomentum().perp();
	//float eta = L1TrackRef -> getMomentum().eta();
	//float phi = L1TrackRef -> getMomentum().phi();
	//float dummyMass = 0. ;
	
 	//math::PtEtaPhiMLorentzVector TrackP4( pt, eta, phi, dummyMass );
 	float px = L1TrackRef -> getMomentum().x();
	float py = L1TrackRef -> getMomentum().y();
	float pz = L1TrackRef -> getMomentum().z();
        math::XYZTLorentzVector TrackP4(px,py,pz,pt);
	int ibx = 0;
 	L1TrackEmParticle trkEm( TrackP4, 
				 EGammaRef,
				 ibx );
		// zvertex of the L1Track:
        float z = trackIter->getVertex().z();
	trkEm.setTrkzVtx( z );

		// isolation w.r.t. tracks :
 	float isol = -999;	// dummy
	trkEm.setTrkzVtx( isol );

	result -> push_back( trkEm );
 }  // end loop over EGamma objects

 iEvent.put( result );

}


// ------------ method called once each job just before starting event loop  ------------
void
L1TrackEmParticleProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
L1TrackEmParticleProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
L1TrackEmParticleProducer::beginRun(edm::Run& iRun, edm::EventSetup const& iSetup)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
L1TrackEmParticleProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
L1TrackEmParticleProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
L1TrackEmParticleProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
L1TrackEmParticleProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1TrackEmParticleProducer);



