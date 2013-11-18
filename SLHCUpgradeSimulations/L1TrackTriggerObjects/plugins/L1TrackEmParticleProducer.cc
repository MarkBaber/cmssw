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

#include "DataFormats/L1TrackTrigger/interface/L1TrackEmParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TrackEmParticleFwd.h"

#include "DataFormats/Math/interface/LorentzVector.h"


// for L1Tracks:
#include "SimDataFormats/SLHC/interface/StackedTrackerTypes.h"

#include <string>
#include "TMath.h"


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
	edm::InputTag L1EGammaInputTag;
	edm::InputTag L1TrackInputTag;
	std::string label;

} ;


//
// constructors and destructor
//
L1TrackEmParticleProducer::L1TrackEmParticleProducer(const edm::ParameterSet& iConfig)
{

   L1EGammaInputTag = iConfig.getParameter<edm::InputTag>("L1EGammaInputTag") ;
   L1TrackInputTag = iConfig.getParameter<edm::InputTag>("L1TrackInputTag");
   label = iConfig.getParameter<std::string>("label");
   

   produces<L1TrackEmParticleCollection>(label);
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
 iEvent.getByLabel(L1EGammaInputTag,EGammaHandle);
 std::vector<L1EmParticle>::const_iterator egIter ;

 edm::Handle<L1TkTrackCollectionType> L1TkTrackHandle;
 iEvent.getByLabel(L1TrackInputTag, L1TkTrackHandle);
 L1TkTrackCollectionType::const_iterator trackIter;

 int ieg = 0;
 for (egIter = EGammaHandle->begin();  egIter != EGammaHandle->end(); ++egIter) {

    edm::Ref< L1EmParticleCollection > EGammaRef( EGammaHandle, ieg );
    ieg ++;

    int bx = egIter -> bx() ;

    if (bx == 0) {

	// calculate the isolation of the L1EG object with
	// respect to L1Tracks.
	// here dummy

	float trkisol = -999;

    	const math::XYZTLorentzVector P4 = egIter -> p4() ;
	L1TrackEmParticle trkEm(  P4,
				EGammaRef,
				trkisol );
    
        result -> push_back( trkEm );

     }  // endif bx==0

 }  // end loop over EGamma objects

 iEvent.put( result, label );

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



