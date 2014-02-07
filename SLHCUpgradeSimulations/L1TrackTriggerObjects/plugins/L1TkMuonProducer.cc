// -*- C++ -*-
//
//
// dummy producer for a L1TkMuonParticle
// This is just an interface, taking the muon objects created
// by PierLuigi's code, and putting them into a collection of
// L1TkMuonParticle.
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

#include "DataFormats/L1TrackTrigger/interface/L1TkMuonParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkMuonParticleFwd.h"

#include "DataFormats/L1TrackTrigger/interface/L1TrackPrimaryVertex.h"

#include "DataFormats/Math/interface/LorentzVector.h"


// for L1Tracks:
#include "SimDataFormats/SLHC/interface/StackedTrackerTypes.h"

#include <string>
#include "TMath.h"


using namespace l1extra ;

//
// class declaration
//

class L1TkMuonParticleProducer : public edm::EDProducer {
   public:

  typedef L1TkTrack_PixelDigi_                          L1TkTrackType;
  typedef std::vector< L1TkTrackType >                               L1TkTrackCollectionType;

      explicit L1TkMuonParticleProducer(const edm::ParameterSet&);
      ~L1TkMuonParticleProducer();

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
	
	 edm::InputTag L1PLInputTag;  // inputTag for PierLuigi's objects


} ;


//
// constructors and destructor
//
L1TkMuonParticleProducer::L1TkMuonParticleProducer(const edm::ParameterSet& iConfig)
{

   L1PLInputTag = iConfig.getParameter<edm::InputTag>("L1PLInputTag");

   produces<L1TkMuonParticleCollection>();
}

L1TkMuonParticleProducer::~L1TkMuonParticleProducer() {
}

// ------------ method called to produce the data  ------------
void
L1TkMuonParticleProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

 std::auto_ptr<L1TkMuonParticleCollection> result(new L1TkMuonParticleCollection);

	// PL: the muon objects from PierLuigi
/*
 edm::Handle<XXXCollection> XXXHandle;
 iEvent.getByLabel(L1PLInputTag,XXXHandle);
 std::vector<XXX>::const_iterator muIter ;

 if (!XXXHandle.isValid() ) {
          LogError("L1TkMuonParticleProducer")
            << "\nWarning: L1XXXCollectionType with " << L1PLInputTag
            << "\nrequested in configuration, but not found in the event. Exit."
            << std::endl;
           return;
 }

	// Now loop over the muons of Pierluigi 

 int imu = 0;
 for (muIter = XXXHandle->begin();  muIter != XXXHandle->end(); ++muIter) {

    edm::Ref< XXXCollection > muRef( XXXHandle, imu );
    imu ++;

    // int bx = egIter -> bx() ;	// if PL's objects have a bx method
    int bx = 0;    // else...

    if (bx == 0) {

	edm::Ptr< L1TkTrackType > L1TrackPtr ;

	// PL : get the matched L1Track from PL's object
	// L1TrackPtr  = muRef -> getRefToTheL1Track() ;
	
            float px = L1TrackPtr -> getMomentum().x();
            float py = L1TrackPtr -> getMomentum().y();
            float pz = L1TrackPtr -> getMomentum().z();
            float e = sqrt( px*px + py*py + pz*pz );    // massless particle
            math::XYZTLorentzVector TrackP4(px,py,pz,e);

	// the code may calculate a tracker-based isolation variable,
	// or pick it up from PL's object if it is there.
	// for the while, dummy.
	float trkisol = -999;


	L1TkMuonParticle trkMu(  P4,
				// muRef,  
				L1TrackPtr,	
				trkisol );
    
     }  // endif bx==0

 }  // end loop over Pierluigi's objects
*/

 iEvent.put( result );

}

// --------------------------------------------------------------------------------------


// ------------ method called once each job just before starting event loop  ------------
void
L1TkMuonParticleProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
L1TkMuonParticleProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
L1TkMuonParticleProducer::beginRun(edm::Run& iRun, edm::EventSetup const& iSetup)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
L1TkMuonParticleProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
L1TkMuonParticleProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
L1TkMuonParticleProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
L1TkMuonParticleProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1TkMuonParticleProducer);



