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

#include "DataFormats/L1TrackTrigger/interface/L1TrackElectronParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TrackElectronParticleFwd.h"

#include "DataFormats/Math/interface/LorentzVector.h"


// for L1Tracks:
#include "SimDataFormats/SLHC/interface/StackedTrackerTypes.h"

#include <string>
#include "TMath.h"


using namespace l1extra ;

//
// class declaration
//

class L1TrackElectronParticleProducer : public edm::EDProducer {
   public:

  typedef L1TkTrack_PixelDigi_                          L1TkTrackType;
  typedef std::vector< L1TkTrackType >                               L1TkTrackCollectionType;

      explicit L1TrackElectronParticleProducer(const edm::ParameterSet&);
      ~L1TrackElectronParticleProducer();

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
L1TrackElectronParticleProducer::L1TrackElectronParticleProducer(const edm::ParameterSet& iConfig)
{

   L1EGammaInputTag = iConfig.getParameter<edm::InputTag>("L1EGammaInputTag") ;
   L1TrackInputTag = iConfig.getParameter<edm::InputTag>("L1TrackInputTag");
   label = iConfig.getParameter<std::string>("label");
   

   produces<L1TrackElectronParticleCollection>(label);
}

L1TrackElectronParticleProducer::~L1TrackElectronParticleProducer() {
}

// ------------ method called to produce the data  ------------
void
L1TrackElectronParticleProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

 std::auto_ptr<L1TrackElectronParticleCollection> result(new L1TrackElectronParticleCollection);

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

    int ibx = egIter -> bx();
    if (ibx != 0) continue;

    float eta = egIter -> eta();
    float phi = egIter -> phi();
    
	// match the L1EG object with a L1Track
	// here dummy : I simply take the closest track
	// and require that DR < 0.5

	float drmin = 999;
	int itr = 0;
	int itrack = -1;
	for (trackIter = L1TkTrackHandle->begin(); trackIter != L1TkTrackHandle->end(); ++trackIter) {
		float Eta = trackIter->getMomentum().eta();
		float Phi = trackIter->getMomentum().phi();
		float deta = eta - Eta;
		float dphi = phi - Phi;
		if (dphi < 0) dphi = dphi + 2.*TMath::Pi();
		if (dphi > TMath::Pi()) dphi = 2.*TMath::Pi() - dphi;
		float dR = sqrt( deta*deta + dphi*dphi );
		if (dR < drmin) {
		  drmin = dR;
		  itrack = itr;
		}
		itr ++ ;
	}

	if (drmin < 0.5 ) {	// found a L1Track matched to the L1EG object

	    edm::Ptr< L1TkTrackType > L1TrackPtr( L1TkTrackHandle, itrack) ;
	    
 	    float px = L1TrackPtr -> getMomentum().x();
	    float py = L1TrackPtr -> getMomentum().y();
	    float pz = L1TrackPtr -> getMomentum().z();
	    float e = sqrt( px*px + py*py + pz*pz );	// massless particle
            math::XYZTLorentzVector TrackP4(px,py,pz,e);

	    float trkisol = -999; 	// dummy

 	    L1TrackElectronParticle trkEm( TrackP4, 
				 EGammaRef,
				 L1TrackPtr, 
			         trkisol );

	/*
	Note : to create an L1TrackElectronParticle without a reference to a L1Track 
	(e.g. when matching the L1EGamma with stubs and not tracks :

	   The 4-momentum of the electron is then probably that of the L1EGamma
	   object - see L1TrackEmProducer.cc to retrieve P4
	   One could then do :

	   edm::Ptr< L1TkTrackType > L1TrackPtrNull;     //  null pointer
	   L1TrackElectronParticle trkEm( P4,
                                 EGammaRef,
                                 L1TrackPtrNull, 
                                 trkisol,
                                 ibx );

	   and one can set the "z" of the electron, as determined by the 
	   algorithm, via :
	   
	   trkEm.setTrkzVtx( z );
	*/

	    result -> push_back( trkEm );

	}  // endif drmin < 0.5

 }  // end loop over EGamma objects

 iEvent.put( result, label );

}


// ------------ method called once each job just before starting event loop  ------------
void
L1TrackElectronParticleProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
L1TrackElectronParticleProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
L1TrackElectronParticleProducer::beginRun(edm::Run& iRun, edm::EventSetup const& iSetup)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
L1TrackElectronParticleProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
L1TrackElectronParticleProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
L1TrackElectronParticleProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
L1TrackElectronParticleProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1TrackElectronParticleProducer);



