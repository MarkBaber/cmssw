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

class L1TrackEmParticleProducer : public edm::EDProducer {
   public:

  typedef L1TkTrack_PixelDigi_                          L1TkTrackType;
  typedef std::vector< L1TkTrackType >                               L1TkTrackCollectionType;

      explicit L1TrackEmParticleProducer(const edm::ParameterSet&);
      ~L1TrackEmParticleProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

      float DeltaPhi(float phi1, float phi2) ;
      float deltaR(float eta1, float eta2, float phi1, float phi2) ;
      float CorrectedEta(float eta, float zv);


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
        edm::InputTag L1VertexInputTag; // used only when VtxConstrain = True.

	std::string label;

	float ZMAX;		// |z_track| < ZMAX in cm
	float CHI2MAX;		
	float DRmin;
	float DRmax;
	float PTmin;
	bool VtxConstrain;	// use the primary vertex (default = false)
	float DeltaZMax;	// | z_track - z_primaryvtx | < DeltaZMax in cm. 
				// Used only when VtxConstrain = True.
	float RelIsoCut;
} ;


//
// constructors and destructor
//
L1TrackEmParticleProducer::L1TrackEmParticleProducer(const edm::ParameterSet& iConfig)
{

   label = iConfig.getParameter<std::string>("label");  // label of the collection produced
							// e.g. EG or IsoEG if all objects are kept
							// EGIsoTrk or IsoEGIsoTrk if only the EG or IsoEG
							// objects that pass a cut RelIso < RelIsoCut are written
							// in the new collection.

   L1EGammaInputTag = iConfig.getParameter<edm::InputTag>("L1EGammaInputTag") ;
   L1TrackInputTag = iConfig.getParameter<edm::InputTag>("L1TrackInputTag");
   L1VertexInputTag = iConfig.getParameter<edm::InputTag>("L1VertexInputTag");	// only used when VtxConstrain = true

	// parameters for the calculation of the isolation :
   ZMAX = (float)iConfig.getParameter<double>("ZMAX");
   CHI2MAX = (float)iConfig.getParameter<double>("CHI2MAX");
   DRmin = (float)iConfig.getParameter<double>("DRmin");
   DRmax = (float)iConfig.getParameter<double>("DRmax");
   VtxConstrain = iConfig.getParameter<bool>("VtxConstrain");
   DeltaZMax = (float)iConfig.getParameter<double>("DeltaZMax");
	// cut applied on the relative isolation (if this number is <= 0, no cut is applied)
   RelIsoCut = (float)iConfig.getParameter<double>("RelIsoCut");

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

	// the L1EGamma objects
 edm::Handle<L1EmParticleCollection> EGammaHandle;
 iEvent.getByLabel(L1EGammaInputTag,EGammaHandle);
 std::vector<L1EmParticle>::const_iterator egIter ;

	// the L1Tracks
 edm::Handle<L1TkTrackCollectionType> L1TkTrackHandle;
 iEvent.getByLabel(L1TrackInputTag, L1TkTrackHandle);
 L1TkTrackCollectionType::const_iterator trackIter;

	// the primary vertex (used only if VtxConstrain = true)
 float zvtxL1tk = -999;
 if (VtxConstrain) {
 	edm::Handle<L1TrackPrimaryVertexCollection> L1VertexHandle;
	iEvent.getByLabel(L1VertexInputTag,L1VertexHandle);
	std::vector<L1TrackPrimaryVertex>::const_iterator vtxIter = L1VertexHandle->begin();
	   // by convention, the first vertex in the collection is the one that should
	   // be used by default
	zvtxL1tk = vtxIter -> getZvertex();
 }


	// Now loop over the L1EGamma objects
 int ieg = 0;
 for (egIter = EGammaHandle->begin();  egIter != EGammaHandle->end(); ++egIter) {

    edm::Ref< L1EmParticleCollection > EGammaRef( EGammaHandle, ieg );
    ieg ++;

    int bx = egIter -> bx() ;

    if (bx == 0) {

	float eta = egIter -> eta();
	if (VtxConstrain) {
		// The eta of the L1EG object is seen from (0,0,0).
		// if VtxConstrain = true, use the zvtxL1tk to correct eta.
	   eta = CorrectedEta( (float)eta, zvtxL1tk);
	}
	float phi = egIter -> phi();
	float et = egIter -> et();

	if (et < 10.) continue;

	// calculate the isolation of the L1EG object with
	// respect to L1Tracks.

	float trkisol = -999;
	float sumPt = 0;

	std::cout << " here an EG w et = " << et << std::endl;

	for (trackIter = L1TkTrackHandle->begin(); trackIter != L1TkTrackHandle->end(); ++trackIter) {

	   float Pt = trackIter->getMomentum().perp();
	   float Eta = trackIter->getMomentum().eta();
	   float Phi = trackIter->getMomentum().phi();
	   float z  = trackIter->getVertex().z();
	   if (fabs(z) > ZMAX) continue;
    	   if (Pt < PTmin) continue;
	   float chi2 = trackIter->getChi2();
	   if (chi2 > CHI2MAX) continue;

	   if (VtxConstrain) {
	     if ( zvtxL1tk > -999 && fabs( z - zvtxL1tk) >= DeltaZMax) continue;
	   }

	   float dr = deltaR(Eta, eta, Phi,phi);
	   if (dr < DRmax && dr >= DRmin)  {
		std::cout << " a track in the cone, z Pt = " << z << " " << Pt << std::endl;
		sumPt += Pt;
	   }

	}  // end loop over tracks

	if (et > 0) trkisol = sumPt / et;	// relative isolation

    	const math::XYZTLorentzVector P4 = egIter -> p4() ;
	L1TrackEmParticle trkEm(  P4,
				EGammaRef,
				trkisol );
    
        if (RelIsoCut <= 0) {
		// write the L1TrackEm particle to the collection, 
		// irrespective of its relative isolation
		result -> push_back( trkEm );
	}
	else {
		// the object is written to the collection only
		// if it passes the isolation cut
		if (trkisol <= RelIsoCut) result -> push_back( trkEm );
	}

     }  // endif bx==0

 }  // end loop over EGamma objects

 iEvent.put( result, label );

}

// --------------------------------------------------------------------------------------

float L1TrackEmParticleProducer::CorrectedEta(float eta, float zv)  {

// Correct the eta of teh L1EG object once we know the zvertex

bool IsBarrel = ( fabs(eta) < 1.479 );
float REcal = 129. ;
float ZEcal = 315.4 ;

float theta = 2. * TMath::ATan( TMath::Exp( - eta ) );
if (theta < 0) theta = theta + TMath::Pi();
float tantheta = TMath::Tan( theta );

float delta;
if (IsBarrel) {
        delta = REcal / tantheta ;
}
else {
        if (theta > 0) delta =  ZEcal;
        if (theta < 0) delta = -ZEcal;
}

float tanthetaprime = delta * tantheta / (delta - zv );

float thetaprime = TMath::ATan( tanthetaprime );
if (thetaprime < 0) thetaprime = thetaprime + TMath::Pi();

float etaprime = -TMath::Log( TMath::Tan( thetaprime / 2.) );
return etaprime;

}

// --------------------------------------------------------------------------------------

float L1TrackEmParticleProducer::DeltaPhi(float phi1, float phi2) {
// dPhi between 0 and Pi
   float dphi = phi1 - phi2;
   if (dphi < 0) dphi = dphi + 2.*TMath::Pi();
   if (dphi > TMath::Pi()) dphi = 2.*TMath::Pi() - dphi;
  return dphi;
}

// --------------------------------------------------------------------------------------

float L1TrackEmParticleProducer::deltaR(float eta1, float eta2, float phi1, float phi2) {
    float deta = eta1 - eta2;
    float dphi = DeltaPhi(phi1, phi2);
    float DR= sqrt( deta*deta + dphi*dphi );
    return DR;
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



