// -*- C++ -*-
//
// Package:    L1TrackTriggerObjectsAnalyzer
// Class:      L1TrackTriggerObjectsAnalyzer
// 
/**\class L1TrackTriggerObjectsAnalyzer L1TrackTriggerObjectsAnalyzer.cc SLHCUpgradeSimulations/L1TrackTriggerObjectsAnalyzer/src/L1TrackTriggerObjectsAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Emmanuelle Perez,40 1-A28,+41227671915,
//         Created:  Thu Nov 14 11:22:13 CET 2013
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"


// Gen-level stuff:
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/L1TrackTrigger/interface/L1TrackPrimaryVertex.h"
#include "DataFormats/L1TrackTrigger/interface/L1TrackEtMissParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TrackEtMissParticleFwd.h"
#include "DataFormats/L1TrackTrigger/interface/L1TrackEmParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TrackEmParticleFwd.h"
#include "DataFormats/L1TrackTrigger/interface/L1TrackElectronParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TrackElectronParticleFwd.h"




#include "TFile.h"
#include "TH1F.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

using namespace l1extra;



//
// class declaration
//

class L1TrackTriggerObjectsAnalyzer : public edm::EDAnalyzer {
   public:

  typedef L1TkTrack_PixelDigi_                          L1TkTrackType;
  typedef std::vector< L1TkTrackType >                               L1TkTrackCollectionType;

      explicit L1TrackTriggerObjectsAnalyzer(const edm::ParameterSet&);
      ~L1TrackTriggerObjectsAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      //virtual void endRun(edm::Run const&, edm::EventSetup const&);
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------

	// to test the L1TrackPrimaryVertex :
	edm::InputTag L1VtxInputTag;
	TH1F* h_zgen;
	TH1F* h_dz1;
	TH1F* h_dz2;

	// for L1TrackEtmiss:
	edm::InputTag L1EtMissInputTag;

	// for L1TrackEmParticles
        edm::InputTag L1TrackPhotonsInputTag;
	edm::InputTag L1TrackElectronsInputTag;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
L1TrackTriggerObjectsAnalyzer::L1TrackTriggerObjectsAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

  edm::Service<TFileService> fs;
  int nbins = 25;
  float x1=-25.;
  float x2 = 25.;
   h_zgen = fs -> make<TH1F>("h_zgen",";generated z_{vtx} (cm); Evts",nbins, x1,x2);
   nbins=100;
   x1 = -2;
   x2 = 2;
   h_dz1 = fs -> make<TH1F>("h_dz1",";z_{L1} - z_{gen} (cm); Evts",nbins,x1,x2);
   h_dz2 = fs -> make<TH1F>("h_dz2",";z_{L1} - z_{gen} (cm); Evts",nbins, x1, x2);

  L1VtxInputTag = iConfig.getParameter<edm::InputTag>("L1VtxInputTag") ;
  L1EtMissInputTag = iConfig.getParameter<edm::InputTag>("L1EtMissInputTag");
  L1TrackElectronsInputTag = iConfig.getParameter<edm::InputTag>("L1TrackElectronsInputTag");
  L1TrackPhotonsInputTag = iConfig.getParameter<edm::InputTag>("L1TrackPhotonsInputTag");
}


L1TrackTriggerObjectsAnalyzer::~L1TrackTriggerObjectsAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
L1TrackTriggerObjectsAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   std::cout << " ----  a new event ----- " << std::endl;

	// First, retrieve the generated primary vertex


  edm::Handle<edm::HepMCProduct> HepMCEvt;
  iEvent.getByLabel("generator",HepMCEvt);

     const HepMC::GenEvent* MCEvt = HepMCEvt->GetEvent();
     const double mm=0.1;

     float zvtx_gen = -999;
     for ( HepMC::GenEvent::vertex_const_iterator ivertex = MCEvt->vertices_begin(); ivertex != MCEvt->vertices_end(); ++ivertex )
         {
             bool hasParentVertex = false;
 
             // Loop over the parents looking to see if they are coming from a production vertex
             for (
                 HepMC::GenVertex::particle_iterator iparent = (*ivertex)->particles_begin(HepMC::parents);
                 iparent != (*ivertex)->particles_end(HepMC::parents);
                 ++iparent
             )
                 if ( (*iparent)->production_vertex() )
                 {
                     hasParentVertex = true;
                     break;
                }
 
             // Reject those vertices with parent vertices
             if (hasParentVertex) continue;
 
             // Get the position of the vertex
             HepMC::FourVector pos = (*ivertex)->position();
	     zvtx_gen = pos.z()*mm; 

	     break;  // there should be one single primary vertex

          }  // end loop over gen vertices

     h_zgen -> Fill( zvtx_gen );
     std::cout << " Generated zvertex : " << zvtx_gen << std::endl;



	// ----------------------------------------------------------------------
	//
        // retrieve the L1 vertices
        
 edm::Handle<L1TrackPrimaryVertexCollection> L1VertexHandle;
 iEvent.getByLabel(L1VtxInputTag,L1VertexHandle);
 std::vector<L1TrackPrimaryVertex>::const_iterator vtxIter;
 
 if ( L1VertexHandle.isValid() ) {
     std::cout << " -----  L1TrackPrimaryVertex objects   ----- " << std::endl;
     int ivtx = 0;
	// several algorithms have been run in the L1TrackPrimaryVertexProducer
	// hence there is a collection of L1 primary vertices.
	// the first one is probably the most reliable.

     for (vtxIter = L1VertexHandle->begin(); vtxIter != L1VertexHandle->end(); ++vtxIter) {
        float z = vtxIter -> getZvertex();
        float sum = vtxIter -> getSum();
        std::cout << " a vertex with  z = sum " << z << " " << sum << std::endl;
        ivtx ++;
        if (ivtx == 1) h_dz1 -> Fill( z - zvtx_gen) ;
        if (ivtx == 2) h_dz2 -> Fill( z - zvtx_gen);
     }  
 }

	//
        // ----------------------------------------------------------------------
	// retrieve the EtMiss objects
	//

 edm::Handle<L1TrackEtMissParticleCollection> L1TrackEtMissHandle;
 iEvent.getByLabel(L1EtMissInputTag, L1TrackEtMissHandle);
 std::vector<L1TrackEtMissParticle>::const_iterator etmIter;

 if (L1TrackEtMissHandle.isValid() ) {
    std::cout << " -----  L1TrackEtMiss objects  -----  " << std::endl; 
    for (etmIter = L1TrackEtMissHandle -> begin(); etmIter != L1TrackEtMissHandle->end(); ++etmIter) {
	float etmis = etmIter -> et();
	const edm::Ref< L1TrackPrimaryVertexCollection > vtxRef = etmIter -> getVtxRef();
	float zvtx = vtxRef -> getZvertex();
        float etMissPU = etmIter -> etMissPU();
	std::cout << " ETmiss = " << etmis << " for zvtx = " << zvtx << " and ETmiss from PU = " << etMissPU << std::endl;
    }
 }

        //
        // ----------------------------------------------------------------------
        // retrieve the L1TrackEmParticle objects
	//

 edm::Handle<L1TrackEmParticleCollection> L1TrackPhotonsHandle;
 iEvent.getByLabel(L1TrackPhotonsInputTag, L1TrackPhotonsHandle);
 std::vector<L1TrackEmParticle>::const_iterator phoIter ;

 if ( L1TrackPhotonsHandle.isValid() ) {
    std::cout << " -----   L1TrackEmParticle  objects -----  " << std::endl;
    for (phoIter = L1TrackPhotonsHandle -> begin(); phoIter != L1TrackPhotonsHandle->end(); ++phoIter) {
	float et = phoIter -> pt();
	float phi = phoIter -> phi();
	float eta = phoIter -> eta();
        float trkisol = phoIter -> getTrkIsol() ;
	const edm::Ref< L1EmParticleCollection > EGref = phoIter -> getEGRef();
	float et_L1Calo = EGref -> et();
	float eta_calo = EGref -> eta();
	float phi_calo = EGref -> phi();

	std::cout << " a photon candidate ET eta phi trkisol " << et << " " << eta << " " << phi << " " << trkisol << std::endl;
	std::cout << "                Calo  ET eta phi " << et_L1Calo << " " << eta_calo << " " << phi_calo << std::endl; 
    }
 }


        //
        // ----------------------------------------------------------------------
        // retrieve the L1TrackElectronParticle objects
        //

 edm::Handle<L1TrackElectronParticleCollection> L1TrackElectronsHandle;
 iEvent.getByLabel(L1TrackElectronsInputTag, L1TrackElectronsHandle);
 std::vector<L1TrackElectronParticle>::const_iterator eleIter ;

 if ( L1TrackElectronsHandle.isValid() ) {
    std::cout << " -----   L1TrackElectronParticle  objects -----  " << std::endl;
    for (eleIter = L1TrackElectronsHandle -> begin(); eleIter != L1TrackElectronsHandle->end(); ++eleIter) {
        float et = eleIter -> pt();
        float phi = eleIter -> phi();
        float eta = eleIter -> eta();
    	float trkisol = eleIter -> getTrkIsol() ;
	float ztr = eleIter -> getTrkzVtx() ;
        std::cout << "an electron candidate ET eta phi trkisol ztr " << et << " " << eta << " " << phi << " " << trkisol << " " << ztr << std::endl;

        const edm::Ref< L1EmParticleCollection > EGref = eleIter -> getEGRef();
        if ( EGref.isNonnull() ) {
           float et_L1Calo = EGref -> et();
           float eta_calo = EGref -> eta();
           float phi_calo = EGref -> phi();
           std::cout << "                Calo  ET eta phi " << et_L1Calo << " " << eta_calo << " " << phi_calo << std::endl;
	}
	else {
	    std::cout << " .... edm::Ref to EGamma is unvalid !! ?  " << std::endl;
	}

        const edm::Ptr< L1TkTrackType > TrkRef = eleIter -> getTrkPtr();
	if ( TrkRef.isNonnull() ) {
            float pt_track = TrkRef -> getMomentum().perp();
            float phi_track = TrkRef -> getMomentum().phi();
            float eta_track = TrkRef -> getMomentum().eta();
            float ztrack = TrkRef -> getVertex().z() ;
            std::cout << "                Track PT eta phi ztr " << pt_track << " " << eta_track << " " << phi_track << " " << ztrack << std::endl;
	}
	else {
	    std::cout << " ... edm::Ptr to L1Tracks is unvalid (e.g. electron was matched to stubs) " << std::endl;
	}
    }
 }


}


// ------------ method called once each job just before starting event loop  ------------
void 
L1TrackTriggerObjectsAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
L1TrackTriggerObjectsAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
L1TrackTriggerObjectsAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
L1TrackTriggerObjectsAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
L1TrackTriggerObjectsAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
L1TrackTriggerObjectsAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
L1TrackTriggerObjectsAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1TrackTriggerObjectsAnalyzer);
