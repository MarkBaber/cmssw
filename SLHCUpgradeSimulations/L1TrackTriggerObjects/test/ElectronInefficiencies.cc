// -*- C++ -*-
//
// Package:    ElectronInefficiencies
// Class:      ElectronInefficiencies
// 
/**\class ElectronInefficiencies ElectronInefficiencies.cc SLHCUpgradeSimulations/ElectronInefficiencies/src/ElectronInefficiencies.cc

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
#include "DataFormats/L1TrackTrigger/interface/L1TkEtMissParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkEtMissParticleFwd.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkEmParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkEmParticleFwd.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkElectronParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkElectronParticleFwd.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkJetParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkJetParticleFwd.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkHTMissParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkHTMissParticleFwd.h"

#include "SLHCUpgradeSimulations/L1TrackTriggerObjects/interface/L1TkElectronTrackMatchAlgo.h"


// for inefficiencies

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/Common/interface/DetSetVector.h"

#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"



#include "TFile.h"
#include "TH1F.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

using namespace l1extra;

using namespace L1TkElectronTrackMatchAlgo;

using namespace std;
using namespace edm;


//
// class declaration
//

class ElectronInefficiencies : public edm::EDAnalyzer {
   public:

  typedef L1TkTrack_PixelDigi_                          L1TkTrackType;
  typedef std::vector< L1TkTrackType >                               L1TkTrackCollectionType;

      explicit ElectronInefficiencies(const edm::ParameterSet&);
      ~ElectronInefficiencies();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

   double deltaR(double eta1, double eta2, double phi1, double phi2);
   float deltaR(float eta1, float eta2, float phi1, float phi2);
   double DeltaPhi(double phi1, double phi2);

      void GetHistory(const edm::Event& iEvent, int verbose ) ;


      //virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      //virtual void endRun(edm::Run const&, edm::EventSetup const&);
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      float MyCorrectedEta(float eta, float zv)  ;

      // ----------member data ---------------------------

	int nevts ;	// nb of events where the generated electron is in a given range
	int nevts_Ele ;	// nb of events where the generated electron is in a given range and for which an L1TkElectron is found

	float p_gen;
	float pt_gen;
	float eta_gen;
	float phi_gen;

	int nb_brems_pixel;
	int nb_brems_outer;
	int ilayer_brem_max;
	float brem_max;

        vector<unsigned int> Electron_SimTrackIds;
        vector<int> BremLayers;
        vector<float> BremF;

	bool BREM_IN_PIXELS;
	bool BREM_IN_OUTER ;


	// to test the L1TrackPrimaryVertex :
	TH1F* h_zgen;
	TH1F* h_dz1;
	TH1F* h_dz2;

	TH1F* h_nbrems_pixel;
	TH1F* h_nbrems_outer;
	TH1F* h_ilayer_brem_max;  // layer with max brem, in outer tracker
	TH1F* h_bremmax;	// max brem, in outer tracker

        TH1F* hpass_nbrems_pixel;
        TH1F* hpass_nbrems_outer;
        TH1F* hpass_ilayer_brem_max;  // layer with max brem, in outer tracker
        TH1F* hpass_bremmax;       // max brem, in outer tracker

        TH1F* effi_nbrems_pixel;
        TH1F* effi_nbrems_outer;
        TH1F* effi_ilayer_brem_max;  // layer with max brem, in outer tracker
        TH1F* effi_bremmax;        // max brem, in outer tracker

	edm::InputTag L1TkElectronsInputTag;
        edm::InputTag L1EGammaInputTag;
        edm::InputTag L1TrackInputTag;

	float  BREM_CUT ;

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
ElectronInefficiencies::ElectronInefficiencies(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

  nevts = 0;
  nevts_Ele = 0;

  Electron_SimTrackIds.clear();
  BremLayers.clear();
  BremF.clear();


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

   h_nbrems_pixel = fs -> make<TH1F>("h_nbrems_pixel",";Nb of brems in pixels; Evts",5,-0.5,4.5);
   h_nbrems_outer = fs -> make<TH1F>("h_nbrems_outer",";Nb of brems in outerT; Evts",5,-0.5,4.5);
   h_ilayer_brem_max = fs -> make<TH1F>("h_ilayer_brem_max",";Layer in OT with largest brem; Evts",6,4.5,10.5);
   h_bremmax = fs -> make<TH1F>("h_bremmax",";f_{brem, max} in OT; Evts",10,0.,1.);

   hpass_nbrems_pixel = fs -> make<TH1F>("hpass_nbrems_pixel",";Nb of brems in pixels; Evts",5,-0.5,4.5);
   hpass_nbrems_outer = fs -> make<TH1F>("hpass_nbrems_outer",";Nb of brems in outerT; Evts",5,-0.5,4.5);
   hpass_ilayer_brem_max = fs -> make<TH1F>("hpass_ilayer_brem_max",";Layer in OT with largest brem; Evts",6,4.5,10.5);
   hpass_bremmax = fs -> make<TH1F>("hpass_bremmax",";f_{brem, max} in OT; Evts",10,0.,1.);

   effi_nbrems_pixel = fs -> make<TH1F>("effi_nbrems_pixel",";Nb of brems in pixels; Efficiency",5,-0.5,4.5);
   effi_nbrems_outer = fs -> make<TH1F>("effi_nbrems_outer",";Nb of brems in outerT; Efficiency",5,-0.5,4.5);
   effi_ilayer_brem_max = fs -> make<TH1F>("effi_ilayer_brem_max",";Layer in OT with largest brem; Efficiency",6,4.5,10.5);
   effi_bremmax = fs -> make<TH1F>("effi_bremmax",";f_{brem, max} in OT; Efficiency",10,0.,1.);


 h_nbrems_pixel -> Sumw2();
 h_nbrems_outer -> Sumw2();
 h_ilayer_brem_max -> Sumw2();
 h_bremmax -> Sumw2();
 hpass_nbrems_pixel -> Sumw2();
 hpass_nbrems_outer -> Sumw2();
 hpass_ilayer_brem_max -> Sumw2();
 hpass_bremmax -> Sumw2();

  L1TkElectronsInputTag = iConfig.getParameter<edm::InputTag>("L1TkElectronsInputTag");

   L1EGammaInputTag = iConfig.getParameter<edm::InputTag>("L1EGammaInputTag") ;
   L1TrackInputTag = iConfig.getParameter<edm::InputTag>("L1TrackInputTag");

   BREM_CUT = (float)iConfig.getParameter<double>("BREM_CUT");

}


ElectronInefficiencies::~ElectronInefficiencies()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ElectronInefficiencies::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   //std::cout << " ----  a new event ----- " << std::endl;

  Electron_SimTrackIds.clear();
  BremLayers.clear();
  BremF.clear();

 nb_brems_pixel = 0;
 nb_brems_outer = 0;
 ilayer_brem_max = -1;
 brem_max = 0;


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
     //std::cout << " Generated zvertex : " << zvtx_gen << std::endl;



	// 
	// ----------------------------------------------------------------------
	// look for the generated electron
	//

 eta_gen = -999;
 pt_gen = -999;
 phi_gen = -999;
 p_gen = -999;
 int nElegen = 0;

 for ( HepMC::GenEvent::particle_const_iterator p = MCEvt->particles_begin();
        p != MCEvt->particles_end();
        ++p ) {

        int id = (*p)->pdg_id();
        if (id == 11 || id == -11) {
                HepMC::FourVector momentum = (*p)->momentum();
                //p_gen = momentum.e();
                eta_gen = momentum.eta();
                phi_gen = momentum.phi();
                pt_gen = momentum.perp();
		p_gen = momentum.e();
                //std::cout << " an electron with pt eta phi " << pt_gen << " " << eta_gen << " " << phi_gen << std::endl;
		nElegen ++;
        }
 }

 if (nElegen != 1) std::cout << " ..... nElegen != 1 ?? " << nElegen << std::endl;

 if (pt_gen < 20.) return;
 if ( fabs(eta_gen) > 1. )  return ;

  //std::cout << " an electron with pt eta phi " << pt_gen << " " << eta_gen << " " << phi_gen << std::endl;


  edm::Handle<L1EmParticleCollection> EGammaHandle;
  iEvent.getByLabel(L1EGammaInputTag,EGammaHandle);
  l1extra::L1EmParticleCollection eGammaCollection = (*EGammaHandle.product());
  //  sort(eGammaCollection.begin(), eGammaCollection.end(), L1TkElectron::EtComparator());
  l1extra::L1EmParticleCollection::const_iterator egIter;
  edm::Handle<L1TkTrackCollectionType> L1TkTrackHandle;
  iEvent.getByLabel(L1TrackInputTag, L1TkTrackHandle);
  L1TkTrackCollectionType::const_iterator trackIter;

 if( !EGammaHandle.isValid() )
        {
          edm::LogError("L1TkElectronTrackProducer")
            << "\nWarning: L1EmParticleCollection with " << L1EGammaInputTag
            << "\nrequested in configuration, but not found in the event. Exit"
            << std::endl;
           return;
        }

 if (!L1TkTrackHandle.isValid() ) {
          edm::LogError("L1TkEmParticleProducer")
            << "\nWarning: L1TkTrackCollectionType with " << L1TrackInputTag
            << "\nrequested in configuration, but not found in the event. Exit."
            << std::endl;
           return;
 }


	// loop over the L1EG objects:

  float dr_gen_EG_min = 9999;
  l1extra::L1EmParticleCollection::const_iterator matched_EG   ;

  for (egIter = eGammaCollection.begin();  egIter != eGammaCollection.end(); ++egIter) {
    
    int ibx = egIter -> bx();
    if (ibx != 0) continue;
    
    float e_ele   = egIter->energy();
    float eta_ele = egIter->eta();
    float phi_ele = egIter->phi();
    float et_ele = 0;
    
    if (cosh(eta_ele) > 0.0) et_ele = e_ele/cosh(eta_ele);
    else et_ele = -1.0;
        
        if (et_ele < 10.) continue;
        if (fabs(eta_ele) > 2.5) continue;

    float dr_gen_EG = deltaR( eta_ele, eta_gen, phi_ele, phi_gen );
    if (dr_gen_EG < dr_gen_EG_min) {
	dr_gen_EG_min = dr_gen_EG ;
	if (dr_gen_EG_min < 0.5) matched_EG = egIter;
    }

  }

  if ( dr_gen_EG_min > 0.5 )  return ;	   // no L1EG found that matched the generated electron

 nevts ++;


        //
        // ----------------------------------------------------------------------
        // retrieve the L1TkElectronParticle objects
        //


 edm::Handle<L1TkElectronParticleCollection> L1TkElectronsHandle;
 iEvent.getByLabel(L1TkElectronsInputTag, L1TkElectronsHandle);
 std::vector<L1TkElectronParticle>::const_iterator eleIter ;

 
 int nL1TkEle = 0;

 int nL1TkEle_match = 0;

 if ( L1TkElectronsHandle.isValid() ) {
    //std::cout << " -----   L1TkElectronParticle  objects -----  " << std::endl;
    for (eleIter = L1TkElectronsHandle -> begin(); eleIter != L1TkElectronsHandle->end(); ++eleIter) {
	nL1TkEle ++;
        float phi = eleIter -> phi();
        float eta = eleIter -> eta();
        float dr = deltaR( eta, eta_gen, phi, phi_gen );
	if (dr < 0.5) nL1TkEle_match ++ ;

/*
        float et = eleIter -> pt();
        float phi = eleIter -> phi();
        float eta = eleIter -> eta();
	int bx = eleIter -> bx();
    	float trkisol = eleIter -> getTrkIsol() ;
	float ztr = eleIter -> getTrkzVtx() ;
        std::cout << "an electron candidate ET eta phi bx trkisol ztr " << et << " " << eta << " " << phi << " " << bx << " " << trkisol << " " << ztr << std::endl;

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
*/
    }
 }

 nL1TkEle = nL1TkEle_match ;

 if (nL1TkEle > 0) nevts_Ele  ++;



GetHistory( iEvent,0 );
 
 h_nbrems_pixel -> Fill( nb_brems_pixel );
 h_nbrems_outer -> Fill( nb_brems_outer );
 h_ilayer_brem_max -> Fill( ilayer_brem_max );
 h_bremmax -> Fill( brem_max );

 if (nL1TkEle > 0) {
    hpass_nbrems_pixel -> Fill( nb_brems_pixel );
    hpass_nbrems_outer -> Fill( nb_brems_outer );
    hpass_ilayer_brem_max -> Fill( ilayer_brem_max );
    hpass_bremmax -> Fill( brem_max );

 }

 if (nL1TkEle > 0) { 

/*
cout << "       --- FOUND MATCH.   Print MC history : " << endl;
  std::cout << "  generated electron p pt eta phi " << p_gen << " " << pt_gen << " " << eta_gen << " " << phi_gen << std::endl;
        
  std::cout << " Generated zvertex : " << zvtx_gen << std::endl;

GetHistory( iEvent, 0 );
cout << endl;
*/

     return;
 }


	// now look at those events for which no L1TkElectron is found

  std::cout << " -----    An event with no L1TkElectron      ----- " << std::endl;
  std::cout << "  generated electron p pt eta phi " << p_gen << " " << pt_gen << " " << eta_gen << " " << phi_gen << std::endl;

  std::cout << " Generated zvertex : " << zvtx_gen << std::endl;


 std::cout << " loop over L1EG objects " << std::endl;

  // for (egIter = eGammaCollection.begin();  egIter != eGammaCollection.end(); ++egIter) {

  egIter = matched_EG ;
    
    float e_ele   = egIter->energy();
    float eta_ele = egIter->eta();
    float et_ele = 0;

    if (cosh(eta_ele) > 0.0) et_ele = e_ele/cosh(eta_ele);
    else et_ele = -1.0; 

        //if (et_ele < 10.) continue;

        std::cout << " a L1EG object with ET > 15: ET, eta, phi " << et_ele << " " << eta_ele << " " << egIter->phi() << std::
endl;

    //if (fabs(eta_ele) > 2.5) continue;
    // match the L1EG object with a L1Track
    //float drmin = 999;
    int itr = 0;
    //int itrack = -1;

 float PTMINTRA = 5. ;
 float CHI2MAX = 100. ;

    for (trackIter = L1TkTrackHandle->begin(); trackIter != L1TkTrackHandle->end(); ++trackIter) {
      edm::Ptr< L1TkTrackType > L1TrackPtr( L1TkTrackHandle, itr) ;
      if ( L1TrackPtr->getMomentum().perp() > PTMINTRA && L1TrackPtr->getChi2() < CHI2MAX) {

        double dPhi = 99.;
        float dR = 99.;
        float dEta = 99.;
	float dphiPrime = 99.;
        L1TkElectronTrackMatchAlgo::doMatch(egIter, trackIter, dPhi, dR, dEta, dphiPrime);

	float track_pt = L1TrackPtr->getMomentum().perp()  ;
	float track_eta = L1TrackPtr->getMomentum().eta() ;
	float track_phi = L1TrackPtr->getMomentum().phi() ;
	float track_z = L1TrackPtr->getVertex().z() ;

        //float eta_cor = MyCorrectedEta( eta_ele,  zvtx_gen );
	float eta_cor = MyCorrectedEta( eta_ele,  track_z);
	float eta_atanu = dEta + track_eta;
	std::cout << " eta L1EG w/o and w correction : " << eta_ele << " " << eta_cor << " Atanu: " << eta_atanu << std::endl;

	std::cout << " a track :  dPhi, dR, dEta " << dPhi << " " << dR << " " << dEta << " --- pt eta phi z " << track_pt << " " << track_eta << " " << track_phi << " " << track_z << std::endl;



/*
        if (abs(dPhi) < dPhiCutoff && dR < dRCutoff && abs(dEta) < dEtaCutoff && dR < drmin) {
          drmin = dR;
          itrack = itr;
        }
*/
      }
      itr++;
    }    // end loop over tracks

 // } // end loop over EG objects

// MC history :

//cout << "   	--- Print MC history : " << endl;
//GetHistory( iEvent,1 );
//cout << endl;

}


// --------------------------------------------------------------------------------------

float ElectronInefficiencies::MyCorrectedEta(float eta, float zv)  {

// Correct the eta of the L1EG object once we know the zvertex

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

void ElectronInefficiencies::GetHistory(const edm::Event& iEvent, int verbose ) {


  ///////////////////
  // GET SIMHITS  //
  edm::Handle<edm::PSimHitContainer> simHitslowHandle;
  edm::Handle<edm::PSimHitContainer> simHitshighHandle;
  iEvent.getByLabel( "g4SimHits","TrackerHitsPixelBarrelLowTof",simHitslowHandle);
  iEvent.getByLabel( "g4SimHits","TrackerHitsPixelBarrelHighTof",simHitshighHandle);


 BREM_IN_PIXELS = false;
 BREM_IN_OUTER = false;

PSimHitContainer::const_iterator iterSimHits;

 bool stopHere = false;
 unsigned int iprevious = -1;
 int itrack = 0;

 float p_first_layer = -999;	// momentum of the track after crossing the first layer of the pixel detector

 for ( iterSimHits = simHitslowHandle->begin();
        iterSimHits != simHitslowHandle->end();
        ++iterSimHits) {

   unsigned short processType = iterSimHits -> processType();
   int type = iterSimHits -> particleType();
   //if (type<0) type = -type;

   if (type == 11 && (processType == 2 || processType == 13)) {
         unsigned int trackId = iterSimHits -> trackId();

	if (p_first_layer < 0) {
	    p_first_layer = iterSimHits -> pabs();
	}

        if (verbose > 10) {
   	   float pabs = iterSimHits -> pabs();
   	   float eloss = iterSimHits -> energyLoss();
   	   float tof = iterSimHits -> timeOfFlight();
   	   unsigned int detid = iterSimHits -> detUnitId();
   	   DetId tkId( detid );

   	   int layer = -1;
   	   if (tkId.subdetId() == 1 ) {
        	   PXBDetId pxbdetid(tkId);
        	   layer = pxbdetid.layer();
   	   }
	     cout << " a hit : trackId partType processType " << trackId << " " << type << " " << processType << " p eloss tof " << pabs << " " << eloss << " " << tof << " layer " << layer << endl;
	}   //endif verbose > 10

        if (! stopHere) {
           if (trackId != iprevious) Electron_SimTrackIds.push_back( trackId );
           iprevious = trackId;
           itrack ++;
        }
   }
   else {
        stopHere = true;
   }
 }  // end loop over simhits

   int ntracks = itrack;

 if (verbose > 1) {
 	cout << " Tracks associated to the primary : " ;
 	for (unsigned int l=0; l < Electron_SimTrackIds.size(); l++) {
    	cout << Electron_SimTrackIds[l] << " " ;
 	}
 	cout << endl;
 }  // endif verbose > 1

 float dfirst = fabs( p_first_layer - p_gen) / p_gen ;
 if (dfirst > BREM_CUT) {
	if (verbose > 0) {
	   cout << " ....... interaction before crossing the first tracker layer ..... " << endl;
	}
 }

        // ----------------------------------------------------------------------------
        // Determine if there has been a hard brem
        
 //cout << " second loop over simhits " << endl;
 float pprevious = p_first_layer ;
 itrack = 0;

 for ( iterSimHits = simHitslowHandle->begin();
        iterSimHits != simHitslowHandle->end();
        ++iterSimHits) {

   itrack ++;
   if (itrack > ntracks) break;

   unsigned short processType = iterSimHits -> processType();
   int type = iterSimHits -> particleType();
   float pabs = iterSimHits -> pabs();
   // unsigned int trackId = iterSimHits -> trackId();

   if (type == 11 && processType == 13) {

        unsigned int detid = iterSimHits -> detUnitId();
        DetId tkId( detid );

        int layer = -1;
        if (tkId.subdetId() == 1 ) {
                PXBDetId pxbdetid(tkId);
                layer = pxbdetid.layer();
        }
        float f = ( pabs - pprevious ) / pprevious;
        float eloss = iterSimHits -> energyLoss();
	if (verbose > 1) cout << " .... layer " << layer << " previous " << pprevious << " eloss " << eloss << " pabs: " << pabs << endl;

        if ( fabs( f ) > BREM_CUT) {   	// e.g. 10%
		BremLayers.push_back( layer );
        	BremF.push_back( f ) ;

		if (layer >=1 && layer <=4) {
			BREM_IN_PIXELS = true;
			nb_brems_pixel ++;
		}
		if (layer >=5 && layer <=10) {
			BREM_IN_OUTER = true;
			nb_brems_outer ++;
		}
	}

	if (layer >=5 && layer <=10) {
	    if ( fabs( f ) > brem_max) {
		brem_max = fabs( f ) ;
		ilayer_brem_max = layer;
	    }
	}
	

   } // end brem

   pprevious = pabs;

 }  // loop over simHits

 if (verbose > 0) {

    if (BREM_IN_PIXELS) cout << " ..... BREM IN PIXELS .... " << endl;

    if (BREM_IN_OUTER) {
	cout << " ..... BREM IN OUTER ..... " << endl;
 	cout << " Brem layers "  << endl ;
 	for (unsigned int l=0; l < BremLayers.size(); l++) {
    	cout << " layer : " << BremLayers[l] << " " << " dp/p " << BremF[l] << endl;
 	}
 	cout << endl;
    }
 }  //endif verbose > 0


}



// ------------ method called once each job just before starting event loop  ------------
void 
ElectronInefficiencies::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ElectronInefficiencies::endJob() 
{
 std::cout << " Number of events studied : " << nevts << std::endl;
 std::cout << " with L1TkEle : " << nevts_Ele << std::endl;
 std::cout << " efficiency : " << (float)nevts_Ele/(float)nevts << std::endl;

 effi_nbrems_pixel -> Divide( hpass_nbrems_pixel, h_nbrems_pixel, 1.0, 1.0, "B");
 effi_nbrems_outer -> Divide( hpass_nbrems_outer, h_nbrems_outer, 1.0, 1.0, "B");
 effi_ilayer_brem_max  -> Divide(hpass_ilayer_brem_max, h_ilayer_brem_max, 1.0, 1.0, "B");
 effi_bremmax -> Divide(hpass_bremmax, h_bremmax, 1.0, 1.0, "B");

}


double ElectronInefficiencies::deltaR(double eta1, double eta2, double phi1, double phi2) {
    double deta = eta1 - eta2;
    //double dphi = phi1 - phi2;
    //if (dphi < 0) dphi = dphi + 2.*TMath::Pi();
    //if (dphi > TMath::Pi()) dphi = 2.*TMath::Pi() - dphi;
    double dphi = DeltaPhi(phi1, phi2);
    double DR= sqrt( deta*deta + dphi*dphi );
    return DR;
}   
    
float ElectronInefficiencies::deltaR(float eta1, float eta2, float phi1, float phi2) {
        double DR = deltaR((double)eta1, (double)eta2, (double)phi1, (double)phi2);
        float DRf = DR; 
        return DRf;
}

double ElectronInefficiencies::DeltaPhi(double phi1, double phi2) {
  double dphi = phi1 - phi2;
    if (dphi < 0) dphi = dphi + 2.*TMath::Pi();
    if (dphi > TMath::Pi()) dphi = 2.*TMath::Pi() - dphi;
  return dphi;
}


// ------------ method called when starting to processes a run  ------------
/*
void 
ElectronInefficiencies::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
ElectronInefficiencies::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
ElectronInefficiencies::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
ElectronInefficiencies::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ElectronInefficiencies::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ElectronInefficiencies);
