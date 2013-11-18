#include "Rtypes.h" 
#include "Math/Cartesian3D.h" 
#include "Math/Polar3D.h" 
#include "Math/CylindricalEta3D.h" 
#include "Math/PxPyPzE4D.h" 
#include <boost/cstdint.hpp> 

#include "DataFormats/L1TrackTrigger/interface/L1TrackPrimaryVertex.h"
#include "DataFormats/L1TrackTrigger/interface/L1TrackEtMissParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TrackEtMissParticleFwd.h"

#include "DataFormats/L1TrackTrigger/interface/L1TrackEmParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TrackEmParticleFwd.h"

#include "DataFormats/L1TrackTrigger/interface/L1TrackElectronParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TrackElectronParticleFwd.h"

#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/Common/interface/RefProd.h"

namespace {
  struct dictionary {


	// L1 Primary Vertex
     L1TrackPrimaryVertex trzvtx;
     edm::Wrapper<L1TrackPrimaryVertexCollection> trzvtxColl;
     edm::Ref< L1TrackPrimaryVertexCollection > trkvtxRef ;


	// L1TrackEtMiss... following L1EtMiss...
     l1extra::L1TrackEtMissParticle TketMiss ;
     l1extra::L1TrackEtMissParticleCollection TketMissColl ;
     edm::Wrapper<l1extra::L1TrackEtMissParticle> w_TketMiss;
     edm::Wrapper<l1extra::L1TrackEtMissParticleCollection> w_TketMissColl;
     //l1extra::L1TrackEtMissParticleRef refTkEtMiss ;
     //l1extra::L1TrackEtMissParticleRefVector refTkVecEtMiss ;
     //l1extra::L1TrackEtMissParticleVectorRef vecTkRefEtMissColl ;
     //l1extra::L1TrackEtMissParticleRefProd refTkProdEtMiss ;
     //edm::reftobase::Holder<reco::Candidate, l1extra::L1TrackEtMissParticleRef> rtbTkm1;
     //edm::reftobase::Holder<reco::Candidate, l1extra::L1TrackEtMissParticleRefProd> rtbTkm2;

	// L1TrackEmParticle
     l1extra::L1TrackEmParticleCollection trkemColl ;
     edm::Wrapper<l1extra::L1TrackEmParticleCollection> w_trkemColl;
     l1extra::L1TrackEmParticleRef reftrkEm ;
     //l1extra::L1TrackEmParticleRefVector refVectrkEmColl ;
     //l1extra::L1TrackEmParticleVectorRef vecReftrkEmColl ;
     //edm::reftobase::Holder<reco::Candidate, l1extra::L1TrackEmParticleRef> rtbtrke;

        // L1TrackElectronParticle
     l1extra::L1TrackElectronParticleCollection trkeleColl ;
     edm::Wrapper<l1extra::L1TrackElectronParticleCollection> w_trkeleColl;
     l1extra::L1TrackElectronParticleRef reftrkEle ;



  };
}
