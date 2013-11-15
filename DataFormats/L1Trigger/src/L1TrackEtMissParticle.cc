#include "DataFormats/L1Trigger/interface/L1TrackEtMissParticle.h"

using namespace l1extra ;

L1TrackEtMissParticle::L1TrackEtMissParticle()
{
}

L1TrackEtMissParticle::L1TrackEtMissParticle(
        const LorentzVector& p4,
        EtMissType type,
        const double& etTotal,
        const edm::Ref< L1TrackPrimaryVertexCollection >& avtxRef,
        int bx )
   : LeafCandidate( ( char ) 0, p4 ),
     type_( type ),
     etTot_( etTotal ),
     vtxRef_( avtxRef ),
     bx_( bx )
{
}


