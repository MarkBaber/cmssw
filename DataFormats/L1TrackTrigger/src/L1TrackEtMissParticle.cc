#include "DataFormats/L1TrackTrigger/interface/L1TrackEtMissParticle.h"

using namespace l1extra ;

L1TrackEtMissParticle::L1TrackEtMissParticle()
{
}

L1TrackEtMissParticle::L1TrackEtMissParticle(
        const LorentzVector& p4,
        EtMissType type,
        const double& etTotal,
        const double& etMissPU,
        const double& etTotalPU,
        const edm::Ref< L1TrackPrimaryVertexCollection >& avtxRef,
        int bx )
   : LeafCandidate( ( char ) 0, p4 ),
     type_( type ),
     etTot_( etTotal ),
     etMissPU_ ( etMissPU ),
     etTotalPU_ ( etTotalPU ),
     vtxRef_( avtxRef ),
     bx_( bx )
{
}


