import FWCore.ParameterSet.Config as cms

process = cms.Process("ALL")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20) )

from SLHCUpgradeSimulations.L1TrackTriggerObjects.singleElectronFiles_cfi import *

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
	'file:example_w_Tracks_and_vertex.root'
    )
)



# ---- Global Tag and geometry :
#      (needed e.g. when running raw2digi below)

process.load("Configuration.Geometry.GeometryIdeal_cff")

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'POSTLS261_V3::All', '')


# ---------------------------------------------------------------------------
#
# --- Run the L1TrackEtMiss producer :

# needs the L1TrackPrimaryVertex
# TrkMET calculated from all tracks that have |z - z_vertex | < DeltaZ,
# pt > Ptmin and at least nStubsmin stubs.
# To check the latter condition, the file must contain the stubs.

process.L1TrackEtMiss = cms.EDProducer('L1TrackEtMissProducer',
     L1VtxLabel = cms.InputTag("L1TrackPrimaryVertex"),
     ZMAX = cms.double ( 25. ) ,        # in cm
     CHI2MAX = cms.double( 100. ),
     DeltaZ = cms.double( 0.05 ),       # in cm
     Ptmin = cms.double( 2. ),
     nStubsmin = cms.int32( 4 )
)

process.pEtMiss = cms.Path( process.L1TrackEtMiss )

#
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
#
# --- test L1TrackEmParticle

# --- First, create l1extra objects for L1EGamma 
#     here I use the Run1 trigger. The produced objects dont really
#     make sense, but OK to check the technicalities.

        # raw2digi to get the gct digis
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.p0 = cms.Path( process.RawToDigi )
        # run L1Reco
process.load('Configuration.StandardSequences.L1Reco_cff')
process.L1Reco = cms.Path( process.l1extraParticles )

# --- Now run the L1TrackEmParticleProducer 

# "photons" :

process.L1TrackPhotons = cms.EDProducer("L1TrackEmParticleProducer",
        L1TrackInputTag = cms.InputTag("L1Tracks","Level1TkTracks"),
        L1EGammaInputTag = cms.InputTag("l1extraParticles","NonIsolated"),
        label = cms.string("NonIsolated")
)
process.pPhotons = cms.Path( process.L1TrackPhotons )

# "electrons" :

process.L1TrackElectrons = cms.EDProducer("L1TrackElectronParticleProducer",
        L1TrackInputTag = cms.InputTag("L1Tracks","Level1TkTracks"),
        L1EGammaInputTag = cms.InputTag("l1extraParticles","NonIsolated"),
        label = cms.string("NonIsolated")
)
process.pElectrons = cms.Path( process.L1TrackElectrons )


#
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
#
# --- Run the analyzer

process.ana = cms.EDAnalyzer( 'L1TrackTriggerObjectsAnalyzer' ,
    L1VtxInputTag = cms.InputTag("L1TrackPrimaryVertex"),
    L1EtMissInputTag = cms.InputTag("L1TrackEtMiss","MET"),
    L1TrackElectronsInputTag = cms.InputTag("L1TrackElectrons","NonIsolated"),
    L1TrackPhotonsInputTag = cms.InputTag("L1TrackPhotons","NonIsolated")
)


# root file with histograms produced by the analyzer
# (mostly empty currently)

filename = "ana.root"
process.TFileService = cms.Service("TFileService", fileName = cms.string(filename), closeFileFast = cms.untracked.bool(True)
)

process.pAna = cms.Path( process.ana )

#
# ---------------------------------------------------------------------------






process.Out = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "example_all.root" ),
    fastCloning = cms.untracked.bool( False ),
    outputCommands = cms.untracked.vstring( 'drop *')
)

process.Out.outputCommands.append( 'keep *_*_*_ALL' )
process.Out.outputCommands.append( 'keep *_*_*_VTX' )
process.Out.outputCommands.append('keep *_generator_*_*')
process.Out.outputCommands.append('keep *_*gen*_*_*')
process.Out.outputCommands.append('keep *_*Gen*_*_*')
process.Out.outputCommands.append('keep *_rawDataCollector_*_*')
process.Out.outputCommands.append('keep *_L1TkStubsFromPixelDigis_StubsPass_*')

process.FEVToutput_step = cms.EndPath(process.Out)




