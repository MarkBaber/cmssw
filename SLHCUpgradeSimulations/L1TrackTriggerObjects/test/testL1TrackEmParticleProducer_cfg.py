import FWCore.ParameterSet.Config as cms

process = cms.Process("Ele")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    'file:example_w_Tracks_and_vertex.root'
    #'/store/cmst3/user/eperez/L1TrackTrigger/612_SLHC6/muDST/TTbar/BE5D/m1_TTbar_BE5D.root'
    )
)


process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'POSTLS261_V3::All', '')


# --- creates l1extra objects for L1EGamma (here for the Run1 trigger !)

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


process.Out = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "example.root" ),
    fastCloning = cms.untracked.bool( False ),
    outputCommands = cms.untracked.vstring( 'drop *')
)

#process.Out.outputCommands.append( 'keep *_*_*_Ele' )
#process.Out.outputCommands.append('keep *_generator_*_*')
process.Out.outputCommands.append('keep *')

process.FEVToutput_step = cms.EndPath(process.Out)




