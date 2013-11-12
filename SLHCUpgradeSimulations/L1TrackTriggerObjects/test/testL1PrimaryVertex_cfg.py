import FWCore.ParameterSet.Config as cms

process = cms.Process("VTX")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#'/store/mc/Summer13/PYTHIA6_Tauola_TTbar_TuneZ2star_14TeV/GEN-SIM-DIGI-RAW/UpgradePhase2BE_2013_DR61SLHCx_PU140Bx25_POSTLS261_V2-v1/10000/F82A0812-97C4-E211-A076-00304867926C.root'
    'file:/tmp/eperez/file_with_L1Tracks.root'
    )
)


# ---- Global Tag :
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS1', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'POSTLS261_V3::All', '')

process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DReco_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5D_cff')

process.load('Geometry.TrackerGeometryBuilder.StackedTrackerGeometry_cfi')


process.vtx = cms.EDProducer('L1TrackPrimaryVertexProducer',
     ZMAX = cms.double ( 25. ) ,
     ZSTEP = cms.double( 0.05 ),
     CHI2MAX = cms.double( 100. )
)

process.p = cms.Path( process.vtx )

process.Out = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "example.root" ),
    fastCloning = cms.untracked.bool( False ),
    outputCommands = cms.untracked.vstring( 'drop *')
)

process.Out.outputCommands.append( 'keep *_*_*_VTX' )
process.FEVToutput_step = cms.EndPath(process.Out)




