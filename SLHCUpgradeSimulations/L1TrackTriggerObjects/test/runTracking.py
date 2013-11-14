import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("essai")

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(1),
    limit = cms.untracked.int32(10000000)
)

process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')

	# -- for LB_6PS geometry :
#process.load('Configuration.Geometry.GeometryExtendedPhase2TkLB6PSReco_cff')
#process.load('Configuration.Geometry.GeometryExtendedPhase2TkLB6PS_cff')
	# -- for the Barrel-Endcap geometry:
#process.load('Configuration.Geometry.GeometryExtendedPhase2TkBEReco_cff')
#process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE_cff')
	# -- for BE5D:
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DReco_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5D_cff')

process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedGauss_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring(
    '/store/mc/UpgFall13d/PYTHIA6_Tauola_TTbar_TuneZ2star_14TeV/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/FC5E79F1-F838-E311-B9A0-002618FDA237.root'
    )
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(4)
)


# ---- Global Tag :
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'POSTLS261_V3::All', '')


# ---- L1Tracking :

process.load('Geometry.TrackerGeometryBuilder.StackedTrackerGeometry_cfi')
process.BeamSpotFromSim =cms.EDProducer("BeamSpotFromSimProducer")

process.load("SLHCUpgradeSimulations.L1TrackTrigger.L1TTrack_cfi")
	# -- for LB_6PS :
#process.L1Tracks.geometry = cms.untracked.string('LB_6PS')
	# -- for BE7D :
#process.L1Tracks.geometry = cms.untracked.string('BE')
	# -- for BE5D:
process.L1Tracks.geometry = cms.untracked.string('BE5D')


process.pL1Tracks = cms.Path( process.BeamSpotFromSim*process.L1Tracks )


process.Out = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "file:/tmp/eperez/file_with_L1Tracks.root" ),
    fastCloning = cms.untracked.bool( False ),
	outputCommands = cms.untracked.vstring( 'keep *')
)        

        # generator
process.Out.outputCommands.append('keep *_generator_*_*')

        # L1Tracks
process.Out.outputCommands.append('keep *_L1Tracks_*_*')
process.Out.outputCommands.append('drop *_L1TkClustersFromPixelDigis_*_*')
process.Out.outputCommands.append('keep *_L1TkStubsFromPixelDigis_*_*')
process.Out.outputCommands.append('drop *_L1TkStubsFromPixelDigis_StubsFail_*')


process.FEVToutput_step = cms.EndPath(process.Out)

process.schedule = cms.Schedule(process.pL1Tracks, process.FEVToutput_step)

