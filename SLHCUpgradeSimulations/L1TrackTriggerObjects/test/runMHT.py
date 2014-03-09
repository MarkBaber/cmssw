# Import configurations
import FWCore.ParameterSet.Config as cms

# set up process
process = cms.Process("L1EG")

#from SLHCUpgradeSimulations.L1TrackTriggerObjects.minBiasFiles_p1_cfi import *


process.source = cms.Source("PoolSource",
   #fileNames = minBiasFiles_p1
   fileNames = cms.untracked.vstring(
     #"root://eoscms///store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/008E2E98-0A39-E311-833F-0025905938D4.root"
    #'/store/mc/UpgFall13d/HToTauTau_125_14TeV_powheg_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/FAF4B2D5-0539-E311-B5C3-002618FDA279.root'
    #'file:example_w_Tracks_and_vertex.root'
   '/store/cmst3/user/eperez/L1TrackTrigger/612_SLHC6/muDST/TTbar/BE5D/zmatchingOff/m1_TTbar_BE5D.root'
   )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )


# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport = cms.untracked.PSet(
    #reportEvery = cms.untracked.int32(500),
    reportEvery = cms.untracked.int32(10),
    limit = cms.untracked.int32(10000000)
)      
       
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

# Load geometry
#process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DReco_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5D_cff')
process.load('Geometry.TrackerGeometryBuilder.StackedTrackerGeometry_cfi')

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
                            
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'POSTLS261_V3::All', '')


process.load("Configuration.StandardSequences.Services_cff")
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/L1HwVal_cff')
process.load("Configuration.StandardSequences.RawToDigi_Data_cff") ###check this for MC!
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
#process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
process.load("JetMETCorrections.Configuration.DefaultJEC_cff")


# ---------------------------------------------------------------------------
#
# ---  Run the SLHCCaloSequence

process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load("SLHCUpgradeSimulations.L1CaloTrigger.SLHCCaloTrigger_cff")


process.p = cms.Path( 
    process.RawToDigi+
    process.SLHCCaloTrigger
    )

# bug fix for missing HCAL TPs in MC RAW
process.p.insert(1, process.valHcalTriggerPrimitiveDigis)
from SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff import HcalTPGCoderULUT
HcalTPGCoderULUT.LUTGenerationMode = cms.bool(True)
process.valRctDigis.hcalDigis             = cms.VInputTag(cms.InputTag('valHcalTriggerPrimitiveDigis'))
process.L1CaloTowerProducer.HCALDigis =  cms.InputTag("valHcalTriggerPrimitiveDigis")

	# run L1Reco to produce the L1EG objects corresponding
	# to the current trigger
#process.load('Configuration.StandardSequences.L1Reco_cff')
#process.L1Reco = cms.Path( process.l1extraParticles )

#
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
#
# -- The HLT HI jets :

# --- Run the calo local reconstruction
process.towerMaker.hbheInput = cms.InputTag("hbheprereco")
process.towerMakerWithHO.hbheInput = cms.InputTag("hbheprereco")
process.reconstruction_step = cms.Path( process.calolocalreco )


# --- Produce the  HLT HeavyIon jets :
process.load("RecoHI.HiJetAlgos.HiRecoJets_TTI_cff")
process.hireco = cms.Path( process.hiRecoJets )

# --- Put them into "L1Jets"
process.L1JetsFromHIHLTJets = cms.EDProducer("L1JetsFromHIHLTJets",
        ETAMIN = cms.double(0),
        ETAMAX = cms.double(3.),
        HIJetsInputTag = cms.InputTag("iterativeConePu5CaloJets")
)
process.pL1Jets = cms.Path( process.L1JetsFromHIHLTJets )


# --- Produce L1TkJets from the HeavyIon jets
process.L1TkJets = cms.EDProducer("L1TkJetProducer",
        L1CentralJetInputTag = cms.InputTag("L1JetsFromHIHLTJets"), 
        L1TrackInputTag = cms.InputTag("L1Tracks","Level1TkTracks"),
)
process.pJets = cms.Path( process.L1TkJets )

# --- Produce L1TkJets from the L1Jets of the SLHCCalo sequence:
process.L1TkJetsL1 = cms.EDProducer("L1TkJetProducer",
        L1CentralJetInputTag = cms.InputTag("L1TowerJetPUSubtractedProducer","PUSubCen8x8"),
        L1TrackInputTag = cms.InputTag("L1Tracks","Level1TkTracks"),
)
process.pJetsL1 = cms.Path( process.L1TkJetsL1 )


#
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
#
# -- Now produce MHT (and HT)

	# from the HLT Heavy Ion jets :
process.L1TkHTMissHI = cms.EDProducer("L1TkHTMissProducer",
        L1TkJetInputTag = cms.InputTag("L1TkJets","Central"),
        DeltaZ = cms.double( 999 ),   #  in mm. Here dummy cut, since I dont have the zvtx of the jets 
        PrimaryVtxConstrain = cms.bool( True ),
        #L1VertexInputTag = cms.InputTag("L1TkPrimaryVertex")
        L1VertexInputTag = cms.InputTag("L1TrackPrimaryVertex")
)
process.pHTM = cms.Path( process.L1TkHTMissHI )


	# from the L1jets created by the SLHCCaloSequence
process.L1TkHTMissL1 = cms.EDProducer("L1TkHTMissProducer",
        L1TkJetInputTag = cms.InputTag("L1TkJetsL1","Central"),
        DeltaZ = cms.double( 999 ),   #  in mm. Here dummy cut, since I dont have the zvtx of the jets 
        PrimaryVtxConstrain = cms.bool( True ),
        #L1VertexInputTag = cms.InputTag("L1TkPrimaryVertex")
        L1VertexInputTag = cms.InputTag("L1TrackPrimaryVertex")
)
process.pHTM2 = cms.Path( process.L1TkHTMissL1 )







process.Out = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "MHT.root" ),
    fastCloning = cms.untracked.bool( False ),
    outputCommands = cms.untracked.vstring( 'drop *')
)

#process.Out.outputCommands.append( 'keep *_*_*_L1EG' )
process.Out.outputCommands.append('keep *_generator_*_*')
process.Out.outputCommands.append('keep *_SLHCL1ExtraParticles_*_*')
process.Out.outputCommands.append('keep *_SLHCL1ExtraParticlesNewClustering_*_*')
process.Out.outputCommands.append('keep *_l1extraParticles_*_*')
process.Out.outputCommands.append('keep *_L1TkHTMiss*_*_*')
process.Out.outputCommands.append('keep *_L1TkJets_*_*')
process.Out.outputCommands.append('keep *_L1TkJetsL1_*_*')
process.Out.outputCommands.append('keep *_iterativeConePu5CaloJets_*_*')


#process.Out.outputCommands.append('keep *')

process.FEVToutput_step = cms.EndPath(process.Out)



