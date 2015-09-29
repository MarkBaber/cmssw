import FWCore.ParameterSet.Config as cms

process = cms.Process('TESTING')

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContentCosmics_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EDMtoMEAtRunEnd_cff')
process.load('Configuration.StandardSequences.HarvestingCosmics_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load("DQMServices.Components.MEtoEDMConverter_cfi")

process.load('HLTriggerOffline.SUSYBSM.SusyExoValidation_cff')
process.load('HLTriggerOffline.SUSYBSM.SUSYBSM_alphaT_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    # fileNames = cms.untracked.vstring(
    #     '/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/AODSIM/AVE30BX50_tsg_PHYS14_ST_V1-v1/00000/4E9B4A18-2F8F-E411-93FB-0025905AA9CC.root',
    #     '/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/AODSIM/AVE30BX50_tsg_PHYS14_ST_V1-v1/00000/BC7F361B-2F8F-E411-9DC0-0025905B860E.root',
    #     ),
    # secondaryFileNames = cms.untracked.vstring(
    #     'file:/afs/cern.ch/user/m/mbaber/WORK/private/DQM/DQM_090915/CMSSW_7_4_10_patch1/src/outputDQM.root'),
 fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/user/m/mbaber/WORK/private/DQM/DQM_090915/CMSSW_7_4_10_patch1/src/outputDQM.root'),
)




from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup', '')

process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_MEtoEDMConverter_*_*'
    ),
    fileName = cms.untracked.string('file:DQM.root'),
)


process.HLTAlphaTSeq = cms.Sequence(process.SUSY_HLT_HT200_alphaT0p63+
                                    process.SUSY_HLT_HT250_alphaT0p58+
                                    process.SUSY_HLT_HT300_alphaT0p53+ 
                                    process.SUSY_HLT_HT350_alphaT0p52+ 
                                    process.SUSY_HLT_HT400_alphaT0p51+ 
                                    process.SUSY_HLT_HT200_alphaT0p57+
                                    process.SUSY_HLT_HT250_alphaT0p55+
                                    process.SUSY_HLT_HT300_alphaT0p54+ 
                                    process.SUSY_HLT_HT350_alphaT0p53+ 
                                    process.SUSY_HLT_HT400_alphaT0p52+ 
                                    process.SUSY_HLT_HT200_alphaT0p51
                                    )

process.run_module = cms.Path(process.HLTAlphaTSeq+process.MEtoEDMConverter)
process.outpath    = cms.EndPath(process.out)
process.schedule   = cms.Schedule(process.run_module, process.outpath)



