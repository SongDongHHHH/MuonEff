import FWCore.ParameterSet.Config as cms

process = cms.Process("MuonVali")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

from Configuration.EventContent.EventContent_cff import *
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.fileNames = ['root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/tt8888tt/muon/nopu/0AACA53D-C33D-E611-B180-0025905A612C.root',]

process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('drop *', "keep *_MEtoEDMConverter_*_"+processName),
    fileName = cms.untracked.string('validationEDM.root')
)
process.outpath = cms.EndPath(process.out)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.load('Configuration/StandardSequences/RawToDigi_cff')
process.raw2digi_step = cms.Path(process.RawToDigi)

process.load("Configuration/StandardSequences/SimulationRandomNumberGeneratorSeeds_cff")

process.load("DQMServices.Components.MEtoEDMConverter_cfi")
process.MEtoEDMConverter_step = cms.Path(process.MEtoEDMConverter)

process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration.Geometry.GeometryExtended2023MuonReco_cff')
process.load('Configuration.Geometry.GeometryExtended2023Muon_cff')
#process.load('Configuration.Geometry.GeometryExtended2019Reco_cff')
#process.load('Configuration.Geometry.GeometryExtended2019_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgrade2019', '')

#---- Validation stuffs ----#
## Default validation modules
process.load("Configuration.StandardSequences.Validation_cff")
process.validation_step = cms.Path(process.validation)
## Load muon validation modules
#process.recoMuonVMuAssoc.outputFileName = 'validationME.root'
process.muonValidation_step = cms.Path(process.recoMuonValidation)

process.load("SimMuon.MCTruth.MuonAssociatorByHits_cfi")
process.muonAssociatorByHitsCommonParameters.useGEMs = cms.bool(True)

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.MuonVali = cms.EDAnalyzer("MuonVali",
    maxEventsToPrint = cms.untracked.int32(50),
    printVertex = cms.untracked.bool(False),
    printOnlyHardInteraction = cms.untracked.bool(False), # Print only status=3 particles. This will not work for Pythia8, which does not have any such particles.
    src = cms.InputTag("genParticles")
)
process.p = cms.Path(process.MuonVali)

