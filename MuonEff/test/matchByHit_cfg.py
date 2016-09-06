import FWCore.ParameterSet.Config as cms

process = cms.Process("matchByHit")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

from Configuration.EventContent.EventContent_cff import *
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.fileNames = ['root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/tt8888tt/muon/nopu/0AACA53D-C33D-E611-B180-0025905A612C.root',]

process.MessageLogger.categories = cms.untracked.vstring('testAssociatorRecoMuon', 'muonAssociatorByHitsHelper')
process.MessageLogger.cout = cms.untracked.PSet(
    noTimeStamps = cms.untracked.bool(True),
    threshold = cms.untracked.string('INFO'),
    INFO = cms.untracked.PSet(limit = cms.untracked.int32(0)),
    default = cms.untracked.PSet(limit = cms.untracked.int32(0)),
    testAssociatorRecoMuon = cms.untracked.PSet(limit = cms.untracked.int32(10000000))
)
process.MessageLogger.cerr = cms.untracked.PSet(placeholder = cms.untracked.bool(True))

# Mixing Module
process.load("SimTracker.TrackerMaterialAnalysis.randomNumberGeneratorService_cfi")
process.load("SimGeneral.MixingModule.mixNoPU_cfi")

process.load("CondCore.CondDB.CondDB_cfi")
# Standard Sequences
process.load("Configuration.StandardSequences.Validation_cff")
process.load("Configuration.StandardSequences.GeometryConf")
#process.load("SimGeneral.TrackingAnalysis.trackingParticles_cfi")
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
process.load("CalibCalorimetry.HcalPlugins.Hcal_FrontierConditions_cff")
process.load("Configuration.Geometry.GeometryECALHCAL_cff") 
process.load("Geometry.HcalEventSetup.HcalGeometry_cfi")
process.load("Geometry.HcalEventSetup.HcalTopology_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
process.load("SimGeneral.TrackingAnalysis.trackingParticleNumberOfLayersProducer_cff")
process.load("SimMuon.MCTruth.NewMuonAssociatorByHits_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

#process.GlobalTag.globaltag = cms.string("START37_V3::All")
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.matchByHit = cms.EDAnalyzer("matchByHit",
    vertexs = cms.InputTag("offlinePrimaryVertices"),
    genParticles = cms.InputTag("genParticles"),
    muons = cms.InputTag("muons"),
    muonsTag = cms.InputTag("muons"),
    trackType = cms.string("segments"),  # or 'inner','outer','global'
    tpTag    = cms.InputTag("mix"),
    associatorLabel = cms.string("NewMuonAssociatorByHits"),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("out.root"
))

#process.p = cms.Path(process.matchByHit)
process.skim = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("muons"), minNumber = cms.uint32(1))
process.test = cms.Path(process.skim+process.mix * process.trackingParticleNumberOfLayersProducer* process.matchByHit)
#process.test = cms.Path(process.skim+process.mix * process.trackingParticles* process.matchByHit)
