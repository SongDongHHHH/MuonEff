import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("MuonEff")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())

#process.source.fileNames = ['root://cms-xrdr.sdfarm.kr:1094///xrd/store/relval/CMSSW_8_1_0_pre10/RelValZMM_13/GEN-SIM-RECO/81X_mcRun2_asymptotic_v5_2023D1-v1/00000/2251A418-9868-E611-B8F6-0CC47A4C8E26.root']
dir = os.environ["CMSSW_BASE"]+'/src/MuonEff/MuonEff/doc/'
filelst = open(dir+"pu140.txt", "r")
process.source.fileNames = filelst.readlines()

process.MuonEff = cms.EDAnalyzer("MuonEff",
    vertexs = cms.InputTag("offlinePrimaryVertices"),
    genParticles = cms.InputTag("genParticles"),
    muons = cms.InputTag("muons"),
    MuonGEMHits = cms.InputTag("g4SimHits", "MuonGEMHits")
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("out.root"
))

process.p = cms.Path(process.MuonEff)
