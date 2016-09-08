import FWCore.ParameterSet.Config as cms

process = cms.Process("MuonDetail")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("SimMuon.MCTruth.MuonAssociatorByHits_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.fileNames = ['root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/tt8888tt/muon/nopu/0AACA53D-C33D-E611-B180-0025905A612C.root',
                            #'root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/tt8888tt/muon/nopu/34742466-C03D-E611-AA9B-0025905B85A2.root',
                            #'root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/tt8888tt/muon/nopu/40245918-C33D-E611-A87C-0025905A6066.root',
                            #'root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/tt8888tt/muon/nopu/4640E9AE-B63D-E611-88D4-0CC47A78A4A6.root',
                            'root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/tt8888tt/muon/nopu/D243AE89-B93D-E611-96AC-0CC47A4D7626.root',]
"""
process.source.fileNames = ['root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/tt8888tt/muon/pu35/063C4686-F73E-E611-BC5A-0025905B85BA.root',
                            'root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/tt8888tt/muon/pu35/0EC6E088-1D3F-E611-8BEE-0CC47A4C8EC8.root',
                            'root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/tt8888tt/muon/pu35/1095AD36-F03E-E611-960C-0025905B85B8.root',
                            'root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/tt8888tt/muon/pu35/2A07DF30-F73E-E611-A74F-0025905B8576.root',
                            'root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/tt8888tt/muon/pu35/301357F1-E83E-E611-B59A-0CC47A78A42C.root',
                            'root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/tt8888tt/muon/pu35/301720A3-FB3E-E611-A2D2-0CC47A78A41C.root',
                            'root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/tt8888tt/muon/pu35/44731022-EF3E-E611-9021-0025905A611E.root',
                            'root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/tt8888tt/muon/pu35/7E91C8CA-F93E-E611-8406-0CC47A4C8E34.root',
                            'root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/tt8888tt/muon/pu35/8AF3EB07-F93E-E611-BB46-0025905B860E.root',
                            'root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/tt8888tt/muon/pu35/963BBFB9-F73E-E611-8A80-0025905B85BC.root',
                            'root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/tt8888tt/muon/pu35/AADC239E-F93E-E611-A7D6-0CC47A745250.root',
                            'root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/tt8888tt/muon/pu35/B0A81A83-EC3E-E611-8B4F-0CC47A4D76B8.root',
                            'root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/tt8888tt/muon/pu35/BCCB7CB7-F93E-E611-8BDF-0CC47A4D7600.root',
                            'root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/tt8888tt/muon/pu35/CA8F1FDC-F13E-E611-B5E1-0CC47A4C8EE2.root',
                            'root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/tt8888tt/muon/pu35/E24B1756-ED3E-E611-913E-0025905A60DA.root',
                            'root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/tt8888tt/muon/pu35/E2B08273-F03E-E611-BD0C-0CC47A4C8F30.root',
                            'root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/tt8888tt/muon/pu35/E8483BF0-EF3E-E611-A7E9-0CC47A78A30E.root',
                            'root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/tt8888tt/muon/pu35/ECC3668F-EC3E-E611-A407-0025905B85DE.root',
                            'root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/tt8888tt/muon/pu35/F6544920-E83E-E611-B935-002618FDA287.root',]
process.source.fileNames = ['root://cms-xrdr.sdfarm.kr:1094///xrd/store/relval/CMSSW_8_1_0_pre8/RelValZMM_13/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v1-v1/10000/04B1EEE2-9744-E611-B590-0025905B856E.root',
                            'root://cms-xrdr.sdfarm.kr:1094///xrd/store/relval/CMSSW_8_1_0_pre8/RelValZMM_13/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v1-v1/10000/0AA6BD15-8C44-E611-8EA2-0CC47A4C8EA8.root',
                            'root://cms-xrdr.sdfarm.kr:1094///xrd/store/relval/CMSSW_8_1_0_pre8/RelValZMM_13/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v1-v1/10000/4A0C22A2-8B44-E611-9BF6-0025905A612A.root',
                            'root://cms-xrdr.sdfarm.kr:1094///xrd/store/relval/CMSSW_8_1_0_pre8/RelValZMM_13/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v1-v1/10000/64E5DDCA-8B44-E611-A17F-0025905A6138.root',
                            'root://cms-xrdr.sdfarm.kr:1094///xrd/store/relval/CMSSW_8_1_0_pre8/RelValZMM_13/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v1-v1/10000/66687402-9844-E611-B599-0CC47A4D76BE.root',
                            'root://cms-xrdr.sdfarm.kr:1094///xrd/store/relval/CMSSW_8_1_0_pre8/RelValZMM_13/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v1-v1/10000/7AF34711-8A44-E611-96C0-0CC47A74525A.root',
                            'root://cms-xrdr.sdfarm.kr:1094///xrd/store/relval/CMSSW_8_1_0_pre8/RelValZMM_13/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v1-v1/10000/DE591F56-8A44-E611-BCEA-0025905A60B4.root',
                            'root://cms-xrdr.sdfarm.kr:1094///xrd/store/relval/CMSSW_8_1_0_pre8/RelValZMM_13/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v1-v1/10000/F27D7C1B-8A44-E611-9929-0025905A612A.root']
"""


process.MuonDetail = cms.EDAnalyzer("MuonDetail",
    vertexs = cms.InputTag("offlinePrimaryVertices"),
    genParticles = cms.InputTag("genParticles"),
    muons = cms.InputTag("muons"),
    isFake = cms.bool(True),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("out.root"
))


process.p = cms.Path(process.MuonDetail)
