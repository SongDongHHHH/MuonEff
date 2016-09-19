import FWCore.ParameterSet.Config as cms

process = cms.Process("MuonTMVA")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("SimMuon.MCTruth.MuonAssociatorByHits_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.fileNames = ['root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/tt8888tt/muon/nopu/0AACA53D-C33D-E611-B180-0025905A612C.root',
                            'root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/tt8888tt/muon/nopu/34742466-C03D-E611-AA9B-0025905B85A2.root',
                            'root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/tt8888tt/muon/nopu/40245918-C33D-E611-A87C-0025905A6066.root',
                            'root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/tt8888tt/muon/nopu/4640E9AE-B63D-E611-88D4-0CC47A78A4A6.root',
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
"""


process.MuonTMVA = cms.EDAnalyzer("MuonTMVA",
    vertexs = cms.InputTag("offlinePrimaryVertices"),
    genParticles = cms.InputTag("genParticles"),
    muons = cms.InputTag("muons"),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("out.root"
))


process.p = cms.Path(process.MuonTMVA)
