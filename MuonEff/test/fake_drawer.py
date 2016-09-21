import ROOT, copy, os

def getTH1(tree, title, binning, plotvar, cut):
    hist = ROOT.TH1D("name", title, binning[0], binning[1], binning[2])
    tree.Project("name", plotvar, cut)
    hist.SetStats(0)
    #hist.Sumw2()
    hist.SetLineWidth(3)
    hist.Scale(1/hist.GetEntries())
    return copy.deepcopy(hist)

plotvarlst = ['isGlobalMuon',
              'isPFMuon',
              'normalizedChi2',
              'chi2LocalPosition',
              'trkKink',
              'segmentCompatibility',
              'numberOfValidMuonHits',
              'numberOfMatchedStations',
              'pv0pos_dxy',
              'pv0pos_dz',
              'numberOfValidPixelHits',
              'trackerLayersWithMeasurement']
binninglst = [[2,0,2], [2,0,2], [40,0,10], [40,0,20], [50,0,50], [50,0,1], [60,0,60], [6,0,6], [50,0,1], [50,0,5], [15,0,15], [15,0,15]]
cut = [1,1,0,0,0,10.,0,1,0.2,0.5,0,5]

names=["nopu","pu35","pu140","pu200"]
datadir = os.environ["CMSSW_BASE"]+'/src/MuonEff/MuonEff/test/formva/'
for name in names:
    for i, plotvar in enumerate(plotvarlst):
    #for i, plotvar in enumerate(['MuonPt','MuonEta','MuonPhi']):
        title = "%s_%s"%(name, plotvar)
        tfile = ROOT.TFile(datadir+"%s.root"%name)
        tree = tfile.Get("MuonDetail/fake")
        h_fake = getTH1(tree, title, binninglst[i], plotvar, "1")
        h_fake.SetLineColor(2)

        tree = tfile.Get("MuonDetail/true")
        h_true = getTH1(tree, title, binninglst[i], plotvar, "1")
        h_true.SetLineColor(4)

        cn = ROOT.TCanvas()
        h_fake.SetMaximum(max(h_true.GetMaximum(), h_fake.GetMaximum())*1.2)
        h_fake.Draw()
        h_true.Draw("same")
        leg = ROOT.TLegend(0.65,0.8,0.85,0.6)
        leg.SetTextSize(0.05)
        leg.AddEntry(h_fake,"Fake","l")
        leg.AddEntry(h_true,"True","l")
        leg.Draw()
        cn.Print("%s.png"%title)

