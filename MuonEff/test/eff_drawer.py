import ROOT, copy

def rebin(h, binning):
    xaxis = h.GetXaxis()
    binSize1 = h.GetBinWidth(1)
    binSize2 = abs(binning[2]-binning[1])/float(binning[0])

    h.Rebin(int(binSize2/binSize1))
    xaxis.SetRangeUser(binning[1],binning[2])

def histinit(h, binning):
    h.SetStats(0)
    rebin(h, binning)

def setHistStyle(h,var,i):
    if "eff" in var:
        h.SetMaximum(1.05)
        h.SetMinimum(0.4)
        h.SetFillColor(i+2)
        h.SetLineColor(0)
        h.SetMarkerStyle(0)
        drawoption="e5same"
    if "fake" in var:
        h.SetMaximum(max(h.GetMaximum() for h in hlist)*2)
        h.SetLineColor(i+2)
        h.SetLineWidth(2)
        drawoption="histsame"
    return drawoption
    

filenames=["pu0","pu140","pu200"]
plotvar=["efficPt","effic","fakeratePt","fakerate"]
ids = ["Tight","Loose"]

for i, p in enumerate(plotvar):
    for filename in filenames:
        hlist = []
        for id in ids:
            tfile = ROOT.TFile(filename+".root")
            tdir="DQMData/Run 1/Muons/Run summary/RecoMuonV/MultiTrack/bestMuon%s5_tpTo%sSel05MuonAssociation/"%(id,id)
            h = tfile.Get(tdir+p)
            h.SetStats(0)
            hlist.append(copy.deepcopy(h))

        cn = ROOT.TCanvas()
        leg = ROOT.TLegend(0.55,0.4,0.8,0.58)

        for i,h in enumerate(hlist):
            drawoption = setHistStyle(h,p,i)
            h.Draw(drawoption)
            leg.AddEntry(h, h.GetTitle(), "f")

        leg.SetBorderSize(0)
        leg.Draw()
        cn.Print("%s_%s.png"%(filename,p))

