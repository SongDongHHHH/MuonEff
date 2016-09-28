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

#binninglst = [[30,0,60],[12,-2.4,2.4],[6,-3,3]]
binninglst = [[20,5,105],[24,-2.4,2.4],[30,-3,3]]
filenames=["run2","nopu","pu140","pu200"]
for filename in filenames:
    tfile = ROOT.TFile("../"+filename+".root")
    nevents = tfile.Get("MuonEff/nevents").GetEntries()
    for i, plotvar in enumerate(['pt','eta','phi']):
        gen = tfile.Get("MuonEff/gen_%s"%(plotvar))
        histinit(gen,binninglst[i])
        gen.Sumw2()
        eff=gen.Clone()
        lst=[]
        fakelst=[]
        for id in ['tight','medium','loose']:
            reco = tfile.Get("MuonEff/gen%s_%s"%(id[0].capitalize(),plotvar))
            histinit(reco,binninglst[i])
            reco.Sumw2()
            eff.Divide(reco, gen, 1, 1, "B")
            eff.SetTitle(id)
            lst.append(copy.deepcopy(eff))

            fake = tfile.Get("MuonEff/fake%s_%s"%(id[0].capitalize(),plotvar))
            fake.Scale(1/nevents)
            histinit(fake,binninglst[i])
            fake.SetTitle(id)
            fakelst.append(copy.deepcopy(fake))

        eff.Reset()
        eff.SetMaximum(1.05)
        eff.SetMinimum(0.4)
        eff.SetLineColor(0)
        eff.SetTitle("%s_%s"%(filename,plotvar))

        drawoption = "e5same"
        if 'pu' not in filename: drawoption = "csame"
        cn = ROOT.TCanvas()
        eff.Draw()
        leg = ROOT.TLegend(0.58,0.14,0.88,0.4)
        for i, h in enumerate(lst):
            h.SetFillColor(i+2)
            h.SetLineColor(i+2)
            h.SetLineWidth(2)
            h.SetFillStyle(3001)
            h.SetMarkerStyle(0)
            h.Draw(drawoption)
            leg.AddEntry(h,h.GetTitle(),"f")

        leg.SetBorderSize(0)
        leg.Draw()

        tex = ROOT.TLatex()
        tex.SetNDC()
        tex.SetTextFont(42)
        tex.SetTextSize(0.03)
        tex.DrawLatex(0.2, 0.3, "#frac{Number of matched gen muon}{Number of total gen muon}")

        cn.Print("%s_%s.png"%(filename,plotvar))

        cn2 = ROOT.TCanvas()
        leg = ROOT.TLegend(0.58,0.58,0.88,0.88)
        #leg.SetNColumns(2)
        fake = copy.deepcopy(fakelst[0])
        fake.SetLineColor(0)
        fake.SetTitle("fake %s %s"%(filename,plotvar))
        fake.SetMaximum(max(h.GetMaximum() for h in fakelst)*2.2)
        fake.Draw()
        for i, h in enumerate(fakelst):
            h.SetLineColor(i+2)
            h.SetLineWidth(2)
            h.Draw("same")
            leg.AddEntry(h,h.GetTitle(),"f")

        leg.SetBorderSize(0)
        leg.Draw()

        tex = ROOT.TLatex()
        tex.SetNDC()
        tex.SetTextFont(42)
        tex.SetTextSize(0.03)
        tex.DrawLatex(0.15, 0.75, "#frac{Number of non matched reco muon}{Number of total reco muon}")

        cn2.Print("fake_%s_%s.png"%(filename,plotvar))


