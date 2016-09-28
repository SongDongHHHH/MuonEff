import ROOT, copy, os, sys
import MuonEff.MuonEff.CMS_lumi as CMS_lumi
import MuonEff.MuonEff.tdrstyle as tdrstyle
ROOT.gROOT.SetBatch(True)
tdrstyle.setTDRStyle()

def rebin(h,binning):
    xaxis = h.GetXaxis()
    binSize1 = h.GetBinWidth(1)
    binSize2 = abs(binning[2]-binning[1])/float(binning[0])
    h.Rebin(int(binSize2/binSize1))
    xaxis.SetRangeUser(binning[1],binning[2])

def setMarkerStyle(h,color,style):
    h.SetMarkerColor(color)
    h.SetMarkerStyle(style)
    h.SetMarkerSize(1.5)
    h.SetLineColor(color)
    h.SetLineWidth(2)

def getHist(filename,title,id):
    tfile = ROOT.TFile(filename+".root")
    tag = "bestMuon%s5_tpTo%sSel05MuonAssociation/"%(id,id)
    tdir = "DQMData/Run 1/Muons/Run summary/RecoMuonV/MultiTrack/"+tag
    h = tfile.Get(tdir+plotvar)
    h.SetStats(0)
    h.SetTitle(title)
    return copy.deepcopy(h)

def setCanvas(canv,W_ref,H_ref):
    W = W_ref
    H = H_ref
    T = 0.08*H_ref
    B = 0.12*H_ref 
    L = 0.12*W_ref
    R = 0.04*W_ref
    canv.SetFillColor(0)
    canv.SetBorderMode(0)
    canv.SetFrameFillStyle(0)
    canv.SetFrameBorderMode(0)
    canv.SetLeftMargin( L/W )
    canv.SetRightMargin( R/W )
    canv.SetTopMargin( T/H )
    canv.SetBottomMargin( B/H )
    canv.SetTickx(0)
    canv.SetTicky(0)
    return canv

id = sys.argv[1]
datadir = os.environ["CMSSW_BASE"]+'/src/MuonEff/MuonEff/test/'

for plotvar in ["effic","efficPt","fakerate","fakeratePt","invptres_vs_eta_Sigma","invptres_vs_pt_Sigma"]:
    #Get histos
    h_ph2pu0 = getHist("pu0", "PhaseII PU0", id)
    h_ph2pu140 = getHist("pu140", "PhaseII PU140", id)
    h_ph2pu200 = getHist("pu200", "PhaseII PU200", id)
    hlist = [h_ph2pu0, h_ph2pu140, h_ph2pu200]

    #Set init histo
    nbins = h_ph2pu0.GetNbinsX()
    h_init = ROOT.TH1F("","",nbins,h_ph2pu0.GetBinLowEdge(1),h_ph2pu0.GetBinLowEdge(nbins+1))

    #Set axis
    x_name = "Muon "
    if ("Pt" in plotvar) or ("_pt" in plotvar):
        x_name = x_name+"p_{T}"
    else:
        x_name = x_name+"|#eta|"

    y_name = id+" Muon "
    if "eff" in plotvar:
        h_init.SetMaximum(1.2)
        h_init.SetMinimum(0.6)
        y_name = y_name+"Efficiency"
    elif "fake" in plotvar:
        h_init.SetMaximum(max(h.GetMaximum() for h in hlist)*1.8)
        y_name = y_name+"Fake Rate"
    else:
        y_name = "#sigma(p_{T})/p_{T}"
    h_init.GetXaxis().SetTitle(x_name)
    h_init.GetYaxis().SetTitle(y_name)
    h_init.GetYaxis().SetTitleOffset(0.98)


    ############ Plot design ##############
    name = "%s_%s"%(plotvar,id)

    iPos = 0
    iPeriod = 0
    if( iPos==0 ): CMS_lumi.relPosX = 0.12

    #Set canvas
    H_ref = 600; 
    W_ref = 800; 
    W = W_ref
    H = H_ref
    canv = ROOT.TCanvas(name,name,50,50,W,H)
    setCanvas(canv,W_ref,H_ref)
    h_init.Draw()

    #Plot style
    setMarkerStyle(h_ph2pu0, 4, 20) #blue, circle
    setMarkerStyle(h_ph2pu140, 1, 34) #black, cross
    setMarkerStyle(h_ph2pu200, 2, 21) #red, square

    #Legend and drawing
    leg = ROOT.TLegend(0.6,0.72,0.85,0.88)
    for h in hlist:
        h.Draw("e1same")
        leg.AddEntry(h,h.GetTitle(),"p")
    leg.SetTextFont(61)
    leg.SetTextSize(0.04)
    leg.SetBorderSize(0)
    leg.Draw()

    #Sample name text
    tex2 = ROOT.TLatex()
    tex2.SetNDC()
    tex2.SetTextFont(42)
    tex2.SetTextSize(0.04)
    tex2.DrawLatex(0.18, 0.8, "Z/#gamma^{*}#rightarrow#font[12]{#mu#mu}, p_{T} > 5 GeV")

    #CMS_lumi setting
    CMS_lumi.extraText = "Simulation"
    CMS_lumi.lumi_sqrtS = "14 TeV"
    CMS_lumi.CMS_lumi(canv, iPeriod, iPos)

    if "Sigma" in plotvar: canv.SetLogy()
    canv.Modified()
    canv.Update()
    canv.SaveAs(name+".png")

