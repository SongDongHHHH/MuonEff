#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"


#include <TTree.h>
#include <TFile.h>
#include <TLorentzVector.h>
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class MuonDetail : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit MuonDetail(const edm::ParameterSet&);
      ~MuonDetail();
      TLorentzVector ToTLorentzVector(const reco::Muon& t) { return TLorentzVector(t.px(), t.py(), t.pz(), t.energy()); }
      TLorentzVector ToTLorentzVector(const reco::GenParticle& t) { return TLorentzVector(t.px(), t.py(), t.pz(), t.energy()); }

      void collectGenMuons(const std::vector<reco::GenParticle> &genParticles, std::vector<reco::GenParticle>& genMuons, std::vector<reco::GenParticle>& genMuonsFromZ) const;
      void treereset();
      void fillHitMap(const reco::Muon& mu, TH1D* h_hitmap) const;

   private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

      TTree * tr_true;
      TTree * tr_fake;
      bool b_global; bool b_pf;
      float b_chi2pos; float b_trkKink; float b_segcompati;
      float b_chi2; int b_nglobalhits; int b_nstations;
      float b_trackdxy; float b_trackdz;
      int b_ninnerhits; float b_trackerlayers;
      float b_mu_pt; float b_mu_eta; float b_mu_phi;

      edm::EDGetTokenT<std::vector<reco::GenParticle> > genParticleToken_;
      edm::EDGetTokenT<std::vector<reco::Muon> > muonToken_;
      edm::EDGetTokenT<std::vector<reco::Vertex> > vtxToken_;
};

using namespace std;
using namespace reco;
using namespace edm;

MuonDetail::MuonDetail(const edm::ParameterSet& iConfig)
{
  genParticleToken_ = consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticles"));
  muonToken_ = consumes<std::vector<reco::Muon> >(iConfig.getParameter<edm::InputTag>("muons"));
  vtxToken_ = consumes<std::vector<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("vertexs"));

  usesResource("TFileService");
  edm::Service<TFileService> fs;
  tr_true = fs->make<TTree>("true", "true");
  tr_fake = fs->make<TTree>("fake", "fake");

  tr_true->Branch("muonPt", &b_mu_pt, "muonPt/F");
  tr_true->Branch("muonEta", &b_mu_eta, "muonEta/F");
  tr_true->Branch("muonPhi", &b_mu_phi, "muonPhi/F");
  tr_true->Branch("isGlobalMuon", &b_global, "isGlobalMuon/O");
  tr_true->Branch("isPFMuon", &b_pf, "isPFMuon/O");
  tr_true->Branch("chi2LocalPosition", &b_chi2pos, "chi2LocalPosition/F");
  tr_true->Branch("trkKink", &b_trkKink, "trkKink/F");
  tr_true->Branch("segmentCompatibility", &b_segcompati, "segmentCompatibility/F");
  tr_true->Branch("normalizedChi2", &b_chi2, "normalizedChi2/F");
  tr_true->Branch("numberOfValidMuonHits", &b_nglobalhits, "numberOfValidMuonHits/I");
  tr_true->Branch("numberOfMatchedStations", &b_nstations, "numberOfMatchedStations/I");
  tr_true->Branch("pv0pos_dxy", &b_trackdxy, "pv0pos_dxy/F");
  tr_true->Branch("pv0pos_dz", &b_trackdz, "pv0pos_dz/F");
  tr_true->Branch("numberOfValidPixelHits", &b_ninnerhits, "numberOfValidPixelHits/I");
  tr_true->Branch("trackerLayersWithMeasurement", &b_trackerlayers, "trackerLayersWithMeasurement/F");

  tr_fake->Branch("muonPt", &b_mu_pt, "muonPt/F");
  tr_fake->Branch("muonEta", &b_mu_eta, "muonEta/F");
  tr_fake->Branch("muonPhi", &b_mu_phi, "muonPhi/F");
  tr_fake->Branch("isGlobalMuon", &b_global, "isGlobalMuon/O");
  tr_fake->Branch("isPFMuon", &b_pf, "isPFMuon/O");
  tr_fake->Branch("chi2LocalPosition", &b_chi2pos, "chi2LocalPosition/F");
  tr_fake->Branch("trkKink", &b_trkKink, "trkKink/F");
  tr_fake->Branch("segmentCompatibility", &b_segcompati, "segmentCompatibility/F");
  tr_fake->Branch("normalizedChi2", &b_chi2, "normalizedChi2/F");
  tr_fake->Branch("numberOfValidMuonHits", &b_nglobalhits, "numberOfValidMuonHits/I");
  tr_fake->Branch("numberOfMatchedStations", &b_nstations, "numberOfMatchedStations/I");
  tr_fake->Branch("pv0pos_dxy", &b_trackdxy, "pv0pos_dxy/F");
  tr_fake->Branch("pv0pos_dz", &b_trackdz, "pv0pos_dz/F");
  tr_fake->Branch("numberOfValidPixelHits", &b_ninnerhits, "numberOfValidPixelHits/I");
  tr_fake->Branch("trackerLayersWithMeasurement", &b_trackerlayers, "trackerLayersWithMeasurement/F");
}

MuonDetail::~MuonDetail(){}

void MuonDetail::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices); 
  if (vertices->empty()) { cout << "noPV" << endl; return; }
  auto pv0 = vertices->front();

  ///// gen muon select
  std::vector<reco::GenParticle> genMuons;
  std::vector<reco::GenParticle> genMuonsFromZ;
  edm::Handle<std::vector<reco::GenParticle>> genParticles;
  iEvent.getByToken(genParticleToken_, genParticles);
  collectGenMuons(*genParticles, genMuons, genMuonsFromZ);

  edm::Handle<std::vector<reco::Muon> > muons;
  iEvent.getByToken(muonToken_, muons);
  double rcut = 0.3;
  for (auto& mu : *muons) {
    bool fake = true;
    for (auto& gen : genMuons){
      double dR = reco::deltaR(gen, mu);
      if (dR > rcut) continue;
      fake = false;
    }

    treereset();
    b_mu_pt = mu.pt(); b_mu_eta = mu.eta(); b_mu_phi = mu.phi();

    //loose
    b_global = mu.isGlobalMuon();
    b_pf = mu.isPFMuon();

    //medium
    if ( mu.globalTrack().isNonnull() ){
        b_chi2 = mu.globalTrack()->normalizedChi2();
    }
    b_chi2pos = mu.combinedQuality().chi2LocalPosition;
    b_trkKink = mu.combinedQuality().trkKink;
    b_segcompati = muon::segmentCompatibility(mu);

    //tight
    b_nstations = mu.numberOfMatchedStations();
    if ( mu.globalTrack().isNonnull() ){
        b_nglobalhits = mu.globalTrack()->hitPattern().numberOfValidMuonHits();
    }
    if ( mu.muonBestTrack().isNonnull() ){
        b_trackdxy = fabs(mu.muonBestTrack()->dxy(pv0.position()));
        b_trackdz = fabs(mu.muonBestTrack()->dz(pv0.position()));
    }
    if ( mu.innerTrack().isNonnull() ){
        b_ninnerhits = mu.innerTrack()->hitPattern().numberOfValidPixelHits();
        b_trackerlayers = mu.innerTrack()->hitPattern().trackerLayersWithMeasurement();
    }

    if (!fake) {tr_true->Fill();}
    else {tr_fake->Fill();}

  }

}

void MuonDetail::collectGenMuons(const std::vector<reco::GenParticle>& genParticles, std::vector<reco::GenParticle>& genMuons, std::vector<reco::GenParticle>& genMuonsFromZ) const
{
  for (auto& gen : genParticles) {
    if (std::abs(gen.pdgId())!=13) continue;
    genMuons.push_back(gen);
    if (gen.pt()<5) continue;
    if (std::abs(gen.eta())>2.4) continue;
    bool isMotherZ = false;
    for ( size_t i=0, n=gen.numberOfMothers(); i<n; ++i ) {
        if (std::abs(gen.mother(i)->pdgId()) == 23) { isMotherZ = true; break; }
    }
    if (!isMotherZ) continue;
    genMuonsFromZ.push_back(gen);
  }
}

void MuonDetail::treereset()
{
  b_mu_pt = -9; b_mu_eta = -9; b_mu_phi = -9;
  b_global = -9; b_pf = -9;
  b_chi2pos = -9; b_trkKink = -9; b_segcompati = -9;
  b_chi2 = -9; b_nglobalhits = -9; b_nstations = -9;
  b_trackdxy = -9; b_trackdz = -9;
  b_ninnerhits = -9; b_trackerlayers = -9;
}

void MuonDetail::fillHitMap(const reco::Muon& mu, TH1D* h_hitmap) const
{
  int map = mu.stationMask();
  for (int i=7; i>=0; i--){
     if (map-pow(2,i)>=0) {map=map-pow(2,i); h_hitmap->Fill(i);}
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonDetail);


