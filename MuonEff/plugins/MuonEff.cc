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
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"

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

class MuonEff : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit MuonEff(const edm::ParameterSet&);
      ~MuonEff();
      TLorentzVector ToTLorentzVector(const reco::Muon& t) { return TLorentzVector(t.px(), t.py(), t.pz(), t.energy()); }
      TLorentzVector ToTLorentzVector(const reco::GenParticle& t) { return TLorentzVector(t.px(), t.py(), t.pz(), t.energy()); }
      bool isGEMrange(const reco::Muon& m) { return (abs(m.eta()) >1.6 && abs(m.eta()) <2.4); }
      bool isGEMrange(const reco::GenParticle& m) { return (abs(m.eta()) >1.6 && abs(m.eta()) <2.4); }

      void collectGenMuons(const std::vector<reco::GenParticle> &genParticles, std::vector<reco::GenParticle>& genMuons, std::vector<reco::GenParticle>& genMuonsFromZ) const;
      bool tightByStep(const reco::Muon& mu, reco::Vertex pv0) const;
      void fillHitMap(const reco::Muon& mu, TH1D* h_hitmap) const;

   private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

      TH1D* h_nevents;
      TH1D* h_genMuons_pt; TH1D* h_genMuons_eta; TH1D* h_genMuons_phi;
      TH1D* h_tight_pt; TH1D* h_tight_eta; TH1D* h_tight_phi;
      TH1D* h_tightByStep_pt; TH1D* h_tightByStep_eta; TH1D* h_tightByStep_phi;
      TH1D* h_medium_pt; TH1D* h_medium_eta; TH1D* h_medium_phi;
      TH1D* h_loose_pt; TH1D* h_loose_eta; TH1D* h_loose_phi;
      TH1D* h_tightFake_pt; TH1D* h_tightFake_eta; TH1D* h_tightFake_phi;
      TH1D* h_tightByStepFake_pt; TH1D* h_tightByStepFake_eta; TH1D* h_tightByStepFake_phi;
      TH1D* h_mediumFake_pt; TH1D* h_mediumFake_eta; TH1D* h_mediumFake_phi;
      TH1D* h_looseFake_pt; TH1D* h_looseFake_eta; TH1D* h_looseFake_phi;
      TH1D* h_tightdetail;
      TH1D* h_hitmap; TH1D* h_hitmap_inrange;
      TH1D* h_tighttest; TH1D* h_tightBySteptest; TH1D* h_mediumtest; TH1D* h_loosetest;

      TH1D* h_nGemhits_true; TH1D* h_nGemhits_fake;

      edm::EDGetTokenT<std::vector<reco::GenParticle> > genParticleToken_;
      edm::EDGetTokenT<std::vector<reco::Muon> > muonToken_;
      edm::EDGetTokenT<std::vector<reco::Vertex> > vtxToken_;
      edm::EDGetTokenT<std::vector<PSimHit> > MuonGEMHitsToken_;

};

using namespace std;
using namespace reco;
using namespace edm;

MuonEff::MuonEff(const edm::ParameterSet& iConfig)
{
  genParticleToken_ = consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticles"));
  muonToken_ = consumes<std::vector<reco::Muon> >(iConfig.getParameter<edm::InputTag>("muons"));
  vtxToken_ = consumes<std::vector<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("vertexs"));
  MuonGEMHitsToken_ = consumes<std::vector<PSimHit> >(iConfig.getParameter<edm::InputTag>("MuonGEMHits"));

  usesResource("TFileService");
  edm::Service<TFileService> fs;

  h_nevents  = fs->make<TH1D>("nevents", "nevents", 1, 0, 1);

  h_genMuons_pt  = fs->make<TH1D>("gen_pt", "gen_pt", 200, 5, 205);
  h_genMuons_eta = fs->make<TH1D>("gen_eta", "gen_eta", 56, 0, 2.8);
  h_genMuons_phi = fs->make<TH1D>("gen_phi", "gen_phi", 90, -3, 3);

  h_tight_pt   = fs->make<TH1D>("genT_pt", "genMatched tight pt", 200, 5, 205);
  h_tight_eta  = fs->make<TH1D>("genT_eta", "genMatched tight eta", 56, 0, 2.8);
  h_tight_phi  = fs->make<TH1D>("genT_phi", "genMatched tight phi", 90, -3, 3);
  h_tightByStep_pt   = fs->make<TH1D>("genTs_pt", "genMatched tightByStep pt", 200, 5, 205);
  h_tightByStep_eta  = fs->make<TH1D>("genTs_eta", "genMatched tightByStep eta", 56, 0, 2.8);
  h_tightByStep_phi  = fs->make<TH1D>("genTs_phi", "genMatched tightByStep phi", 90, -3, 3);
  h_medium_pt  = fs->make<TH1D>("genM_pt", "genMatched medium pt", 200, 5, 205);
  h_medium_eta = fs->make<TH1D>("genM_eta", "genMatched medium eta", 56, 0, 2.8);
  h_medium_phi = fs->make<TH1D>("genM_phi", "genMatched medium phi", 90, -3, 3);
  h_loose_pt   = fs->make<TH1D>("genL_pt", "genMatched loose pt", 200, 5, 205);
  h_loose_eta  = fs->make<TH1D>("genL_eta", "genMatched loose eta", 56, 0, 2.8);
  h_loose_phi  = fs->make<TH1D>("genL_phi", "genMatched loose phi", 90, -3, 3);

  h_tightFake_pt   = fs->make<TH1D>("fakeT_pt", "reco tightfake pt", 200, 0, 200);
  h_tightFake_eta  = fs->make<TH1D>("fakeT_eta", "reco tightfake eta", 56, 0, 2.8);
  h_tightFake_phi  = fs->make<TH1D>("fakeT_phi", "reco tightfake phi", 90, -3, 3);
  h_tightByStepFake_pt   = fs->make<TH1D>("fakeTs_pt", "reco tightByStepfake pt", 200, 0, 200);
  h_tightByStepFake_eta  = fs->make<TH1D>("fakeTs_eta", "reco tightByStepfake eta", 56, 0, 2.8);
  h_tightByStepFake_phi  = fs->make<TH1D>("fakeTs_phi", "reco tightByStepfake phi", 90, -3, 3);
  h_mediumFake_pt  = fs->make<TH1D>("fakeM_pt", "reco mediumfake pt", 200, 0, 200);
  h_mediumFake_eta = fs->make<TH1D>("fakeM_eta", "reco mediumfake eta", 56, 0, 2.8);
  h_mediumFake_phi = fs->make<TH1D>("fakeM_phi", "reco mediumfake phi", 90, -3, 3);
  h_looseFake_pt   = fs->make<TH1D>("fakeL_pt", "reco loosefake pt", 200, 0, 200);
  h_looseFake_eta  = fs->make<TH1D>("fakeL_eta", "reco loosefake eta", 56, 0, 2.8);
  h_looseFake_phi  = fs->make<TH1D>("fakeL_phi", "reco loosefake phi", 90, -3, 3);

  h_tightdetail  = fs->make<TH1D>("tightDetail", "tightDetail", 8, 0, 8);
  h_hitmap  = fs->make<TH1D>("hitmap", "hitmap", 8, 0, 8);
  h_hitmap_inrange  = fs->make<TH1D>("hitmap_inrange", "hitmap inrange", 8, 0, 8);

  h_tighttest  = fs->make<TH1D>("genMatched tight test", "genMatched tight test", 15, -3, 3);
  h_tightBySteptest  = fs->make<TH1D>("genMatched tightByStepbystep test", "genMatched tightByStep test", 15, -3, 3);
  h_mediumtest  = fs->make<TH1D>("genMatched medium test", "genMatched medium test", 15, -3, 3);
  h_loosetest  = fs->make<TH1D>("genMatched loose test", "genMatched loose test", 15, -3, 3);

  h_nGemhits_true = fs->make<TH1D>("nGemhits_true", "nGemhits_true", 30,0,30);
  h_nGemhits_fake = fs->make<TH1D>("nGemhits_fake", "nGemhits_fake", 30,0,30);
}

MuonEff::~MuonEff(){}

void MuonEff::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  h_nevents->Fill(0);

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

  edm::Handle<std::vector<PSimHit>> MuonGEMHits;
  iEvent.getByToken(MuonGEMHitsToken_, MuonGEMHits);

  ///// Gen-Reco matching
  edm::Handle<std::vector<reco::Muon> > muons;
  iEvent.getByToken(muonToken_, muons);
  for (auto& gen : genMuonsFromZ){
    h_genMuons_pt->Fill(gen.pt());
    h_genMuons_eta->Fill(abs(gen.eta()));
    h_genMuons_phi->Fill(gen.phi());

    double dR = 0.1;
    const reco::Muon* muMatched = nullptr;
    for (auto& mu : *muons) {
      double dR_tmp = reco::deltaR(gen, mu);
      if (dR < dR_tmp) continue;
      dR = dR_tmp;
      muMatched = &mu;
    }
    if (muMatched == nullptr) continue;
    auto& mu = *muMatched;
      
    if ( muon::isTightMuon(mu, pv0) ){
      h_tight_pt->Fill(gen.pt());
      h_tight_eta->Fill(abs(gen.eta()));
      h_tight_phi->Fill(gen.phi());
    }
    if ( tightByStep(mu,pv0) ){
      h_tightByStep_pt->Fill(gen.pt());
      h_tightByStep_eta->Fill(abs(gen.eta()));
      h_tightByStep_phi->Fill(gen.phi());
    }
    if ( muon::isMediumMuon(mu) ){
      h_medium_pt->Fill(gen.pt());
      h_medium_eta->Fill(abs(gen.eta()));
      h_medium_phi->Fill(gen.phi());
    }
    if ( muon::isLooseMuon(mu) ){
      h_loose_pt->Fill(gen.pt());
      h_loose_eta->Fill(abs(gen.eta()));
      h_loose_phi->Fill(gen.phi());
    }

    // hitmap
    fillHitMap(mu, h_hitmap);
    if (gen.eta()>0.2 && gen.eta()<0.4 && gen.phi()>1.4 && gen.phi()<1.8 && muon::isMediumMuon(mu) && !(muon::isTightMuon(mu, pv0)) ){
      h_tightdetail->Fill(0);
      fillHitMap(mu, h_hitmap_inrange);
      if ( (mu.isGlobalMuon() && mu.globalTrack()->normalizedChi2()<10.) ) h_tightdetail->Fill(1);
      if ( (mu.isGlobalMuon() && mu.globalTrack()->hitPattern().numberOfValidMuonHits() > 0) ) h_tightdetail->Fill(2);
      if ( (mu.numberOfMatchedStations() > 1) ) h_tightdetail->Fill(3);
      if ( (fabs(mu.muonBestTrack()->dxy(pv0.position())) < 0.2) ) h_tightdetail->Fill(4);
      if ( (fabs(mu.muonBestTrack()->dz(pv0.position())) < 0.5) ) h_tightdetail->Fill(5);
      if ( (mu.innerTrack()->hitPattern().numberOfValidPixelHits() > 0) ) h_tightdetail->Fill(6);
      if ( (mu.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5) ) h_tightdetail->Fill(7);
    }

  }

  //// fake rate
  double rcut = 0.3;
  for (auto& mu : *muons) {
    bool fake = true;
    for (auto& gen : genMuons){
      double dR = reco::deltaR(gen, mu);
      if (dR > rcut) continue;
      fake = false;
    }
    if (fake){
      if (isGEMrange(mu)) h_nGemhits_fake->Fill((*MuonGEMHits).size());
      if ( muon::isTightMuon(mu, pv0) ){
        h_tightFake_pt->Fill(mu.pt());
        h_tightFake_eta->Fill(abs(mu.eta()));
        h_tightFake_phi->Fill(mu.phi());
      }
      if ( tightByStep(mu, pv0) ){
        h_tightByStepFake_pt->Fill(mu.pt());
        h_tightByStepFake_eta->Fill(abs(mu.eta()));
        h_tightByStepFake_phi->Fill(mu.phi());
      }
      if ( muon::isMediumMuon(mu) ){
        h_mediumFake_pt->Fill(mu.pt());
        h_mediumFake_eta->Fill(abs(mu.eta()));
        h_mediumFake_phi->Fill(mu.phi());
      }
      if ( muon::isLooseMuon(mu) ){
        h_looseFake_pt->Fill(mu.pt());
        h_looseFake_eta->Fill(abs(mu.eta()));
        h_looseFake_phi->Fill(mu.phi());
      }
    }
    else {
      if (isGEMrange(mu)) h_nGemhits_true->Fill((*MuonGEMHits).size());
    }
  }

}

void MuonEff::collectGenMuons(const std::vector<reco::GenParticle>& genParticles, std::vector<reco::GenParticle>& genMuons, std::vector<reco::GenParticle>& genMuonsFromZ) const
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

bool MuonEff::tightByStep(const reco::Muon& mu, reco::Vertex pv0) const
{
  if ( !(mu.isGlobalMuon()) ) return false;
  if ( !(mu.isPFMuon()) ) return false;
  if ( !(mu.globalTrack().isNonnull()) )  return false;
  if ( !(mu.muonBestTrack().isNonnull()) ) return false;
  if ( !(mu.innerTrack().isNonnull()) ) return false;
  if ( !(mu.globalTrack()->normalizedChi2()<10.) ) return false;
  if ( !(mu.globalTrack()->hitPattern().numberOfValidMuonHits() > 0) ) return false;
  if ( !(mu.numberOfMatchedStations() > 1) ) return false;
  if ( !(fabs(mu.muonBestTrack()->dxy(pv0.position())) < 0.2) ) return false;
  //if ( !(fabs(mu.muonBestTrack()->dz(pv0.position())) < 0.5) ) return false;
  if ( !(mu.innerTrack()->hitPattern().numberOfValidPixelHits() > 0) ) return false;
  if ( !(mu.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5) ) return false;
  return true;
}

void MuonEff::fillHitMap(const reco::Muon& mu, TH1D* h_hitmap) const
{
  int map = mu.stationMask();
  for (int i=7; i>=0; i--){
     if (map-pow(2,i)>=0) {map=map-pow(2,i); h_hitmap->Fill(i);}
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonEff);


