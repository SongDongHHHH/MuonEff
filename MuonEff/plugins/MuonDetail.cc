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
      bool tightByStep(const reco::Muon& mu, reco::Vertex pv0) const;
      void fillHitMap(const reco::Muon& mu, TH1D* h_hitmap) const;

   private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

      TTree * tr;
      bool b_global; bool b_pf;
      float b_chi2; int b_nglobalhits; int b_nstations;
      float b_trackdxy; float b_trackdz;
      int b_ninnerhits; float b_trackerlayers;

      edm::EDGetTokenT<std::vector<reco::GenParticle> > genParticleToken_;
      edm::EDGetTokenT<std::vector<reco::Muon> > muonToken_;
      edm::EDGetTokenT<std::vector<reco::Vertex> > vtxToken_;
      bool isFake;
};

using namespace std;
using namespace reco;
using namespace edm;

MuonDetail::MuonDetail(const edm::ParameterSet& iConfig)
{
  genParticleToken_ = consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticles"));
  muonToken_ = consumes<std::vector<reco::Muon> >(iConfig.getParameter<edm::InputTag>("muons"));
  vtxToken_ = consumes<std::vector<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("vertexs"));
  isFake = iConfig.getParameter<bool>("isFake");

  usesResource("TFileService");
  edm::Service<TFileService> fs;
  tr = fs->make<TTree>("tree", "tree");

  tr->Branch("isGlobalMuon", &b_global, "isGlobalMuon/O");
  tr->Branch("isPFMuon", &b_pf, "isPFMuon/O");
  tr->Branch("normalizedChi2", &b_chi2, "normalizedChi2/F");
  tr->Branch("numberOfValidMuonHits", &b_nglobalhits, "numberOfValidMuonHits/I");
  tr->Branch("numberOfMatchedStations", &b_nstations, "numberOfMatchedStations/I");
  tr->Branch("pv0pos_dxy", &b_trackdxy, "pv0pos_dxy/F");
  tr->Branch("pv0pos_dz", &b_trackdz, "pv0pos_dz/F");
  tr->Branch("numberOfValidPixelHits", &b_ninnerhits, "numberOfValidPixelHits/I");
  tr->Branch("trackerLayersWithMeasurement", &b_trackerlayers, "trackerLayersWithMeasurement/F");
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

  ///// id
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

    if (fake == isFake){
      b_global = mu.isGlobalMuon();
      b_pf = mu.isPFMuon();
      b_nstations = mu.numberOfMatchedStations();
      if ( mu.globalTrack().isNonnull() ){
          b_chi2 = mu.globalTrack()->normalizedChi2();
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
    }

    tr->Fill();
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

bool MuonDetail::tightByStep(const reco::Muon& mu, reco::Vertex pv0) const
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

void MuonDetail::fillHitMap(const reco::Muon& mu, TH1D* h_hitmap) const
{
  int map = mu.stationMask();
  for (int i=7; i>=0; i--){
     if (map-pow(2,i)>=0) {map=map-pow(2,i); h_hitmap->Fill(i);}
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonDetail);


