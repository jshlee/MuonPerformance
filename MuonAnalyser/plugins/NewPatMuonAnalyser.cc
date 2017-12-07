#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "PhysicsTools/PatUtils/interface/MiniIsolation.h"
#include "DataFormats/MuonReco/interface/MuonSimInfo.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/GEMGeometry/interface/ME0Geometry.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

#include <vector>
#include <TTree.h>
#include <TFile.h>
#include <TLorentzVector.h>

using namespace std;

class NewPatMuonAnalyser : public edm::one::EDAnalyzer<edm::one::SharedResources>  {

public:
  explicit NewPatMuonAnalyser(const edm::ParameterSet&);
  ~NewPatMuonAnalyser(){};

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  void setMvaValues(string tmvaWeight_);
  void getMvaValues(const pat::Muon& mu, reco::Vertex pv0);
  bool isTightMuonCustom(const reco::Muon& mu, reco::Vertex pv0) const;

  float ptRel(const reco::Candidate::LorentzVector& muP4, const reco::Candidate::LorentzVector& jetP4, bool subtractMuon);
  void initTreeValues();

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

  void setBranches(TTree *tree);
  void fillBranches(TTree *tree, TLorentzVector &tlv, edm::RefToBase<pat::Muon> muref);

  bool isFromZ(const reco::GenParticle &gen);

  edm::EDGetTokenT<std::vector<reco::Vertex>> verticesToken_;
  edm::EDGetTokenT<edm::View<pat::Muon>> muonsToken_;

  edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo>> putoken_;

  edm::EDGetTokenT<edm::ValueMap<float> > puppiNewIso_ch_;
  edm::EDGetTokenT<edm::ValueMap<float> > puppiNewIso_nh_;
  edm::EDGetTokenT<edm::ValueMap<float> > puppiNewIso_ph_;
  edm::Handle<edm::ValueMap<float>> puppiNewIso_ch;
  edm::Handle<edm::ValueMap<float>> puppiNewIso_nh;
  edm::Handle<edm::ValueMap<float>> puppiNewIso_ph;
  edm::EDGetTokenT<edm::ValueMap<float> > trkNewIso_;
  edm::Handle<edm::ValueMap<float>> trkNewIso;
  edm::EDGetTokenT<edm::ValueMap<float> > pfNewIso_ch_;
  edm::EDGetTokenT<edm::ValueMap<float> > pfNewIso_nh_;
  edm::EDGetTokenT<edm::ValueMap<float> > pfNewIso_ph_;
  edm::Handle<edm::ValueMap<float>> pfNewIso_ch;
  edm::Handle<edm::ValueMap<float>> pfNewIso_nh;
  edm::Handle<edm::ValueMap<float>> pfNewIso_ph;

  reco::Vertex priVertex_;
  
  const ME0Geometry* ME0Geometry_;
  
  TTree* genttree_;
  TTree* recottree_;
  TH1D* h_nevents;

  int b_pu_density, b_pu_numInteractions;
  int b_nvertex;
  
  TLorentzVector b_muon;
  int b_muon_pdgId;
  int b_muon_mpdgId;
  int b_muon_no;
  int b_muon_simType;
  bool b_muon_isSignalMuon;

  float b_muon_embeddedMva;
  float b_muon_bdtMva;

  float b_muon_istracker, b_muon_isglobal, b_muon_ispf;
  float b_muon_chi2, b_muon_chi2pos, b_muon_trkKink, b_muon_segmentCompatibility, b_muon_nstations, b_muon_nglobalhits;
  float b_muon_trackdxy, b_muon_trackdz, b_muon_ninnerhits, b_muon_trackerlayers, b_muon_innerquality, b_muon_caloCompatibility;

  TMVA::Reader* bdt_;

  float b_muon_miniIso_ch;
  float b_muon_miniIso_nh;
  float b_muon_miniIso_ph;
  float b_muon_miniIso_pu;

  float b_muon_poszPV0, b_muon_poszMuon;
  bool b_muon_isTight, b_muon_isMedium, b_muon_isLoose;
  bool b_muon_isTightCustom;

  float b_muon_PFIso03, b_muon_PFIso04;
  float b_muon_PFIso03ChargedHadronPt, b_muon_PFIso03NeutralHadronEt, b_muon_PFIso03PhotonEt, b_muon_PFIso03PUPt;
  float b_muon_PFIso04ChargedHadronPt, b_muon_PFIso04NeutralHadronEt, b_muon_PFIso04PhotonEt, b_muon_PFIso04PUPt;
  float b_muon_TrkIso03, b_muon_TrkIso05;
  float b_muon_puppiIso, b_muon_puppiIsoNoLep;
  float b_muon_puppiIso_ChargedHadron, b_muon_puppiIso_NeutralHadron, b_muon_puppiIso_Photon;  
  float b_muon_puppiIsoNoLep_ChargedHadron, b_muon_puppiIsoNoLep_NeutralHadron, b_muon_puppiIsoNoLep_Photon;  

  float b_muon_puppiNewIso_ch, b_muon_puppiNewIso_nh, b_muon_puppiNewIso_ph, b_muon_puppiNewIso;
  float b_muon_trkNewIso;
  float b_muon_pfNewIso_ch, b_muon_pfNewIso_nh, b_muon_pfNewIso_ph, b_muon_pfNewIso;

  string tmvaWeight_;
};

NewPatMuonAnalyser::NewPatMuonAnalyser(const edm::ParameterSet& iConfig):
  verticesToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
  muonsToken_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
  prunedGenToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned")))
{
  putoken_ = consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("addPileupInfo"));
  tmvaWeight_ = edm::FileInPath(iConfig.getParameter<std::string>("tmvaWeightLabel")).fullPath();
  
  puppiNewIso_ch_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("puppiNewIso_ch"));
  puppiNewIso_nh_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("puppiNewIso_nh"));
  puppiNewIso_ph_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("puppiNewIso_ph"));
  trkNewIso_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("trkNewIso"));
  pfNewIso_ch_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("pfNewIso_ch"));
  pfNewIso_nh_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("pfNewIso_nh"));
  pfNewIso_ph_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("pfNewIso_ph"));
  
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  h_nevents = fs->make<TH1D>("nevents", "nevents", 1, 0, 1);
  
  genttree_ = fs->make<TTree>("gen", "gen");
  setBranches(genttree_);
  recottree_ = fs->make<TTree>("reco", "reco");
  setBranches(recottree_);
}

void NewPatMuonAnalyser::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  h_nevents->Fill(0.5);
  
  edm::ESHandle<ME0Geometry> hGeom;
  iSetup.get<MuonGeometryRecord>().get(hGeom);
  ME0Geometry_ =( &*hGeom);
  
  using namespace edm;
  Handle<std::vector<reco::Vertex>> vertices;
  iEvent.getByToken(verticesToken_, vertices);
  // Vertices
  int prVtx = -1;
  for (size_t i = 0; i < vertices->size(); i++) {
    if (vertices->at(i).isFake()) continue;
    if (vertices->at(i).ndof() <= 4) continue;
    if (prVtx < 0) prVtx = i;
  }
  if (prVtx < 0) return;
  priVertex_ = vertices->at(prVtx);
  b_nvertex = vertices->size();

  Handle<View<pat::Muon> > muons;
  iEvent.getByToken(muonsToken_, muons);

  Handle<edm::View<reco::GenParticle> > pruned;
  iEvent.getByToken(prunedGenToken_,pruned);

  iEvent.getByToken(puppiNewIso_ch_, puppiNewIso_ch);
  iEvent.getByToken(puppiNewIso_nh_, puppiNewIso_nh);
  iEvent.getByToken(puppiNewIso_ph_, puppiNewIso_ph);  
  iEvent.getByToken(trkNewIso_, trkNewIso);
  iEvent.getByToken(pfNewIso_ch_, pfNewIso_ch);
  iEvent.getByToken(pfNewIso_nh_, pfNewIso_nh);
  iEvent.getByToken(pfNewIso_ph_, pfNewIso_ph);  

  edm::Handle<std::vector <PileupSummaryInfo> > PupInfo;
  iEvent.getByToken(putoken_, PupInfo);
  b_pu_density = 0; b_pu_numInteractions = 0;
  std::vector<PileupSummaryInfo>::const_iterator ipu;
  for (ipu = PupInfo->begin(); ipu != PupInfo->end(); ++ipu) {
    if ( ipu->getBunchCrossing() != 0 ) continue; // storing detailed PU info only for BX=0
    for (unsigned int i=0; i<ipu->getPU_zpositions().size(); ++i) {
      if ( abs((ipu->getPU_zpositions())[i] - priVertex_.position().z()) < 0.1 )
      ++b_pu_density;
      ++b_pu_numInteractions;
    }
  }

  // gen muon loop
  b_muon_no = 0;
  for (const reco::GenParticle &gen : *pruned) {
    if (abs(gen.pdgId()) != 13) continue;
    //if ( !isFromZ(gen) ) continue;

    TLorentzVector gentlv(gen.momentum().x(), gen.momentum().y(), gen.momentum().z(), gen.energy() );

    // getting matched reco muon
    edm::RefToBase<pat::Muon> recomuRef;
    for (size_t i = 0; i < muons->size(); i++) {
      auto muon = muons->at(i);

      //matching GEN muon & SIM muon using deltaR
      TLorentzVector simtlv;
      simtlv.SetPtEtaPhiM(muon.simPt(), muon.simEta(), muon.simPhi(), 0.10566);
      if (gentlv.DeltaR(simtlv) < 0.1){
            recomuRef = muons->refAt(i);
            break;
      }
    }

    //if ((recomuRef.get())->simBX() != 0) { continue; }
    fillBranches(genttree_, gentlv, recomuRef);
  }

  // reco muons loop
  b_muon_no = 0;
  for (size_t i = 0; i < muons->size(); i++) {
    edm::RefToBase<pat::Muon> recomuRef = muons->refAt(i);
    auto muon = muons->at(i);

    TLorentzVector recotlv(muon.momentum().x(), muon.momentum().y(), muon.momentum().z(), muon.energy() );
    fillBranches(recottree_, recotlv, recomuRef);
  }

  return;
}

void NewPatMuonAnalyser::fillBranches(TTree *tree, TLorentzVector &tlv, edm::RefToBase<pat::Muon> muref)
{
  initTreeValues();
  b_muon = tlv;
  ++b_muon_no;

  if (muref.isNonnull()){
    auto muon = muref.get();
  
    b_muon_poszPV0  = priVertex_.position().z();
    b_muon_poszMuon = muon->vz();
    
    // pdg id
    b_muon_pdgId = muon->pdgId();
    b_muon_mpdgId = muon->simMotherPdgId();

    // muon id
    b_muon_isTightCustom = isTightMuonCustom(*muon, priVertex_);
    b_muon_isTight = muon::isTightMuon(*muon, priVertex_);
    b_muon_isMedium = muon::isMediumMuon(*muon);
    b_muon_isLoose = muon::isLooseMuon(*muon);
    
    // matched with any gen?
    b_muon_isSignalMuon = (muon->simType() == reco::MuonSimType::MatchedPrimaryMuon);
    b_muon_simType = muon->simType();

    // mva values
    b_muon_embeddedMva = muon->mvaValue();
    b_muon_bdtMva = muon->mvaValue();

    // miniIsolation
    b_muon_miniIso_ch = muon->miniPFIsolation().chargedHadronIso();
    b_muon_miniIso_nh = muon->miniPFIsolation().neutralHadronIso();
    b_muon_miniIso_ph = muon->miniPFIsolation().photonIso();
    b_muon_miniIso_pu = muon->miniPFIsolation().puChargedHadronIso(); 
    
    // old isolations
    b_muon_TrkIso03 = muon->isolationR03().sumPt/muon->pt();
    b_muon_TrkIso05 = muon->isolationR05().sumPt/muon->pt();
    
    b_muon_PFIso03ChargedHadronPt = muon->pfIsolationR03().sumChargedHadronPt;
    b_muon_PFIso03NeutralHadronEt = muon->pfIsolationR03().sumNeutralHadronEt;
    b_muon_PFIso03PhotonEt        = muon->pfIsolationR03().sumPhotonEt;
    b_muon_PFIso03PUPt            = muon->pfIsolationR03().sumPUPt;

    b_muon_PFIso04ChargedHadronPt = muon->pfIsolationR04().sumChargedHadronPt;
    b_muon_PFIso04NeutralHadronEt = muon->pfIsolationR04().sumNeutralHadronEt;
    b_muon_PFIso04PhotonEt        = muon->pfIsolationR04().sumPhotonEt;
    b_muon_PFIso04PUPt            = muon->pfIsolationR04().sumPUPt;   
    
    b_muon_PFIso04 = (muon->pfIsolationR04().sumChargedHadronPt + TMath::Max(0.,muon->pfIsolationR04().sumNeutralHadronEt + muon->pfIsolationR04().sumPhotonEt - 0.5*muon->pfIsolationR04().sumPUPt))/muon->pt();
    b_muon_PFIso03 = (muon->pfIsolationR03().sumChargedHadronPt + TMath::Max(0.,muon->pfIsolationR03().sumNeutralHadronEt + muon->pfIsolationR03().sumPhotonEt - 0.5*muon->pfIsolationR03().sumPUPt))/muon->pt();
    
    b_muon_puppiIso_ChargedHadron = muon->puppiChargedHadronIso();
    b_muon_puppiIso_NeutralHadron = muon->puppiNeutralHadronIso();
    b_muon_puppiIso_Photon = muon->puppiPhotonIso();
    b_muon_puppiIso = (b_muon_puppiIso_ChargedHadron+b_muon_puppiIso_NeutralHadron+b_muon_puppiIso_Photon)/muon->pt();
    b_muon_puppiIsoNoLep_ChargedHadron = muon->puppiNoLeptonsChargedHadronIso();
    b_muon_puppiIsoNoLep_NeutralHadron = muon->puppiNoLeptonsNeutralHadronIso();
    b_muon_puppiIsoNoLep_Photon = muon->puppiNoLeptonsPhotonIso();
    b_muon_puppiIsoNoLep = (b_muon_puppiIsoNoLep_ChargedHadron+b_muon_puppiIsoNoLep_NeutralHadron+b_muon_puppiIsoNoLep_Photon)/muon->pt(); 

    b_muon_puppiNewIso_ch = (*puppiNewIso_ch)[muref];
    b_muon_puppiNewIso_nh = (*puppiNewIso_nh)[muref];
    b_muon_puppiNewIso_ph = (*puppiNewIso_ph)[muref];
    b_muon_puppiNewIso    = ( b_muon_puppiNewIso_ch + b_muon_puppiNewIso_nh + b_muon_puppiNewIso_ph )/muon->pt();
    b_muon_trkNewIso = (*trkNewIso)[muref] / muon->pt();
    b_muon_pfNewIso_ch = (*pfNewIso_ch)[muref];
    b_muon_pfNewIso_nh = (*pfNewIso_nh)[muref];
    b_muon_pfNewIso_ph = (*pfNewIso_ph)[muref];
    b_muon_pfNewIso    = ( b_muon_pfNewIso_ch + max(0.0, b_muon_pfNewIso_nh + b_muon_pfNewIso_ph - 0.5 ) )/ muon->pt();
    
  }
  tree->Fill();
}

void NewPatMuonAnalyser::setMvaValues(string tmvaWeight_)
{
  bdt_ = new TMVA::Reader();
  bdt_->AddVariable("muon_istracker",&b_muon_istracker);
  bdt_->AddVariable("muon_isglobal",&b_muon_isglobal);
  bdt_->AddVariable("muon_ispf",&b_muon_ispf);
  bdt_->AddVariable("muon_normalizedChi2",&b_muon_chi2);
  bdt_->AddVariable("muon_chi2LocalPosition",&b_muon_chi2pos);
  bdt_->AddVariable("muon_trkKink",&b_muon_trkKink);
  bdt_->AddVariable("muon_segmentCompatibility",&b_muon_segmentCompatibility);
  bdt_->AddVariable("muon_numberOfMatchedStations",&b_muon_nstations);
  bdt_->AddVariable("muon_numberOfValidMuonHits",&b_muon_nglobalhits);
  bdt_->AddVariable("muon_pv0pos_dxy",&b_muon_trackdxy);
  bdt_->AddVariable("muon_numberOfValidPixelHits",&b_muon_ninnerhits);
  bdt_->AddVariable("muon_trackerLayersWithMeasurement",&b_muon_trackerlayers);
  bdt_->AddVariable("muon_innerquality",&b_muon_innerquality);
  bdt_->AddVariable("muon_caloCompatibility",&b_muon_caloCompatibility);
  bdt_->BookMVA("BDT", tmvaWeight_);
}

void NewPatMuonAnalyser::getMvaValues(const pat::Muon& mu, reco::Vertex pv0)
{
  float dummyVal = -999;

  b_muon_istracker = mu.isTrackerMuon();
  b_muon_isglobal = mu.isGlobalMuon();
  b_muon_ispf = mu.isPFMuon();
  if ( mu.globalTrack().isNonnull() ){ b_muon_chi2 = mu.globalTrack()->normalizedChi2(); }
  else { b_muon_chi2 = dummyVal; }
  b_muon_chi2pos = mu.combinedQuality().chi2LocalPosition;
  b_muon_trkKink = mu.combinedQuality().trkKink;
  b_muon_segmentCompatibility = muon::segmentCompatibility(mu);
  b_muon_nstations = mu.numberOfMatchedStations();
  if ( mu.globalTrack().isNonnull() ){ b_muon_nglobalhits = mu.globalTrack()->hitPattern().numberOfValidMuonHits(); }
  else { b_muon_nglobalhits = dummyVal; }
  if ( mu.muonBestTrack().isNonnull() ){ b_muon_trackdxy = abs(mu.muonBestTrack()->dxy(pv0.position())); }
  else { b_muon_trackdxy = dummyVal; }
  if ( mu.innerTrack().isNonnull() ){
      b_muon_ninnerhits = mu.innerTrack()->hitPattern().numberOfValidPixelHits();
      b_muon_trackerlayers = mu.innerTrack()->hitPattern().trackerLayersWithMeasurement();
      b_muon_innerquality = mu.innerTrack()->quality(reco::Track::highPurity);
  }
  else {
      b_muon_ninnerhits = dummyVal;
      b_muon_trackerlayers = dummyVal;
      b_muon_innerquality = dummyVal;
  }
  b_muon_caloCompatibility = mu.caloCompatibility();
}

bool NewPatMuonAnalyser::isFromZ(const reco::GenParticle &gen)
{
  for (unsigned int i = 0; i < gen.numberOfMothers(); ++i){
    if (gen.mother(i)->pdgId() == 23 || gen.mother(i)->pdgId() == 25){
      return true;
    }
  }
  return false;
}

bool NewPatMuonAnalyser::isTightMuonCustom(const reco::Muon& muon, reco::Vertex vtx) const
{
  if ( !muon.isPFMuon() || !muon.isGlobalMuon() ) return false;
  
  bool muID = muon::isGoodMuon(muon,muon::GlobalMuonPromptTight) && (muon.numberOfMatchedStations() > 1);
      
  bool hits = muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 &&
    muon.innerTrack()->hitPattern().numberOfValidPixelHits() > 0; 
  
  bool ip = fabs(muon.muonBestTrack()->dxy(vtx.position())) < 0.2;
  //&& fabs(muon.muonBestTrack()->dz(vtx.position())) < 0.5;
  
  return muID && hits && ip;
}

void NewPatMuonAnalyser::setBranches(TTree *tree)
{
  tree->Branch("nvertex", &b_nvertex, "nvertex/I");
  tree->Branch("pu_density", &b_pu_density, "pu_density/I");
  tree->Branch("pu_numInteractions", &b_pu_numInteractions, "pu_numInteractions/I");
  tree->Branch("muon", "TLorentzVector", &b_muon);  
  tree->Branch("muon_no", &b_muon_no, "muon_no/I");
  tree->Branch("muon_pdgId", &b_muon_pdgId, "muon_pdgId/I");
  tree->Branch("muon_mpdgId", &b_muon_mpdgId, "muon_mpdgId/I");

  tree->Branch("muon_poszPV0",&b_muon_poszPV0,"muon_poszPV0/F");
  tree->Branch("muon_poszMuon",&b_muon_poszMuon,"muon_poszMuon/F");
  tree->Branch("muon_isSignalMuon", &b_muon_isSignalMuon, "muon_isSignalMuon/O");
  tree->Branch("muon_simType", &b_muon_simType, "muon_simType/I");

  tree->Branch("muon_miniIso_ch", &b_muon_miniIso_ch, "muon_miniIso_ch/F");
  tree->Branch("muon_miniIso_nh", &b_muon_miniIso_nh, "muon_miniIso_nh/F");
  tree->Branch("muon_miniIso_ph", &b_muon_miniIso_ph, "muon_miniIso_ph/F");
  tree->Branch("muon_miniIso_pu", &b_muon_miniIso_pu, "muon_miniIso_pu/F");

  tree->Branch("muon_embeddedMva", &b_muon_embeddedMva, "muon_embeddedMva/F");
  tree->Branch("muon_bdtMva", &b_muon_bdtMva, "muon_bdtMva/F");

  tree->Branch("muon_isTightCustom", &b_muon_isTightCustom, "muon_isTightCustom/O");
  tree->Branch("muon_isTight", &b_muon_isTight, "muon_isTight/O");
  tree->Branch("muon_isMedium", &b_muon_isMedium, "muon_isMedium/O");
  tree->Branch("muon_isLoose", &b_muon_isLoose, "muon_isLoose/O");

  tree->Branch("muon_TrkIsolation03",&b_muon_TrkIso03,"muon_TrkIsolation03/F");
  tree->Branch("muon_TrkIsolation05",&b_muon_TrkIso05,"muon_TrkIsolation05/F");
  tree->Branch("muon_PFIsolation04",&b_muon_PFIso04,"muon_PFIsolation04/F");
  tree->Branch("muon_PFIsolation03",&b_muon_PFIso03,"muon_PFIsolation03/F");
  tree->Branch("muon_PFIso03ChargedHadronPt",&b_muon_PFIso03ChargedHadronPt,"muon_PFIso03ChargedHadronPt/F");
  tree->Branch("muon_PFIso03NeutralHadronEt",&b_muon_PFIso03NeutralHadronEt,"muon_PFIso03NeutralHadronEt/F");
  tree->Branch("muon_PFIso03PhotonEt",&b_muon_PFIso03PhotonEt,"muon_PFIso03PhotonEt/F");
  tree->Branch("muon_PFIso03PUPt",&b_muon_PFIso03PUPt,"muon_PFIso03PUPt/F");
  tree->Branch("muon_PFIso04ChargedHadronPt",&b_muon_PFIso04ChargedHadronPt,"muon_PFIso04ChargedHadronPt/F");
  tree->Branch("muon_PFIso04NeutralHadronEt",&b_muon_PFIso04NeutralHadronEt,"muon_PFIso04NeutralHadronEt/F");
  tree->Branch("muon_PFIso04PhotonEt",&b_muon_PFIso04PhotonEt,"muon_PFIso04PhotonEt/F");
  tree->Branch("muon_PFIso04PUPt",&b_muon_PFIso04PUPt,"muon_PFIso04PUPt/F");
  tree->Branch("muon_puppiIso",&b_muon_puppiIso,"muon_puppiIso/F");
  tree->Branch("muon_puppiIso_ChargedHadron",&b_muon_puppiIso_ChargedHadron,"muon_puppiIso_ChargedHadron/F");
  tree->Branch("muon_puppiIso_NeutralHadron",&b_muon_puppiIso_NeutralHadron,"muon_puppiIso_NeutralHadron/F");
  tree->Branch("muon_puppiIso_Photon",&b_muon_puppiIso_Photon,"muon_puppiIso_Photon/F");
  tree->Branch("muon_puppiIsoNoLep",&b_muon_puppiIsoNoLep,"muon_puppiIsoNoLep/F");
  tree->Branch("muon_puppiIsoNoLep_ChargedHadron",&b_muon_puppiIsoNoLep_ChargedHadron,"muon_puppiIsoNoLep_ChargedHadron/F");
  tree->Branch("muon_puppiIsoNoLep_NeutralHadron",&b_muon_puppiIsoNoLep_NeutralHadron,"muon_puppiIsoNoLep_NeutralHadron/F");
  tree->Branch("muon_puppiIsoNoLep_Photon",&b_muon_puppiIsoNoLep_Photon,"muon_puppiIsoNoLep_Photon/F");
  
  tree->Branch("muon_puppiNewIso_ch",&b_muon_puppiNewIso_ch,"muon_puppiNewIso_ch/F");
  tree->Branch("muon_puppiNewIso_nh",&b_muon_puppiNewIso_nh,"muon_puppiNewIso_nh/F");
  tree->Branch("muon_puppiNewIso_ph",&b_muon_puppiNewIso_ph,"muon_puppiNewIso_ph/F");
  tree->Branch("muon_puppiNewIso",&b_muon_puppiNewIso,"muon_puppiNewIso/F");
  tree->Branch("muon_trkNewIso",&b_muon_trkNewIso,"muon_trkNewIso/F");
  tree->Branch("muon_pfNewIso_ch",&b_muon_pfNewIso_ch,"muon_pfNewIso_ch/F");
  tree->Branch("muon_pfNewIso_nh",&b_muon_pfNewIso_nh,"muon_pfNewIso_nh/F");
  tree->Branch("muon_pfNewIso_ph",&b_muon_pfNewIso_ph,"muon_pfNewIso_ph/F");
  tree->Branch("muon_pfNewIso",&b_muon_pfNewIso,"muon_pfNewIso/F");
}

void NewPatMuonAnalyser::initTreeValues(){
  b_muon.Clear();

  b_muon_no = 0;
  b_muon_pdgId = 0; b_muon_mpdgId = 0;

  b_nvertex = 0; b_pu_density = 0; b_pu_numInteractions = 0;
  b_muon_poszPV0  = 0; b_muon_poszMuon = 0;
    
  b_muon_isSignalMuon = 0; b_muon_simType = -99;

  b_muon_isTight = 0; b_muon_isMedium = 0; b_muon_isLoose = 0;
  b_muon_isTightCustom = 0;
  
  b_muon_embeddedMva = -99;
  b_muon_bdtMva = -99;

  b_muon_miniIso_ch = -999;
  b_muon_miniIso_nh= -999;
  b_muon_miniIso_ph = -999;
  b_muon_miniIso_pu = -999;

  b_muon_PFIso04 = 0;  b_muon_PFIso03 = 0;
  b_muon_PFIso03ChargedHadronPt = 0; b_muon_PFIso03NeutralHadronEt = 0;
  b_muon_PFIso03PhotonEt = 0; b_muon_PFIso03PUPt = 0;
  b_muon_PFIso04ChargedHadronPt = 0; b_muon_PFIso04NeutralHadronEt = 0;
  b_muon_PFIso04PhotonEt = 0; b_muon_PFIso04PUPt = 0;
  b_muon_TrkIso05 = 0;  b_muon_TrkIso03 = 0;
  b_muon_puppiIso = 0; b_muon_puppiIso_ChargedHadron = 0; b_muon_puppiIso_NeutralHadron = 0; b_muon_puppiIso_Photon = 0;
  b_muon_puppiIsoNoLep = 0; b_muon_puppiIsoNoLep_ChargedHadron = 0; b_muon_puppiIsoNoLep_NeutralHadron = 0; b_muon_puppiIsoNoLep_Photon = 0;  

  b_muon_puppiNewIso_ch = -1; b_muon_puppiNewIso_nh = -1; b_muon_puppiNewIso_ph = -1; b_muon_puppiNewIso = -1;
  b_muon_trkNewIso = -1;
  b_muon_pfNewIso_ch = -1; b_muon_pfNewIso_nh = -1; b_muon_pfNewIso_ph = -1; b_muon_pfNewIso = -1;
}

void NewPatMuonAnalyser::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(NewPatMuonAnalyser);
