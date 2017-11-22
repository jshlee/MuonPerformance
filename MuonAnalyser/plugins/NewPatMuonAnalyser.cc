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

private:
  void initTreeValues();
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

  void setBranches(TTree *tree);
  void fillBranches(TTree *tree, TLorentzVector &tlv, edm::RefToBase<pat::Muon> muref,bool isSignal, int pdgId, pat::PFIsolation miniiso);

  
  edm::EDGetTokenT<std::vector<reco::Vertex>> verticesToken_;
  edm::EDGetTokenT<edm::View<pat::Muon>> muonsToken_;

  edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo>> putoken_;
  edm::EDGetTokenT<std::vector<pat::PackedCandidate>> pfCandsToken_; 
  std::vector<double> miniIsoParams_ ;

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
  //edm::EDGetTokenT<edm::ValueMap<float> > pfNewIso_pu_;
  edm::Handle<edm::ValueMap<float>> pfNewIso_ch;
  edm::Handle<edm::ValueMap<float>> pfNewIso_nh;
  edm::Handle<edm::ValueMap<float>> pfNewIso_ph;
  //edm::Handle<edm::ValueMap<float>> pfNewIso_pu;
  edm::EDGetTokenT<edm::ValueMap<float> > minipuppiNewIso_ch_;
  edm::EDGetTokenT<edm::ValueMap<float> > minipuppiNewIso_nh_;
  edm::EDGetTokenT<edm::ValueMap<float> > minipuppiNewIso_ph_;
  edm::Handle<edm::ValueMap<float>> minipuppiNewIso_ch;
  edm::Handle<edm::ValueMap<float>> minipuppiNewIso_nh;
  edm::Handle<edm::ValueMap<float>> minipuppiNewIso_ph;
  edm::EDGetTokenT<edm::ValueMap<float> > minitrkNewIso_;
  edm::Handle<edm::ValueMap<float>> minitrkNewIso;
  edm::EDGetTokenT<edm::ValueMap<float> > minipfNewIso_ch_;
  edm::EDGetTokenT<edm::ValueMap<float> > minipfNewIso_nh_;
  edm::EDGetTokenT<edm::ValueMap<float> > minipfNewIso_ph_;
  //edm::EDGetTokenT<edm::ValueMap<float> > minipfNewIso_pu_;
  edm::Handle<edm::ValueMap<float>> minipfNewIso_ch;
  edm::Handle<edm::ValueMap<float>> minipfNewIso_nh;
  edm::Handle<edm::ValueMap<float>> minipfNewIso_ph;
  //edm::Handle<edm::ValueMap<float>> minipfNewIso_pu;

  reco::Vertex priVertex_;
  
  const ME0Geometry* ME0Geometry_;
  
  TTree* genttree_;
  TTree* recottree_;
  TH1D* h_nevents;

  int b_pu_density, b_pu_numInteractions;
  int b_nvertex;
  
  TLorentzVector b_muon;
  bool b_muon_signal;
  int b_muon_pdgId;
  int b_muon_no;

  float b_muon_miniIso_ch;
  float b_muon_miniIso_nh;
  float b_muon_miniIso_ph;
  float b_muon_miniIso_pu;

  float b_muon_poszPV0, b_muon_poszMuon;
  bool b_muon_isTight, b_muon_isMedium, b_muon_isLoose;

  float b_muon_PFIso04; float b_muon_PFIso03;
  float b_muon_PFIso03ChargedHadronPt, b_muon_PFIso03NeutralHadronEt;
  float b_muon_PFIso03PhotonEt, b_muon_PFIso03PUPt;
  float b_muon_PFIso04ChargedHadronPt, b_muon_PFIso04NeutralHadronEt;
  float b_muon_PFIso04PhotonEt, b_muon_PFIso04PUPt;
  float b_muon_TrkIso05; float b_muon_TrkIso03;
  float b_muon_puppiIso, b_muon_puppiIsoNoLep;
  float b_muon_puppiIso_ChargedHadron, b_muon_puppiIso_NeutralHadron, b_muon_puppiIso_Photon;  
  float b_muon_puppiIsoNoLep_ChargedHadron, b_muon_puppiIsoNoLep_NeutralHadron, b_muon_puppiIsoNoLep_Photon;  

  float b_muon_puppiNewIso_ch, b_muon_puppiNewIso_nh, b_muon_puppiNewIso_ph, b_muon_puppiNewIso;
  float b_muon_trkNewIso;
  float b_muon_pfNewIso_ch, b_muon_pfNewIso_nh, b_muon_pfNewIso_ph, /*b_muon_pfNewIso_pu, */b_muon_pfNewIso;
  float b_muon_minipuppiNewIso_ch, b_muon_minipuppiNewIso_nh, b_muon_minipuppiNewIso_ph, b_muon_minipuppiNewIso_pu, b_muon_minipuppiNewIso;
  float b_muon_minitrkNewIso;
  float b_muon_minipfNewIso_ch, b_muon_minipfNewIso_nh, b_muon_minipfNewIso_ph, /*b_muon_minipfNewIso_pu, */b_muon_minipfNewIso;

  float b_muon_PFIsoFixOnlyCH;
  float b_muon_puppiIsoFixOnlyCH;
  
  float b_muon_PFIsoRepTrk;
  float b_muon_puppiIsoRepTrk;

};
NewPatMuonAnalyser::NewPatMuonAnalyser(const edm::ParameterSet& iConfig):
  verticesToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
  muonsToken_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
  prunedGenToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned")))
{
  pfCandsToken_ = consumes<std::vector<pat::PackedCandidate> >(iConfig.getParameter<edm::InputTag>("pfCands"));
  miniIsoParams_ = iConfig.getParameter<std::vector<double> >("miniIsoParams");

  putoken_ = consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("addPileupInfo"));
  
  puppiNewIso_ch_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("puppiNewIso_ch"));
  puppiNewIso_nh_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("puppiNewIso_nh"));
  puppiNewIso_ph_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("puppiNewIso_ph"));
  trkNewIso_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("trkNewIso"));
  pfNewIso_ch_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("pfNewIso_ch"));
  pfNewIso_nh_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("pfNewIso_nh"));
  pfNewIso_ph_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("pfNewIso_ph"));
  //pfNewIso_pu_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("pfNewIso_pu"));
  minipuppiNewIso_ch_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("minipuppiNewIso_ch"));
  minipuppiNewIso_nh_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("minipuppiNewIso_nh"));
  minipuppiNewIso_ph_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("minipuppiNewIso_ph"));
  minitrkNewIso_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("minitrkNewIso"));
  minipfNewIso_ch_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("minipfNewIso_ch"));
  minipfNewIso_nh_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("minipfNewIso_nh"));
  minipfNewIso_ph_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("minipfNewIso_ph"));
  //minipfNewIso_pu_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("minipfNewIso_pu"));

  
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
  //iEvent.getByToken(pfNewIso_pu_, pfNewIso_pu);  
  iEvent.getByToken(minipuppiNewIso_ch_, minipuppiNewIso_ch);
  iEvent.getByToken(minipuppiNewIso_nh_, minipuppiNewIso_nh);
  iEvent.getByToken(minipuppiNewIso_ph_, minipuppiNewIso_ph);  
  iEvent.getByToken(minitrkNewIso_, minitrkNewIso);
  iEvent.getByToken(minipfNewIso_ch_, minipfNewIso_ch);
  iEvent.getByToken(minipfNewIso_nh_, minipfNewIso_nh);
  iEvent.getByToken(minipfNewIso_ph_, minipfNewIso_ph);  
  //iEvent.getByToken(minipfNewIso_pu_, minipfNewIso_pu);  
    

  Handle<std::vector<pat::PackedCandidate>> pfCands;
  iEvent.getByToken(pfCandsToken_, pfCands);

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

  // gen muons
  b_muon_no = 0;
  for (const reco::GenParticle &gen : *pruned) {
     TLorentzVector gentlv(gen.momentum().x(), gen.momentum().y(), gen.momentum().z(), gen.energy() );
     edm::RefToBase<pat::Muon> muref;
     pat::PFIsolation miniiso;
     for (size_t i = 0; i < muons->size(); i++) {
       auto muon = muons->at(i);

       //!isSignal
       if (muon.simExtType() == reco::ExtendedMuonSimType::ExtNotMatched ||
           muon.simExtType() == reco::ExtendedMuonSimType::ExtUnknown) { continue; }

       if (gen.pt() == muon.simPt() && gen.eta() == muon.simEta() && gen.phi() == muon.simPhi()){
	      muref = muons->refAt(i);
          TLorentzVector recotlv(muon.momentum().x(), muon.momentum().y(), muon.momentum().z(), muon.energy() );
          miniiso = pat::getMiniPFIsolation(&(*pfCands), muon.p4(),
                                              miniIsoParams_[0], miniIsoParams_[1], miniIsoParams_[2],
                                              miniIsoParams_[3], miniIsoParams_[4], miniIsoParams_[5],
                                              miniIsoParams_[6], miniIsoParams_[7], miniIsoParams_[8]);
          break;
       }
     }
     fillBranches(genttree_, gentlv, muref, true, gen.pdgId(), miniiso);
  }

  // reco muons
  b_muon_no = 0;
  for (size_t i = 0; i < muons->size(); i++) {
    edm::RefToBase<pat::Muon> muref = muons->refAt(i);
    auto muon = muons->at(i);
    if (abs(muon.eta()) > 2.8) continue;

    bool isMatchedMuon = true;
    //Check using siminfo
    if (muon.simExtType() == reco::ExtendedMuonSimType::ExtNotMatched ||
        muon.simExtType() == reco::ExtendedMuonSimType::ExtUnknown) {
          isMatchedMuon = false;
    }

    bool isSignalMuon = false;
    // is it from pion, keon, or non primary?
    if (muon.simExtType() == reco::ExtendedMuonSimType::MatchedMuonFromPiKppMuX ||
        muon.simExtType() == reco::ExtendedMuonSimType::MatchedMuonFromPiKNotppMuX ||
        muon.simExtType() == reco::ExtendedMuonSimType::MatchedMuonFromNonPrimaryParticle) {
        if ( isMatchedMuon ) isSignalMuon = true;
    }

    bool isGhostMuon = false;
    // is it ghost muon?
    if ( muon.simExtType() < 0){ isGhostMuon = true; }

    cout << isSignalMuon << isMatchedMuon << isGhostMuon << endl;

    TLorentzVector recotlv(muon.momentum().x(), muon.momentum().y(), muon.momentum().z(), muon.energy() );
    pat::PFIsolation miniiso = pat::getMiniPFIsolation(&(*pfCands), muon.p4(),
                                                        miniIsoParams_[0], miniIsoParams_[1], miniIsoParams_[2],
                                                        miniIsoParams_[3], miniIsoParams_[4], miniIsoParams_[5],
                                                        miniIsoParams_[6], miniIsoParams_[7], miniIsoParams_[8]);
    fillBranches(recottree_, recotlv, muref, isMatchedMuon, muon.pdgId(), miniiso);
  }

  return;
}


void NewPatMuonAnalyser::fillBranches(TTree *tree, TLorentzVector &tlv, edm::RefToBase<pat::Muon> muref, bool isSignal, int pdgId, pat::PFIsolation miniiso)
{
  initTreeValues();
  b_muon = tlv;
  b_muon_signal = isSignal;
  b_muon_pdgId = pdgId;
  ++b_muon_no;

  b_muon_miniIso_ch = miniiso.chargedHadronIso();
  b_muon_miniIso_nh = miniiso.neutralHadronIso();
  b_muon_miniIso_ph = miniiso.photonIso();
  b_muon_miniIso_pu = miniiso.puChargedHadronIso(); 
  
  if (muref.isNonnull()){
    auto muon = muref;
    
    b_muon_poszPV0  = priVertex_.position().z();
    b_muon_poszMuon = muon->vz();
    
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

    b_muon_isTight = muon::isTightMuon(*muon, priVertex_);
    b_muon_isMedium = muon::isMediumMuon(*muon);
    b_muon_isLoose = muon::isLooseMuon(*muon);
    
    b_muon_puppiNewIso_ch = (*puppiNewIso_ch)[muref];
    b_muon_puppiNewIso_nh = (*puppiNewIso_nh)[muref];
    b_muon_puppiNewIso_ph = (*puppiNewIso_ph)[muref];
    b_muon_puppiNewIso    = ( b_muon_puppiNewIso_ch + b_muon_puppiNewIso_nh + b_muon_puppiNewIso_ph )/muon->pt();
    b_muon_trkNewIso = (*trkNewIso)[muref] / muon->pt();
    b_muon_pfNewIso_ch = (*pfNewIso_ch)[muref];
    b_muon_pfNewIso_nh = (*pfNewIso_nh)[muref];
    b_muon_pfNewIso_ph = (*pfNewIso_ph)[muref];
    //b_muon_pfNewIso_pu = muon->pfIsolationR03().sumPUPt;
    //b_muon_pfNewIso_pu = (*pfNewIso_pu)[muref];
    b_muon_pfNewIso    = ( b_muon_pfNewIso_ch + max(0.0, b_muon_pfNewIso_nh + b_muon_pfNewIso_ph - 0.5 ) )/ muon->pt();
    //b_muon_pfNewIso    = ( b_muon_pfNewIso_ch + max(0.0, b_muon_pfNewIso_nh + b_muon_pfNewIso_ph - 0.5 * b_muon_pfNewIso_pu) ) / muon->pt();
    b_muon_minipuppiNewIso_ch = (*minipuppiNewIso_ch)[muref];
    b_muon_minipuppiNewIso_nh = (*minipuppiNewIso_nh)[muref];
    b_muon_minipuppiNewIso_ph = (*minipuppiNewIso_ph)[muref];
    //b_muon_minipuppiNewIso_pu = muon->pfIsolationR03().sumPUPt;
    b_muon_minipuppiNewIso    = ( b_muon_minipuppiNewIso_ch + b_muon_minipuppiNewIso_nh + b_muon_minipuppiNewIso_ph ) / muon->pt();
    b_muon_minitrkNewIso = (*minitrkNewIso)[muref] / muon->pt();
    b_muon_minipfNewIso_ch = (*minipfNewIso_ch)[muref];
    b_muon_minipfNewIso_nh = (*minipfNewIso_nh)[muref];
    b_muon_minipfNewIso_ph = (*minipfNewIso_ph)[muref];
    //b_muon_minipfNewIso_pu = muon->pfIsolationR03().sumPUPt;
    //b_muon_minipfNewIso_pu = (*minipfNewIso_pu)[muref];
    b_muon_minipfNewIso    = ( b_muon_minipfNewIso_ch + max(0.0, b_muon_minipfNewIso_nh + b_muon_minipfNewIso_ph - 0.5 ) )/ muon->pt();
    //b_muon_minipfNewIso    = ( b_muon_minipfNewIso_ch + max(0.0, b_muon_minipfNewIso_nh + b_muon_minipfNewIso_ph - 0.5 * b_muon_minipfNewIso_pu) ) / muon->pt();
    
  }
  tree->Fill();
}

void NewPatMuonAnalyser::setBranches(TTree *tree)
{
  tree->Branch("nvertex", &b_nvertex, "nvertex/I");
  tree->Branch("pu_density", &b_pu_density, "pu_density/I");
  tree->Branch("pu_numInteractions", &b_pu_numInteractions, "pu_numInteractions/I");
  tree->Branch("muon", "TLorentzVector", &b_muon);  
  tree->Branch("muon_no", &b_muon_no, "muon_no/I");
  tree->Branch("muon_pdgId", &b_muon_pdgId, "muon_pdgId/I");
  tree->Branch("muon_poszPV0",&b_muon_poszPV0,"muon_poszPV0/F");
  tree->Branch("muon_poszMuon",&b_muon_poszMuon,"muon_poszMuon/F");
  tree->Branch("muon_signal", &b_muon_signal, "muon_signal/O");

  tree->Branch("muon_miniIso_ch", &b_muon_miniIso_ch, "muon_miniIso_ch/F");
  tree->Branch("muon_miniIso_nh", &b_muon_miniIso_nh, "muon_miniIso_nh/F");
  tree->Branch("muon_miniIso_ph", &b_muon_miniIso_ph, "muon_miniIso_ph/F");
  tree->Branch("muon_miniIso_pu", &b_muon_miniIso_pu, "muon_miniIso_pu/F");

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
  //tree->Branch("muon_pfNewIso_pu",&b_muon_pfNewIso_pu,"muon_pfNewIso_pu/F");
  tree->Branch("muon_pfNewIso",&b_muon_pfNewIso,"muon_pfNewIso/F");
  tree->Branch("muon_minipuppiNewIso_ch",&b_muon_minipuppiNewIso_ch,"muon_minipuppiNewIso_ch/F");
  tree->Branch("muon_minipuppiNewIso_nh",&b_muon_minipuppiNewIso_nh,"muon_minipuppiNewIso_nh/F");
  tree->Branch("muon_minipuppiNewIso_ph",&b_muon_minipuppiNewIso_ph,"muon_minipuppiNewIso_ph/F");
  tree->Branch("muon_minipuppiNewIso",&b_muon_minipuppiNewIso,"muon_minipuppiNewIso/F");
  tree->Branch("muon_minitrkNewIso",&b_muon_minitrkNewIso,"muon_minitrkNewIso/F");
  tree->Branch("muon_minipfNewIso_ch",&b_muon_minipfNewIso_ch,"muon_minipfNewIso_ch/F");
  tree->Branch("muon_minipfNewIso_nh",&b_muon_minipfNewIso_nh,"muon_minipfNewIso_nh/F");
  tree->Branch("muon_minipfNewIso_ph",&b_muon_minipfNewIso_ph,"muon_minipfNewIso_ph/F");
  //tree->Branch("muon_minipfNewIso_pu",&b_muon_minipfNewIso_pu,"muon_minipfNewIso_pu/F");
  tree->Branch("muon_minipfNewIso",&b_muon_minipfNewIso,"muon_minipfNewIso/F");

  tree->Branch("muon_PFIsoFixOnlyCH",&b_muon_PFIsoFixOnlyCH,"muon_PFIsoFixOnlyCH/F");
  tree->Branch("muon_puppiIsoFixOnlyCH",&b_muon_puppiIsoFixOnlyCH,"muon_puppiIsoFixOnlyCH/F");
  tree->Branch("muon_PFIsoRepTrk",&b_muon_PFIsoRepTrk,"muon_PFIsoRepTrk/F");
  tree->Branch("muon_puppiIsoRepTrk",&b_muon_puppiIsoRepTrk,"muon_puppiIsoRepTrk/F");
}

void NewPatMuonAnalyser::initTreeValues(){
  b_muon_miniIso_ch = -999;
  b_muon_miniIso_nh= -999;
  b_muon_miniIso_ph = -999;
  b_muon_miniIso_pu = -999;

  b_muon_poszPV0  = 0;
  b_muon_poszMuon = 0;
    
  b_muon_isTight = 0; b_muon_isMedium = 0; b_muon_isLoose = 0;
  
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
  b_muon_pfNewIso_ch = -1; b_muon_pfNewIso_nh = -1; b_muon_pfNewIso_ph = -1; /*b_muon_pfNewIso_pu = -1;*/ b_muon_pfNewIso = -1;
  b_muon_minipuppiNewIso_ch = -1; b_muon_minipuppiNewIso_nh = -1; b_muon_minipuppiNewIso_ph = -1; b_muon_minipuppiNewIso_pu = -1; b_muon_minipuppiNewIso = -1;
  b_muon_minitrkNewIso = -1;
  b_muon_minipfNewIso_ch = -1; b_muon_minipfNewIso_nh = -1; b_muon_minipfNewIso_ph = -1; /*b_muon_minipfNewIso_pu = -1; */b_muon_minipfNewIso = -1;

  b_muon_PFIsoFixOnlyCH = 0;
  b_muon_puppiIsoFixOnlyCH = 0;
  
  b_muon_PFIsoRepTrk = 0;
  b_muon_puppiIsoRepTrk = 0;
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
