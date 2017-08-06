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

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/GEMGeometry/interface/ME0Geometry.h"

#include <vector>
#include <TTree.h>
#include <TFile.h>
#include <TLorentzVector.h>

using namespace std;

class PatMuonAnalyser : public edm::one::EDAnalyzer<edm::one::SharedResources>  {

public:
  explicit PatMuonAnalyser(const edm::ParameterSet&);
  ~PatMuonAnalyser(){};

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

  void setBranches(TTree *tree);
  void fillBranches(TTree *tree, TLorentzVector &tlv, edm::RefToBase<pat::Muon> muref,bool isSignal, int pdgId);

  bool isSignalMuon(const reco::GenParticle &gen);
  
  bool isME0MuonSelNew(reco::Muon muon, double dEtaCut, double dPhiCut, double dPhiBendCut);
  
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
  edm::EDGetTokenT<edm::ValueMap<float> > pfNewIso_ch_;
  edm::EDGetTokenT<edm::ValueMap<float> > pfNewIso_nh_;
  edm::EDGetTokenT<edm::ValueMap<float> > pfNewIso_ph_;
  edm::Handle<edm::ValueMap<float>> pfNewIso_ch;
  edm::Handle<edm::ValueMap<float>> pfNewIso_nh;
  edm::Handle<edm::ValueMap<float>> pfNewIso_ph;

  edm::EDGetTokenT<edm::ValueMap<float> > puppiNewIsoPt02_ch_;
  edm::EDGetTokenT<edm::ValueMap<float> > puppiNewIsoPt02_nh_;
  edm::EDGetTokenT<edm::ValueMap<float> > puppiNewIsoPt02_ph_;
  edm::Handle<edm::ValueMap<float>> puppiNewIsoPt02_ch;
  edm::Handle<edm::ValueMap<float>> puppiNewIsoPt02_nh;
  edm::Handle<edm::ValueMap<float>> puppiNewIsoPt02_ph;
  edm::EDGetTokenT<edm::ValueMap<float> > pfNewIsoPt02_ch_;
  edm::EDGetTokenT<edm::ValueMap<float> > pfNewIsoPt02_nh_;
  edm::EDGetTokenT<edm::ValueMap<float> > pfNewIsoPt02_ph_;
  edm::Handle<edm::ValueMap<float>> pfNewIsoPt02_ch;
  edm::Handle<edm::ValueMap<float>> pfNewIsoPt02_nh;
  edm::Handle<edm::ValueMap<float>> pfNewIsoPt02_ph;
  
  edm::EDGetTokenT<edm::ValueMap<float> > puppiNewIsoPt04_ch_;
  edm::EDGetTokenT<edm::ValueMap<float> > puppiNewIsoPt04_nh_;
  edm::EDGetTokenT<edm::ValueMap<float> > puppiNewIsoPt04_ph_;
  edm::Handle<edm::ValueMap<float>> puppiNewIsoPt04_ch;
  edm::Handle<edm::ValueMap<float>> puppiNewIsoPt04_nh;
  edm::Handle<edm::ValueMap<float>> puppiNewIsoPt04_ph;
  edm::EDGetTokenT<edm::ValueMap<float> > pfNewIsoPt04_ch_;
  edm::EDGetTokenT<edm::ValueMap<float> > pfNewIsoPt04_nh_;
  edm::EDGetTokenT<edm::ValueMap<float> > pfNewIsoPt04_ph_;
  edm::Handle<edm::ValueMap<float>> pfNewIsoPt04_ch;
  edm::Handle<edm::ValueMap<float>> pfNewIsoPt04_nh;
  edm::Handle<edm::ValueMap<float>> pfNewIsoPt04_ph;

  edm::EDGetTokenT<edm::ValueMap<float> > puppiNewIsoPt06_ch_;
  edm::EDGetTokenT<edm::ValueMap<float> > puppiNewIsoPt06_nh_;
  edm::EDGetTokenT<edm::ValueMap<float> > puppiNewIsoPt06_ph_;
  edm::Handle<edm::ValueMap<float>> puppiNewIsoPt06_ch;
  edm::Handle<edm::ValueMap<float>> puppiNewIsoPt06_nh;
  edm::Handle<edm::ValueMap<float>> puppiNewIsoPt06_ph;
  edm::EDGetTokenT<edm::ValueMap<float> > pfNewIsoPt06_ch_;
  edm::EDGetTokenT<edm::ValueMap<float> > pfNewIsoPt06_nh_;
  edm::EDGetTokenT<edm::ValueMap<float> > pfNewIsoPt06_ph_;
  edm::Handle<edm::ValueMap<float>> pfNewIsoPt06_ch;
  edm::Handle<edm::ValueMap<float>> pfNewIsoPt06_nh;
  edm::Handle<edm::ValueMap<float>> pfNewIsoPt06_ph;
  
  edm::EDGetTokenT<edm::ValueMap<float> > puppiNewIsoPt08_ch_;
  edm::EDGetTokenT<edm::ValueMap<float> > puppiNewIsoPt08_nh_;
  edm::EDGetTokenT<edm::ValueMap<float> > puppiNewIsoPt08_ph_;
  edm::Handle<edm::ValueMap<float>> puppiNewIsoPt08_ch;
  edm::Handle<edm::ValueMap<float>> puppiNewIsoPt08_nh;
  edm::Handle<edm::ValueMap<float>> puppiNewIsoPt08_ph;
  edm::EDGetTokenT<edm::ValueMap<float> > pfNewIsoPt08_ch_;
  edm::EDGetTokenT<edm::ValueMap<float> > pfNewIsoPt08_nh_;
  edm::EDGetTokenT<edm::ValueMap<float> > pfNewIsoPt08_ph_;
  edm::Handle<edm::ValueMap<float>> pfNewIsoPt08_ch;
  edm::Handle<edm::ValueMap<float>> pfNewIsoPt08_nh;
  edm::Handle<edm::ValueMap<float>> pfNewIsoPt08_ph;
  
  edm::EDGetTokenT<edm::ValueMap<float> > puppiNewIsoPt10_ch_;
  edm::EDGetTokenT<edm::ValueMap<float> > puppiNewIsoPt10_nh_;
  edm::EDGetTokenT<edm::ValueMap<float> > puppiNewIsoPt10_ph_;
  edm::Handle<edm::ValueMap<float>> puppiNewIsoPt10_ch;
  edm::Handle<edm::ValueMap<float>> puppiNewIsoPt10_nh;
  edm::Handle<edm::ValueMap<float>> puppiNewIsoPt10_ph;
  edm::EDGetTokenT<edm::ValueMap<float> > pfNewIsoPt10_ch_;
  edm::EDGetTokenT<edm::ValueMap<float> > pfNewIsoPt10_nh_;
  edm::EDGetTokenT<edm::ValueMap<float> > pfNewIsoPt10_ph_;
  edm::Handle<edm::ValueMap<float>> pfNewIsoPt10_ch;
  edm::Handle<edm::ValueMap<float>> pfNewIsoPt10_nh;
  edm::Handle<edm::ValueMap<float>> pfNewIsoPt10_ph;

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
  float b_muon_poszPV0, b_muon_poszMuon;
  bool b_muon_isTight, b_muon_isMedium, b_muon_isLoose;
  bool b_muon_isME0MuonTight, b_muon_isME0MuonMedium, b_muon_isME0MuonLoose;

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
  float b_muon_pfNewIso_ch, b_muon_pfNewIso_nh, b_muon_pfNewIso_ph, b_muon_pfNewIso_pu, b_muon_pfNewIso;

  float b_muon_puppiNewIsoPt02_ch, b_muon_puppiNewIsoPt02_nh, b_muon_puppiNewIsoPt02_ph, b_muon_puppiNewIsoPt02;
  float b_muon_pfNewIsoPt02_ch, b_muon_pfNewIsoPt02_nh, b_muon_pfNewIsoPt02_ph, b_muon_pfNewIsoPt02_pu, b_muon_pfNewIsoPt02;
  
  float b_muon_puppiNewIsoPt04_ch, b_muon_puppiNewIsoPt04_nh, b_muon_puppiNewIsoPt04_ph, b_muon_puppiNewIsoPt04;
  float b_muon_pfNewIsoPt04_ch, b_muon_pfNewIsoPt04_nh, b_muon_pfNewIsoPt04_ph, b_muon_pfNewIsoPt04_pu, b_muon_pfNewIsoPt04;

  float b_muon_puppiNewIsoPt06_ch, b_muon_puppiNewIsoPt06_nh, b_muon_puppiNewIsoPt06_ph, b_muon_puppiNewIsoPt06;
  float b_muon_pfNewIsoPt06_ch, b_muon_pfNewIsoPt06_nh, b_muon_pfNewIsoPt06_ph, b_muon_pfNewIsoPt06_pu, b_muon_pfNewIsoPt06;

  float b_muon_puppiNewIsoPt08_ch, b_muon_puppiNewIsoPt08_nh, b_muon_puppiNewIsoPt08_ph, b_muon_puppiNewIsoPt08;
  float b_muon_pfNewIsoPt08_ch, b_muon_pfNewIsoPt08_nh, b_muon_pfNewIsoPt08_ph, b_muon_pfNewIsoPt08_pu, b_muon_pfNewIsoPt08;

  float b_muon_puppiNewIsoPt10_ch, b_muon_puppiNewIsoPt10_nh, b_muon_puppiNewIsoPt10_ph, b_muon_puppiNewIsoPt10;
  float b_muon_pfNewIsoPt10_ch, b_muon_pfNewIsoPt10_nh, b_muon_pfNewIsoPt10_ph, b_muon_pfNewIsoPt10_pu, b_muon_pfNewIsoPt10;

  float b_muon_PFIsoFixOnlyCH;
  float b_muon_puppiIsoFixOnlyCH;
  
  float b_muon_PFIsoRepTrk;
  float b_muon_puppiIsoRepTrk;

};
PatMuonAnalyser::PatMuonAnalyser(const edm::ParameterSet& iConfig):
  verticesToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
  muonsToken_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
  prunedGenToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned")))
{
  putoken_ = consumes<std::vector<PileupSummaryInfo>>(edm::InputTag("addPileupInfo"));
  
  puppiNewIso_ch_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("puppiNewIso_ch"));
  puppiNewIso_nh_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("puppiNewIso_nh"));
  puppiNewIso_ph_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("puppiNewIso_ph"));
  pfNewIso_ch_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("pfNewIso_ch"));
  pfNewIso_nh_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("pfNewIso_nh"));
  pfNewIso_ph_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("pfNewIso_ph"));

  puppiNewIsoPt02_ch_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("puppiNewIsoPt02_ch"));
  puppiNewIsoPt02_nh_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("puppiNewIsoPt02_nh"));
  puppiNewIsoPt02_ph_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("puppiNewIsoPt02_ph"));
  pfNewIsoPt02_ch_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("pfNewIsoPt02_ch"));
  pfNewIsoPt02_nh_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("pfNewIsoPt02_nh"));
  pfNewIsoPt02_ph_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("pfNewIsoPt02_ph"));
  
  puppiNewIsoPt04_ch_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("puppiNewIsoPt04_ch"));
  puppiNewIsoPt04_nh_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("puppiNewIsoPt04_nh"));
  puppiNewIsoPt04_ph_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("puppiNewIsoPt04_ph"));
  pfNewIsoPt04_ch_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("pfNewIsoPt04_ch"));
  pfNewIsoPt04_nh_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("pfNewIsoPt04_nh"));
  pfNewIsoPt04_ph_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("pfNewIsoPt04_ph"));
  
  puppiNewIsoPt06_ch_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("puppiNewIsoPt06_ch"));
  puppiNewIsoPt06_nh_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("puppiNewIsoPt06_nh"));
  puppiNewIsoPt06_ph_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("puppiNewIsoPt06_ph"));
  pfNewIsoPt06_ch_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("pfNewIsoPt06_ch"));
  pfNewIsoPt06_nh_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("pfNewIsoPt06_nh"));
  pfNewIsoPt06_ph_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("pfNewIsoPt06_ph"));
  
  puppiNewIsoPt08_ch_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("puppiNewIsoPt08_ch"));
  puppiNewIsoPt08_nh_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("puppiNewIsoPt08_nh"));
  puppiNewIsoPt08_ph_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("puppiNewIsoPt08_ph"));
  pfNewIsoPt08_ch_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("pfNewIsoPt08_ch"));
  pfNewIsoPt08_nh_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("pfNewIsoPt08_nh"));
  pfNewIsoPt08_ph_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("pfNewIsoPt08_ph"));
  
  puppiNewIsoPt10_ch_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("puppiNewIsoPt10_ch"));
  puppiNewIsoPt10_nh_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("puppiNewIsoPt10_nh"));
  puppiNewIsoPt10_ph_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("puppiNewIsoPt10_ph"));
  pfNewIsoPt10_ch_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("pfNewIsoPt10_ch"));
  pfNewIsoPt10_nh_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("pfNewIsoPt10_nh"));
  pfNewIsoPt10_ph_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("pfNewIsoPt10_ph"));
  
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  h_nevents = fs->make<TH1D>("nevents", "nevents", 1, 0, 1);
  
  genttree_ = fs->make<TTree>("gen", "gen");
  setBranches(genttree_);
  recottree_ = fs->make<TTree>("reco", "reco");
  setBranches(recottree_);
}

void PatMuonAnalyser::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
  iEvent.getByToken(pfNewIso_ch_, pfNewIso_ch);
  iEvent.getByToken(pfNewIso_nh_, pfNewIso_nh);
  iEvent.getByToken(pfNewIso_ph_, pfNewIso_ph);  
    
  iEvent.getByToken(puppiNewIsoPt02_ch_, puppiNewIsoPt02_ch);
  iEvent.getByToken(puppiNewIsoPt02_nh_, puppiNewIsoPt02_nh);
  iEvent.getByToken(puppiNewIsoPt02_ph_, puppiNewIsoPt02_ph);  
  iEvent.getByToken(pfNewIsoPt02_ch_, pfNewIsoPt02_ch);
  iEvent.getByToken(pfNewIsoPt02_nh_, pfNewIsoPt02_nh);
  iEvent.getByToken(pfNewIsoPt02_ph_, pfNewIsoPt02_ph);

  iEvent.getByToken(puppiNewIsoPt04_ch_, puppiNewIsoPt04_ch);
  iEvent.getByToken(puppiNewIsoPt04_nh_, puppiNewIsoPt04_nh);
  iEvent.getByToken(puppiNewIsoPt04_ph_, puppiNewIsoPt04_ph);  
  iEvent.getByToken(pfNewIsoPt04_ch_, pfNewIsoPt04_ch);
  iEvent.getByToken(pfNewIsoPt04_nh_, pfNewIsoPt04_nh);
  iEvent.getByToken(pfNewIsoPt04_ph_, pfNewIsoPt04_ph);  
  
  iEvent.getByToken(puppiNewIsoPt06_ch_, puppiNewIsoPt06_ch);
  iEvent.getByToken(puppiNewIsoPt06_nh_, puppiNewIsoPt06_nh);
  iEvent.getByToken(puppiNewIsoPt06_ph_, puppiNewIsoPt06_ph);  
  iEvent.getByToken(pfNewIsoPt06_ch_, pfNewIsoPt06_ch);
  iEvent.getByToken(pfNewIsoPt06_nh_, pfNewIsoPt06_nh);
  iEvent.getByToken(pfNewIsoPt06_ph_, pfNewIsoPt06_ph);  

  iEvent.getByToken(puppiNewIsoPt08_ch_, puppiNewIsoPt08_ch);
  iEvent.getByToken(puppiNewIsoPt08_nh_, puppiNewIsoPt08_nh);
  iEvent.getByToken(puppiNewIsoPt08_ph_, puppiNewIsoPt08_ph);  
  iEvent.getByToken(pfNewIsoPt08_ch_, pfNewIsoPt08_ch);
  iEvent.getByToken(pfNewIsoPt08_nh_, pfNewIsoPt08_nh);
  iEvent.getByToken(pfNewIsoPt08_ph_, pfNewIsoPt08_ph);  

  iEvent.getByToken(puppiNewIsoPt10_ch_, puppiNewIsoPt10_ch);
  iEvent.getByToken(puppiNewIsoPt10_nh_, puppiNewIsoPt10_nh);
  iEvent.getByToken(puppiNewIsoPt10_ph_, puppiNewIsoPt10_ph);  
  iEvent.getByToken(pfNewIsoPt10_ch_, pfNewIsoPt10_ch);
  iEvent.getByToken(pfNewIsoPt10_nh_, pfNewIsoPt10_nh);
  iEvent.getByToken(pfNewIsoPt10_ph_, pfNewIsoPt10_ph);  

  // edm::Handle<std::vector <PileupSummaryInfo> > PupInfo;
  // iEvent.getByToken(putoken, PupInfo);
  // b_pu_density = 0; b_pu_numInteractions = 0;
  // std::vector<PileupSummaryInfo>::const_iterator ipu;
  // for (ipu = PupInfo->begin(); ipu != PupInfo->end(); ++ipu) {
  //   if ( ipu->getBunchCrossing() != 0 ) continue; // storing detailed PU info only for BX=0
  //   for (unsigned int i=0; i<ipu->getPU_zpositions().size(); ++i) {
  //     if ( abs((ipu->getPU_zpositions())[i] - simVertex_.position().z()) < 0.1 )
  // 	++b_pu_density;
  //     ++b_pu_numInteractions;
  //   }
  // }

  b_muon_no = 0;
  for (const reco::GenParticle &gen : *pruned) {
    if (!isSignalMuon(gen)) continue;
    TLorentzVector gentlv(gen.momentum().x(), gen.momentum().y(), gen.momentum().z(), gen.energy() );

    
    edm::RefToBase<pat::Muon> muref;

    for (size_t i = 0; i < muons->size(); i++) {
      auto muon = muons->at(i);
      
      TLorentzVector recotlv(muon.momentum().x(), muon.momentum().y(), muon.momentum().z(), muon.energy() );
      if (gentlv.DeltaR(recotlv) < 0.1){
	muref = muons->refAt(i);
	break;
      }
    }
    
    fillBranches(genttree_, gentlv, muref, true, gen.pdgId());
  }

  
  b_muon_no = 0;
  for (size_t i = 0; i < muons->size(); i++) {
    edm::RefToBase<pat::Muon> muref = muons->refAt(i);
    auto muon = muons->at(i);
  
    if (abs(muon.eta()) > 2.8) continue;

    bool isSignal = false;
    const reco::GenParticle *gen = muon.genLepton();
    int pdgId = 0;
    if (gen){
      pdgId = gen->pdgId();
      isSignal = isSignalMuon(*gen);
    }
  
    TLorentzVector recotlv(muon.momentum().x(), muon.momentum().y(), muon.momentum().z(), muon.energy() );
    fillBranches(recottree_, recotlv, muref, isSignal, pdgId);
  }
  
  return;
}

bool PatMuonAnalyser::isSignalMuon(const reco::GenParticle &gen)
{
  if (abs(gen.pdgId()) != 13) return false;
  
  for (unsigned int i = 0; i < gen.numberOfMothers(); ++i){
    //In case of pdgId() = 23, indicate Z-boson. if it's 25, that becomes higgs.
    if (gen.mother(i)->pdgId() == 23 || gen.mother(i)->pdgId() == 25){
      return true;
    }
  }
  return false;
  
}

void PatMuonAnalyser::fillBranches(TTree *tree, TLorentzVector &tlv, edm::RefToBase<pat::Muon> muref, bool isSignal, int pdgId)
{
  b_muon = tlv;
  b_muon_signal = isSignal;
  b_muon_pdgId = pdgId;
  ++b_muon_no;
  
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
  b_muon_pfNewIso_ch = -1; b_muon_pfNewIso_nh = -1; b_muon_pfNewIso_ph = -1; b_muon_pfNewIso_pu = -1; b_muon_pfNewIso = -1;
    
  b_muon_puppiNewIsoPt02_ch = -1; b_muon_puppiNewIsoPt02_nh = -1; b_muon_puppiNewIsoPt02_ph = -1; b_muon_puppiNewIsoPt02 = -1;
  b_muon_pfNewIsoPt02_ch = -1; b_muon_pfNewIsoPt02_nh = -1; b_muon_pfNewIsoPt02_ph = -1; b_muon_pfNewIsoPt02_pu = -1; b_muon_pfNewIsoPt02 = -1;

  b_muon_puppiNewIsoPt04_ch = -1; b_muon_puppiNewIsoPt04_nh = -1; b_muon_puppiNewIsoPt04_ph = -1; b_muon_puppiNewIsoPt04 = -1;
  b_muon_pfNewIsoPt04_ch = -1; b_muon_pfNewIsoPt04_nh = -1; b_muon_pfNewIsoPt04_ph = -1; b_muon_pfNewIsoPt04_pu = -1; b_muon_pfNewIsoPt04 = -1;

  b_muon_puppiNewIsoPt06_ch = -1; b_muon_puppiNewIsoPt06_nh = -1; b_muon_puppiNewIsoPt06_ph = -1; b_muon_puppiNewIsoPt06 = -1;
  b_muon_pfNewIsoPt06_ch = -1; b_muon_pfNewIsoPt06_nh = -1; b_muon_pfNewIsoPt06_ph = -1; b_muon_pfNewIsoPt06_pu = -1; b_muon_pfNewIsoPt06 = -1;

  b_muon_puppiNewIsoPt08_ch = -1; b_muon_puppiNewIsoPt08_nh = -1; b_muon_puppiNewIsoPt08_ph = -1; b_muon_puppiNewIsoPt08 = -1;
  b_muon_pfNewIsoPt08_ch = -1; b_muon_pfNewIsoPt08_nh = -1; b_muon_pfNewIsoPt08_ph = -1; b_muon_pfNewIsoPt08_pu = -1; b_muon_pfNewIsoPt08 = -1;

  b_muon_puppiNewIsoPt10_ch = -1; b_muon_puppiNewIsoPt10_nh = -1; b_muon_puppiNewIsoPt10_ph = -1; b_muon_puppiNewIsoPt10 = -1;
  b_muon_pfNewIsoPt10_ch = -1; b_muon_pfNewIsoPt10_nh = -1; b_muon_pfNewIsoPt10_ph = -1; b_muon_pfNewIsoPt10_pu = -1; b_muon_pfNewIsoPt10 = -1;

  b_muon_PFIsoFixOnlyCH = 0;
  b_muon_puppiIsoFixOnlyCH = 0;
  
  b_muon_PFIsoRepTrk = 0;
  b_muon_puppiIsoRepTrk = 0;

  b_muon_isME0MuonTight = 0; b_muon_isME0MuonMedium = 0; b_muon_isME0MuonLoose = 0;
  
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
    
    double mom = muon->p();
    double dPhiCut_ = std::min(std::max(1.2/mom,1.2/100),0.056);
    double dPhiBendCut_ = std::min(std::max(0.2/mom,0.2/100),0.0096);

    b_muon_isME0MuonLoose = isME0MuonSelNew(*muon, 0.077, dPhiCut_, dPhiBendCut_);
    
    b_muon_puppiNewIso_ch = (*puppiNewIso_ch)[muref];
    b_muon_puppiNewIso_nh = (*puppiNewIso_nh)[muref];
    b_muon_puppiNewIso_ph = (*puppiNewIso_ph)[muref];
    b_muon_puppiNewIso    = ( b_muon_puppiNewIso_ch + b_muon_puppiNewIso_nh + b_muon_puppiNewIso_ph )/muon->pt();
    b_muon_pfNewIso_ch = (*pfNewIso_ch)[muref];
    b_muon_pfNewIso_nh = (*pfNewIso_nh)[muref];
    b_muon_pfNewIso_ph = (*pfNewIso_ph)[muref];
    b_muon_pfNewIso_pu = muon->pfIsolationR03().sumPUPt;
    b_muon_pfNewIso    = ( b_muon_pfNewIso_ch + max(0.0, b_muon_pfNewIso_nh + b_muon_pfNewIso_ph - 0.5 * b_muon_pfNewIso_pu) ) / muon->pt();

    b_muon_puppiNewIsoPt02_ch = (*puppiNewIsoPt02_ch)[muref];
    b_muon_puppiNewIsoPt02_nh = (*puppiNewIsoPt02_nh)[muref];
    b_muon_puppiNewIsoPt02_ph = (*puppiNewIsoPt02_ph)[muref];
    b_muon_puppiNewIsoPt02    = ( b_muon_puppiNewIsoPt02_ch + b_muon_puppiNewIsoPt02_nh + b_muon_puppiNewIsoPt02_ph )/muon->pt();    
    b_muon_pfNewIsoPt02_ch = (*pfNewIsoPt02_ch)[muref];
    b_muon_pfNewIsoPt02_nh = (*pfNewIsoPt02_nh)[muref];
    b_muon_pfNewIsoPt02_ph = (*pfNewIsoPt02_ph)[muref];
    b_muon_pfNewIsoPt02_pu = muon->pfIsolationR03().sumPUPt;
    b_muon_pfNewIsoPt02    = ( b_muon_pfNewIsoPt02_ch + max(0.0, b_muon_pfNewIsoPt02_nh + b_muon_pfNewIsoPt02_ph - 0.5 * b_muon_pfNewIsoPt02_pu) ) / muon->pt();
    
    b_muon_puppiNewIsoPt04_ch = (*puppiNewIsoPt04_ch)[muref];
    b_muon_puppiNewIsoPt04_nh = (*puppiNewIsoPt04_nh)[muref];
    b_muon_puppiNewIsoPt04_ph = (*puppiNewIsoPt04_ph)[muref];
    b_muon_puppiNewIsoPt04    = ( b_muon_puppiNewIsoPt04_ch + b_muon_puppiNewIsoPt04_nh + b_muon_puppiNewIsoPt04_ph )/muon->pt();    
    b_muon_pfNewIsoPt04_ch = (*pfNewIsoPt04_ch)[muref];
    b_muon_pfNewIsoPt04_nh = (*pfNewIsoPt04_nh)[muref];
    b_muon_pfNewIsoPt04_ph = (*pfNewIsoPt04_ph)[muref];
    b_muon_pfNewIsoPt04_pu = muon->pfIsolationR03().sumPUPt;
    b_muon_pfNewIsoPt04    = ( b_muon_pfNewIsoPt04_ch + max(0.0, b_muon_pfNewIsoPt04_nh + b_muon_pfNewIsoPt04_ph - 0.5 * b_muon_pfNewIsoPt04_pu) ) / muon->pt();

    b_muon_puppiNewIsoPt06_ch = (*puppiNewIsoPt06_ch)[muref];
    b_muon_puppiNewIsoPt06_nh = (*puppiNewIsoPt06_nh)[muref];
    b_muon_puppiNewIsoPt06_ph = (*puppiNewIsoPt06_ph)[muref];
    b_muon_puppiNewIsoPt06    = ( b_muon_puppiNewIsoPt06_ch + b_muon_puppiNewIsoPt06_nh + b_muon_puppiNewIsoPt06_ph )/muon->pt();    
    b_muon_pfNewIsoPt06_ch = (*pfNewIsoPt06_ch)[muref];
    b_muon_pfNewIsoPt06_nh = (*pfNewIsoPt06_nh)[muref];
    b_muon_pfNewIsoPt06_ph = (*pfNewIsoPt06_ph)[muref];
    b_muon_pfNewIsoPt06_pu = muon->pfIsolationR03().sumPUPt;
    b_muon_pfNewIsoPt06    = ( b_muon_pfNewIsoPt06_ch + max(0.0, b_muon_pfNewIsoPt06_nh + b_muon_pfNewIsoPt06_ph - 0.5 * b_muon_pfNewIsoPt06_pu) ) / muon->pt();
    
    b_muon_puppiNewIsoPt08_ch = (*puppiNewIsoPt08_ch)[muref];
    b_muon_puppiNewIsoPt08_nh = (*puppiNewIsoPt08_nh)[muref];
    b_muon_puppiNewIsoPt08_ph = (*puppiNewIsoPt08_ph)[muref];
    b_muon_puppiNewIsoPt08    = ( b_muon_puppiNewIsoPt08_ch + b_muon_puppiNewIsoPt08_nh + b_muon_puppiNewIsoPt08_ph )/muon->pt();    
    b_muon_pfNewIsoPt08_ch = (*pfNewIsoPt08_ch)[muref];
    b_muon_pfNewIsoPt08_nh = (*pfNewIsoPt08_nh)[muref];
    b_muon_pfNewIsoPt08_ph = (*pfNewIsoPt08_ph)[muref];
    b_muon_pfNewIsoPt08_pu = muon->pfIsolationR03().sumPUPt;
    b_muon_pfNewIsoPt08    = ( b_muon_pfNewIsoPt08_ch + max(0.0, b_muon_pfNewIsoPt08_nh + b_muon_pfNewIsoPt08_ph - 0.5 * b_muon_pfNewIsoPt08_pu) ) / muon->pt();
    
    b_muon_puppiNewIsoPt10_ch = (*puppiNewIsoPt10_ch)[muref];
    b_muon_puppiNewIsoPt10_nh = (*puppiNewIsoPt10_nh)[muref];
    b_muon_puppiNewIsoPt10_ph = (*puppiNewIsoPt10_ph)[muref];
    b_muon_puppiNewIsoPt10    = ( b_muon_puppiNewIsoPt10_ch + b_muon_puppiNewIsoPt10_nh + b_muon_puppiNewIsoPt10_ph )/muon->pt();    
    b_muon_pfNewIsoPt10_ch = (*pfNewIsoPt10_ch)[muref];
    b_muon_pfNewIsoPt10_nh = (*pfNewIsoPt10_nh)[muref];
    b_muon_pfNewIsoPt10_ph = (*pfNewIsoPt10_ph)[muref];
    b_muon_pfNewIsoPt10_pu = muon->pfIsolationR03().sumPUPt;
    b_muon_pfNewIsoPt10    = ( b_muon_pfNewIsoPt10_ch + max(0.0, b_muon_pfNewIsoPt10_nh + b_muon_pfNewIsoPt10_ph - 0.5 * b_muon_pfNewIsoPt10_pu) ) / muon->pt();
    
    bool ipxy = false, ipz = false, validPxlHit = false, highPurity = false;
    if (muon->innerTrack().isNonnull()){
      ipxy = std::abs(muon->muonBestTrack()->dxy(priVertex_.position())) < 0.2;
      ipz = std::abs(muon->muonBestTrack()->dz((priVertex_.position()))) < 0.5;
      validPxlHit = muon->innerTrack()->hitPattern().numberOfValidPixelHits() > 0;
      highPurity = muon->innerTrack()->quality(reco::Track::highPurity);
    }
    // isMediumME0 - just loose with track requirements for now, this needs to be updated
    b_muon_isME0MuonMedium = isME0MuonSelNew(*muon, 0.077, dPhiCut_, dPhiBendCut_) && ipxy && validPxlHit && highPurity;

    // tighter cuts for tight ME0
    dPhiCut_ = std::min(std::max(1.2/mom,1.2/100),0.032);
    dPhiBendCut_ = std::min(std::max(0.2/mom,0.2/100),0.0041);
    b_muon_isME0MuonTight = isME0MuonSelNew(*muon, 0.048, dPhiCut_, dPhiBendCut_) && ipxy && ipz && validPxlHit && highPurity;
  }
  tree->Fill();
}

bool PatMuonAnalyser::isME0MuonSelNew(reco::Muon muon, double dEtaCut, double dPhiCut, double dPhiBendCut)
{
  bool result = false;
  bool isME0 = muon.isME0Muon();
    
  if(isME0){
      
    double deltaEta = 999;
    double deltaPhi = 999;
    double deltaPhiBend = 999;

    const std::vector<reco::MuonChamberMatch>& chambers = muon.matches();
    for( std::vector<reco::MuonChamberMatch>::const_iterator chamber = chambers.begin(); chamber != chambers.end(); ++chamber ){
        
      if (chamber->detector() == 5){
          
	for ( std::vector<reco::MuonSegmentMatch>::const_iterator segment = chamber->me0Matches.begin(); segment != chamber->me0Matches.end(); ++segment ){

	  LocalPoint trk_loc_coord(chamber->x, chamber->y, 0);
	  LocalPoint seg_loc_coord(segment->x, segment->y, 0);
	  LocalVector trk_loc_vec(chamber->dXdZ, chamber->dYdZ, 1);
	  LocalVector seg_loc_vec(segment->dXdZ, segment->dYdZ, 1);
            
	  const ME0Chamber * me0chamber = ME0Geometry_->chamber(chamber->id);
            
	  GlobalPoint trk_glb_coord = me0chamber->toGlobal(trk_loc_coord);
	  GlobalPoint seg_glb_coord = me0chamber->toGlobal(seg_loc_coord);
            
	  //double segDPhi = segment->me0SegmentRef->deltaPhi();
	  // need to check if this works
	  double segDPhi = me0chamber->computeDeltaPhi(seg_loc_coord, seg_loc_vec);
	  double trackDPhi = me0chamber->computeDeltaPhi(trk_loc_coord, trk_loc_vec);
            
	  deltaEta = std::abs(trk_glb_coord.eta() - seg_glb_coord.eta() );
	  deltaPhi = std::abs(trk_glb_coord.phi() - seg_glb_coord.phi() );
	  deltaPhiBend = std::abs(segDPhi - trackDPhi);
            
	  if (deltaEta < dEtaCut && deltaPhi < dPhiCut && deltaPhiBend < dPhiBendCut) result = true;
            
	}
      }
    }
      
  }
    
  return result;
    
}

void PatMuonAnalyser::setBranches(TTree *tree)
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
  tree->Branch("muon_isTight", &b_muon_isTight, "muon_isTight/O");
  tree->Branch("muon_isMedium", &b_muon_isMedium, "muon_isMedium/O");
  tree->Branch("muon_isLoose", &b_muon_isLoose, "muon_isLoose/O");
  tree->Branch("muon_isME0MuonTight", &b_muon_isME0MuonTight, "muon_isME0MuonTight/O");
  tree->Branch("muon_isME0MuonMedium", &b_muon_isME0MuonMedium, "muon_isME0MuonMedium/O");
  tree->Branch("muon_isME0MuonLoose", &b_muon_isME0MuonLoose, "muon_isME0MuonLoose/O");

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
  tree->Branch("muon_pfNewIso_ch",&b_muon_pfNewIso_ch,"muon_pfNewIso_ch/F");
  tree->Branch("muon_pfNewIso_nh",&b_muon_pfNewIso_nh,"muon_pfNewIso_nh/F");
  tree->Branch("muon_pfNewIso_ph",&b_muon_pfNewIso_ph,"muon_pfNewIso_ph/F");
  tree->Branch("muon_pfNewIso_pu",&b_muon_pfNewIso_pu,"muon_pfNewIso_pu/F");
  tree->Branch("muon_pfNewIso",&b_muon_pfNewIso,"muon_pfNewIso/F");

  tree->Branch("muon_puppiNewIsoPt02_ch",&b_muon_puppiNewIsoPt02_ch,"muon_puppiNewIsoPt02_ch/F");
  tree->Branch("muon_puppiNewIsoPt02_nh",&b_muon_puppiNewIsoPt02_nh,"muon_puppiNewIsoPt02_nh/F");
  tree->Branch("muon_puppiNewIsoPt02_ph",&b_muon_puppiNewIsoPt02_ph,"muon_puppiNewIsoPt02_ph/F");
  tree->Branch("muon_puppiNewIsoPt02",&b_muon_puppiNewIsoPt02,"muon_puppiNewIsoPt02/F");
  tree->Branch("muon_pfNewIsoPt02_ch",&b_muon_pfNewIsoPt02_ch,"muon_pfNewIsoPt02_ch/F");
  tree->Branch("muon_pfNewIsoPt02_nh",&b_muon_pfNewIsoPt02_nh,"muon_pfNewIsoPt02_nh/F");
  tree->Branch("muon_pfNewIsoPt02_ph",&b_muon_pfNewIsoPt02_ph,"muon_pfNewIsoPt02_ph/F");
  tree->Branch("muon_pfNewIsoPt02_pu",&b_muon_pfNewIsoPt02_pu,"muon_pfNewIsoPt02_pu/F");
  tree->Branch("muon_pfNewIsoPt02",&b_muon_pfNewIsoPt02,"muon_pfNewIsoPt02/F");
  
  tree->Branch("muon_puppiNewIsoPt04_ch",&b_muon_puppiNewIsoPt04_ch,"muon_puppiNewIsoPt04_ch/F");
  tree->Branch("muon_puppiNewIsoPt04_nh",&b_muon_puppiNewIsoPt04_nh,"muon_puppiNewIsoPt04_nh/F");
  tree->Branch("muon_puppiNewIsoPt04_ph",&b_muon_puppiNewIsoPt04_ph,"muon_puppiNewIsoPt04_ph/F");
  tree->Branch("muon_puppiNewIsoPt04",&b_muon_puppiNewIsoPt04,"muon_puppiNewIsoPt04/F");
  tree->Branch("muon_pfNewIsoPt04_ch",&b_muon_pfNewIsoPt04_ch,"muon_pfNewIsoPt04_ch/F");
  tree->Branch("muon_pfNewIsoPt04_nh",&b_muon_pfNewIsoPt04_nh,"muon_pfNewIsoPt04_nh/F");
  tree->Branch("muon_pfNewIsoPt04_ph",&b_muon_pfNewIsoPt04_ph,"muon_pfNewIsoPt04_ph/F");
  tree->Branch("muon_pfNewIsoPt04_pu",&b_muon_pfNewIsoPt04_pu,"muon_pfNewIsoPt04_pu/F");
  tree->Branch("muon_pfNewIsoPt04",&b_muon_pfNewIsoPt04,"muon_pfNewIsoPt04/F");
  
  tree->Branch("muon_puppiNewIsoPt06_ch",&b_muon_puppiNewIsoPt06_ch,"muon_puppiNewIsoPt06_ch/F");
  tree->Branch("muon_puppiNewIsoPt06_nh",&b_muon_puppiNewIsoPt06_nh,"muon_puppiNewIsoPt06_nh/F");
  tree->Branch("muon_puppiNewIsoPt06_ph",&b_muon_puppiNewIsoPt06_ph,"muon_puppiNewIsoPt06_ph/F");
  tree->Branch("muon_puppiNewIsoPt06",&b_muon_puppiNewIsoPt06,"muon_puppiNewIsoPt06/F");
  tree->Branch("muon_pfNewIsoPt06_ch",&b_muon_pfNewIsoPt06_ch,"muon_pfNewIsoPt06_ch/F");
  tree->Branch("muon_pfNewIsoPt06_nh",&b_muon_pfNewIsoPt06_nh,"muon_pfNewIsoPt06_nh/F");
  tree->Branch("muon_pfNewIsoPt06_ph",&b_muon_pfNewIsoPt06_ph,"muon_pfNewIsoPt06_ph/F");
  tree->Branch("muon_pfNewIsoPt06_pu",&b_muon_pfNewIsoPt06_pu,"muon_pfNewIsoPt06_pu/F");
  tree->Branch("muon_pfNewIsoPt06",&b_muon_pfNewIsoPt06,"muon_pfNewIsoPt06/F");
  
  tree->Branch("muon_puppiNewIsoPt08_ch",&b_muon_puppiNewIsoPt08_ch,"muon_puppiNewIsoPt08_ch/F");
  tree->Branch("muon_puppiNewIsoPt08_nh",&b_muon_puppiNewIsoPt08_nh,"muon_puppiNewIsoPt08_nh/F");
  tree->Branch("muon_puppiNewIsoPt08_ph",&b_muon_puppiNewIsoPt08_ph,"muon_puppiNewIsoPt08_ph/F");
  tree->Branch("muon_puppiNewIsoPt08",&b_muon_puppiNewIsoPt08,"muon_puppiNewIsoPt08/F");
  tree->Branch("muon_pfNewIsoPt08_ch",&b_muon_pfNewIsoPt08_ch,"muon_pfNewIsoPt08_ch/F");
  tree->Branch("muon_pfNewIsoPt08_nh",&b_muon_pfNewIsoPt08_nh,"muon_pfNewIsoPt08_nh/F");
  tree->Branch("muon_pfNewIsoPt08_ph",&b_muon_pfNewIsoPt08_ph,"muon_pfNewIsoPt08_ph/F");
  tree->Branch("muon_pfNewIsoPt08_pu",&b_muon_pfNewIsoPt08_pu,"muon_pfNewIsoPt08_pu/F");
  tree->Branch("muon_pfNewIsoPt08",&b_muon_pfNewIsoPt08,"muon_pfNewIsoPt08/F");
  
  tree->Branch("muon_puppiNewIsoPt10_ch",&b_muon_puppiNewIsoPt10_ch,"muon_puppiNewIsoPt10_ch/F");
  tree->Branch("muon_puppiNewIsoPt10_nh",&b_muon_puppiNewIsoPt10_nh,"muon_puppiNewIsoPt10_nh/F");
  tree->Branch("muon_puppiNewIsoPt10_ph",&b_muon_puppiNewIsoPt10_ph,"muon_puppiNewIsoPt10_ph/F");
  tree->Branch("muon_puppiNewIsoPt10",&b_muon_puppiNewIsoPt10,"muon_puppiNewIsoPt10/F");
  tree->Branch("muon_pfNewIsoPt10_ch",&b_muon_pfNewIsoPt10_ch,"muon_pfNewIsoPt10_ch/F");
  tree->Branch("muon_pfNewIsoPt10_nh",&b_muon_pfNewIsoPt10_nh,"muon_pfNewIsoPt10_nh/F");
  tree->Branch("muon_pfNewIsoPt10_ph",&b_muon_pfNewIsoPt10_ph,"muon_pfNewIsoPt10_ph/F");
  tree->Branch("muon_pfNewIsoPt10_pu",&b_muon_pfNewIsoPt10_pu,"muon_pfNewIsoPt10_pu/F");
  tree->Branch("muon_pfNewIsoPt10",&b_muon_pfNewIsoPt10,"muon_pfNewIsoPt10/F");
  
  tree->Branch("muon_PFIsoFixOnlyCH",&b_muon_PFIsoFixOnlyCH,"muon_PFIsoFixOnlyCH/F");
  tree->Branch("muon_puppiIsoFixOnlyCH",&b_muon_puppiIsoFixOnlyCH,"muon_puppiIsoFixOnlyCH/F");
  tree->Branch("muon_PFIsoRepTrk",&b_muon_PFIsoRepTrk,"muon_PFIsoRepTrk/F");
  tree->Branch("muon_puppiIsoRepTrk",&b_muon_puppiIsoRepTrk,"muon_puppiIsoRepTrk/F");
}

void
PatMuonAnalyser::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PatMuonAnalyser);
