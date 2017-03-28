#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonTime.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "RecoParticleFlow/PFProducer/interface/PFMuonAlgo.h"

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/Associations/interface/MuonToTrackingParticleAssociator.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimMuon/MCTruth/interface/MuonToSimAssociatorByHits.h"
#include "SimTracker/Common/interface/TrackingParticleSelector.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"

#include "MuonPerformance/MuonAnalyser/src/TMVAClassification_BDT.class.C"
#include "MuonPerformance/MuonAnalyser/src/TMVAClassification_MLP.class.C"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

#include <TTree.h>
#include <TFile.h>
#include <TLorentzVector.h>

using namespace std;
using namespace reco;
using namespace edm;

enum particleType{
  CH = 0, 
  NH = 1,
  PH = 2,
  OTHER = 100000
};

class MuonAnalyser : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  struct PUPPIIsolation {
    float charged_hadrons;
    float neutral_hadrons;
    float photons;
    float NoLep_charged_hadrons;
    float NoLep_neutral_hadrons;
    float NoLep_photons;
  };
  
  struct puppiIso {
    int nNumCands;
    int nNumCandsInR05;
    int nNumCandsOutR05;
    
    int nNumCandsInR05OT;
    
    int nNumCandsInR05CH;
    int nNumCandsInR05NH;
    int nNumCandsInR05PH;
    
    int nNumCandsInR05CHApp;
    int nNumCandsInR05NHApp;
    int nNumCandsInR05PHApp;
    
    double combined;
    double withLep;
    double withoutlep;
    
    double combined03;
    double withLep03;
    double withoutlep03;
    
    double combined05;
    double withLep05;
    double withoutlep05;
    
    double withLep03CH;
    double withLep03NH;
    double withLep03PH;
    
    double withLep04CH;
    double withLep04NH;
    double withLep04PH;
    
    double withLep05CH;
    double withLep05NH;
    double withLep05PH;
    
    double withoutlep03CH;
    double withoutlep03NH;
    double withoutlep03PH;
    
    double withoutlep04CH;
    double withoutlep04NH;
    double withoutlep04PH;
    
    double withoutlep05CH;
    double withoutlep05NH;
    double withoutlep05PH;
  };

public:
  explicit MuonAnalyser(const edm::ParameterSet&);
  ~MuonAnalyser();
  
  bool isGlobalTightMuon( const Muon* muonRef );
  bool isTrackerTightMuon( const reco::Muon *muonRef );
  bool isIsolatedMuon( const reco::Muon *muonRef );

  bool isLooseMuonCustom(const reco::Muon& mu) const;
  bool isTightMuonCustom(const reco::Muon& mu, reco::Vertex pv0) const;
  bool isTightMuonCustomOptimized(const reco::Muon& mu, reco::Vertex pv0) const;
  
  bool isLooseMod(const reco::Muon *muon);
  bool isTightMod(const reco::VertexCollection* vertices, const SimVertex &simPVh, const reco::Muon *muon, bool useIPxy, bool useIPz);
  
  std::vector<double> collectTMVAvalues(const reco::Muon& mu, reco::Vertex pv0) const;
  int nGEMhit(const reco::Muon * mu) const;
  int nME0hit(const reco::Muon * mu) const;
  
  puppiIso getPuppiIso(const reco::Muon *mu, const vector< pat::PackedCandidate> *pcs, int nIsReco) const;
  bool isNH( long pdgid ) const;
  bool isCH( long pdgid ) const;
  bool isPH( long pdgid ) const;

  void setBranches(TTree *tree);
  void fillBranches(TTree *tree, TLorentzVector tlv, edm::RefToBase<reco::Muon> muref, bool isSignal, int pdgId);
  
private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

  TH1D* h_nevents;
  TH1D* h_vertex;
  
  TH1D* h_vertexBS;
  TH1D* h_vertex1D;
  TH1D* h_vertex1DBS;
  TH1D* h_vertex4D;
  TH1D* h_vertex4DBS;
  
  int b_nvertex;

  TTree* genttree_;
  TTree* recottree_;

  TLorentzVector b_muon;
  bool b_muon_signal;
  int b_muon_pdgId;
  float b_muon_pTresolution, b_muon_pTinvresolution;
  bool b_muon_isTightOptimized, b_muon_isTightCustom, b_muon_isTight, b_muon_isMedium, b_muon_isLoose;
  bool b_muon_isME0Muon, b_muon_isGEMMuon, b_muon_isRPCMuon, b_muon_isCaloMuon, b_muon_isTrackerMuon;
  bool b_muon_isGlobalMuon, b_muon_isStandAloneMuon, b_muon_isPFMuon;
  bool b_muon_isLooseMod;
  bool b_muon_isTightModNoIP, b_muon_isTightModIPxy, b_muon_isTightModIPz, b_muon_isTightModIPxyz;

  bool b_muon_global; bool b_muon_pf;
  float b_muon_chi2pos; float b_muon_trkKink; float b_muon_segcompati;
  float b_muon_chi2; int b_muon_nglobalhits; int b_muon_nstations;
  float b_muon_trackdxy; float b_muon_trackdz;
  int b_muon_ninnerhits; float b_muon_trackerlayers;
  float b_muon_poszPV0, b_muon_poszSimPV, b_muon_poszMuon;

  float b_muon_ME0deltaX, b_muon_ME0deltaY, b_muon_ME0deltaDXDZ, b_muon_ME0deltaDYDZ, b_muon_ME0pullX, b_muon_ME0pullY, b_muon_ME0dPhi; int b_muon_ME0noRecHit;
  float b_muon_GE11deltaX, b_muon_GE11deltaY, b_muon_GE11deltaDXDZ, b_muon_GE11deltaDYDZ, b_muon_GE11pullX, b_muon_GE11pullY, b_muon_GE11dPhi; int b_muon_GE11noRecHit;
  float b_muon_GE21deltaX, b_muon_GE21deltaY, b_muon_GE21deltaDXDZ, b_muon_GE21deltaDYDZ, b_muon_GE21pullX, b_muon_GE21pullY, b_muon_GE21dPhi; int b_muon_GE21noRecHit;
  
  float b_muon_PFIso04; float b_muon_PFIso03;
  float b_muon_PFIso03ChargedHadronPt, b_muon_PFIso03NeutralHadronEt;
  float b_muon_PFIso03PhotonEt, b_muon_PFIso03PUPt;
  float b_muon_PFIso04ChargedHadronPt, b_muon_PFIso04NeutralHadronEt;
  float b_muon_PFIso04PhotonEt, b_muon_PFIso04PUPt;
  float b_muon_TrkIso05; float b_muon_TrkIso03;
  float b_muon_puppiIso, b_muon_puppiIsoNoLep;
  float b_muon_puppiIso_ChargedHadron, b_muon_puppiIso_NeutralHadron, b_muon_puppiIso_Photon;  
  float b_muon_puppiIsoNoLep_ChargedHadron, b_muon_puppiIsoNoLep_NeutralHadron, b_muon_puppiIsoNoLep_Photon;  
  float b_muon_puppiIsoWithLep, b_muon_puppiIsoWithoutLep, b_muon_puppiIsoCombined;
  float b_muon_puppiIsoWithLep03, b_muon_puppiIsoWithoutLep03, b_muon_puppiIsoCombined03;
  float b_muon_puppiIsoWithLep05, b_muon_puppiIsoWithoutLep05, b_muon_puppiIsoCombined05;
  float b_muon_puppiIsoWithLep03ChargedHadron, b_muon_puppiIsoWithLep03NeutralHadron, b_muon_puppiIsoWithLep03Photon;
  float b_muon_puppiIsoWithLep04ChargedHadron, b_muon_puppiIsoWithLep04NeutralHadron, b_muon_puppiIsoWithLep04Photon;
  float b_muon_puppiIsoWithLep05ChargedHadron, b_muon_puppiIsoWithLep05NeutralHadron, b_muon_puppiIsoWithLep05Photon;
  float b_muon_puppiIsoWithoutLep03ChargedHadron, b_muon_puppiIsoWithoutLep03NeutralHadron, b_muon_puppiIsoWithoutLep03Photon;
  float b_muon_puppiIsoWithoutLep04ChargedHadron, b_muon_puppiIsoWithoutLep04NeutralHadron, b_muon_puppiIsoWithoutLep04Photon;
  float b_muon_puppiIsoWithoutLep05ChargedHadron, b_muon_puppiIsoWithoutLep05NeutralHadron, b_muon_puppiIsoWithoutLep05Photon;
  int b_muon_puppiIsoNumOfCands;
  int b_muon_puppiIsoNumOfCandsInR05;
  int b_muon_puppiIsoNumOfCandsOutR05;
  int b_muon_puppiIsoNumOfCandsInR05OT;
  int b_muon_puppiIsoNumOfCandsInR05CH, b_muon_puppiIsoNumOfCandsInR05CHApp;
  int b_muon_puppiIsoNumOfCandsInR05NH, b_muon_puppiIsoNumOfCandsInR05NHApp;
  int b_muon_puppiIsoNumOfCandsInR05PH, b_muon_puppiIsoNumOfCandsInR05PHApp;
  bool b_muon_isMuon;
  int b_muon_numberOfValidMuonGEMHits, b_muon_numberOfValidMuonME0Hits;

  float b_muon_tmva_bdt, b_muon_tmva_mlp;
  
  float b_tmva_bdt; float b_tmva_mlp;
  ReadBDT* bdt_;
  ReadMLP* mlp_;

  edm::Handle<edm::ValueMap<float>> PUPPIIsolation_charged_hadrons;
  edm::Handle<edm::ValueMap<float>> PUPPIIsolation_neutral_hadrons;
  edm::Handle<edm::ValueMap<float>> PUPPIIsolation_photons;
  edm::Handle<edm::ValueMap<float>> PUPPINoLeptonsIsolation_charged_hadrons;
  edm::Handle<edm::ValueMap<float>> PUPPINoLeptonsIsolation_neutral_hadrons;
  edm::Handle<edm::ValueMap<float>> PUPPINoLeptonsIsolation_photons;
  const reco::VertexCollection* vertexes_;
  SimVertex simVertex_;
  const std::vector<pat::PackedCandidate> * candidates_;

  edm::EDGetTokenT<std::vector<reco::Vertex> > vtxToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex> > vtx1DToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex> > vtx1DBSToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex> > vtx4DToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex> > vtx4DBSToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex> > vtxBSToken_;
  edm::EDGetTokenT<TrackingParticleCollection> simToken_;
  edm::EDGetTokenT<std::vector<SimVertex> > simVertexToken_;
  edm::EDGetTokenT<edm::View<reco::Muon> > muonToken_;
  edm::EDGetTokenT<reco::MuonToTrackingParticleAssociator> muAssocToken_;
  edm::EDGetTokenT <std::vector< pat::PackedCandidate> > tokenPackedCandidate ;
  edm::EDGetTokenT<edm::ValueMap<float> > PUPPIIsolation_charged_hadrons_;
  edm::EDGetTokenT<edm::ValueMap<float> > PUPPIIsolation_neutral_hadrons_;
  edm::EDGetTokenT<edm::ValueMap<float> > PUPPIIsolation_photons_;
  edm::EDGetTokenT<edm::ValueMap<float> > PUPPINoLeptonsIsolation_charged_hadrons_;
  edm::EDGetTokenT<edm::ValueMap<float> > PUPPINoLeptonsIsolation_neutral_hadrons_;
  edm::EDGetTokenT<edm::ValueMap<float> > PUPPINoLeptonsIsolation_photons_;

  TrackingParticleSelector tpSelector_;
};

MuonAnalyser::MuonAnalyser(const edm::ParameterSet& pset)
{
  vtxToken_     = consumes<vector<Vertex> >(pset.getParameter<edm::InputTag>("primaryVertex"));
  vtx1DToken_   = consumes<vector<Vertex> >(pset.getParameter<edm::InputTag>("primaryVertex1D"));
  vtx1DBSToken_ = consumes<vector<Vertex> >(pset.getParameter<edm::InputTag>("primaryVertex1DBS"));
  vtx4DToken_   = consumes<vector<Vertex> >(pset.getParameter<edm::InputTag>("primaryVertex4D"));
  vtx4DBSToken_ = consumes<vector<Vertex> >(pset.getParameter<edm::InputTag>("primaryVertex4DBS"));
  vtxBSToken_   = consumes<vector<Vertex> >(pset.getParameter<edm::InputTag>("primaryVertexBS"));
  
  simToken_ = consumes<TrackingParticleCollection>(pset.getParameter<InputTag>("simLabel"));
  simVertexToken_ = consumes<std::vector<SimVertex> >(pset.getParameter<edm::InputTag> ("simVertexCollection"));  
  muonToken_ = consumes<View<Muon> >(pset.getParameter<InputTag>("muonLabel"));
  muAssocToken_ = consumes<reco::MuonToTrackingParticleAssociator>(pset.getParameter<InputTag>("muAssocLabel"));
  tokenPackedCandidate = consumes <std::vector< pat::PackedCandidate> > ( edm::InputTag( std::string("packedPFCandidates"), std::string(""),std::string("") ) );

  PUPPIIsolation_charged_hadrons_ = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("puppiIsolationChargedHadrons"));
  PUPPIIsolation_neutral_hadrons_ = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("puppiIsolationNeutralHadrons"));
  PUPPIIsolation_photons_ = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("puppiIsolationPhotons"));
  PUPPINoLeptonsIsolation_charged_hadrons_ = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("puppiNoLepIsolationChargedHadrons"));
  PUPPINoLeptonsIsolation_neutral_hadrons_ = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("puppiNoLepIsolationNeutralHadrons"));
  PUPPINoLeptonsIsolation_photons_ = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("puppiNoLepIsolationPhotons"));
  
  ParameterSet tpset = pset.getParameter<ParameterSet>("tpSelector");
  tpSelector_ = TrackingParticleSelector(tpset.getParameter<double>("ptMin"),
                                         tpset.getParameter<double>("minRapidity"),
                                         tpset.getParameter<double>("maxRapidity"),
                                         tpset.getParameter<double>("tip"),
                                         tpset.getParameter<double>("lip"),
                                         tpset.getParameter<int>("minHit"),
                                         tpset.getParameter<bool>("signalOnly"),
                                         tpset.getParameter<bool>("intimeOnly"),
                                         tpset.getParameter<bool>("chargedOnly"),
                                         tpset.getParameter<bool>("stableOnly"),
                                         tpset.getParameter<std::vector<int> >("pdgId"));
  
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  h_nevents = fs->make<TH1D>("nevents", "nevents", 1, 0, 1);
  h_vertex = fs->make<TH1D>("vertex reco vs sim", "vertex reco vs sim", 100, -1, 1);
  
  h_vertexBS   = fs->make<TH1D>("vertex reco with BS vs sim", "vertex reco with BS vs sim", 100, -1, 1);
  h_vertex1D   = fs->make<TH1D>("vertex reco 1D vs sim", "vertex reco 1D vs sim", 100, -1, 1);
  h_vertex1DBS = fs->make<TH1D>("vertex reco 1D with BS vs sim", "vertex reco 1D with BS vs sim", 100, -1, 1);
  h_vertex4D   = fs->make<TH1D>("vertex reco 4D vs sim", "vertex reco 4D vs sim", 100, -1, 1);
  h_vertex4DBS = fs->make<TH1D>("vertex reco 4D with BS vs sim", "vertex reco 4D with BS vs sim", 100, -1, 1);

  genttree_ = fs->make<TTree>("gen", "gen");
  setBranches(genttree_);
  recottree_ = fs->make<TTree>("reco", "reco");
  setBranches(recottree_);

  string dummy[] = { "recoMuon_isGlobalMuon", "recoMuon_isPFMuon", "recoMuon_normalizedChi2", "recoMuon_chi2LocalPosition", "recoMuon_trkKink", "recoMuon_segmentCompatibility", "recoMuon_numberOfMatchedStations", "recoMuon_numberOfValidMuonHits", "recoMuon_pv0pos_dxy", "recoMuon_pv0pos_dz", "recoMuon_numberOfValidPixelHits", "recoMuon_trackerLayersWithMeasurement" };
  vector< string > dummy_label;
  dummy_label.assign(dummy, dummy+12);
  bdt_ = new ReadBDT(dummy_label);
  mlp_ = new ReadMLP(dummy_label);

}
MuonAnalyser::~MuonAnalyser(){}
void MuonAnalyser::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  h_nevents->Fill(0.5);

  Handle<VertexCollection> vertices;
  Handle<VertexCollection> vertices1D, vertices1DBS, vertices4D, vertices4DBS, verticesBS;
  
  iEvent.getByToken(vtxToken_, vertices); 
  iEvent.getByToken(vtx1DToken_,   vertices1D); 
  iEvent.getByToken(vtx1DBSToken_, vertices1DBS); 
  iEvent.getByToken(vtx4DToken_,   vertices4D); 
  iEvent.getByToken(vtx4DBSToken_, vertices4DBS); 
  iEvent.getByToken(vtxBSToken_,   verticesBS); 
  
  vertexes_ = vertices.product();
  if (vertexes_->empty()) { cout << "no PV" << endl; return; }

  Handle<std::vector<SimVertex> > simVertexCollection;
  iEvent.getByToken(simVertexToken_, simVertexCollection);  
  simVertex_ = simVertexCollection->at(0);
  
  h_vertex->Fill(vertexes_->at(0).position().z() - simVertex_.position().z());  
  h_vertexBS->Fill(  verticesBS->front().position().z()   - simVertex_.position().z());
  h_vertex1D->Fill(  vertices1D->front().position().z()   - simVertex_.position().z());
  h_vertex1DBS->Fill(vertices1DBS->front().position().z() - simVertex_.position().z());
  h_vertex4D->Fill(  vertices4D->front().position().z()   - simVertex_.position().z());
  h_vertex4DBS->Fill(vertices4DBS->front().position().z() - simVertex_.position().z());
  
  Handle<TrackingParticleCollection> simHandle;
  iEvent.getByToken(simToken_, simHandle);

  Handle<View<Muon> > muonHandle;
  iEvent.getByToken(muonToken_, muonHandle);

  edm::Handle<std::vector< pat::PackedCandidate> > Candidates_Collection ;
  iEvent.getByToken( tokenPackedCandidate , Candidates_Collection ) ; 
  candidates_ = Candidates_Collection.product();
  if (candidates_->empty()) { cout << "no PackedCandidate :: something is wrong" << endl;}

  Handle<MuonToTrackingParticleAssociator> associatorBase;
  iEvent.getByToken(muAssocToken_, associatorBase);
  MuonToTrackingParticleAssociator const* assoByHits = associatorBase.product();
  reco::MuonToSimCollection muonToSimColl; reco::SimToMuonCollection simToMuonColl;  
  assoByHits->associateMuons(muonToSimColl, simToMuonColl, muonHandle, reco::GlobalTk, simHandle);
  reco::MuonToSimCollection trkToSimColl; reco::SimToMuonCollection simToTrkColl;  
  assoByHits->associateMuons(trkToSimColl, simToTrkColl, muonHandle, reco::InnerTk, simHandle);

  iEvent.getByToken(PUPPIIsolation_charged_hadrons_, PUPPIIsolation_charged_hadrons);
  iEvent.getByToken(PUPPIIsolation_neutral_hadrons_, PUPPIIsolation_neutral_hadrons);
  iEvent.getByToken(PUPPIIsolation_photons_, PUPPIIsolation_photons);  
  iEvent.getByToken(PUPPINoLeptonsIsolation_charged_hadrons_, PUPPINoLeptonsIsolation_charged_hadrons);
  iEvent.getByToken(PUPPINoLeptonsIsolation_neutral_hadrons_, PUPPINoLeptonsIsolation_neutral_hadrons);
  iEvent.getByToken(PUPPINoLeptonsIsolation_photons_, PUPPINoLeptonsIsolation_photons);  
  
  // gen muon loop for efficinecy 
  for (TrackingParticleCollection::size_type i=0; i<simHandle->size(); i++) {
    TrackingParticleRef simRef(simHandle, i);
    const TrackingParticle* simTP = simRef.get();
    // select good gen muons
    //if ( ! tpSelector_(*simTP) ) continue;
    bool isSignalMuon = abs(simTP->pdgId())==13 && !simTP->genParticles().empty() && (simTP->eventId().event() == 0) && (simTP->eventId().bunchCrossing() == 0);
    if (!isSignalMuon) continue;
    
    TLorentzVector gentlv(simRef->momentum().x(), simRef->momentum().y(), simRef->momentum().z(), simRef->energy() );
    edm::RefToBase<reco::Muon> muonRef;
    if ( simToMuonColl.find(simRef) != simToMuonColl.end() ) {
      vector<pair<RefToBase<Muon>, double> > MuRefV = simToMuonColl[simRef];      
      if ( !MuRefV.empty()) {
    	muonRef = MuRefV.begin()->first;
      }
    }
    else if ( simToTrkColl.find(simRef) != simToTrkColl.end() ) {
      vector<pair<RefToBase<Muon>, double> > trRefV = simToTrkColl[simRef];      
      if ( !trRefV.empty()) {
	edm::RefToBase<reco::Muon> trkMu = trRefV.begin()->first;
	const Muon* mu = trkMu.get();
	if (!mu->isGlobalMuon())
	  muonRef = trkMu;
      }
    }
    
    fillBranches(genttree_, gentlv, muonRef, true, simTP->pdgId());
  }

  
  // reco muon loop - for fake rate
  for (size_t i = 0; i < muonHandle->size(); ++i) {    
    edm::RefToBase<reco::Muon> muRef = muonHandle->refAt(i);
    const Muon* mu = muRef.get();
    TLorentzVector recotlv(mu->momentum().x(), mu->momentum().y(), mu->momentum().z(), mu->energy() );
    
    int pdgId = 0;
    const TrackingParticle* simTP = NULL;
    if (mu->isGlobalMuon()){
      if ( muonToSimColl.find(muRef) != muonToSimColl.end() ) {
	auto trkRefV = muonToSimColl[muRef];
	if ( !trkRefV.empty()) {
	  simTP = trkRefV.begin()->first.get();
	  pdgId = simTP->pdgId();	
	}
      }
    }
    else {
      if ( trkToSimColl.find(muRef) != trkToSimColl.end() ) {
	auto trkRefV = trkToSimColl[muRef];
	if ( !trkRefV.empty()) {
	  simTP = trkRefV.begin()->first.get();
	  pdgId = simTP->pdgId();	
	}
      }
    }
    bool isSignalMuon = false;
    if (simTP)
      isSignalMuon = abs(simTP->pdgId())==13 && !simTP->genParticles().empty() && (simTP->eventId().event() == 0) && (simTP->eventId().bunchCrossing() == 0);
    fillBranches(recottree_, recotlv, muRef, isSignalMuon, pdgId);
  }

}

void MuonAnalyser::fillBranches(TTree *tree, TLorentzVector tlv, edm::RefToBase<reco::Muon> muref, bool isSignal, int pdgId)
{
  b_muon = tlv;
  b_muon_signal = isSignal;
  b_muon_pdgId = pdgId;
  reco::Vertex pv0 = vertexes_->at(0);
    
  b_muon_pTresolution = 0; b_muon_pTinvresolution = 0;
  
  b_muon_isTightOptimized = 0; b_muon_isTightCustom = 0; b_muon_isTight = 0; b_muon_isMedium = 0; b_muon_isLoose = 0;
  b_muon_isME0Muon = 0; b_muon_isGEMMuon = 0; b_muon_isRPCMuon = 0; b_muon_isCaloMuon = 0; b_muon_isTrackerMuon = 0;
  b_muon_isGlobalMuon = 0; b_muon_isStandAloneMuon = 0; b_muon_isPFMuon = 0;
  b_muon_isLooseMod = 0;
  b_muon_isTightModNoIP = 0; b_muon_isTightModIPxy = 0; b_muon_isTightModIPz = 0; b_muon_isTightModIPxyz = 0;

  b_muon_ME0deltaX = 100; b_muon_ME0deltaY = 0; b_muon_ME0deltaDXDZ = 0; b_muon_ME0deltaDYDZ = 0; b_muon_ME0noRecHit = 0; b_muon_ME0pullX = 0; b_muon_ME0pullY = 0; b_muon_ME0dPhi = 0;
  b_muon_GE11deltaX = 100; b_muon_GE11deltaY = 0; b_muon_GE11deltaDXDZ = 0; b_muon_GE11deltaDYDZ = 0; b_muon_GE11noRecHit = 0; b_muon_GE11pullX = 0; b_muon_GE11pullY = 0; b_muon_GE11dPhi = 0;
  b_muon_GE21deltaX = 100; b_muon_GE21deltaY = 0; b_muon_GE21deltaDXDZ = 0; b_muon_GE21deltaDYDZ = 0; b_muon_GE21noRecHit = 0; b_muon_GE21pullX = 0; b_muon_GE21pullY = 0; b_muon_GE21dPhi = 0;

  b_muon_global = 0;  b_muon_pf = 0;
  b_muon_chi2pos = 0;  b_muon_trkKink = 0;  b_muon_segcompati = 0;
  b_muon_chi2 = 0;  b_muon_nglobalhits = 0;  b_muon_nstations = 0;
  b_muon_trackdxy = 0;  b_muon_trackdz = 0;
  b_muon_ninnerhits = 0;  b_muon_trackerlayers = 0;
  b_muon_poszPV0 = 0; b_muon_poszSimPV = 0; b_muon_poszMuon = 0;
  b_muon_PFIso04 = 0;  b_muon_PFIso03 = 0;
  b_muon_PFIso03ChargedHadronPt = 0; b_muon_PFIso03NeutralHadronEt = 0;
  b_muon_PFIso03PhotonEt = 0; b_muon_PFIso03PUPt = 0;
  b_muon_PFIso04ChargedHadronPt = 0; b_muon_PFIso04NeutralHadronEt = 0;
  b_muon_PFIso04PhotonEt = 0; b_muon_PFIso04PUPt = 0;
  b_muon_TrkIso05 = 0;  b_muon_TrkIso03 = 0;
  b_muon_puppiIso = 0; b_muon_puppiIso_ChargedHadron = 0; b_muon_puppiIso_NeutralHadron = 0; b_muon_puppiIso_Photon = 0;
  b_muon_puppiIsoNoLep = 0; b_muon_puppiIsoNoLep_ChargedHadron = 0; b_muon_puppiIsoNoLep_NeutralHadron = 0; b_muon_puppiIsoNoLep_Photon = 0;  
  b_muon_puppiIsoWithLep = 0; b_muon_puppiIsoWithoutLep = 0; b_muon_puppiIsoCombined = 0;
  b_muon_puppiIsoWithLep03 = 0; b_muon_puppiIsoWithoutLep03 = 0; b_muon_puppiIsoCombined03 = 0;
  b_muon_puppiIsoWithLep05 = 0; b_muon_puppiIsoWithoutLep05 = 0; b_muon_puppiIsoCombined05 = 0;
  b_muon_puppiIsoWithLep03ChargedHadron = 0; b_muon_puppiIsoWithLep03NeutralHadron = 0; b_muon_puppiIsoWithLep03Photon = 0;
  b_muon_puppiIsoWithLep04ChargedHadron = 0; b_muon_puppiIsoWithLep04NeutralHadron = 0; b_muon_puppiIsoWithLep04Photon = 0;
  b_muon_puppiIsoWithLep05ChargedHadron = 0; b_muon_puppiIsoWithLep05NeutralHadron = 0; b_muon_puppiIsoWithLep05Photon = 0;
  b_muon_puppiIsoWithoutLep03ChargedHadron = 0; b_muon_puppiIsoWithoutLep03NeutralHadron = 0; b_muon_puppiIsoWithoutLep03Photon = 0;
  b_muon_puppiIsoWithoutLep04ChargedHadron = 0; b_muon_puppiIsoWithoutLep04NeutralHadron = 0; b_muon_puppiIsoWithoutLep04Photon = 0;
  b_muon_puppiIsoWithoutLep05ChargedHadron = 0; b_muon_puppiIsoWithoutLep05NeutralHadron = 0; b_muon_puppiIsoWithoutLep05Photon = 0;
  b_muon_puppiIsoNumOfCands = 0;
  b_muon_puppiIsoNumOfCandsInR05 = 0;
  b_muon_puppiIsoNumOfCandsOutR05 = 0;
  b_muon_puppiIsoNumOfCandsInR05OT = 0;
  b_muon_puppiIsoNumOfCandsInR05CH = 0; b_muon_puppiIsoNumOfCandsInR05CHApp = 0;
  b_muon_puppiIsoNumOfCandsInR05NH = 0; b_muon_puppiIsoNumOfCandsInR05NHApp = 0;
  b_muon_puppiIsoNumOfCandsInR05PH = 0; b_muon_puppiIsoNumOfCandsInR05PHApp = 0;
  b_muon_isMuon = 0;
  b_muon_numberOfValidMuonGEMHits = 0; b_muon_numberOfValidMuonME0Hits = 0;

  b_muon_tmva_bdt = 0; b_muon_tmva_mlp = 0;
  b_tmva_bdt = 0;  b_tmva_mlp = 0;

  const Muon* mu = muref.get();
  if (mu){
    b_muon_pTresolution = (b_muon.Pt()-mu->pt())/b_muon.Pt();
    b_muon_pTinvresolution = (1/b_muon.Pt() - 1/mu->pt())/(1/b_muon.Pt());
    b_muon_poszPV0       = pv0.position().z();
    b_muon_poszSimPV     = simVertex_.position().z();
    b_muon_poszMuon      = mu->muonBestTrack()->vz();
    
    b_muon_TrkIso03 = mu->isolationR03().sumPt/mu->pt();
    b_muon_TrkIso05 = mu->isolationR05().sumPt/mu->pt();
    b_muon_PFIso04 = (mu->pfIsolationR04().sumChargedHadronPt + TMath::Max(0.,mu->pfIsolationR04().sumNeutralHadronEt + mu->pfIsolationR04().sumPhotonEt - 0.5*mu->pfIsolationR04().sumPUPt))/mu->pt();
    b_muon_PFIso03 = (mu->pfIsolationR03().sumChargedHadronPt + TMath::Max(0.,mu->pfIsolationR03().sumNeutralHadronEt + mu->pfIsolationR03().sumPhotonEt - 0.5*mu->pfIsolationR03().sumPUPt))/mu->pt();
    
    b_muon_PFIso03ChargedHadronPt = mu->pfIsolationR03().sumChargedHadronPt;
    b_muon_PFIso03NeutralHadronEt = mu->pfIsolationR03().sumNeutralHadronEt;
    b_muon_PFIso03PhotonEt        = mu->pfIsolationR03().sumPhotonEt;
    b_muon_PFIso03PUPt            = mu->pfIsolationR03().sumPUPt;

    b_muon_PFIso04ChargedHadronPt = mu->pfIsolationR04().sumChargedHadronPt;
    b_muon_PFIso04NeutralHadronEt = mu->pfIsolationR04().sumNeutralHadronEt;
    b_muon_PFIso04PhotonEt        = mu->pfIsolationR04().sumPhotonEt;
    b_muon_PFIso04PUPt            = mu->pfIsolationR04().sumPUPt;   

    b_muon_puppiIso_ChargedHadron = (*PUPPIIsolation_charged_hadrons)[muref];
    b_muon_puppiIso_NeutralHadron = (*PUPPIIsolation_neutral_hadrons)[muref];
    b_muon_puppiIso_Photon = (*PUPPIIsolation_photons)[muref];
    b_muon_puppiIso = (b_muon_puppiIso_ChargedHadron+b_muon_puppiIso_NeutralHadron+b_muon_puppiIso_Photon)/mu->pt();
    b_muon_puppiIsoNoLep_ChargedHadron = (*PUPPINoLeptonsIsolation_charged_hadrons)[muref];
    b_muon_puppiIsoNoLep_NeutralHadron = (*PUPPINoLeptonsIsolation_neutral_hadrons)[muref];
    b_muon_puppiIsoNoLep_Photon = (*PUPPINoLeptonsIsolation_photons)[muref];
    b_muon_puppiIsoNoLep = (b_muon_puppiIsoNoLep_ChargedHadron+b_muon_puppiIsoNoLep_NeutralHadron+b_muon_puppiIsoNoLep_Photon)/mu->pt();
    
    puppiIso puppiIsoValues = getPuppiIso( mu, candidates_, 0);
    //cout << "pIso.combined "<< pIso.combined  <<endl;

    b_muon_puppiIsoWithLep    = puppiIsoValues.withLep;
    b_muon_puppiIsoWithoutLep = puppiIsoValues.withoutlep;
    b_muon_puppiIsoCombined   = puppiIsoValues.combined;

    b_muon_puppiIsoWithLep03    = puppiIsoValues.withLep03;
    b_muon_puppiIsoWithoutLep03 = puppiIsoValues.withoutlep03;
    b_muon_puppiIsoCombined03   = puppiIsoValues.combined03;

    b_muon_puppiIsoWithLep05    = puppiIsoValues.withLep05;
    b_muon_puppiIsoWithoutLep05 = puppiIsoValues.withoutlep05;
    b_muon_puppiIsoCombined05   = puppiIsoValues.combined05;
    
    b_muon_puppiIsoWithLep03ChargedHadron = puppiIsoValues.withLep03CH;
    b_muon_puppiIsoWithLep03NeutralHadron = puppiIsoValues.withLep03NH;
    b_muon_puppiIsoWithLep03Photon        = puppiIsoValues.withLep03PH;
    
    b_muon_puppiIsoWithLep04ChargedHadron = puppiIsoValues.withLep04CH;
    b_muon_puppiIsoWithLep04NeutralHadron = puppiIsoValues.withLep04NH;
    b_muon_puppiIsoWithLep04Photon        = puppiIsoValues.withLep04PH;
    
    b_muon_puppiIsoWithLep05ChargedHadron = puppiIsoValues.withLep05CH;
    b_muon_puppiIsoWithLep05NeutralHadron = puppiIsoValues.withLep05NH;
    b_muon_puppiIsoWithLep05Photon        = puppiIsoValues.withLep05PH;
    
    b_muon_puppiIsoWithoutLep03ChargedHadron = puppiIsoValues.withoutlep03CH;
    b_muon_puppiIsoWithoutLep03NeutralHadron = puppiIsoValues.withoutlep03NH;
    b_muon_puppiIsoWithoutLep03Photon        = puppiIsoValues.withoutlep03PH;
    
    b_muon_puppiIsoWithoutLep04ChargedHadron = puppiIsoValues.withoutlep04CH;
    b_muon_puppiIsoWithoutLep04NeutralHadron = puppiIsoValues.withoutlep04NH;
    b_muon_puppiIsoWithoutLep04Photon        = puppiIsoValues.withoutlep04PH;
    
    b_muon_puppiIsoWithoutLep05ChargedHadron = puppiIsoValues.withoutlep05CH;
    b_muon_puppiIsoWithoutLep05NeutralHadron = puppiIsoValues.withoutlep05NH;
    b_muon_puppiIsoWithoutLep05Photon        = puppiIsoValues.withoutlep05PH;
    
    b_muon_puppiIsoNumOfCands       = puppiIsoValues.nNumCands;
    b_muon_puppiIsoNumOfCandsInR05  = puppiIsoValues.nNumCandsInR05;
    b_muon_puppiIsoNumOfCandsOutR05 = puppiIsoValues.nNumCandsOutR05;
    
    b_muon_puppiIsoNumOfCandsInR05OT = puppiIsoValues.nNumCandsInR05OT;
    
    b_muon_puppiIsoNumOfCandsInR05CH = puppiIsoValues.nNumCandsInR05CH;
    b_muon_puppiIsoNumOfCandsInR05NH = puppiIsoValues.nNumCandsInR05NH;
    b_muon_puppiIsoNumOfCandsInR05PH = puppiIsoValues.nNumCandsInR05PH;
    
    b_muon_puppiIsoNumOfCandsInR05CHApp = puppiIsoValues.nNumCandsInR05CHApp;
    b_muon_puppiIsoNumOfCandsInR05NHApp = puppiIsoValues.nNumCandsInR05NHApp;
    b_muon_puppiIsoNumOfCandsInR05PHApp = puppiIsoValues.nNumCandsInR05PHApp;
    
    b_muon_isTightOptimized = isTightMuonCustomOptimized(*mu, pv0);
    b_muon_isTightCustom = isTightMuonCustom(*mu, pv0);
    b_muon_isTight = muon::isTightMuon(*mu, pv0);
    b_muon_isMedium = muon::isMediumMuon(*mu);
    b_muon_isLoose = muon::isLooseMuon(*mu);
    b_muon_isME0Muon = mu->isME0Muon();
    b_muon_isGEMMuon = mu->isGEMMuon();
    b_muon_isRPCMuon = mu->isRPCMuon();
    b_muon_isCaloMuon = mu->isCaloMuon();
    b_muon_isTrackerMuon = mu->isTrackerMuon();
    b_muon_isMuon = mu->isMuon();
    b_muon_isGlobalMuon = mu->isGlobalMuon();
    b_muon_isStandAloneMuon = mu->isStandAloneMuon();
    b_muon_isPFMuon = mu->isPFMuon();
    
    b_muon_isLooseMod = isLooseMod(mu);
    b_muon_isTightModNoIP  = isTightMod(vertexes_, simVertex_, mu, false, false);
    b_muon_isTightModIPxy  = isTightMod(vertexes_, simVertex_, mu, true,  false);
    b_muon_isTightModIPz   = isTightMod(vertexes_, simVertex_, mu, false, true);
    b_muon_isTightModIPxyz = isTightMod(vertexes_, simVertex_, mu, true,  true);
    
    float me0SegX = 100;
    float ge11SegX = 100;
    float ge21SegX = 100;
    for (auto chamber : mu->matches()){
      for (auto segment : chamber.me0Matches){
	if (chamber.detector() == 5){
	  auto me0Segment = (*segment.me0SegmentRef);
	  me0SegX = abs( chamber.x - segment.x );	  
	  if (me0SegX < abs(b_muon_ME0deltaX)){
	    b_muon_ME0deltaX    = ( chamber.x - segment.x );
	    b_muon_ME0deltaY    = ( chamber.y - segment.y );
	    b_muon_ME0pullX  = (chamber.x - segment.x) / std::sqrt(chamber.xErr + segment.xErr);
	    b_muon_ME0pullY  = (chamber.y - segment.y) / std::sqrt(chamber.yErr + segment.yErr);
	    b_muon_ME0dPhi = atan(chamber.dXdZ) - atan(segment.dXdZ);	      
	    b_muon_ME0deltaDXDZ = ( chamber.dXdZ - segment.dXdZ );
	    b_muon_ME0deltaDYDZ = ( chamber.dYdZ - segment.dYdZ );
	    b_muon_ME0noRecHit  = me0Segment.nRecHits();
	  }
	}
      }
      for (auto segment : chamber.gemMatches){
	if (chamber.detector() == 4){
	  auto gemSegment = (*segment.gemSegmentRef);
	  if (gemSegment.gemDetId().station() == 1){
	    ge11SegX = abs( chamber.x - segment.x );
	    if (ge11SegX < abs(b_muon_GE11deltaX)){
	      b_muon_GE11deltaX    = ( chamber.x - segment.x );
	      b_muon_GE11deltaY    = ( chamber.y - segment.y );
	      b_muon_GE11pullX  = (chamber.x - segment.x) / std::sqrt(chamber.xErr + segment.xErr);
	      b_muon_GE11pullY  = (chamber.y - segment.y) / std::sqrt(chamber.yErr + segment.yErr);
	      b_muon_GE11dPhi = atan(chamber.dXdZ) - atan(segment.dXdZ);	      
	      b_muon_GE11deltaDXDZ = ( chamber.dXdZ - segment.dXdZ );
	      b_muon_GE11deltaDYDZ = ( chamber.dYdZ - segment.dYdZ );
	      b_muon_GE11noRecHit  = gemSegment.nRecHits();	      
	    }
	  }
	  if (gemSegment.gemDetId().station() == 2){
	    ge21SegX = abs( chamber.x - segment.x );	  
	    if (ge21SegX < abs(b_muon_GE21deltaX)){
	      b_muon_GE21deltaX    = ( chamber.x - segment.x );
	      b_muon_GE21deltaY    = ( chamber.y - segment.y );
	      b_muon_GE21pullX  = (chamber.x - segment.x) / std::sqrt(chamber.xErr + segment.xErr);
	      b_muon_GE21pullY  = (chamber.y - segment.y) / std::sqrt(chamber.yErr + segment.yErr);
	      b_muon_GE21dPhi = atan(chamber.dXdZ) - atan(segment.dXdZ);	      
	      b_muon_GE21deltaDXDZ = ( chamber.dXdZ - segment.dXdZ );
	      b_muon_GE21deltaDYDZ = ( chamber.dYdZ - segment.dYdZ );
	      b_muon_GE21noRecHit  = gemSegment.nRecHits();	      
	    }
	  }
	}
      }
    }
    
    const reco::Track* muonTrack = 0;  
    if ( mu->globalTrack().isNonnull() ) muonTrack = mu->globalTrack().get();
    else if ( mu->outerTrack().isNonnull()  ) muonTrack = mu->outerTrack().get();
    if (muonTrack){
      b_muon_numberOfValidMuonGEMHits = muonTrack->hitPattern().numberOfValidMuonGEMHits();
      b_muon_numberOfValidMuonME0Hits = muonTrack->hitPattern().numberOfValidMuonME0Hits();
    }

    std::vector<double> tmvaValues = collectTMVAvalues(*mu, pv0);
    b_muon_tmva_bdt = bdt_->GetMvaValue(tmvaValues);
    b_muon_tmva_mlp = mlp_->GetMvaValue(tmvaValues);

    b_muon_global = tmvaValues[0];
    b_muon_pf = tmvaValues[1];
    b_muon_chi2 = tmvaValues[2];
    b_muon_chi2pos = tmvaValues[3];
    b_muon_trkKink = tmvaValues[4];
    b_muon_segcompati = tmvaValues[5];
    b_muon_nstations = tmvaValues[6];
    b_muon_nglobalhits = tmvaValues[7];
    b_muon_trackdxy = tmvaValues[8];
    b_muon_trackdz = tmvaValues[9];
    b_muon_ninnerhits = tmvaValues[10];
    b_muon_trackerlayers =tmvaValues[11];    
  }
  tree->Fill();
}

std::vector<double> MuonAnalyser::collectTMVAvalues(const reco::Muon& mu, reco::Vertex pv0) const
{
  std::vector<double> values;
  int dummyVal = -9;

  //Loose
  values.push_back(mu.isGlobalMuon());
  values.push_back(mu.isPFMuon());
  //Medium
  if ( mu.globalTrack().isNonnull() ){ values.push_back(mu.globalTrack()->normalizedChi2()); }
  else { values.push_back(dummyVal); }
  values.push_back(mu.combinedQuality().chi2LocalPosition);
  values.push_back(mu.combinedQuality().trkKink);
  values.push_back(muon::segmentCompatibility(mu));
  //Tight
  values.push_back(mu.numberOfMatchedStations());
  if ( mu.globalTrack().isNonnull() ){ values.push_back(mu.globalTrack()->hitPattern().numberOfValidMuonHits()); }
  else { values.push_back(dummyVal); }
  if ( mu.muonBestTrack().isNonnull() ){
    values.push_back(abs(mu.muonBestTrack()->dxy(pv0.position())));
    values.push_back(abs(mu.muonBestTrack()->dz(pv0.position())));
  }
  else {
    values.push_back(dummyVal);
    values.push_back(dummyVal);
  }
  if ( mu.innerTrack().isNonnull() ){
    values.push_back(mu.innerTrack()->hitPattern().numberOfValidPixelHits());
    values.push_back(mu.innerTrack()->hitPattern().trackerLayersWithMeasurement());
  }
  else {
    values.push_back(dummyVal);
    values.push_back(dummyVal);
  }

  return values;
}

bool MuonAnalyser::isLooseMuonCustom(const reco::Muon& mu) const
{
  if ( !(mu.isPFMuon()) ) return false;
  if ( !(mu.isGlobalMuon() || mu.isTrackerMuon()) ) return false;
  return true;
}

bool MuonAnalyser::isTightMuonCustom(const reco::Muon& mu, reco::Vertex pv0) const
{
  if ( !(mu.isGlobalMuon()) ) return false;
  if ( !(mu.isPFMuon()) ) return false;
  if ( !(mu.globalTrack().isNonnull()) )  return false;
  if ( !(mu.muonBestTrack().isNonnull()) ) return false;
  if ( !(mu.innerTrack().isNonnull()) ) return false;
  if ( !(mu.globalTrack()->normalizedChi2()<10.) ) return false;
  if ( !(mu.globalTrack()->hitPattern().numberOfValidMuonHits() > 0) ) return false;
  if ( !(mu.numberOfMatchedStations() > 1) ) return false;
  if ( !(abs(mu.muonBestTrack()->dxy(pv0.position())) < 0.2) ) return false;
  //if ( !(abs(mu.muonBestTrack()->dz(pv0.position())) < 0.5) ) return false;
  if ( !(mu.innerTrack()->hitPattern().numberOfValidPixelHits() > 0) ) return false;
  if ( !(mu.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5) ) return false;
  return true;
}

bool MuonAnalyser::isTightMuonCustomOptimized(const reco::Muon& mu, reco::Vertex pv0) const
{
  if ( !(mu.isGlobalMuon()) ) return false;
  if ( !(mu.isPFMuon()) ) return false;
  if ( !(mu.globalTrack().isNonnull()) )  return false;
  if ( !(mu.muonBestTrack().isNonnull()) ) return false;
  if ( !(mu.innerTrack().isNonnull()) ) return false;
  if ( !(mu.globalTrack()->normalizedChi2() < 2.) ) return false; // < 10.
  if ( !(mu.globalTrack()->hitPattern().numberOfValidMuonHits() > 10) ) return false; // > 0
  if ( !(mu.numberOfMatchedStations() > 1) ) return false;
  if ( !(abs(mu.muonBestTrack()->dxy(pv0.position())) < 0.02) ) return false; // < 0.2
  //if ( !(abs(mu.muonBestTrack()->dz(pv0.position())) < 0.5) ) return false;
  if ( !(mu.innerTrack()->hitPattern().numberOfValidPixelHits() > 3) ) return false; // > 0
  if ( !(mu.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 7) ) return false; // > 5
  return true;
}

int MuonAnalyser::nGEMhit(const reco::Muon* muon) const
{
  int noRecHitGEM = 0;
  int noRecHitMuon = 0;
  int noRecHit = 0;
  const reco::Track* muonTrack = 0;  
  if ( muon->globalTrack().isNonnull() ) muonTrack = muon->globalTrack().get();
  else if ( muon->outerTrack().isNonnull()  ) muonTrack = muon->outerTrack().get();
  if (muonTrack){
    // cout << " numberOfValidHits        "<< muonTrack->hitPattern().numberOfValidHits() <<endl;
    // cout << " numberOfValidTrackerHits "<< muonTrack->hitPattern().numberOfValidTrackerHits() <<endl;
    // cout << " numberOfMuonHits         "<< muonTrack->hitPattern().numberOfMuonHits() <<endl;
    // cout << " numberOfValidMuonHits    "<< muonTrack->hitPattern().numberOfValidMuonHits() <<endl;
    // cout << " numberOfValidMuonDTHits  "<< muonTrack->hitPattern().numberOfValidMuonDTHits() <<endl;
    // cout << " numberOfValidMuonCSCHits "<< muonTrack->hitPattern().numberOfValidMuonCSCHits() <<endl;
    // cout << " numberOfValidMuonRPCHits "<< muonTrack->hitPattern().numberOfValidMuonRPCHits() <<endl;
    // cout << " numberOfValidMuonGEMHits "<< muonTrack->hitPattern().numberOfValidMuonGEMHits() <<endl;
    // cout << " numberOfLostMuonGEMHits  "<< muonTrack->hitPattern().numberOfLostMuonGEMHits() <<endl;
    // cout << " numberOfBadMuonGEMHits   "<< muonTrack->hitPattern().numberOfBadMuonGEMHits() <<endl;
    // cout << " numberOfValidMuonME0Hits "<< muonTrack->hitPattern().numberOfValidMuonME0Hits() <<endl;
    for(auto i=muonTrack->recHitsBegin(); i!=muonTrack->recHitsEnd(); i++) {
      DetId hitId = (*i)->geographicalId();
      if (!(*i)->isValid() ) continue;
      if ((*i)->recHits().size()) noRecHit +=(*i)->recHits().size();
      else noRecHit++;
      
      if (hitId.det()!=DetId::Muon) continue;
      if (hitId.subdetId() == MuonSubdetId::GEM) ++noRecHitGEM;

      if ((*i)->recHits().size()) noRecHitMuon +=(*i)->recHits().size();
      else noRecHitMuon++;
      //      cout << "(*i)->size()="<< (*i)->recHits().size()<< " det="<< hitId.det() << " subdet=" << hitId.subdetId() <<endl;
    }
    
  }
  // cout << " noRecHit     "<< noRecHit <<endl;
  // cout << " noRecHitMuon "<< noRecHitMuon <<endl;
  // cout << " noRecHitGEM  "<< noRecHitGEM <<endl;
  return noRecHitGEM;
}

MuonAnalyser::puppiIso MuonAnalyser::getPuppiIso(const reco::Muon *mu, const vector< pat::PackedCandidate> *pcs, int nIsReco) const
{
  puppiIso puppivalues;
  
  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  // 
  // 
  // ATENTION : in the source, this value is taken via configuration; 
  //     "mixFraction" and "dR".
  //     They should be get via configuration.
  //     (But perhaps dR_threshold will be not changed since the def. value is from 
  //     original usage of PUPPI, Jet algorithm.)
  // 
  // 
  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  const double _mix_fraction_ = 0.5;
  //const double dR_threshold = 0.4;
  
  //double dR2_threshold = dR_threshold * dR_threshold;
  double dR2_threshold03 = 0.3 * 0.3;
  double dR2_threshold04 = 0.4 * 0.4;
  double dR2_threshold05 = 0.5 * 0.5;

  double val_PuppiWithLep03    [3]= {0,0,0} ;
  double val_PuppiWithoutLep03 [3]= {0,0,0} ;
  double val_PuppiWithLep04    [3]= {0,0,0} ;
  double val_PuppiWithoutLep04 [3]= {0,0,0} ;
  double val_PuppiWithLep05    [3]= {0,0,0} ;
  double val_PuppiWithoutLep05 [3]= {0,0,0} ;
  
  int nNumCands       = 0;
  int nNumCandsInR05  = 0;
  int nNumCandsOutR05 = 0;
  
  int nNumCandsInR05OT = 0;
  
  int nNumCandsInR05CH = 0;
  int nNumCandsInR05NH = 0;
  int nNumCandsInR05PH = 0;
  
  int nNumCandsInR05CHApp = 0;
  int nNumCandsInR05NHApp = 0;
  int nNumCandsInR05PHApp = 0;
  
  double dTrkSumPT();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  // loop ever all the candidates, and accumulate PT deposit around the lepton.
  for( std::vector<pat::PackedCandidate>::const_iterator cand = pcs -> begin();
       cand != pcs->end();
       cand ++ )
    {
      // calc DR
      nNumCands++;

      double d_eta = abs( cand->eta() - mu->eta() ) ;
      double d_phi = abs( cand->phi() - mu->phi() ) ; 
      d_phi = ( d_phi < M_PI ) ? d_phi : 2 * M_PI - d_phi ; 
      double dR2 = d_eta * d_eta  + d_phi * d_phi ;
    
      if( dR2 > dR2_threshold05 ) {
	nNumCandsOutR05++;
	continue ;
      }
    
      nNumCandsInR05++;
    
      long nIDAbs = abs(cand -> pdgId());

      // check particleTyple (CH/NH/PH or other). remove 'other'.
      const particleType pType =
	isCH( nIDAbs ) ? CH :
	isNH( nIDAbs ) ? NH :
	isPH( nIDAbs ) ? PH : OTHER ;
    
      if( pType == OTHER ) {
	if( cand -> pdgId() != 1 && cand -> pdgId() != 2 // d quark and u quark
	    && nIDAbs != 11  // electron
	    && nIDAbs != 13) // muon
	  {
	    std::cout <<"candidate with PDGID = " << cand -> pdgId() << " is not CH/NH/PH/e/mu or 1/2 "
	      "(and this is removed from isolation calculation)"  << std::endl ; 
	  }
      
	nNumCandsInR05OT++;
      
	continue ;
      }
    
      if( pType == CH ) {nNumCandsInR05CH++; }
      if( pType == NH ) nNumCandsInR05NH++;
      if( pType == PH ) nNumCandsInR05PH++;
    
      // Check particleType dependent DR cut (remove overlapped candiadte)
      // The threshold values were taken from 'MuonPFIsolationSequence_cff.py'.
      if( pType == CH && dR2 < 0.0001*0.0001 ) {
	continue ;
      } else {
	nNumCandsInR05CHApp++;
      }
    
      if( pType == NH && dR2 < 0.01  *0.01   ) {
	continue ;
      } else {
	nNumCandsInR05NHApp++;
      }
    
      if( pType == PH && dR2 < 0.01  *0.01   ) {
	continue ;
      } else {
	nNumCandsInR05PHApp++;
      }

      // The candidate passed all the selection.
      // Now, add its PT to the variable with weight.

      val_PuppiWithLep05   [ pType ] += cand -> pt() * cand -> puppiWeight() ;
      val_PuppiWithoutLep05[ pType ] += cand -> pt() * cand -> puppiWeightNoLep();
      
      if ( dR2 <= dR2_threshold04 ) {
        val_PuppiWithLep04   [ pType ] += cand -> pt() * cand -> puppiWeight() ;
        val_PuppiWithoutLep04[ pType ] += cand -> pt() * cand -> puppiWeightNoLep();
        
        if ( dR2 <= dR2_threshold03 ) {
          val_PuppiWithLep03   [ pType ] += cand -> pt() * cand -> puppiWeight() ;
          val_PuppiWithoutLep03[ pType ] += cand -> pt() * cand -> puppiWeightNoLep();
        }
      }

   
    }// end of candidate LOOP.

  const double reliso_Puppi_withLep03    = ( val_PuppiWithLep03   [CH] + val_PuppiWithLep03   [NH] + val_PuppiWithLep03   [PH] ) / mu->pt() ;
  const double reliso_Puppi_withoutlep03 = ( val_PuppiWithoutLep03[CH] + val_PuppiWithoutLep03[NH] + val_PuppiWithoutLep03[PH] ) / mu->pt() ;

  const double reliso_Puppi_withLep04    = ( val_PuppiWithLep04   [CH] + val_PuppiWithLep04   [NH] + val_PuppiWithLep04   [PH] ) / mu->pt() ;
  const double reliso_Puppi_withoutlep04 = ( val_PuppiWithoutLep04[CH] + val_PuppiWithoutLep04[NH] + val_PuppiWithoutLep04[PH] ) / mu->pt() ;

  const double reliso_Puppi_withLep05    = ( val_PuppiWithLep05   [CH] + val_PuppiWithLep05   [NH] + val_PuppiWithLep05   [PH] ) / mu->pt() ;
  const double reliso_Puppi_withoutlep05 = ( val_PuppiWithoutLep05[CH] + val_PuppiWithoutLep05[NH] + val_PuppiWithoutLep05[PH] ) / mu->pt() ;

  puppivalues.withLep03    = reliso_Puppi_withLep03;
  puppivalues.withoutlep03 = reliso_Puppi_withoutlep03;
  puppivalues.combined03   = _mix_fraction_ * reliso_Puppi_withLep03 + ( 1.0 - _mix_fraction_) * reliso_Puppi_withoutlep03;

  puppivalues.withLep      = reliso_Puppi_withLep04;
  puppivalues.withoutlep   = reliso_Puppi_withoutlep04;
  puppivalues.combined     = _mix_fraction_ * reliso_Puppi_withLep04 + ( 1.0 - _mix_fraction_) * reliso_Puppi_withoutlep04;

  puppivalues.withLep05    = reliso_Puppi_withLep05;
  puppivalues.withoutlep05 = reliso_Puppi_withoutlep05;
  puppivalues.combined05   = _mix_fraction_ * reliso_Puppi_withLep05 + ( 1.0 - _mix_fraction_) * reliso_Puppi_withoutlep05;
  
  puppivalues.withLep03CH = val_PuppiWithLep03[CH];
  puppivalues.withLep03NH = val_PuppiWithLep03[NH];
  puppivalues.withLep03PH = val_PuppiWithLep03[PH];
  
  puppivalues.withLep04CH = val_PuppiWithLep04[CH];
  puppivalues.withLep04NH = val_PuppiWithLep04[NH];
  puppivalues.withLep04PH = val_PuppiWithLep04[PH];
  
  puppivalues.withLep05CH = val_PuppiWithLep05[CH];
  puppivalues.withLep05NH = val_PuppiWithLep05[NH];
  puppivalues.withLep05PH = val_PuppiWithLep05[PH];
  
  puppivalues.withoutlep03CH = val_PuppiWithoutLep03[CH];
  puppivalues.withoutlep03NH = val_PuppiWithoutLep03[NH];
  puppivalues.withoutlep03PH = val_PuppiWithoutLep03[PH];
  
  puppivalues.withoutlep04CH = val_PuppiWithoutLep04[CH];
  puppivalues.withoutlep04NH = val_PuppiWithoutLep04[NH];
  puppivalues.withoutlep04PH = val_PuppiWithoutLep04[PH];
  
  puppivalues.withoutlep05CH = val_PuppiWithoutLep05[CH];
  puppivalues.withoutlep05NH = val_PuppiWithoutLep05[NH];
  puppivalues.withoutlep05PH = val_PuppiWithoutLep05[PH];
  
  puppivalues.nNumCands       = nNumCands;
  puppivalues.nNumCandsInR05  = nNumCandsInR05;
  puppivalues.nNumCandsOutR05 = nNumCandsOutR05;
  
  puppivalues.nNumCandsInR05OT = nNumCandsInR05OT;
  
  puppivalues.nNumCandsInR05CH = nNumCandsInR05CH;
  puppivalues.nNumCandsInR05NH = nNumCandsInR05NH;
  puppivalues.nNumCandsInR05PH = nNumCandsInR05PH;
  
  puppivalues.nNumCandsInR05CHApp = nNumCandsInR05CHApp;
  puppivalues.nNumCandsInR05NHApp = nNumCandsInR05NHApp;
  puppivalues.nNumCandsInR05PHApp = nNumCandsInR05PHApp;
  
  return puppivalues;
}

bool MuonAnalyser::isNH( long pdgidAbs ) const{
  //     pdgId = cms.vint32(111,130,310,2112),
  if( pdgidAbs == 111 )  return true ; // pion0
  if( pdgidAbs == 130 )  return true ; // Kaon0L
  if( pdgidAbs == 310 )  return true ; // Kaon0S
  if( pdgidAbs == 2112 ) return true ; // Neutron
  return false;
}
bool MuonAnalyser::isCH( long pdgidAbs ) const{
  //  pdgId = cms.vint32(211,-211,321,-321,999211,2212,-2212),
  if( pdgidAbs == 211    ) return true ; // Pion+
  if( pdgidAbs == 321    ) return true ; // Kaon+
  if( pdgidAbs == 999211 ) return true ; // ???
  if( pdgidAbs == 2212   ) return true ; // Proton
  return false;
}
bool MuonAnalyser::isPH( long pdgidAbs ) const{
  if( pdgidAbs == 22 ) return true ; // Photon
  return false;
}

bool MuonAnalyser::isGlobalTightMuon( const reco::Muon *muonRef ) {

  //if ( !muonRef.isNonnull() ) return false;

  if ( !muonRef->isGlobalMuon() ) return false;
  if ( !muonRef->isStandAloneMuon() ) return false;
 
 
  if ( muonRef->isTrackerMuon() ) { 
   
    bool result = muon::isGoodMuon(*muonRef,muon::GlobalMuonPromptTight);
   
    bool isTM2DCompatibilityTight =  muon::isGoodMuon(*muonRef,muon::TM2DCompatibilityTight);   
    int nMatches = muonRef->numberOfMatches();
    bool quality = nMatches > 2 || isTM2DCompatibilityTight;
   
    return result && quality;
   
  } else {
 
    reco::TrackRef standAloneMu = muonRef->standAloneMuon();
   
    // No tracker muon -> Request a perfect stand-alone muon, or an even better global muon
    bool result = false;
      
    // Check the quality of the stand-alone muon : 
    // good chi**2 and large number of hits and good pt error
    if ( ( standAloneMu->hitPattern().numberOfValidMuonDTHits() < 22 &&
	   standAloneMu->hitPattern().numberOfValidMuonCSCHits() < 15 ) ||
	 standAloneMu->normalizedChi2() > 10. || 
	 standAloneMu->ptError()/standAloneMu->pt() > 0.20 ) {
      result = false;
    } else { 
      
      reco::TrackRef combinedMu = muonRef->combinedMuon();
      reco::TrackRef trackerMu = muonRef->track();
            
      // If the stand-alone muon is good, check the global muon
      if ( combinedMu->normalizedChi2() > standAloneMu->normalizedChi2() ) {
	// If the combined muon is worse than the stand-alone, it 
	// means that either the corresponding tracker track was not 
	// reconstructed, or that the sta muon comes from a late 
	// pion decay (hence with a momentum smaller than the track)
	// Take the stand-alone muon only if its momentum is larger
	// than that of the track
	result = standAloneMu->pt() > trackerMu->pt() ;
      } else { 
	// If the combined muon is better (and good enough), take the 
	// global muon
	result = 
	  combinedMu->ptError()/combinedMu->pt() < 
	    std::min(0.20,standAloneMu->ptError()/standAloneMu->pt());
      }
    }      

    return result;    
  }

  return false;

}

bool MuonAnalyser::isTrackerTightMuon( const reco::Muon *muonRef ) {

  //if ( !muonRef.isNonnull() ) return false;
    
  if( !muonRef->isTrackerMuon() ) return false;
  
  reco::TrackRef trackerMu = muonRef->track();
  const reco::Track& track = *trackerMu;
  
  unsigned nTrackerHits =  track.hitPattern().numberOfValidTrackerHits();
  
  if(nTrackerHits<=12) return false;
  
  bool isAllArbitrated = muon::isGoodMuon(*muonRef,muon::AllArbitrated);
  
  bool isTM2DCompatibilityTight = muon::isGoodMuon(*muonRef,muon::TM2DCompatibilityTight);
  
  if(!isAllArbitrated || !isTM2DCompatibilityTight)  return false;

  if((trackerMu->ptError()/trackerMu->pt() > 0.10)){
    //std::cout<<" PT ERROR > 10 % "<< trackerMu->pt() <<std::endl;
    return false;
  }
  return true;
  
}

bool MuonAnalyser::isIsolatedMuon( const reco::Muon *muonRef ){


  //if ( !muonRef.isNonnull() ) return false;
  if ( !muonRef->isIsolationValid() ) return false;
  
  // Isolated Muons which are missed by standard cuts are nearly always global+tracker
  if ( !muonRef->isGlobalMuon() ) return false;

  // If it's not a tracker muon, only take it if there are valid muon hits

  reco::TrackRef standAloneMu = muonRef->standAloneMuon();

  if ( !muonRef->isTrackerMuon() ){
    if(standAloneMu->hitPattern().numberOfValidMuonDTHits() == 0 &&
       standAloneMu->hitPattern().numberOfValidMuonCSCHits() ==0) return false;
  }
  
  // for isolation, take the smallest pt available to reject fakes

  reco::TrackRef combinedMu = muonRef->combinedMuon();
  double smallestMuPt = combinedMu->pt();
  
  if(standAloneMu->pt()<smallestMuPt) smallestMuPt = standAloneMu->pt();
  
  if(muonRef->isTrackerMuon())
    {
      reco::TrackRef trackerMu = muonRef->track();
      if(trackerMu->pt() < smallestMuPt) smallestMuPt= trackerMu->pt();
    }
     
  double sumPtR03 = muonRef->isolationR03().sumPt;
  double emEtR03 = muonRef->isolationR03().emEt;
  double hadEtR03 = muonRef->isolationR03().hadEt;
  
  double relIso = (sumPtR03 + emEtR03 + hadEtR03)/smallestMuPt;

  if(relIso<0.1) return true;
  else return false;
}

bool MuonAnalyser::isLooseMod(const reco::Muon *muon)
{
  bool isPF= isGlobalTightMuon(muon) || isTrackerTightMuon(muon) || isIsolatedMuon(muon);
  bool isGLB = muon->isGlobalMuon();
  bool isTrk = muon->isTrackerMuon();
    
  return ( isPF && (isGLB || isTrk) );
}

bool MuonAnalyser::isTightMod(const reco::VertexCollection* vertices, const SimVertex &simPVh, const reco::Muon *muon, bool useIPxy, bool useIPz)
{
  bool result = false;
    
  if (muon->muonBestTrack().isNonnull() && muon->innerTrack().isNonnull() && muon->globalTrack().isNonnull()){
        
    //std::vector<double> vtxCoord = findSimVtx(iEvent);
    std::vector<double> vtxCoord;
    vtxCoord.push_back(1.0);
        
    //GlobalPoint point(vtxCoord[1],vtxCoord[2],vtxCoord[3]);
    //GlobalPoint pointDY(vtxCoord[4],vtxCoord[5],vtxCoord[6]);
    GlobalPoint point(simPVh.position().x(), simPVh.position().y(), simPVh.position().z());
    GlobalPoint pointDY(simPVh.position().x(), simPVh.position().y(), simPVh.position().z());
        
    //double muonX = muon->vx();
    //double muonY = muon->vy();
    //double muonZ = muon->vz();
        
    double muonZ = pointDY.z();
        
    //edm::Handle<reco::VertexCollection> vertexHandle; // quark2 : It is given by argument
    //iEvent.getByToken(vtx_Token,vertexHandle);
        
    double distInit = 24;
    int indexFinal = 0;
    for(int i = 0; i < (int)vertices->size(); i++){
            
      //double vtxX = (*vertices)[i].x();
      //double vtxY = (*vertices)[i].y();
      double vtxZ = (*vertices)[i].z();
            
      double dist = abs(muonZ - vtxZ);
      //std::cout<<"dist "<<dist<<std::endl;
      if(dist < distInit){
                
	distInit = dist;
	indexFinal = i;
                
      }
            
    }
    //std::cout<<distInit<<" "<<indexFinal<<std::endl;
        
    double ipxySim = 999;
    double ipzSim = 999;
        
    if(vtxCoord[0] > 1.5 && vtxCoord[0] < 3.5){//Mu and nu gun samples
            
      ipxySim = abs(muon->muonBestTrack()->dxy(math::XYZPoint(point.x(),point.y(),point.z())));
      ipzSim = abs(muon->muonBestTrack()->dz(math::XYZPoint(point.x(),point.y(),point.z())));
            
    }
    else if(vtxCoord[0] > 0.5 && vtxCoord[0] < 1.5){//DY samples
            
      ipxySim = abs(muon->muonBestTrack()->dxy(math::XYZPoint(pointDY.x(),pointDY.y(),pointDY.z())));
      ipzSim = abs(muon->muonBestTrack()->dz(math::XYZPoint(pointDY.x(),pointDY.y(),pointDY.z())));
            
    }
    bool ipxySimBool = ipxySim < 0.2;
    bool ipzSimBool = ipzSim < 0.5;
    //std::cout<<"vx: "<<point.x()<<" vy: "<<point.y()<<" vz: "<<point.z()<<" |Dxy|: "<<ipxySim<<" "<<ipxySimBool<<" |Dz|: "<<ipzSim<<" "<<ipzSimBool<<std::endl;
    //std::cout<<"vx: "<<pointDY.x()<<" vy: "<<pointDY.y()<<" vz: "<<pointDY.z()<<" |Dxy|: "<<ipxySim<<" "<<ipxySimBool<<" |Dz|: "<<ipzSim<<" "<<ipzSimBool<<std::endl;
        
    bool trkLayMeas = muon->innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5;
    bool isGlb = muon->isGlobalMuon();
    bool isPF= isGlobalTightMuon(muon) || isTrackerTightMuon(muon) || isIsolatedMuon(muon);
    bool chi2 = muon->globalTrack()->normalizedChi2() < 10.;
    bool validHits = muon->globalTrack()->hitPattern().numberOfValidMuonHits() > 0;
    bool matchedSt = muon->numberOfMatchedStations() > 1;
        
    bool ipxy = false;
    bool ipz = false;
    if(useIPxy == true){
            
      if(vertices->size() !=0 && vtxCoord[0] > 0.5 && vtxCoord[0] < 1.5){
                
	ipxy = abs(muon->muonBestTrack()->dxy((*vertices)[indexFinal].position())) < 0.2;
	//std::cout<<"vx: "<<pointDY.x()<<" vy: "<<pointDY.y()<<" vz: "<<pointDY.z()<<" |Dxy|: "<<ipxy<<std::endl;
                
      }
      else if(vtxCoord[0] > 1.5 && vtxCoord[0] < 3.5){
                
	ipxy = ipxySimBool;
	//std::cout<<"vx: "<<point.x()<<" vy: "<<point.y()<<" vz: "<<point.z()<<" |Dxy|: "<<ipxy<<std::endl;
                
      }
            
    }
    else if(useIPxy == false) ipxy = true;
        
    if(useIPz == true){
            
      if(vertices->size() !=0 && vtxCoord[0] > 0.5 && vtxCoord[0] < 1.5){
                
	ipz = abs(muon->muonBestTrack()->dz((*vertices)[indexFinal].position())) < 0.5;
	//std::cout<<"vx: "<<pointDY.x()<<" vy: "<<pointDY.y()<<" vz: "<<pointDY.z()<<" |Dz|: "<<ipz<<std::endl;
                
      }
      else if(vtxCoord[0] > 1.5 && vtxCoord[0] < 3.5){
                
	ipz = ipzSimBool;
	//std::cout<<"vx: "<<point.x()<<" vy: "<<point.y()<<" vz: "<<point.z()<<" |Dz|: "<<ipz<<std::endl;
                
      }
            
    }
    else if(useIPz == false) ipz = true;
        
    bool validPxlHit = muon->innerTrack()->hitPattern().numberOfValidPixelHits() > 0;
    //bool validPxlHit = muon->innerTrack()->hitPattern().pixelLayersWithMeasurement(3,2) > 0;
    //bool validPxlHit = muon->innerTrack()->hitPattern().pixelLayersWithMeasurement(4,3) > 0;
        
    //std::cout<<trkLayMeas<<" "<<isGlb<<" "<<isPF<<" "<<chi2<<" "<<validHits<<" "<<matchedSt<<" "<<ipxy<<" "<<ipz<<" "<<validPxlHit<<std::endl;
        
    if(trkLayMeas && isGlb && isPF && chi2 && validHits && matchedSt && ipxy && ipz && validPxlHit) result = true;
        
  }
    
  return result;
}
void MuonAnalyser::setBranches(TTree *tree)
{
  tree->Branch("nvertex", &b_nvertex, "nvertex/I");
  tree->Branch("muon", "TLorentzVector", &b_muon);  
  tree->Branch("muon_pdgId", &b_muon_pdgId, "muon_pdgId/I");
  tree->Branch("muon_signal", &b_muon_signal, "muon_signal/O");
  tree->Branch("muon_pTresolution",&b_muon_pTresolution,"muon_pTresolution/F");
  tree->Branch("muon_pTinvresolution",&b_muon_pTinvresolution,"muon_pTinvresolution/F");
  tree->Branch("muon_isTight", &b_muon_isTight, "muon_isTight/O");
  tree->Branch("muon_isMedium", &b_muon_isMedium, "muon_isMedium/O");
  tree->Branch("muon_isLoose", &b_muon_isLoose, "muon_isLoose/O");
  tree->Branch("muon_isMuon", &b_muon_isMuon, "muon_isMuon/O");
  tree->Branch("muon_isTrackerMuon", &b_muon_isTrackerMuon, "muon_isTrackerMuon/O");
  tree->Branch("muon_isGlobalMuon", &b_muon_isGlobalMuon, "muon_isGlobalMuon/O");
  tree->Branch("muon_isStandAloneMuon", &b_muon_isStandAloneMuon, "muon_isStandAloneMuon/O");
  tree->Branch("muon_isPFMuon", &b_muon_isPFMuon, "muon_isPFMuon/O");
  tree->Branch("muon_isME0Muon", &b_muon_isME0Muon, "muon_isME0Muon/O");
  tree->Branch("muon_isGEMMuon", &b_muon_isGEMMuon, "muon_isGEMMuon/O");
  tree->Branch("muon_isRPCMuon", &b_muon_isRPCMuon, "muon_isRPCMuon/O");
  tree->Branch("muon_isCaloMuon", &b_muon_isCaloMuon, "muon_isCaloMuon/O");

  tree->Branch("muon_isLooseMod", &b_muon_isLooseMod, "muon_isLooseMod/O");
  tree->Branch("muon_isTightOptimized", &b_muon_isTightOptimized, "muon_isTightOptimized/O");
  tree->Branch("muon_isTightCustom", &b_muon_isTightCustom, "muon_isTightCustom/O");
  tree->Branch("muon_isTightModNoIP", &b_muon_isTightModNoIP, "muon_isTightModNoIP/O");
  tree->Branch("muon_isTightModIPxy", &b_muon_isTightModIPxy, "muon_isTightModIPxy/O");
  tree->Branch("muon_isTightModIPz", &b_muon_isTightModIPz, "muon_isTightModIPz/O");
  tree->Branch("muon_isTightModIPxyz", &b_muon_isTightModIPxyz, "muon_isTightModIPxyz/O");
  
  tree->Branch("muon_isGlobalMuon", &b_muon_global, "muon_isGlobalMuon/O");
  tree->Branch("muon_isPFMuon", &b_muon_pf, "muon_isPFMuon/O");
  tree->Branch("muon_normalizedChi2", &b_muon_chi2, "muon_normalizedChi2/F");
  tree->Branch("muon_chi2LocalPosition", &b_muon_chi2pos, "muon_chi2LocalPosition/F");
  tree->Branch("muon_trkKink", &b_muon_trkKink, "muon_trkKink/F");
  tree->Branch("muon_segmentCompatibility", &b_muon_segcompati, "muon_segmentCompatibility/F");
  tree->Branch("muon_numberOfValidMuonHits", &b_muon_nglobalhits, "muon_numberOfValidMuonHits/I");
  tree->Branch("muon_numberOfMatchedStations", &b_muon_nstations, "muon_numberOfMatchedStations/I");
  tree->Branch("muon_pv0pos_dxy", &b_muon_trackdxy, "muon_pv0pos_dxy/F");
  tree->Branch("muon_pv0pos_dz", &b_muon_trackdz, "muon_pv0pos_dz/F");
  tree->Branch("muon_numberOfValidPixelHits", &b_muon_ninnerhits, "muon_numberOfValidPixelHits/I");
  tree->Branch("muon_trackerLayersWithMeasurement", &b_muon_trackerlayers, "muon_trackerLayersWithMeasurement/F");
  tree->Branch("muon_tmva_bdt", &b_muon_tmva_bdt, "muon_tmva_bdt/F");
  tree->Branch("muon_tmva_mlp", &b_muon_tmva_mlp, "muon_tmva_mlp/F");  

  tree->Branch("muon_ME0deltaX", &b_muon_ME0deltaX, "muon_ME0deltaX/F");  
  tree->Branch("muon_ME0deltaY", &b_muon_ME0deltaY, "muon_ME0deltaY/F");  
  tree->Branch("muon_ME0pullX", &b_muon_ME0pullX, "muon_ME0pullX/F");  
  tree->Branch("muon_ME0pullY", &b_muon_ME0pullY, "muon_ME0pullY/F");  
  tree->Branch("muon_ME0dPhi", &b_muon_ME0dPhi, "muon_ME0dPhi/F");  
  tree->Branch("muon_ME0deltaDXDZ", &b_muon_ME0deltaDXDZ, "muon_ME0deltaDXDZ/F");  
  tree->Branch("muon_ME0deltaDYDZ", &b_muon_ME0deltaDYDZ, "muon_ME0deltaDYDZ/F");  
  tree->Branch("muon_ME0noRecHit", &b_muon_ME0noRecHit, "muon_ME0noRecHit/I");  

  tree->Branch("muon_GE11deltaX", &b_muon_GE11deltaX, "muon_GE11deltaX/F");  
  tree->Branch("muon_GE11deltaY", &b_muon_GE11deltaY, "muon_GE11deltaY/F");  
  tree->Branch("muon_GE11pullX", &b_muon_GE11pullX, "muon_GE11pullX/F");  
  tree->Branch("muon_GE11pullY", &b_muon_GE11pullY, "muon_GE11pullY/F");  
  tree->Branch("muon_GE11dPhi", &b_muon_GE11dPhi, "muon_GE11dPhi/F");  
  tree->Branch("muon_GE11deltaDXDZ", &b_muon_GE11deltaDXDZ, "muon_GE11deltaDXDZ/F");  
  tree->Branch("muon_GE11deltaDYDZ", &b_muon_GE11deltaDYDZ, "muon_GE11deltaDYDZ/F");  
  tree->Branch("muon_GE11noRecHit", &b_muon_GE11noRecHit, "muon_GE11noRecHit/I");  
  
  tree->Branch("muon_GE21deltaX", &b_muon_GE21deltaX, "muon_GE21deltaX/F");  
  tree->Branch("muon_GE21deltaY", &b_muon_GE21deltaY, "muon_GE21deltaY/F");  
  tree->Branch("muon_GE21pullX", &b_muon_GE21pullX, "muon_GE21pullX/F");  
  tree->Branch("muon_GE21pullY", &b_muon_GE21pullY, "muon_GE21pullY/F");  
  tree->Branch("muon_GE21dPhi", &b_muon_GE21dPhi, "muon_GE21dPhi/F");  
  tree->Branch("muon_GE21deltaDXDZ", &b_muon_GE21deltaDXDZ, "muon_GE21deltaDXDZ/F");  
  tree->Branch("muon_GE21deltaDYDZ", &b_muon_GE21deltaDYDZ, "muon_GE21deltaDYDZ/F");  
  tree->Branch("muon_GE21noRecHit", &b_muon_GE21noRecHit, "muon_GE21noRecHit/I");  
  
  tree->Branch("muon_numberOfValidMuonGEMHits",&b_muon_numberOfValidMuonGEMHits,"muon_numberOfValidMuonGEMHits/I");
  tree->Branch("muon_numberOfValidMuonME0Hits",&b_muon_numberOfValidMuonME0Hits,"muon_numberOfValidMuonME0Hits/I");
  tree->Branch("muon_poszPV0",&b_muon_poszPV0,"muon_poszPV0/F");
  tree->Branch("muon_poszSimPV",&b_muon_poszSimPV,"muon_poszSimPV/F");
  tree->Branch("muon_poszMuon",&b_muon_poszMuon,"muon_poszMuon/F");
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
  tree->Branch("muon_puppiIsoWithLep",&b_muon_puppiIsoWithLep,"muon_puppiIsoWithLep/F");
  tree->Branch("muon_puppiIsoWithoutLep",&b_muon_puppiIsoWithoutLep,"muon_puppiIsoWithoutLep/F");
  tree->Branch("muon_puppiIsoCombined",&b_muon_puppiIsoCombined,"muon_puppiIsoCombined/F");
  tree->Branch("muon_puppiIsoWithLep03",&b_muon_puppiIsoWithLep03,"muon_puppiIsoWithLep03/F");
  tree->Branch("muon_puppiIsoWithoutLep03",&b_muon_puppiIsoWithoutLep03,"muon_puppiIsoWithoutLep03/F");
  tree->Branch("muon_puppiIsoCombined03",&b_muon_puppiIsoCombined03,"muon_puppiIsoCombined03/F");
  tree->Branch("muon_puppiIsoWithLep05",&b_muon_puppiIsoWithLep05,"muon_puppiIsoWithLep05/F");
  tree->Branch("muon_puppiIsoWithoutLep05",&b_muon_puppiIsoWithoutLep05,"muon_puppiIsoWithoutLep05/F");
  tree->Branch("muon_puppiIsoCombined05",&b_muon_puppiIsoCombined05,"muon_puppiIsoCombined05/F");
  tree->Branch("muon_puppiIsoWithLep03ChargedHadron",&b_muon_puppiIsoWithLep03ChargedHadron,"muon_puppiIsoWithLep03ChargedHadron/F");
  tree->Branch("muon_puppiIsoWithLep03NeutralHadron",&b_muon_puppiIsoWithLep03NeutralHadron,"muon_puppiIsoWithLep03NeutralHadron/F");
  tree->Branch("muon_puppiIsoWithLep03Photon",&b_muon_puppiIsoWithLep03Photon,"muon_puppiIsoWithLep03Photon/F");
  tree->Branch("muon_puppiIsoWithLep04ChargedHadron",&b_muon_puppiIsoWithLep04ChargedHadron,"muon_puppiIsoWithLep04ChargedHadron/F");
  tree->Branch("muon_puppiIsoWithLep04NeutralHadron",&b_muon_puppiIsoWithLep04NeutralHadron,"muon_puppiIsoWithLep04NeutralHadron/F");
  tree->Branch("muon_puppiIsoWithLep04Photon",&b_muon_puppiIsoWithLep04Photon,"muon_puppiIsoWithLep04Photon/F");
  tree->Branch("muon_puppiIsoWithLep05ChargedHadron",&b_muon_puppiIsoWithLep05ChargedHadron,"muon_puppiIsoWithLep05ChargedHadron/F");
  tree->Branch("muon_puppiIsoWithLep05NeutralHadron",&b_muon_puppiIsoWithLep05NeutralHadron,"muon_puppiIsoWithLep05NeutralHadron/F");
  tree->Branch("muon_puppiIsoWithLep05Photon",&b_muon_puppiIsoWithLep05Photon,"muon_puppiIsoWithLep05Photon/F");
  tree->Branch("muon_puppiIsoWithoutLep03ChargedHadron",&b_muon_puppiIsoWithoutLep03ChargedHadron,"muon_puppiIsoWithoutLep03ChargedHadron/F");
  tree->Branch("muon_puppiIsoWithoutLep03NeutralHadron",&b_muon_puppiIsoWithoutLep03NeutralHadron,"muon_puppiIsoWithoutLep03NeutralHadron/F");
  tree->Branch("muon_puppiIsoWithoutLep03Photon",&b_muon_puppiIsoWithoutLep03Photon,"muon_puppiIsoWithoutLep03Photon/F");
  tree->Branch("muon_puppiIsoWithoutLep04ChargedHadron",&b_muon_puppiIsoWithoutLep04ChargedHadron,"muon_puppiIsoWithoutLep04ChargedHadron/F");
  tree->Branch("muon_puppiIsoWithoutLep04NeutralHadron",&b_muon_puppiIsoWithoutLep04NeutralHadron,"muon_puppiIsoWithoutLep04NeutralHadron/F");
  tree->Branch("muon_puppiIsoWithoutLep04Photon",&b_muon_puppiIsoWithoutLep04Photon,"muon_puppiIsoWithoutLep04Photon/F");
  tree->Branch("muon_puppiIsoWithoutLep05ChargedHadron",&b_muon_puppiIsoWithoutLep05ChargedHadron,"muon_puppiIsoWithoutLep05ChargedHadron/F");
  tree->Branch("muon_puppiIsoWithoutLep05NeutralHadron",&b_muon_puppiIsoWithoutLep05NeutralHadron,"muon_puppiIsoWithoutLep05NeutralHadron/F");
  tree->Branch("muon_puppiIsoWithoutLep05Photon",&b_muon_puppiIsoWithoutLep05Photon,"muon_puppiIsoWithoutLep05Photon/F");
  tree->Branch("muon_puppiIsoNumOfCands",&b_muon_puppiIsoNumOfCands,"muon_puppiIsoNumOfCands/I");
  tree->Branch("muon_puppiIsoNumOfCandsInR05",&b_muon_puppiIsoNumOfCandsInR05,"muon_puppiIsoNumOfCandsInR05/I");
  tree->Branch("muon_puppiIsoNumOfCandsOutR05",&b_muon_puppiIsoNumOfCandsOutR05,"muon_puppiIsoNumOfCandsOutR05/I");
  tree->Branch("muon_puppiIsoNumOfCandsInR05OT",&b_muon_puppiIsoNumOfCandsInR05OT,"muon_puppiIsoNumOfCandsInR05OT/I");
  tree->Branch("muon_puppiIsoNumOfCandsInR05CH",&b_muon_puppiIsoNumOfCandsInR05CH,"muon_puppiIsoNumOfCandsInR05CH/I");
  tree->Branch("muon_puppiIsoNumOfCandsInR05NH",&b_muon_puppiIsoNumOfCandsInR05NH,"muon_puppiIsoNumOfCandsInR05NH/I");
  tree->Branch("muon_puppiIsoNumOfCandsInR05PH",&b_muon_puppiIsoNumOfCandsInR05PH,"muon_puppiIsoNumOfCandsInR05PH/I");
  tree->Branch("muon_puppiIsoNumOfCandsInR05CHApp",&b_muon_puppiIsoNumOfCandsInR05CHApp,"muon_puppiIsoNumOfCandsInR05CHApp/I");
  tree->Branch("muon_puppiIsoNumOfCandsInR05NHApp",&b_muon_puppiIsoNumOfCandsInR05NHApp,"muon_puppiIsoNumOfCandsInR05NHApp/I");
  tree->Branch("muon_puppiIsoNumOfCandsInR05PHApp",&b_muon_puppiIsoNumOfCandsInR05PHApp,"muon_puppiIsoNumOfCandsInR05PHApp/I");
}
DEFINE_FWK_MODULE(MuonAnalyser);
