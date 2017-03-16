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
  bool isTightMod(const Handle<VertexCollection> &vertexHandle, const SimVertex &simPVh, const reco::Muon *muon, bool useIPxy, bool useIPz);
  
  std::vector<double> collectTMVAvalues(const reco::Muon& mu, reco::Vertex pv0) const;
  int nGEMhit(const reco::Muon * mu) const;
  int nME0hit(const reco::Muon * mu) const;
  void treereset();
  
  puppiIso getPuppiIso(const reco::Muon *mu, const vector< pat::PackedCandidate> *pcs, int nIsReco) const;
  bool isNH( long pdgid ) const;
  bool isCH( long pdgid ) const;
  bool isPH( long pdgid ) const;

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
  TLorentzVector b_genMuon;
  bool b_genMuon_isTightOptimized, b_genMuon_isTightCustom, b_genMuon_isTight, b_genMuon_isMedium, b_genMuon_isLoose;
  bool b_genMuon_isME0Muon, b_genMuon_isGEMMuon, b_genMuon_isMuon, b_genMuon_signal;
  bool b_genMuon_isGlobalMuon, b_genMuon_isStandAloneMuon;
  bool b_genMuon_isLooseMod;
  bool b_genMuon_isTightModNoIP, b_genMuon_isTightModIPxy, b_genMuon_isTightModIPz, b_genMuon_isTightModIPxyz;
  float b_genMuon_pTresolution;
  int b_genMuon_noRecHitGEM;
  int b_genMuon_noRecHitME0;
  
  float b_genMuon_poszPV0, b_genMuon_poszSimPV;
  float b_genMuon_poszPV01D, b_genMuon_poszPV01DBS, b_genMuon_poszPV04D, b_genMuon_poszPV04DBS, b_genMuon_poszPV0BS;
  float b_genMuon_poszLPVOfCand, b_genMuon_varposzCand, b_genMuon_poszMuon;
  float b_genMuon_pfIso03; float b_genMuon_pfIso04;
  float b_genMuon_pfIso03ChargedHadronPt, b_genMuon_pfIso03NeutralHadronEt;
  float b_genMuon_pfIso03PhotonEt, b_genMuon_pfIso03PUPt;
  float b_genMuon_pfIso04ChargedHadronPt, b_genMuon_pfIso04NeutralHadronEt;
  float b_genMuon_pfIso04PhotonEt, b_genMuon_pfIso04PUPt;
  float b_genMuon_TrkIso05; float b_genMuon_TrkIso03;
  float b_genMuon_puppiIsoWithLep, b_genMuon_puppiIsoWithoutLep, b_genMuon_puppiIsoCombined;
  float b_genMuon_puppiIsoWithLep03, b_genMuon_puppiIsoWithoutLep03, b_genMuon_puppiIsoCombined03;
  float b_genMuon_puppiIsoWithLep05, b_genMuon_puppiIsoWithoutLep05, b_genMuon_puppiIsoCombined05;
  float b_genMuon_puppiIsoWithLep03ChargedHadron, b_genMuon_puppiIsoWithLep03NeutralHadron, b_genMuon_puppiIsoWithLep03Photon;
  float b_genMuon_puppiIsoWithLep04ChargedHadron, b_genMuon_puppiIsoWithLep04NeutralHadron, b_genMuon_puppiIsoWithLep04Photon;
  float b_genMuon_puppiIsoWithLep05ChargedHadron, b_genMuon_puppiIsoWithLep05NeutralHadron, b_genMuon_puppiIsoWithLep05Photon;
  float b_genMuon_puppiIsoWithoutLep03ChargedHadron, b_genMuon_puppiIsoWithoutLep03NeutralHadron, b_genMuon_puppiIsoWithoutLep03Photon;
  float b_genMuon_puppiIsoWithoutLep04ChargedHadron, b_genMuon_puppiIsoWithoutLep04NeutralHadron, b_genMuon_puppiIsoWithoutLep04Photon;
  float b_genMuon_puppiIsoWithoutLep05ChargedHadron, b_genMuon_puppiIsoWithoutLep05NeutralHadron, b_genMuon_puppiIsoWithoutLep05Photon;
  int b_genMuon_puppiIsoNumOfCands;
  int b_genMuon_puppiIsoNumOfCandsInR05;
  int b_genMuon_puppiIsoNumOfCandsOutR05;
  int b_genMuon_puppiIsoNumOfCandsInR05OT;
  int b_genMuon_puppiIsoNumOfCandsInR05CH, b_genMuon_puppiIsoNumOfCandsInR05CHApp;
  int b_genMuon_puppiIsoNumOfCandsInR05NH, b_genMuon_puppiIsoNumOfCandsInR05NHApp;
  int b_genMuon_puppiIsoNumOfCandsInR05PH, b_genMuon_puppiIsoNumOfCandsInR05PHApp;
  int b_genMuon_numberOfValidMuonGEMHits, b_genMuon_numberOfValidMuonME0Hits;
  float b_genMuon_tmva_bdt, b_genMuon_tmva_mlp;

  float b_genMuon_deltaXME0, b_genMuon_deltaYME0;
  float b_genMuon_deltaXErrME0, b_genMuon_deltaYErrME0;
  float b_genMuon_deltaDXDZME0, b_genMuon_deltaDYDZME0;
  float b_genMuon_deltaDXDZErrME0, b_genMuon_deltaDYDZErrME0;

  TLorentzVector b_recoMuon;
  bool b_recoMuon_signal;
  bool b_recoMuon_isTightOptimized, b_recoMuon_isTightCustom, b_recoMuon_isTight, b_recoMuon_isMedium, b_recoMuon_isLoose;
  bool b_recoMuon_isME0Muon, b_recoMuon_isGEMMuon;
  bool b_recoMuon_isGlobalMuon, b_recoMuon_isStandAloneMuon;
  bool b_recoMuon_isLooseMod;
  bool b_recoMuon_isTightModNoIP, b_recoMuon_isTightModIPxy, b_recoMuon_isTightModIPz, b_recoMuon_isTightModIPxyz;
  int b_recoMuon_noChamberMatch;
  int b_recoMuon_noSegment, b_recoMuon_noSegmentDT, b_recoMuon_noSegmentCSC, b_recoMuon_noSegmentRPC, b_recoMuon_noSegmentGEM, b_recoMuon_noSegmentME0;
  int b_recoMuon_noRecHitGEM, b_recoMuon_noRecHitME0;

  float b_recoMuon_edgeXME0, b_recoMuon_edgeYME0;
  float b_recoMuon_chamberMatchXME0, b_recoMuon_chamberMatchYME0;
  float b_recoMuon_chamberMatchXErrME0, b_recoMuon_chamberMatchYErrME0;
  float b_recoMuon_chamberMatchDXDZME0, b_recoMuon_chamberMatchDYDZME0;
  float b_recoMuon_chamberMatchDXDZErrME0, b_recoMuon_chamberMatchDYDZErrME0;
  float b_recoMuon_distance, b_recoMuon_distErr;
  float b_recoMuon_segmentMatchXME0, b_recoMuon_segmentMatchYME0;
  float b_recoMuon_segmentMatchXErrME0, b_recoMuon_segmentMatchYErrME0;
  float b_recoMuon_segmentMatchDXDZME0, b_recoMuon_segmentMatchDYDZME0;
  float b_recoMuon_segmentMatchDXDZErrME0, b_recoMuon_segmentMatchDYDZErrME0;

  float b_recoMuon_deltaXME0, b_recoMuon_deltaYME0;
  float b_recoMuon_deltaXDivBySegErrME0, b_recoMuon_deltaYDivBySegErrME0, b_recoMuon_deltaXDivByChamErrME0, b_recoMuon_deltaYDivByChamErrME0;
  float b_recoMuon_deltaDXDZME0, b_recoMuon_deltaDYDZME0;
  float b_recoMuon_deltaDXDZDivBySegErrME0, b_recoMuon_deltaDYDZDivBySegErrME0, b_recoMuon_deltaDXDZDivByChamErrME0, b_recoMuon_deltaDYDZDivByChamErrME0;

  bool b_recoMuon_global; bool b_recoMuon_pf;
  float b_recoMuon_chi2pos; float b_recoMuon_trkKink; float b_recoMuon_segcompati;
  float b_recoMuon_chi2; int b_recoMuon_nglobalhits; int b_recoMuon_nstations;
  float b_recoMuon_trackdxy; float b_recoMuon_trackdz;
  int b_recoMuon_ninnerhits; float b_recoMuon_trackerlayers;
  int b_recoMuon_pdgId;
  float b_recoMuon_poszPV0, b_recoMuon_poszSimPV;
  float b_recoMuon_poszPV01D, b_recoMuon_poszPV01DBS, b_recoMuon_poszPV04D, b_recoMuon_poszPV04DBS, b_recoMuon_poszPV0BS;
  float b_recoMuon_poszLPVOfCand, b_recoMuon_poszMuon;
  float b_recoMuon_PFIso04; float b_recoMuon_PFIso03;
  float b_recoMuon_PFIso03ChargedHadronPt, b_recoMuon_PFIso03NeutralHadronEt;
  float b_recoMuon_PFIso03PhotonEt, b_recoMuon_PFIso03PUPt;
  float b_recoMuon_PFIso04ChargedHadronPt, b_recoMuon_PFIso04NeutralHadronEt;
  float b_recoMuon_PFIso04PhotonEt, b_recoMuon_PFIso04PUPt;
  float b_recoMuon_TrkIso05; float b_recoMuon_TrkIso03;
  float b_recoMuon_puppiIsoWithLep, b_recoMuon_puppiIsoWithoutLep, b_recoMuon_puppiIsoCombined;
  float b_recoMuon_puppiIsoWithLep03, b_recoMuon_puppiIsoWithoutLep03, b_recoMuon_puppiIsoCombined03;
  float b_recoMuon_puppiIsoWithLep05, b_recoMuon_puppiIsoWithoutLep05, b_recoMuon_puppiIsoCombined05;
  float b_recoMuon_puppiIsoWithLep03ChargedHadron, b_recoMuon_puppiIsoWithLep03NeutralHadron, b_recoMuon_puppiIsoWithLep03Photon;
  float b_recoMuon_puppiIsoWithLep04ChargedHadron, b_recoMuon_puppiIsoWithLep04NeutralHadron, b_recoMuon_puppiIsoWithLep04Photon;
  float b_recoMuon_puppiIsoWithLep05ChargedHadron, b_recoMuon_puppiIsoWithLep05NeutralHadron, b_recoMuon_puppiIsoWithLep05Photon;
  float b_recoMuon_puppiIsoWithoutLep03ChargedHadron, b_recoMuon_puppiIsoWithoutLep03NeutralHadron, b_recoMuon_puppiIsoWithoutLep03Photon;
  float b_recoMuon_puppiIsoWithoutLep04ChargedHadron, b_recoMuon_puppiIsoWithoutLep04NeutralHadron, b_recoMuon_puppiIsoWithoutLep04Photon;
  float b_recoMuon_puppiIsoWithoutLep05ChargedHadron, b_recoMuon_puppiIsoWithoutLep05NeutralHadron, b_recoMuon_puppiIsoWithoutLep05Photon;
  int b_recoMuon_puppiIsoNumOfCands;
  int b_recoMuon_puppiIsoNumOfCandsInR05;
  int b_recoMuon_puppiIsoNumOfCandsOutR05;
  int b_recoMuon_puppiIsoNumOfCandsInR05OT;
  int b_recoMuon_puppiIsoNumOfCandsInR05CH, b_recoMuon_puppiIsoNumOfCandsInR05CHApp;
  int b_recoMuon_puppiIsoNumOfCandsInR05NH, b_recoMuon_puppiIsoNumOfCandsInR05NHApp;
  int b_recoMuon_puppiIsoNumOfCandsInR05PH, b_recoMuon_puppiIsoNumOfCandsInR05PHApp;
  bool b_recoMuon_isMuon;
  int b_recoMuon_numberOfValidMuonGEMHits, b_recoMuon_numberOfValidMuonME0Hits;

  float b_recoMuon_tmva_bdt, b_recoMuon_tmva_mlp;
   
  std::vector<double> tmvaValues;
  float b_tmva_bdt; float b_tmva_mlp;
  ReadBDT* bdt_;
  ReadMLP* mlp_;
  
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
  genttree_->Branch("nvertex", &b_nvertex, "nvertex/I");
  genttree_->Branch("genMuon", "TLorentzVector", &b_genMuon);
  genttree_->Branch("genMuon_isTightOptimized", &b_genMuon_isTightOptimized, "genMuon_isTightOptimized/O");
  genttree_->Branch("genMuon_isTightCustom", &b_genMuon_isTightCustom, "genMuon_isTightCustom/O");
  genttree_->Branch("genMuon_isTight", &b_genMuon_isTight, "genMuon_isTight/O");
  genttree_->Branch("genMuon_isMedium", &b_genMuon_isMedium, "genMuon_isMedium/O");
  genttree_->Branch("genMuon_isLoose", &b_genMuon_isLoose, "genMuon_isLoose/O");
  genttree_->Branch("genMuon_noRecHitGEM", &b_genMuon_noRecHitGEM, "genMuon_noRecHitGEM/I");
  genttree_->Branch("genMuon_noRecHitME0", &b_genMuon_noRecHitME0, "genMuon_noRecHitME0/I");
  genttree_->Branch("genMuon_isME0Muon", &b_genMuon_isME0Muon, "genMuon_isME0Muon/O");
  genttree_->Branch("genMuon_isGEMMuon", &b_genMuon_isGEMMuon, "genMuon_isGEMMuon/O");
  genttree_->Branch("genMuon_isMuon", &b_genMuon_isMuon, "genMuon_isMuon/O");
  genttree_->Branch("genMuon_isGlobalMuon", &b_genMuon_isGlobalMuon, "genMuon_isGlobalMuon/O");
  genttree_->Branch("genMuon_isStandAloneMuon", &b_genMuon_isStandAloneMuon, "genMuon_isStandAloneMuon/O");
  genttree_->Branch("genMuon_isLooseMod", &b_genMuon_isLooseMod, "genMuon_isLooseMod/O");
  genttree_->Branch("genMuon_isTightModNoIP", &b_genMuon_isTightModNoIP, "genMuon_isTightModNoIP/O");
  genttree_->Branch("genMuon_isTightModIPxy", &b_genMuon_isTightModIPxy, "genMuon_isTightModIPxy/O");
  genttree_->Branch("genMuon_isTightModIPz", &b_genMuon_isTightModIPz, "genMuon_isTightModIPz/O");
  genttree_->Branch("genMuon_isTightModIPxyz", &b_genMuon_isTightModIPxyz, "genMuon_isTightModIPxyz/O");
  genttree_->Branch("genMuon_pTresolution", &b_genMuon_pTresolution, "genMuon_pTresolution/F");
  
  genttree_->Branch("genMuon_poszPV0",&b_genMuon_poszPV0,"genMuon_poszPV0/F");
  genttree_->Branch("genMuon_poszPV01D",&b_genMuon_poszPV01D,"genMuon_poszPV01D/F");
  genttree_->Branch("genMuon_poszPV01DBS",&b_genMuon_poszPV01DBS,"genMuon_poszPV01DBS/F");
  genttree_->Branch("genMuon_poszPV04D",&b_genMuon_poszPV04D,"genMuon_poszPV04D/F");
  genttree_->Branch("genMuon_poszPV04DBS",&b_genMuon_poszPV04DBS,"genMuon_poszPV04DBS/F");
  genttree_->Branch("genMuon_poszPV0BS",&b_genMuon_poszPV0BS,"genMuon_poszPV0BS/F");
  genttree_->Branch("genMuon_poszSimPV",&b_genMuon_poszSimPV,"genMuon_poszSimPV/F");
  genttree_->Branch("genMuon_poszLPVOfCand",&b_genMuon_poszLPVOfCand,"genMuon_poszLPVOfCand/F");
  genttree_->Branch("genMuon_varposzCand",&b_genMuon_varposzCand,"genMuon_varposzCand/F");
  genttree_->Branch("genMuon_poszMuon",&b_genMuon_poszMuon,"genMuon_poszMuon/F");
  genttree_->Branch("genMuon_TrkIsolation03",&b_genMuon_TrkIso03,"genMuon_TrkIsolation03/F");
  genttree_->Branch("genMuon_TrkIsolation05",&b_genMuon_TrkIso05,"genMuon_TrkIsolation05/F");
  genttree_->Branch("genMuon_PFIsolation03",&b_genMuon_pfIso03,"genMuon_PFIsolation03/F");
  genttree_->Branch("genMuon_PFIsolation04",&b_genMuon_pfIso04,"genMuon_PFIsolation04/F");
  genttree_->Branch("genMuon_PFIso03ChargedHadronPt",&b_genMuon_pfIso03ChargedHadronPt,"genMuon_PFIso03ChargedHadronPt/F");
  genttree_->Branch("genMuon_PFIso03NeutralHadronEt",&b_genMuon_pfIso03NeutralHadronEt,"genMuon_PFIso03NeutralHadronEt/F");
  genttree_->Branch("genMuon_PFIso03PhotonEt",&b_genMuon_pfIso03PhotonEt,"genMuon_PFIso03PhotonEt/F");
  genttree_->Branch("genMuon_PFIso03PUPt",&b_genMuon_pfIso03PUPt,"genMuon_PFIso03PUPt/F");
  genttree_->Branch("genMuon_PFIso04ChargedHadronPt",&b_genMuon_pfIso04ChargedHadronPt,"genMuon_PFIso04ChargedHadronPt/F");
  genttree_->Branch("genMuon_PFIso04NeutralHadronEt",&b_genMuon_pfIso04NeutralHadronEt,"genMuon_PFIso04NeutralHadronEt/F");
  genttree_->Branch("genMuon_PFIso04PhotonEt",&b_genMuon_pfIso04PhotonEt,"genMuon_PFIso04PhotonEt/F");
  genttree_->Branch("genMuon_PFIso04PUPt",&b_genMuon_pfIso04PUPt,"genMuon_PFIso04PUPt/F");
  genttree_->Branch("genMuon_puppiIsoWithLep",&b_genMuon_puppiIsoWithLep,"genMuon_puppiIsoWithLep/F");
  genttree_->Branch("genMuon_puppiIsoWithoutLep",&b_genMuon_puppiIsoWithoutLep,"genMuon_puppiIsoWithoutLep/F");
  genttree_->Branch("genMuon_puppiIsoCombined",&b_genMuon_puppiIsoCombined,"genMuon_puppiIsoCombined/F");
  genttree_->Branch("genMuon_puppiIsoWithLep03",&b_genMuon_puppiIsoWithLep03,"genMuon_puppiIsoWithLep03/F");
  genttree_->Branch("genMuon_puppiIsoWithoutLep03",&b_genMuon_puppiIsoWithoutLep03,"genMuon_puppiIsoWithoutLep03/F");
  genttree_->Branch("genMuon_puppiIsoCombined03",&b_genMuon_puppiIsoCombined03,"genMuon_puppiIsoCombined03/F");
  genttree_->Branch("genMuon_puppiIsoWithLep05",&b_genMuon_puppiIsoWithLep05,"genMuon_puppiIsoWithLep05/F");
  genttree_->Branch("genMuon_puppiIsoWithoutLep05",&b_genMuon_puppiIsoWithoutLep05,"genMuon_puppiIsoWithoutLep05/F");
  genttree_->Branch("genMuon_puppiIsoCombined05",&b_genMuon_puppiIsoCombined05,"genMuon_puppiIsoCombined05/F");
  genttree_->Branch("genMuon_puppiIsoWithLep03ChargedHadron",&b_genMuon_puppiIsoWithLep03ChargedHadron,"genMuon_puppiIsoWithLep03ChargedHadron/F");
  genttree_->Branch("genMuon_puppiIsoWithLep03NeutralHadron",&b_genMuon_puppiIsoWithLep03NeutralHadron,"genMuon_puppiIsoWithLep03NeutralHadron/F");
  genttree_->Branch("genMuon_puppiIsoWithLep03Photon",&b_genMuon_puppiIsoWithLep03Photon,"genMuon_puppiIsoWithLep03Photon/F");
  genttree_->Branch("genMuon_puppiIsoWithLep04ChargedHadron",&b_genMuon_puppiIsoWithLep04ChargedHadron,"genMuon_puppiIsoWithLep04ChargedHadron/F");
  genttree_->Branch("genMuon_puppiIsoWithLep04NeutralHadron",&b_genMuon_puppiIsoWithLep04NeutralHadron,"genMuon_puppiIsoWithLep04NeutralHadron/F");
  genttree_->Branch("genMuon_puppiIsoWithLep04Photon",&b_genMuon_puppiIsoWithLep04Photon,"genMuon_puppiIsoWithLep04Photon/F");
  genttree_->Branch("genMuon_puppiIsoWithLep05ChargedHadron",&b_genMuon_puppiIsoWithLep05ChargedHadron,"genMuon_puppiIsoWithLep05ChargedHadron/F");
  genttree_->Branch("genMuon_puppiIsoWithLep05NeutralHadron",&b_genMuon_puppiIsoWithLep05NeutralHadron,"genMuon_puppiIsoWithLep05NeutralHadron/F");
  genttree_->Branch("genMuon_puppiIsoWithLep05Photon",&b_genMuon_puppiIsoWithLep05Photon,"genMuon_puppiIsoWithLep05Photon/F");
  genttree_->Branch("genMuon_puppiIsoWithoutLep03ChargedHadron",&b_genMuon_puppiIsoWithoutLep03ChargedHadron,"genMuon_puppiIsoWithoutLep03ChargedHadron/F");
  genttree_->Branch("genMuon_puppiIsoWithoutLep03NeutralHadron",&b_genMuon_puppiIsoWithoutLep03NeutralHadron,"genMuon_puppiIsoWithoutLep03NeutralHadron/F");
  genttree_->Branch("genMuon_puppiIsoWithoutLep03Photon",&b_genMuon_puppiIsoWithoutLep03Photon,"genMuon_puppiIsoWithoutLep03Photon/F");
  genttree_->Branch("genMuon_puppiIsoWithoutLep04ChargedHadron",&b_genMuon_puppiIsoWithoutLep04ChargedHadron,"genMuon_puppiIsoWithoutLep04ChargedHadron/F");
  genttree_->Branch("genMuon_puppiIsoWithoutLep04NeutralHadron",&b_genMuon_puppiIsoWithoutLep04NeutralHadron,"genMuon_puppiIsoWithoutLep04NeutralHadron/F");
  genttree_->Branch("genMuon_puppiIsoWithoutLep04Photon",&b_genMuon_puppiIsoWithoutLep04Photon,"genMuon_puppiIsoWithoutLep04Photon/F");
  genttree_->Branch("genMuon_puppiIsoWithoutLep05ChargedHadron",&b_genMuon_puppiIsoWithoutLep05ChargedHadron,"genMuon_puppiIsoWithoutLep05ChargedHadron/F");
  genttree_->Branch("genMuon_puppiIsoWithoutLep05NeutralHadron",&b_genMuon_puppiIsoWithoutLep05NeutralHadron,"genMuon_puppiIsoWithoutLep05NeutralHadron/F");
  genttree_->Branch("genMuon_puppiIsoWithoutLep05Photon",&b_genMuon_puppiIsoWithoutLep05Photon,"genMuon_puppiIsoWithoutLep05Photon/F");
  genttree_->Branch("genMuon_puppiIsoNumOfCands",&b_genMuon_puppiIsoNumOfCands,"genMuon_puppiIsoNumOfCands/I");
  genttree_->Branch("genMuon_puppiIsoNumOfCandsInR05",&b_genMuon_puppiIsoNumOfCandsInR05,"genMuon_puppiIsoNumOfCandsInR05/I");
  genttree_->Branch("genMuon_puppiIsoNumOfCandsOutR05",&b_genMuon_puppiIsoNumOfCandsOutR05,"genMuon_puppiIsoNumOfCandsOutR05/I");
  genttree_->Branch("genMuon_puppiIsoNumOfCandsInR05OT",&b_genMuon_puppiIsoNumOfCandsInR05OT,"genMuon_puppiIsoNumOfCandsInR05OT/I");
  genttree_->Branch("genMuon_puppiIsoNumOfCandsInR05CH",&b_genMuon_puppiIsoNumOfCandsInR05CH,"genMuon_puppiIsoNumOfCandsInR05CH/I");
  genttree_->Branch("genMuon_puppiIsoNumOfCandsInR05NH",&b_genMuon_puppiIsoNumOfCandsInR05NH,"genMuon_puppiIsoNumOfCandsInR05NH/I");
  genttree_->Branch("genMuon_puppiIsoNumOfCandsInR05PH",&b_genMuon_puppiIsoNumOfCandsInR05PH,"genMuon_puppiIsoNumOfCandsInR05PH/I");
  genttree_->Branch("genMuon_puppiIsoNumOfCandsInR05CHApp",&b_genMuon_puppiIsoNumOfCandsInR05CHApp,"genMuon_puppiIsoNumOfCandsInR05CHApp/I");
  genttree_->Branch("genMuon_puppiIsoNumOfCandsInR05NHApp",&b_genMuon_puppiIsoNumOfCandsInR05NHApp,"genMuon_puppiIsoNumOfCandsInR05NHApp/I");
  genttree_->Branch("genMuon_puppiIsoNumOfCandsInR05PHApp",&b_genMuon_puppiIsoNumOfCandsInR05PHApp,"genMuon_puppiIsoNumOfCandsInR05PHApp/I");
  genttree_->Branch("genMuon_numberOfValidMuonGEMHits",&b_genMuon_numberOfValidMuonGEMHits,"genMuon_numberOfValidMuonGEMHits/I");
  genttree_->Branch("genMuon_numberOfValidMuonME0Hits",&b_genMuon_numberOfValidMuonME0Hits,"genMuon_numberOfValidMuonME0Hits/I");
  genttree_->Branch("genMuon_signal", &b_genMuon_signal, "genMuon_signal/O");
  genttree_->Branch("genMuon_tmva_bdt", &b_genMuon_tmva_bdt, "genMuon_tmva_bdt/F");
  genttree_->Branch("genMuon_tmva_mlp", &b_genMuon_tmva_mlp, "genMuon_tmva_mlp/F");

  genttree_->Branch("genMuon_deltaXME0", &b_genMuon_deltaXME0, "genMuon_deltaXME0/F");
  genttree_->Branch("genMuon_deltaYME0", &b_genMuon_deltaYME0, "genMuon_deltaYME0/F");
  genttree_->Branch("genMuon_deltaXErrME0", &b_genMuon_deltaXErrME0, "genMuon_deltaXErrME0/F");
  genttree_->Branch("genMuon_deltaYErrME0", &b_genMuon_deltaYErrME0, "genMuon_deltaYErrME0/F");
  genttree_->Branch("genMuon_deltaDXDZME0", &b_genMuon_deltaDXDZME0, "genMuon_deltaDXDZME0/F");
  genttree_->Branch("genMuon_deltaDYDZME0", &b_genMuon_deltaDYDZME0, "genMuon_deltaDYDZME0/F");
  genttree_->Branch("genMuon_deltaDXDZErrME0", &b_genMuon_deltaDXDZErrME0, "genMuon_deltaDXDZErrME0/F");
  genttree_->Branch("genMuon_deltaDYDZErrME0", &b_genMuon_deltaDYDZErrME0, "genMuon_deltaDYDZErrME0/F");

  recottree_ = fs->make<TTree>("reco", "reco");
  recottree_->Branch("nvertex", &b_nvertex, "nvertex/I");
  recottree_->Branch("recoMuon", "TLorentzVector", &b_recoMuon);  
  recottree_->Branch("recoMuon_pdgId", &b_recoMuon_pdgId, "recoMuon_pdgId/I");
  recottree_->Branch("recoMuon_signal", &b_recoMuon_signal, "recoMuon_signal/O");
  recottree_->Branch("recoMuon_isTightOptimized", &b_recoMuon_isTightOptimized, "recoMuon_isTightOptimized/O");
  recottree_->Branch("recoMuon_isTightCustom", &b_recoMuon_isTightCustom, "recoMuon_isTightCustom/O");
  recottree_->Branch("recoMuon_isTight", &b_recoMuon_isTight, "recoMuon_isTight/O");
  recottree_->Branch("recoMuon_isMedium", &b_recoMuon_isMedium, "recoMuon_isMedium/O");
  recottree_->Branch("recoMuon_isLoose", &b_recoMuon_isLoose, "recoMuon_isLoose/O");
  recottree_->Branch("recoMuon_isME0Muon", &b_recoMuon_isME0Muon, "recoMuon_isME0Muon/O");
  recottree_->Branch("recoMuon_isGEMMuon", &b_recoMuon_isGEMMuon, "recoMuon_isGEMMuon/O");
  recottree_->Branch("recoMuon_isMuon", &b_recoMuon_isMuon, "recoMuon_isMuon/O");
  recottree_->Branch("recoMuon_isGlobalMuon", &b_recoMuon_isGlobalMuon, "recoMuon_isGlobalMuon/O");
  recottree_->Branch("recoMuon_isStandAloneMuon", &b_recoMuon_isStandAloneMuon, "recoMuon_isStandAloneMuon/O");
  recottree_->Branch("recoMuon_isLooseMod", &b_recoMuon_isLooseMod, "recoMuon_isLooseMod/O");
  recottree_->Branch("recoMuon_isTightModNoIP", &b_recoMuon_isTightModNoIP, "recoMuon_isTightModNoIP/O");
  recottree_->Branch("recoMuon_isTightModIPxy", &b_recoMuon_isTightModIPxy, "recoMuon_isTightModIPxy/O");
  recottree_->Branch("recoMuon_isTightModIPz", &b_recoMuon_isTightModIPz, "recoMuon_isTightModIPz/O");
  recottree_->Branch("recoMuon_isTightModIPxyz", &b_recoMuon_isTightModIPxyz, "recoMuon_isTightModIPxyz/O");
  recottree_->Branch("recoMuon_numberOfValidMuonGEMHits",&b_recoMuon_numberOfValidMuonGEMHits,"recoMuon_numberOfValidMuonGEMHits/I");
  recottree_->Branch("recoMuon_numberOfValidMuonME0Hits",&b_recoMuon_numberOfValidMuonME0Hits,"recoMuon_numberOfValidMuonME0Hits/I");

  recottree_->Branch("recoMuon_poszPV0",&b_recoMuon_poszPV0,"recoMuon_poszPV0/F");
  recottree_->Branch("recoMuon_poszPV01D",&b_recoMuon_poszPV01D,"recoMuon_poszPV01D/F");
  recottree_->Branch("recoMuon_poszPV01DBS",&b_recoMuon_poszPV01DBS,"recoMuon_poszPV01DBS/F");
  recottree_->Branch("recoMuon_poszPV04D",&b_recoMuon_poszPV04D,"recoMuon_poszPV04D/F");
  recottree_->Branch("recoMuon_poszPV04DBS",&b_recoMuon_poszPV04DBS,"recoMuon_poszPV04DBS/F");
  recottree_->Branch("recoMuon_poszPV0BS",&b_recoMuon_poszPV0BS,"recoMuon_poszPV0BS/F");
  recottree_->Branch("recoMuon_poszSimPV",&b_recoMuon_poszSimPV,"recoMuon_poszSimPV/F");
  recottree_->Branch("recoMuon_poszLPVOfCand",&b_recoMuon_poszLPVOfCand,"recoMuon_poszLPVOfCand/F");
  recottree_->Branch("recoMuon_poszMuon",&b_recoMuon_poszMuon,"recoMuon_poszMuon/F");
  recottree_->Branch("recoMuon_TrkIsolation03",&b_recoMuon_TrkIso03,"recoMuon_TrkIsolation03/F");
  recottree_->Branch("recoMuon_TrkIsolation05",&b_recoMuon_TrkIso05,"recoMuon_TrkIsolation05/F");
  recottree_->Branch("recoMuon_PFIsolation04",&b_recoMuon_PFIso04,"recoMuon_PFIsolation04/F");
  recottree_->Branch("recoMuon_PFIsolation03",&b_recoMuon_PFIso03,"recoMuon_PFIsolation03/F");
  recottree_->Branch("recoMuon_PFIso03ChargedHadronPt",&b_recoMuon_PFIso03ChargedHadronPt,"recoMuon_PFIso03ChargedHadronPt/F");
  recottree_->Branch("recoMuon_PFIso03NeutralHadronEt",&b_recoMuon_PFIso03NeutralHadronEt,"recoMuon_PFIso03NeutralHadronEt/F");
  recottree_->Branch("recoMuon_PFIso03PhotonEt",&b_recoMuon_PFIso03PhotonEt,"recoMuon_PFIso03PhotonEt/F");
  recottree_->Branch("recoMuon_PFIso03PUPt",&b_recoMuon_PFIso03PUPt,"recoMuon_PFIso03PUPt/F");
  recottree_->Branch("recoMuon_PFIso04ChargedHadronPt",&b_recoMuon_PFIso04ChargedHadronPt,"recoMuon_PFIso04ChargedHadronPt/F");
  recottree_->Branch("recoMuon_PFIso04NeutralHadronEt",&b_recoMuon_PFIso04NeutralHadronEt,"recoMuon_PFIso04NeutralHadronEt/F");
  recottree_->Branch("recoMuon_PFIso04PhotonEt",&b_recoMuon_PFIso04PhotonEt,"recoMuon_PFIso04PhotonEt/F");
  recottree_->Branch("recoMuon_PFIso04PUPt",&b_recoMuon_PFIso04PUPt,"recoMuon_PFIso04PUPt/F");
  recottree_->Branch("recoMuon_puppiIsoWithLep",&b_recoMuon_puppiIsoWithLep,"recoMuon_puppiIsoWithLep/F");
  recottree_->Branch("recoMuon_puppiIsoWithoutLep",&b_recoMuon_puppiIsoWithoutLep,"recoMuon_puppiIsoWithoutLep/F");
  recottree_->Branch("recoMuon_puppiIsoCombined",&b_recoMuon_puppiIsoCombined,"recoMuon_puppiIsoCombined/F");
  recottree_->Branch("recoMuon_puppiIsoWithLep03",&b_recoMuon_puppiIsoWithLep03,"recoMuon_puppiIsoWithLep03/F");
  recottree_->Branch("recoMuon_puppiIsoWithoutLep03",&b_recoMuon_puppiIsoWithoutLep03,"recoMuon_puppiIsoWithoutLep03/F");
  recottree_->Branch("recoMuon_puppiIsoCombined03",&b_recoMuon_puppiIsoCombined03,"recoMuon_puppiIsoCombined03/F");
  recottree_->Branch("recoMuon_puppiIsoWithLep05",&b_recoMuon_puppiIsoWithLep05,"recoMuon_puppiIsoWithLep05/F");
  recottree_->Branch("recoMuon_puppiIsoWithoutLep05",&b_recoMuon_puppiIsoWithoutLep05,"recoMuon_puppiIsoWithoutLep05/F");
  recottree_->Branch("recoMuon_puppiIsoCombined05",&b_recoMuon_puppiIsoCombined05,"recoMuon_puppiIsoCombined05/F");
  recottree_->Branch("recoMuon_puppiIsoWithLep03ChargedHadron",&b_recoMuon_puppiIsoWithLep03ChargedHadron,"recoMuon_puppiIsoWithLep03ChargedHadron/F");
  recottree_->Branch("recoMuon_puppiIsoWithLep03NeutralHadron",&b_recoMuon_puppiIsoWithLep03NeutralHadron,"recoMuon_puppiIsoWithLep03NeutralHadron/F");
  recottree_->Branch("recoMuon_puppiIsoWithLep03Photon",&b_recoMuon_puppiIsoWithLep03Photon,"recoMuon_puppiIsoWithLep03Photon/F");
  recottree_->Branch("recoMuon_puppiIsoWithLep04ChargedHadron",&b_recoMuon_puppiIsoWithLep04ChargedHadron,"recoMuon_puppiIsoWithLep04ChargedHadron/F");
  recottree_->Branch("recoMuon_puppiIsoWithLep04NeutralHadron",&b_recoMuon_puppiIsoWithLep04NeutralHadron,"recoMuon_puppiIsoWithLep04NeutralHadron/F");
  recottree_->Branch("recoMuon_puppiIsoWithLep04Photon",&b_recoMuon_puppiIsoWithLep04Photon,"recoMuon_puppiIsoWithLep04Photon/F");
  recottree_->Branch("recoMuon_puppiIsoWithLep05ChargedHadron",&b_recoMuon_puppiIsoWithLep05ChargedHadron,"recoMuon_puppiIsoWithLep05ChargedHadron/F");
  recottree_->Branch("recoMuon_puppiIsoWithLep05NeutralHadron",&b_recoMuon_puppiIsoWithLep05NeutralHadron,"recoMuon_puppiIsoWithLep05NeutralHadron/F");
  recottree_->Branch("recoMuon_puppiIsoWithLep05Photon",&b_recoMuon_puppiIsoWithLep05Photon,"recoMuon_puppiIsoWithLep05Photon/F");
  recottree_->Branch("recoMuon_puppiIsoWithoutLep03ChargedHadron",&b_recoMuon_puppiIsoWithoutLep03ChargedHadron,"recoMuon_puppiIsoWithoutLep03ChargedHadron/F");
  recottree_->Branch("recoMuon_puppiIsoWithoutLep03NeutralHadron",&b_recoMuon_puppiIsoWithoutLep03NeutralHadron,"recoMuon_puppiIsoWithoutLep03NeutralHadron/F");
  recottree_->Branch("recoMuon_puppiIsoWithoutLep03Photon",&b_recoMuon_puppiIsoWithoutLep03Photon,"recoMuon_puppiIsoWithoutLep03Photon/F");
  recottree_->Branch("recoMuon_puppiIsoWithoutLep04ChargedHadron",&b_recoMuon_puppiIsoWithoutLep04ChargedHadron,"recoMuon_puppiIsoWithoutLep04ChargedHadron/F");
  recottree_->Branch("recoMuon_puppiIsoWithoutLep04NeutralHadron",&b_recoMuon_puppiIsoWithoutLep04NeutralHadron,"recoMuon_puppiIsoWithoutLep04NeutralHadron/F");
  recottree_->Branch("recoMuon_puppiIsoWithoutLep04Photon",&b_recoMuon_puppiIsoWithoutLep04Photon,"recoMuon_puppiIsoWithoutLep04Photon/F");
  recottree_->Branch("recoMuon_puppiIsoWithoutLep05ChargedHadron",&b_recoMuon_puppiIsoWithoutLep05ChargedHadron,"recoMuon_puppiIsoWithoutLep05ChargedHadron/F");
  recottree_->Branch("recoMuon_puppiIsoWithoutLep05NeutralHadron",&b_recoMuon_puppiIsoWithoutLep05NeutralHadron,"recoMuon_puppiIsoWithoutLep05NeutralHadron/F");
  recottree_->Branch("recoMuon_puppiIsoWithoutLep05Photon",&b_recoMuon_puppiIsoWithoutLep05Photon,"recoMuon_puppiIsoWithoutLep05Photon/F");
  recottree_->Branch("recoMuon_puppiIsoNumOfCands",&b_recoMuon_puppiIsoNumOfCands,"recoMuon_puppiIsoNumOfCands/I");
  recottree_->Branch("recoMuon_puppiIsoNumOfCandsInR05",&b_recoMuon_puppiIsoNumOfCandsInR05,"recoMuon_puppiIsoNumOfCandsInR05/I");
  recottree_->Branch("recoMuon_puppiIsoNumOfCandsOutR05",&b_recoMuon_puppiIsoNumOfCandsOutR05,"recoMuon_puppiIsoNumOfCandsOutR05/I");
  recottree_->Branch("recoMuon_puppiIsoNumOfCandsInR05OT",&b_recoMuon_puppiIsoNumOfCandsInR05OT,"recoMuon_puppiIsoNumOfCandsInR05OT/I");
  recottree_->Branch("recoMuon_puppiIsoNumOfCandsInR05CH",&b_recoMuon_puppiIsoNumOfCandsInR05CH,"recoMuon_puppiIsoNumOfCandsInR05CH/I");
  recottree_->Branch("recoMuon_puppiIsoNumOfCandsInR05NH",&b_recoMuon_puppiIsoNumOfCandsInR05NH,"recoMuon_puppiIsoNumOfCandsInR05NH/I");
  recottree_->Branch("recoMuon_puppiIsoNumOfCandsInR05PH",&b_recoMuon_puppiIsoNumOfCandsInR05PH,"recoMuon_puppiIsoNumOfCandsInR05PH/I");
  recottree_->Branch("recoMuon_puppiIsoNumOfCandsInR05CHApp",&b_recoMuon_puppiIsoNumOfCandsInR05CHApp,"recoMuon_puppiIsoNumOfCandsInR05CHApp/I");
  recottree_->Branch("recoMuon_puppiIsoNumOfCandsInR05NHApp",&b_recoMuon_puppiIsoNumOfCandsInR05NHApp,"recoMuon_puppiIsoNumOfCandsInR05NHApp/I");
  recottree_->Branch("recoMuon_puppiIsoNumOfCandsInR05PHApp",&b_recoMuon_puppiIsoNumOfCandsInR05PHApp,"recoMuon_puppiIsoNumOfCandsInR05PHApp/I");
  recottree_->Branch("recoMuon_noChamberMatch", &b_recoMuon_noChamberMatch, "recoMuon_noChamberMatch/I");
  recottree_->Branch("recoMuon_noSegment", &b_recoMuon_noSegment, "recoMuon_noSegment/I");
  recottree_->Branch("recoMuon_noSegmentDT", &b_recoMuon_noSegmentDT, "recoMuon_noSegmentDT/I");
  recottree_->Branch("recoMuon_noSegmentCSC", &b_recoMuon_noSegmentCSC, "recoMuon_noSegmentCSC/I");
  recottree_->Branch("recoMuon_noSegmentRPC", &b_recoMuon_noSegmentRPC, "recoMuon_noSegmentRPC/I");
  recottree_->Branch("recoMuon_noSegmentGEM", &b_recoMuon_noSegmentGEM, "recoMuon_noSegmentGEM/I");
  recottree_->Branch("recoMuon_noSegmentME0", &b_recoMuon_noSegmentME0, "recoMuon_noSegmentME0/I");

  recottree_->Branch("recoMuon_noRecHitGEM", &b_recoMuon_noRecHitGEM, "recoMuon_noRecHitGEM/I");
  recottree_->Branch("recoMuon_noRecHitME0", &b_recoMuon_noRecHitME0, "recoMuon_noRecHitME0/I");

  recottree_->Branch("recoMuon_edgeXME0", &b_recoMuon_edgeXME0, "recoMuon_edgeXME0/F");
  recottree_->Branch("recoMuon_edgeYME0", &b_recoMuon_edgeYME0, "recoMuon_edgeYME0/F");
  recottree_->Branch("recoMuon_chamberMatchXME0", &b_recoMuon_chamberMatchXME0, "recoMuon_chamberMatchXME0/F");
  recottree_->Branch("recoMuon_chamberMatchYME0", &b_recoMuon_chamberMatchYME0, "recoMuon_chamberMatchYME0/F");
  recottree_->Branch("recoMuon_chamberMatchXErrME0", &b_recoMuon_chamberMatchXErrME0, "recoMuon_chamberMatchXErrME0/F");
  recottree_->Branch("recoMuon_chamberMatchYErrME0", &b_recoMuon_chamberMatchYErrME0, "recoMuon_chamberMatchYErrME0/F");
  recottree_->Branch("recoMuon_chamberMatchDXDZME0", &b_recoMuon_chamberMatchDXDZME0, "recoMuon_chamberMatchDXDZME0/F");
  recottree_->Branch("recoMuon_chamberMatchDYDZME0", &b_recoMuon_chamberMatchDYDZME0, "recoMuon_chamberMatchDYDZME0/F");
  recottree_->Branch("recoMuon_chamberMatchDXDZErrME0", &b_recoMuon_chamberMatchDXDZErrME0, "recoMuon_chamberMatchDXDZErrME0/F");
  recottree_->Branch("recoMuon_chamberMatchDYDZErrME0", &b_recoMuon_chamberMatchDYDZErrME0, "recoMuon_chamberMatchDYDZErrME0/F");

  recottree_->Branch("recoMuon_segmentMatchXME0", &b_recoMuon_segmentMatchXME0, "recoMuon_segmentMatchXME0/F");
  recottree_->Branch("recoMuon_segmentMatchYME0", &b_recoMuon_segmentMatchYME0, "recoMuon_segmentMatchYME0/F");
  recottree_->Branch("recoMuon_segmentMatchXErrME0", &b_recoMuon_segmentMatchXErrME0, "recoMuon_segmentMatchXErrME0/F");
  recottree_->Branch("recoMuon_segmentMatchYErrME0", &b_recoMuon_segmentMatchYErrME0, "recoMuon_segmentMatchYErrME0/F");
  recottree_->Branch("recoMuon_segmentMatchDXDZME0", &b_recoMuon_segmentMatchDXDZME0, "recoMuon_segmentMatchDXDZME0/F");
  recottree_->Branch("recoMuon_segmentMatchDYDZME0", &b_recoMuon_segmentMatchDYDZME0, "recoMuon_segmentMatchDYDZME0/F");
  recottree_->Branch("recoMuon_segmentMatchDXDZErrME0", &b_recoMuon_segmentMatchDXDZErrME0, "recoMuon_segmentMatchDXDZErrME0/F");
  recottree_->Branch("recoMuon_segmentMatchDYDZErrME0", &b_recoMuon_segmentMatchDYDZErrME0, "recoMuon_segmentMatchDYDZErrME0/F");

  recottree_->Branch("recoMuon_deltaXME0", &b_recoMuon_deltaXME0, "recoMuon_deltaXME0/F");
  recottree_->Branch("recoMuon_deltaYME0", &b_recoMuon_deltaYME0, "recoMuon_deltaYME0/F");
  recottree_->Branch("recoMuon_deltaXDivBySegErrME0", &b_recoMuon_deltaXDivBySegErrME0, "recoMuon_deltaXDivBySegErrME0/F");
  recottree_->Branch("recoMuon_deltaYDivBySegErrME0", &b_recoMuon_deltaYDivBySegErrME0, "recoMuon_deltaYDivBySegErrME0/F");
  recottree_->Branch("recoMuon_deltaXDivByChamErrME0", &b_recoMuon_deltaXDivByChamErrME0, "recoMuon_deltaXDivByChamErrME0/F");
  recottree_->Branch("recoMuon_deltaYDivByChamErrME0", &b_recoMuon_deltaYDivByChamErrME0, "recoMuon_deltaYDivByChamErrME0/F");

  recottree_->Branch("recoMuon_deltaDXDZME0", &b_recoMuon_deltaDXDZME0, "recoMuon_deltaDXDZME0/F");
  recottree_->Branch("recoMuon_deltaDYDZME0", &b_recoMuon_deltaDYDZME0, "recoMuon_deltaDYDZME0/F");
  recottree_->Branch("recoMuon_deltaDXDZDivBySegErrME0", &b_recoMuon_deltaDXDZDivBySegErrME0, "recoMuon_deltaDXDZDivBySegErrME0/F");
  recottree_->Branch("recoMuon_deltaDYDZDivBySegErrME0", &b_recoMuon_deltaDYDZDivBySegErrME0, "recoMuon_deltaDYDZDivBySegErrME0/F");
  recottree_->Branch("recoMuon_deltaDXDZDivByChamErrME0", &b_recoMuon_deltaDXDZDivByChamErrME0, "recoMuon_deltaDXDZDivByChamErrME0/F");
  recottree_->Branch("recoMuon_deltaDYDZDivByChamErrME0", &b_recoMuon_deltaDYDZDivByChamErrME0, "recoMuon_deltaDYDZDivByChamErrME0/F");

  recottree_->Branch("recoMuon_distance", &b_recoMuon_distance, "recoMuon_distance/F");
  recottree_->Branch("recoMuon_distErr", &b_recoMuon_distErr, "recoMuon_distErr/F");

  recottree_->Branch("recoMuon_isGlobalMuon", &b_recoMuon_global, "recoMuon_isGlobalMuon/O");
  recottree_->Branch("recoMuon_isPFMuon", &b_recoMuon_pf, "recoMuon_isPFMuon/O");
  recottree_->Branch("recoMuon_normalizedChi2", &b_recoMuon_chi2, "recoMuon_normalizedChi2/F");
  recottree_->Branch("recoMuon_chi2LocalPosition", &b_recoMuon_chi2pos, "recoMuon_chi2LocalPosition/F");
  recottree_->Branch("recoMuon_trkKink", &b_recoMuon_trkKink, "recoMuon_trkKink/F");
  recottree_->Branch("recoMuon_segmentCompatibility", &b_recoMuon_segcompati, "recoMuon_segmentCompatibility/F");
  recottree_->Branch("recoMuon_numberOfValidMuonHits", &b_recoMuon_nglobalhits, "recoMuon_numberOfValidMuonHits/I");
  recottree_->Branch("recoMuon_numberOfMatchedStations", &b_recoMuon_nstations, "recoMuon_numberOfMatchedStations/I");
  recottree_->Branch("recoMuon_pv0pos_dxy", &b_recoMuon_trackdxy, "recoMuon_pv0pos_dxy/F");
  recottree_->Branch("recoMuon_pv0pos_dz", &b_recoMuon_trackdz, "recoMuon_pv0pos_dz/F");
  recottree_->Branch("recoMuon_numberOfValidPixelHits", &b_recoMuon_ninnerhits, "recoMuon_numberOfValidPixelHits/I");
  recottree_->Branch("recoMuon_trackerLayersWithMeasurement", &b_recoMuon_trackerlayers, "recoMuon_trackerLayersWithMeasurement/F");
  recottree_->Branch("recoMuon_tmva_bdt", &b_recoMuon_tmva_bdt, "recoMuon_tmva_bdt/F");
  recottree_->Branch("recoMuon_tmva_mlp", &b_recoMuon_tmva_mlp, "recoMuon_tmva_mlp/F");

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
  
  if (vertices->empty()) { cout << "noPV" << endl; return; }
  auto pv0 = vertices->front();
  b_nvertex = vertices->size();

  Handle<std::vector<SimVertex> > simVertexCollection;
  iEvent.getByToken(simVertexToken_, simVertexCollection);
  const SimVertex simPVh = *(simVertexCollection->begin());
  h_vertex->Fill(pv0.position().z() - simPVh.position().z());
  
  h_vertexBS->Fill(  verticesBS->front().position().z()   - simPVh.position().z());
  h_vertex1D->Fill(  vertices1D->front().position().z()   - simPVh.position().z());
  h_vertex1DBS->Fill(vertices1DBS->front().position().z() - simPVh.position().z());
  h_vertex4D->Fill(  vertices4D->front().position().z()   - simPVh.position().z());
  h_vertex4DBS->Fill(vertices4DBS->front().position().z() - simPVh.position().z());
  
  Handle<TrackingParticleCollection> simHandle;
  iEvent.getByToken(simToken_, simHandle);

  Handle<View<Muon> > muonHandle;
  iEvent.getByToken(muonToken_, muonHandle);

  edm::Handle<std::vector< pat::PackedCandidate> > Candidates_Collection ;
  iEvent.getByToken( tokenPackedCandidate , Candidates_Collection ) ; 
  const std::vector<pat::PackedCandidate> * candidates = Candidates_Collection.product();
  
  reco::MuonToSimCollection muonToSimColl; reco::SimToMuonCollection simToMuonColl;  
  Handle<MuonToTrackingParticleAssociator> associatorBase;
  iEvent.getByToken(muAssocToken_, associatorBase);
  MuonToTrackingParticleAssociator const* assoByHits = associatorBase.product();
  //reco::InnerTk or reco::GlobalTk
  assoByHits->associateMuons(muonToSimColl, simToMuonColl, muonHandle, reco::GlobalTk, simHandle);
  //assoByHits->associateMuons(muonToSimColl, simToMuonColl, muonHandle, reco::InnerTk, simHandle); // For Only GW
  
  vector<const Muon*> signalMuons; signalMuons.clear();
  
  // The following few lines are for finding the vertex of candidates with PUPPI weight non-zero
  double dSumPosZVtxCand;
  double dSumSqrPosZVtxCand;
  int nNumCandCH;
  
  double dPosZVtxCand, dVarPosZVtxCand;
  
  dSumPosZVtxCand = 0.0;
  dSumSqrPosZVtxCand = 0.0;
  nNumCandCH = 0;
  
  for ( std::vector<pat::PackedCandidate>::const_iterator cand = candidates-> begin();
    cand != candidates->end();
    cand ++ )
  {
    if ( !isCH( abs(cand->pdgId()) ) ) {
      continue;
    }
    
    if ( cand->puppiWeight() > 0 || cand->puppiWeightNoLep() > 0 ) {
      dSumPosZVtxCand += cand->vz();
      dSumSqrPosZVtxCand += cand->vz() * cand->vz();
      nNumCandCH++;
    }
  }
  
  // The vertex is found
  dPosZVtxCand = dSumPosZVtxCand / nNumCandCH;
  // Will it be needed?
  dVarPosZVtxCand = dSumSqrPosZVtxCand / nNumCandCH - dPosZVtxCand * dPosZVtxCand;
  
  for (TrackingParticleCollection::size_type i=0; i<simHandle->size(); i++) {
    TrackingParticleRef simRef(simHandle, i);
    const TrackingParticle* simTP = simRef.get();
    if ( ! tpSelector_(*simTP) ) continue;

    b_genMuon = TLorentzVector(simRef->momentum().x(), simRef->momentum().y(), simRef->momentum().z(), simRef->energy() );
    b_genMuon_isTightOptimized = false;
    b_genMuon_isTightCustom = false;
    b_genMuon_isTight = false;
    b_genMuon_isMedium = false;
    b_genMuon_isLoose = false;

    b_genMuon_pTresolution = -9;
    b_genMuon_noRecHitGEM = -1;
	b_genMuon_noRecHitME0 = -1;
    b_genMuon_isME0Muon = false;
    b_genMuon_isGEMMuon = false;
    b_genMuon_isMuon = false;
    b_genMuon_isGlobalMuon = false;
    b_genMuon_isStandAloneMuon = false;
    b_genMuon_signal = false;
    b_genMuon_isLooseMod = false;
    b_genMuon_isTightModNoIP = false;
    b_genMuon_isTightModIPxy = false;
    b_genMuon_isTightModIPz = false;
    b_genMuon_isTightModIPxyz = false;
    b_genMuon_numberOfValidMuonGEMHits = -1;
    b_genMuon_numberOfValidMuonME0Hits = -1;

    b_genMuon_deltaXME0 = 0;
    b_genMuon_deltaYME0 = 0;
    b_genMuon_deltaXErrME0 = 0;
    b_genMuon_deltaYErrME0 = 0;
    b_genMuon_deltaDXDZME0 = 0;
    b_genMuon_deltaDYDZME0 = 0;
    b_genMuon_deltaDXDZErrME0 = 0;
    b_genMuon_deltaDYDZErrME0 = 0;

    if ( simToMuonColl.find(simRef) != simToMuonColl.end() ) {
      vector<pair<RefToBase<Muon>, double> > MuRefV = simToMuonColl[simRef];      
      if ( !MuRefV.empty()) {
	const Muon* mu = MuRefV.begin()->first.get();
	signalMuons.push_back(mu);
	b_genMuon_signal = true;
	b_genMuon_pTresolution = (b_genMuon.Pt()-mu->pt())/b_genMuon.Pt();
    
    b_genMuon_poszPV0       = pv0.position().z();
    b_genMuon_poszSimPV     = simPVh.position().z();
    b_genMuon_poszLPVOfCand = dPosZVtxCand;
    b_genMuon_varposzCand   = dVarPosZVtxCand;
    b_genMuon_poszMuon      = mu->muonBestTrack()->vz();
    
    b_genMuon_poszPV01D   = vertices1D->front().position().z();
    b_genMuon_poszPV01DBS = vertices1DBS->front().position().z();
    b_genMuon_poszPV04D   = vertices4D->front().position().z();
    b_genMuon_poszPV04DBS = vertices4DBS->front().position().z();
    b_genMuon_poszPV0BS   = verticesBS->front().position().z();

    b_genMuon_TrkIso03 = mu->isolationR03().sumPt/mu->pt();
    b_genMuon_TrkIso05 = mu->isolationR05().sumPt/mu->pt();
    b_genMuon_pfIso03 = (mu->pfIsolationR03().sumChargedHadronPt + TMath::Max(0.,mu->pfIsolationR03().sumNeutralHadronEt + mu->pfIsolationR03().sumPhotonEt - 0.5*mu->pfIsolationR03().sumPUPt))/mu->pt();
    b_genMuon_pfIso04 = (mu->pfIsolationR04().sumChargedHadronPt + TMath::Max(0.,mu->pfIsolationR04().sumNeutralHadronEt + mu->pfIsolationR04().sumPhotonEt - 0.5*mu->pfIsolationR04().sumPUPt))/mu->pt();
    
    b_genMuon_pfIso03ChargedHadronPt = mu->pfIsolationR03().sumChargedHadronPt;
    b_genMuon_pfIso03NeutralHadronEt = mu->pfIsolationR03().sumNeutralHadronEt;
    b_genMuon_pfIso03PhotonEt        = mu->pfIsolationR03().sumPhotonEt;
    b_genMuon_pfIso03PUPt            = mu->pfIsolationR03().sumPUPt;
    
    b_genMuon_pfIso04ChargedHadronPt = mu->pfIsolationR04().sumChargedHadronPt;
    b_genMuon_pfIso04NeutralHadronEt = mu->pfIsolationR04().sumNeutralHadronEt;
    b_genMuon_pfIso04PhotonEt        = mu->pfIsolationR04().sumPhotonEt;
    b_genMuon_pfIso04PUPt            = mu->pfIsolationR04().sumPUPt;

	puppiIso puppiIsoGen = getPuppiIso( mu, candidates, 0);
	//cout << "pIso.combined "<< pIso.combined  <<endl;

	b_genMuon_puppiIsoWithLep    = puppiIsoGen.withLep;
	b_genMuon_puppiIsoWithoutLep = puppiIsoGen.withoutlep;
	b_genMuon_puppiIsoCombined   = puppiIsoGen.combined;

	b_genMuon_puppiIsoWithLep03    = puppiIsoGen.withLep03;
	b_genMuon_puppiIsoWithoutLep03 = puppiIsoGen.withoutlep03;
	b_genMuon_puppiIsoCombined03   = puppiIsoGen.combined03;

	b_genMuon_puppiIsoWithLep05    = puppiIsoGen.withLep05;
	b_genMuon_puppiIsoWithoutLep05 = puppiIsoGen.withoutlep05;
	b_genMuon_puppiIsoCombined05   = puppiIsoGen.combined05;
    
    b_genMuon_puppiIsoWithLep03ChargedHadron = puppiIsoGen.withLep03CH;
    b_genMuon_puppiIsoWithLep03NeutralHadron = puppiIsoGen.withLep03NH;
    b_genMuon_puppiIsoWithLep03Photon        = puppiIsoGen.withLep03PH;
    
    b_genMuon_puppiIsoWithLep04ChargedHadron = puppiIsoGen.withLep04CH;
    b_genMuon_puppiIsoWithLep04NeutralHadron = puppiIsoGen.withLep04NH;
    b_genMuon_puppiIsoWithLep04Photon        = puppiIsoGen.withLep04PH;
    
    b_genMuon_puppiIsoWithLep05ChargedHadron = puppiIsoGen.withLep05CH;
    b_genMuon_puppiIsoWithLep05NeutralHadron = puppiIsoGen.withLep05NH;
    b_genMuon_puppiIsoWithLep05Photon        = puppiIsoGen.withLep05PH;
    
    b_genMuon_puppiIsoWithoutLep03ChargedHadron = puppiIsoGen.withoutlep03CH;
    b_genMuon_puppiIsoWithoutLep03NeutralHadron = puppiIsoGen.withoutlep03NH;
    b_genMuon_puppiIsoWithoutLep03Photon        = puppiIsoGen.withoutlep03PH;
    
    b_genMuon_puppiIsoWithoutLep04ChargedHadron = puppiIsoGen.withoutlep04CH;
    b_genMuon_puppiIsoWithoutLep04NeutralHadron = puppiIsoGen.withoutlep04NH;
    b_genMuon_puppiIsoWithoutLep04Photon        = puppiIsoGen.withoutlep04PH;
    
    b_genMuon_puppiIsoWithoutLep05ChargedHadron = puppiIsoGen.withoutlep05CH;
    b_genMuon_puppiIsoWithoutLep05NeutralHadron = puppiIsoGen.withoutlep05NH;
    b_genMuon_puppiIsoWithoutLep05Photon        = puppiIsoGen.withoutlep05PH;
    
    b_genMuon_puppiIsoNumOfCands       = puppiIsoGen.nNumCands;
    b_genMuon_puppiIsoNumOfCandsInR05  = puppiIsoGen.nNumCandsInR05;
    b_genMuon_puppiIsoNumOfCandsOutR05 = puppiIsoGen.nNumCandsOutR05;
    
    b_genMuon_puppiIsoNumOfCandsInR05OT = puppiIsoGen.nNumCandsInR05OT;
    
    b_genMuon_puppiIsoNumOfCandsInR05CH = puppiIsoGen.nNumCandsInR05CH;
    b_genMuon_puppiIsoNumOfCandsInR05NH = puppiIsoGen.nNumCandsInR05NH;
    b_genMuon_puppiIsoNumOfCandsInR05PH = puppiIsoGen.nNumCandsInR05PH;
    
    b_genMuon_puppiIsoNumOfCandsInR05CHApp = puppiIsoGen.nNumCandsInR05CHApp;
    b_genMuon_puppiIsoNumOfCandsInR05NHApp = puppiIsoGen.nNumCandsInR05NHApp;
    b_genMuon_puppiIsoNumOfCandsInR05PHApp = puppiIsoGen.nNumCandsInR05PHApp;
    
	b_genMuon_isTightOptimized = isTightMuonCustomOptimized(*mu, pv0);
	b_genMuon_isTightCustom = isTightMuonCustom(*mu, pv0);
	b_genMuon_isTight = muon::isTightMuon(*mu, pv0);
	b_genMuon_isMedium = muon::isMediumMuon(*mu);
	b_genMuon_isLoose = muon::isLooseMuon(*mu);
	b_genMuon_isME0Muon = mu->isME0Muon();
	b_genMuon_isGEMMuon = mu->isGEMMuon();
	b_genMuon_isMuon = mu->isMuon();
    b_genMuon_isGlobalMuon = mu->isGlobalMuon();
    b_genMuon_isStandAloneMuon = mu->isStandAloneMuon();
    
    b_genMuon_isLooseMod = isLooseMod(mu);
    b_genMuon_isTightModNoIP  = isTightMod(vertices, simPVh, mu, false, false);
    b_genMuon_isTightModIPxy  = isTightMod(vertices, simPVh, mu, true,  false);
    b_genMuon_isTightModIPz   = isTightMod(vertices, simPVh, mu, false, true);
    b_genMuon_isTightModIPxyz = isTightMod(vertices, simPVh, mu, true,  true);
    
    const vector<MuonChamberMatch>& chambers = mu->matches();
    for( std::vector<MuonChamberMatch>::const_iterator chamber = chambers.begin(); chamber != chambers.end(); ++chamber ){
      for ( std::vector<reco::MuonSegmentMatch>::const_iterator segment = chamber->me0Matches.begin(); segment != chamber->me0Matches.end(); ++segment ){
        if (chamber->detector() == 5){
	      auto me0Segment = (*(*segment).me0SegmentRef);
	      b_genMuon_noRecHitME0 = me0Segment.nRecHits();

          b_genMuon_deltaXME0 = fabs( chamber->x - segment->x );
          b_genMuon_deltaYME0 = fabs( chamber->y - segment->y );
          b_genMuon_deltaDXDZME0 = fabs( chamber->dXdZ - segment->dXdZ );
          b_genMuon_deltaDYDZME0 = fabs( chamber->dYdZ - segment->dYdZ );
        }
      }
    }

	const reco::Track* muonTrack = 0;  
	if ( mu->globalTrack().isNonnull() ) muonTrack = mu->globalTrack().get();
	else if ( mu->outerTrack().isNonnull()  ) muonTrack = mu->outerTrack().get();
	if (muonTrack){
	  b_genMuon_noRecHitGEM = nGEMhit(mu);
	  b_genMuon_numberOfValidMuonGEMHits = muonTrack->hitPattern().numberOfValidMuonGEMHits();
	  b_genMuon_numberOfValidMuonME0Hits = muonTrack->hitPattern().numberOfValidMuonME0Hits();
	}

	tmvaValues.clear();
	tmvaValues = collectTMVAvalues(*mu, pv0);
	b_genMuon_tmva_bdt = bdt_->GetMvaValue(tmvaValues);
	b_genMuon_tmva_mlp = mlp_->GetMvaValue(tmvaValues);
      }
    }
    
    genttree_->Fill();
  }

  for (size_t i = 0; i < muonHandle->size(); ++i) {
    treereset();
    
    auto muRef = muonHandle->refAt(i);
    const Muon* mu = muRef.get();

    b_recoMuon = TLorentzVector(mu->momentum().x(), mu->momentum().y(), mu->momentum().z(), mu->energy() );

    b_recoMuon_signal = false;
    for (auto signal : signalMuons){
      if (mu == signal){
	b_recoMuon_signal = true;
	break;
      }
    }
    
    b_recoMuon_pdgId = 0;
    if ( muonToSimColl.find(muRef) != muonToSimColl.end() ) {
      auto trkRefV = muonToSimColl[muRef];
      if ( !trkRefV.empty()) {
	const TrackingParticle* trkParticle = trkRefV.begin()->first.get();
	b_recoMuon_pdgId = trkParticle->pdgId();
	// cout << "trkParticle " << trkParticle->pdgId()
	//      << " pt = " << trkParticle->pt()
	//      << " eta = " << trkParticle->eta()
	//      << endl;
      }
    }
    
    b_recoMuon_isLooseMod = false;
    b_recoMuon_isTightModNoIP  = false;
    b_recoMuon_isTightModIPxy  = false;
    b_recoMuon_isTightModIPz   = false;
    b_recoMuon_isTightModIPxyz = false;
    
    b_recoMuon_poszPV0       = pv0.position().z();
    b_recoMuon_poszSimPV     = simPVh.position().z();
    b_recoMuon_poszLPVOfCand = dPosZVtxCand;
    b_recoMuon_poszMuon      = mu->muonBestTrack()->vz();
    
    b_recoMuon_poszPV01D   = vertices1D->front().position().z();
    b_recoMuon_poszPV01DBS = vertices1DBS->front().position().z();
    b_recoMuon_poszPV04D   = vertices4D->front().position().z();
    b_recoMuon_poszPV04DBS = vertices4DBS->front().position().z();
    b_recoMuon_poszPV0BS   = verticesBS->front().position().z();

    b_recoMuon_TrkIso03 = mu->isolationR03().sumPt/mu->pt();
    b_recoMuon_TrkIso05 = mu->isolationR05().sumPt/mu->pt();
    b_recoMuon_PFIso04 = (mu->pfIsolationR04().sumChargedHadronPt + TMath::Max(0.,mu->pfIsolationR04().sumNeutralHadronEt + mu->pfIsolationR04().sumPhotonEt - 0.5*mu->pfIsolationR04().sumPUPt))/mu->pt();
    b_recoMuon_PFIso03 = (mu->pfIsolationR03().sumChargedHadronPt + TMath::Max(0.,mu->pfIsolationR03().sumNeutralHadronEt + mu->pfIsolationR03().sumPhotonEt - 0.5*mu->pfIsolationR03().sumPUPt))/mu->pt();
    
    b_recoMuon_PFIso03ChargedHadronPt = mu->pfIsolationR03().sumChargedHadronPt;
    b_recoMuon_PFIso03NeutralHadronEt = mu->pfIsolationR03().sumNeutralHadronEt;
    b_recoMuon_PFIso03PhotonEt        = mu->pfIsolationR03().sumPhotonEt;
    b_recoMuon_PFIso03PUPt            = mu->pfIsolationR03().sumPUPt;
    
    b_recoMuon_PFIso04ChargedHadronPt = mu->pfIsolationR04().sumChargedHadronPt;
    b_recoMuon_PFIso04NeutralHadronEt = mu->pfIsolationR04().sumNeutralHadronEt;
    b_recoMuon_PFIso04PhotonEt        = mu->pfIsolationR04().sumPhotonEt;
    b_recoMuon_PFIso04PUPt            = mu->pfIsolationR04().sumPUPt;
    
    puppiIso puppiIsoRec = getPuppiIso( mu, candidates, 1);
    
    b_recoMuon_puppiIsoWithLep    = puppiIsoRec.withLep;
    b_recoMuon_puppiIsoWithoutLep = puppiIsoRec.withoutlep;
    b_recoMuon_puppiIsoCombined   = puppiIsoRec.combined;
    
    b_recoMuon_puppiIsoWithLep03    = puppiIsoRec.withLep03;
    b_recoMuon_puppiIsoWithoutLep03 = puppiIsoRec.withoutlep03;
    b_recoMuon_puppiIsoCombined03   = puppiIsoRec.combined03;
    
    b_recoMuon_puppiIsoWithLep05    = puppiIsoRec.withLep05;
    b_recoMuon_puppiIsoWithoutLep05 = puppiIsoRec.withoutlep05;
    b_recoMuon_puppiIsoCombined05   = puppiIsoRec.combined05;
    
    b_recoMuon_puppiIsoWithLep03ChargedHadron = puppiIsoRec.withLep03CH;
    b_recoMuon_puppiIsoWithLep03NeutralHadron = puppiIsoRec.withLep03NH;
    b_recoMuon_puppiIsoWithLep03Photon        = puppiIsoRec.withLep03PH;
    
    b_recoMuon_puppiIsoWithLep04ChargedHadron = puppiIsoRec.withLep04CH;
    b_recoMuon_puppiIsoWithLep04NeutralHadron = puppiIsoRec.withLep04NH;
    b_recoMuon_puppiIsoWithLep04Photon        = puppiIsoRec.withLep04PH;
    
    b_recoMuon_puppiIsoWithLep05ChargedHadron = puppiIsoRec.withLep05CH;
    b_recoMuon_puppiIsoWithLep05NeutralHadron = puppiIsoRec.withLep05NH;
    b_recoMuon_puppiIsoWithLep05Photon        = puppiIsoRec.withLep05PH;
    
    b_recoMuon_puppiIsoWithoutLep03ChargedHadron = puppiIsoRec.withoutlep03CH;
    b_recoMuon_puppiIsoWithoutLep03NeutralHadron = puppiIsoRec.withoutlep03NH;
    b_recoMuon_puppiIsoWithoutLep03Photon        = puppiIsoRec.withoutlep03PH;
    
    b_recoMuon_puppiIsoWithoutLep04ChargedHadron = puppiIsoRec.withoutlep04CH;
    b_recoMuon_puppiIsoWithoutLep04NeutralHadron = puppiIsoRec.withoutlep04NH;
    b_recoMuon_puppiIsoWithoutLep04Photon        = puppiIsoRec.withoutlep04PH;
    
    b_recoMuon_puppiIsoWithoutLep05ChargedHadron = puppiIsoRec.withoutlep05CH;
    b_recoMuon_puppiIsoWithoutLep05NeutralHadron = puppiIsoRec.withoutlep05NH;
    b_recoMuon_puppiIsoWithoutLep05Photon        = puppiIsoRec.withoutlep05PH;
    
    b_recoMuon_puppiIsoNumOfCands       = puppiIsoRec.nNumCands;
    b_recoMuon_puppiIsoNumOfCandsInR05  = puppiIsoRec.nNumCandsInR05;
    b_recoMuon_puppiIsoNumOfCandsOutR05 = puppiIsoRec.nNumCandsOutR05;
    
    b_recoMuon_puppiIsoNumOfCandsInR05OT = puppiIsoRec.nNumCandsInR05OT;
    
    b_recoMuon_puppiIsoNumOfCandsInR05CH = puppiIsoRec.nNumCandsInR05CH;
    b_recoMuon_puppiIsoNumOfCandsInR05NH = puppiIsoRec.nNumCandsInR05NH;
    b_recoMuon_puppiIsoNumOfCandsInR05PH = puppiIsoRec.nNumCandsInR05PH;
    
    b_recoMuon_puppiIsoNumOfCandsInR05CHApp = puppiIsoRec.nNumCandsInR05CHApp;
    b_recoMuon_puppiIsoNumOfCandsInR05NHApp = puppiIsoRec.nNumCandsInR05NHApp;
    b_recoMuon_puppiIsoNumOfCandsInR05PHApp = puppiIsoRec.nNumCandsInR05PHApp;
    
    b_recoMuon_isTightOptimized = isTightMuonCustomOptimized(*mu, pv0);
    b_recoMuon_isTightCustom = isTightMuonCustom(*mu, pv0);
    b_recoMuon_isTight = muon::isTightMuon(*mu, pv0);
    b_recoMuon_isMedium = muon::isMediumMuon(*mu);
    b_recoMuon_isLoose = muon::isLooseMuon(*mu);
    b_recoMuon_isME0Muon = mu->isME0Muon();
    b_recoMuon_isGEMMuon = mu->isGEMMuon();
    b_recoMuon_isMuon = mu->isMuon();
    b_recoMuon_isGlobalMuon = mu->isGlobalMuon();
    b_recoMuon_isStandAloneMuon = mu->isStandAloneMuon();
    
    b_recoMuon_isLooseMod = isLooseMod(mu);
    b_recoMuon_isTightModNoIP  = isTightMod(vertices, simPVh, mu, false, false);
    b_recoMuon_isTightModIPxy  = isTightMod(vertices, simPVh, mu, true,  false);
    b_recoMuon_isTightModIPz   = isTightMod(vertices, simPVh, mu, false, true);
    b_recoMuon_isTightModIPxyz = isTightMod(vertices, simPVh, mu, true,  true);
    
    b_recoMuon_noRecHitGEM = -1;
    b_recoMuon_numberOfValidMuonGEMHits = -1;
    b_recoMuon_numberOfValidMuonME0Hits = -1;

    const reco::Track* muonTrack = 0;  
    if ( mu->globalTrack().isNonnull() ) muonTrack = mu->globalTrack().get();
    else if ( mu->outerTrack().isNonnull()  ) muonTrack = mu->outerTrack().get();
    if (muonTrack){
      b_recoMuon_noRecHitGEM = nGEMhit(mu);
      b_recoMuon_numberOfValidMuonGEMHits = muonTrack->hitPattern().numberOfValidMuonGEMHits();
      b_recoMuon_numberOfValidMuonME0Hits = muonTrack->hitPattern().numberOfValidMuonME0Hits();
    }

    const vector<MuonChamberMatch>& chambers = mu->matches();
    b_recoMuon_noChamberMatch = chambers.size();
    b_recoMuon_noSegment = 0;
    b_recoMuon_noSegmentDT = 0;
    b_recoMuon_noSegmentCSC = 0;
    b_recoMuon_noSegmentRPC = 0;
    b_recoMuon_noSegmentGEM = 0;
    b_recoMuon_noSegmentME0 = 0;

    b_recoMuon_noRecHitME0 = 0;

    b_recoMuon_edgeXME0 = 0;
    b_recoMuon_edgeYME0 = 0;
    b_recoMuon_chamberMatchXME0 = 0;
    b_recoMuon_chamberMatchYME0 = 0;
    b_recoMuon_chamberMatchXErrME0 = 0;
    b_recoMuon_chamberMatchYErrME0 = 0;
    b_recoMuon_chamberMatchDXDZME0 = 0;
    b_recoMuon_chamberMatchDYDZME0 = 0;    
    b_recoMuon_chamberMatchDXDZErrME0 = 0;    
    b_recoMuon_chamberMatchDYDZErrME0 = 0;    

    b_recoMuon_segmentMatchXME0 = 0;
    b_recoMuon_segmentMatchYME0 = 0;
    b_recoMuon_segmentMatchXErrME0 = 0;
    b_recoMuon_segmentMatchYErrME0 = 0;
    b_recoMuon_segmentMatchDXDZME0 = 0;
    b_recoMuon_segmentMatchDYDZME0 = 0;    
    b_recoMuon_segmentMatchDXDZErrME0 = 0;    
    b_recoMuon_segmentMatchDYDZErrME0 = 0;    

    b_recoMuon_distance = 0;
    b_recoMuon_distErr = 0;

    b_recoMuon_deltaXME0 = 0;
    b_recoMuon_deltaYME0 = 0;
    b_recoMuon_deltaXDivBySegErrME0 = 0;
    b_recoMuon_deltaYDivBySegErrME0 = 0;
    b_recoMuon_deltaXDivByChamErrME0 = 0;
    b_recoMuon_deltaYDivByChamErrME0 = 0;
    b_recoMuon_deltaDXDZME0 = 0;
    b_recoMuon_deltaDYDZME0 = 0;
    b_recoMuon_deltaDXDZDivBySegErrME0 = 0;
    b_recoMuon_deltaDYDZDivBySegErrME0 = 0;
    b_recoMuon_deltaDXDZDivByChamErrME0 = 0;
    b_recoMuon_deltaDYDZDivByChamErrME0 = 0;

    for( std::vector<MuonChamberMatch>::const_iterator chamber = chambers.begin(); chamber != chambers.end(); ++chamber ){

      //cout<< "chamber->detector() "<< chamber->detector()<<endl;
      for ( std::vector<reco::MuonSegmentMatch>::const_iterator segment = chamber->segmentMatches.begin(); segment != chamber->segmentMatches.end(); ++segment ){
	++b_recoMuon_noSegment;
	if (chamber->detector() == 1){
	  ++b_recoMuon_noSegmentDT;
	}
	if (chamber->detector() == 2){
	  ++b_recoMuon_noSegmentCSC;
	}
      }
      
      for ( std::vector<reco::MuonRPCHitMatch>::const_iterator segment = chamber->rpcMatches.begin(); segment != chamber->rpcMatches.end(); ++segment ){
	++b_recoMuon_noSegment;      
	if (chamber->detector() == 3){
	  ++b_recoMuon_noSegmentRPC;
	}
      }
      
      for ( std::vector<reco::MuonSegmentMatch>::const_iterator segment = chamber->gemMatches.begin(); segment != chamber->gemMatches.end(); ++segment ){
	++b_recoMuon_noSegment;      
	if (chamber->detector() == 4){ //gem
	  // auto gemSegment = (*(*segment).gemSegmentRef);
	  // b_recoMuon_noRecHitGEM += gemSegment.nRecHits();
	  ++b_recoMuon_noSegmentGEM;
	}
      }

      for ( std::vector<reco::MuonSegmentMatch>::const_iterator segment = chamber->me0Matches.begin(); segment != chamber->me0Matches.end(); ++segment ){
	++b_recoMuon_noSegment;      
	if (chamber->detector() == 5){ //me0, chamber->detector() = muonSubdetId
	  ++b_recoMuon_noSegmentME0;
	  auto me0Segment = (*(*segment).me0SegmentRef);
	  b_recoMuon_noRecHitME0 = me0Segment.nRecHits();

	  b_recoMuon_edgeXME0 = chamber->edgeX;
	  b_recoMuon_edgeYME0 = chamber->edgeY;
	  b_recoMuon_chamberMatchXME0 = chamber->x; // x position of the track
	  b_recoMuon_chamberMatchYME0 = chamber->y;
	  b_recoMuon_chamberMatchXErrME0 = chamber->xErr;
	  b_recoMuon_chamberMatchYErrME0 = chamber->yErr;
	  b_recoMuon_chamberMatchDXDZME0 = chamber->dXdZ;
	  b_recoMuon_chamberMatchDYDZME0 = chamber->dYdZ;
	  b_recoMuon_chamberMatchDXDZErrME0 = chamber->dXdZErr;
	  b_recoMuon_chamberMatchDYDZErrME0 = chamber->dYdZErr;

	  b_recoMuon_segmentMatchXME0 = segment->x; // x position of the matched segment
	  b_recoMuon_segmentMatchYME0 = segment->y;
	  b_recoMuon_segmentMatchXErrME0 = segment->xErr;
	  b_recoMuon_segmentMatchYErrME0 = segment->yErr;
	  b_recoMuon_segmentMatchDXDZME0 = segment->dXdZ;
	  b_recoMuon_segmentMatchDYDZME0 = segment->dYdZ;
	  b_recoMuon_segmentMatchDXDZErrME0 = segment->dXdZErr;
	  b_recoMuon_segmentMatchDYDZErrME0 = segment->dYdZErr;

	  b_recoMuon_distance = chamber->dist();
	  b_recoMuon_distErr = chamber->distErr();    

	  //b_recoMuon_deltaXME0 = fabs( (b_recoMuon_chamberMatchXME0-b_recoMuon_segmentMatchXME0) / sqrt(pow(b_recoMuon_chamberMatchXErrME0, 2)+pow(b_recoMuon_segmentMatchXErrME0,2)) );
	  b_recoMuon_deltaXME0 = fabs( b_recoMuon_chamberMatchXME0 - b_recoMuon_segmentMatchXME0 );
	  b_recoMuon_deltaYME0 = fabs( b_recoMuon_chamberMatchYME0 - b_recoMuon_segmentMatchYME0 );
	  b_recoMuon_deltaXDivBySegErrME0 = fabs( b_recoMuon_deltaXME0 / b_recoMuon_segmentMatchXErrME0 );
	  b_recoMuon_deltaYDivBySegErrME0 = fabs( b_recoMuon_deltaYME0 / b_recoMuon_segmentMatchYErrME0 );
	  b_recoMuon_deltaXDivByChamErrME0 = fabs( b_recoMuon_deltaXME0 / sqrt(pow(b_recoMuon_chamberMatchXErrME0,2)+pow(b_recoMuon_segmentMatchXErrME0,2)) );
	  b_recoMuon_deltaYDivByChamErrME0 = fabs( b_recoMuon_deltaYME0 / sqrt(pow(b_recoMuon_chamberMatchYErrME0,2)+pow(b_recoMuon_segmentMatchYErrME0,2)) );

	  b_recoMuon_deltaDXDZME0 = fabs( b_recoMuon_chamberMatchDXDZME0-b_recoMuon_segmentMatchDXDZME0 );
	  b_recoMuon_deltaDYDZME0 = fabs( b_recoMuon_chamberMatchDYDZME0-b_recoMuon_segmentMatchDYDZME0 );
	  b_recoMuon_deltaDXDZDivBySegErrME0 = fabs( b_recoMuon_deltaDXDZME0 / b_recoMuon_segmentMatchDXDZErrME0 );
	  b_recoMuon_deltaDYDZDivBySegErrME0 = fabs( b_recoMuon_deltaDYDZME0 / b_recoMuon_segmentMatchDYDZErrME0 );
	  b_recoMuon_deltaDXDZDivByChamErrME0 = fabs( b_recoMuon_deltaDXDZME0 / sqrt(pow(b_recoMuon_chamberMatchDXDZErrME0,2)+pow(b_recoMuon_segmentMatchDXDZErrME0,2)) );
	  b_recoMuon_deltaDYDZDivByChamErrME0 = fabs( b_recoMuon_deltaDYDZME0 / sqrt(pow(b_recoMuon_chamberMatchDYDZErrME0,2)+pow(b_recoMuon_segmentMatchDYDZErrME0,2)) );

	  //cout << "station = " << chamber->station() << endl;
	  //cout << "detector = " << chamber->detector() << endl;

	  //b_recoMuon_isMaskBestInME0ByDX = segment->isMask((*segment).BestInChamberByDX);
	  //Mask += segment->mask;
	  //cout << "flag:   " << (*segment).BestInChamberByDX << endl;
	  //cout << "mask:   " << (*segment).mask << endl;
	  //cout << "isMask: " << segment->isMask((*segment).BestInChamberByDX) << endl;

	}
      }
      
    }
    
    tmvaValues.clear();
    tmvaValues = collectTMVAvalues(*mu, pv0);
    b_recoMuon_tmva_bdt = bdt_->GetMvaValue(tmvaValues);
    b_recoMuon_tmva_mlp = mlp_->GetMvaValue(tmvaValues);

    b_recoMuon_global = tmvaValues[0];
    b_recoMuon_pf = tmvaValues[1];
    b_recoMuon_chi2 = tmvaValues[2];
    b_recoMuon_chi2pos = tmvaValues[3];
    b_recoMuon_trkKink = tmvaValues[4];
    b_recoMuon_segcompati = tmvaValues[5];
    b_recoMuon_nstations = tmvaValues[6];
    b_recoMuon_nglobalhits = tmvaValues[7];
    b_recoMuon_trackdxy = tmvaValues[8];
    b_recoMuon_trackdz = tmvaValues[9];
    b_recoMuon_ninnerhits = tmvaValues[10];
    b_recoMuon_trackerlayers =tmvaValues[11];

    recottree_->Fill();
  }

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
    values.push_back(fabs(mu.muonBestTrack()->dxy(pv0.position())));
    values.push_back(fabs(mu.muonBestTrack()->dz(pv0.position())));
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
  if ( !(fabs(mu.muonBestTrack()->dxy(pv0.position())) < 0.2) ) return false;
  //if ( !(fabs(mu.muonBestTrack()->dz(pv0.position())) < 0.5) ) return false;
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
  if ( !(fabs(mu.muonBestTrack()->dxy(pv0.position())) < 0.02) ) return false; // < 0.2
  //if ( !(fabs(mu.muonBestTrack()->dz(pv0.position())) < 0.5) ) return false;
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

void MuonAnalyser::treereset()
{
  b_recoMuon_global = -9; b_recoMuon_pf = -9;
  b_recoMuon_chi2 = -9; b_recoMuon_chi2pos = -9; b_recoMuon_trkKink = -9; b_recoMuon_segcompati = -9;
  b_recoMuon_nglobalhits = -9; b_recoMuon_nstations = -9;
  b_recoMuon_trackdxy = -9; b_recoMuon_trackdz = -9;
  b_recoMuon_ninnerhits = -9; b_recoMuon_trackerlayers = -9;
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

    double d_eta = fabs( cand->eta() - mu->eta() ) ;
    double d_phi = fabs( cand->phi() - mu->phi() ) ; 
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

bool MuonAnalyser::isTightMod(const Handle<VertexCollection> &vertexHandle, const SimVertex &simPVh, const reco::Muon *muon, bool useIPxy, bool useIPz)
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
        const reco::VertexCollection* vertices = vertexHandle.product();
        
        double distInit = 24;
        int indexFinal = 0;
        for(int i = 0; i < (int)vertices->size(); i++){
            
            //double vtxX = (*vertices)[i].x();
            //double vtxY = (*vertices)[i].y();
            double vtxZ = (*vertices)[i].z();
            
            double dist = fabs(muonZ - vtxZ);
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
            
            ipxySim = fabs(muon->muonBestTrack()->dxy(math::XYZPoint(point.x(),point.y(),point.z())));
            ipzSim = fabs(muon->muonBestTrack()->dz(math::XYZPoint(point.x(),point.y(),point.z())));
            
        }
        else if(vtxCoord[0] > 0.5 && vtxCoord[0] < 1.5){//DY samples
            
            ipxySim = fabs(muon->muonBestTrack()->dxy(math::XYZPoint(pointDY.x(),pointDY.y(),pointDY.z())));
            ipzSim = fabs(muon->muonBestTrack()->dz(math::XYZPoint(pointDY.x(),pointDY.y(),pointDY.z())));
            
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
                
                ipxy = fabs(muon->muonBestTrack()->dxy((*vertices)[indexFinal].position())) < 0.2;
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
                
                ipz = fabs(muon->muonBestTrack()->dz((*vertices)[indexFinal].position())) < 0.5;
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

DEFINE_FWK_MODULE(MuonAnalyser);
