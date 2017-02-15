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
    double combined;
    double withLep;
    double withoutlep;
  };

public:
  explicit MuonAnalyser(const edm::ParameterSet&);
  ~MuonAnalyser();

  bool isLooseMuonCustom(const reco::Muon& mu) const;
  bool isTightMuonCustom(const reco::Muon& mu, reco::Vertex pv0) const;
  bool isTightMuonCustomOptimized(const reco::Muon& mu, reco::Vertex pv0) const;
  std::vector<double> collectTMVAvalues(const reco::Muon& mu, reco::Vertex pv0) const;
  int nGEMhit(const reco::Muon * mu) const;
  void treereset();
  
  puppiIso getPuppiIso(const reco::Muon *mu, const vector< pat::PackedCandidate> *pcs) const;
  bool isNH( long pdgid ) const;
  bool isCH( long pdgid ) const;
  bool isPH( long pdgid ) const;

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

  TTree* genttree_;
  TTree* recottree_;
  TLorentzVector b_genMuon;
  bool b_genMuon_isTightOptimized, b_genMuon_isTight, b_genMuon_isMedium, b_genMuon_isLoose;
  bool b_genMuon_isME0Muon, b_genMuon_isGEMMuon, b_genMuon_isMuon, b_genMuon_signal;
  int b_genMuon_noRecHitGEM;
  float b_genMuon_pfIso03; float b_genMuon_pfIso04;
  float b_genMuon_TrkIso05; float b_genMuon_TrkIso03;
  float b_genMuon_puppiIsoWithLep, b_genMuon_puppiIsoWithoutLep, b_genMuon_puppiIsoCombined;
  int b_genMuon_numberOfValidMuonGEMHits, b_genMuon_numberOfValidMuonME0Hits;
  float b_genMuon_tmva_bdt, b_genMuon_tmva_mlp;

  float b_genMuon_deltaXME0, b_genMuon_deltaYME0;
  float b_genMuon_deltaXErrME0, b_genMuon_deltaYErrME0;
  float b_genMuon_deltaDXDZME0, b_genMuon_deltaDYDZME0;
  float b_genMuon_deltaDXDZErrME0, b_genMuon_deltaDYDZErrME0;

  TLorentzVector b_recoMuon;
  bool b_recoMuon_signal;
  bool b_recoMuon_isTightOptimized, b_recoMuon_isTight, b_recoMuon_isMedium, b_recoMuon_isLoose;
  bool b_recoMuon_isME0Muon, b_recoMuon_isGEMMuon;
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
  float b_recoMuon_PFIso04; float b_recoMuon_PFIso03;
  float b_recoMuon_TrkIso05; float b_recoMuon_TrkIso03;
  float b_recoMuon_puppiIsoWithLep, b_recoMuon_puppiIsoWithoutLep, b_recoMuon_puppiIsoCombined;
  bool b_recoMuon_isMuon;
  int b_recoMuon_numberOfValidMuonGEMHits, b_recoMuon_numberOfValidMuonME0Hits;

  float b_recoMuon_tmva_bdt, b_recoMuon_tmva_mlp;
   
  std::vector<double> tmvaValues;
  float b_tmva_bdt; float b_tmva_mlp;
  ReadBDT* bdt_;
  ReadMLP* mlp_;
  
  edm::EDGetTokenT<std::vector<reco::Vertex> > vtxToken_;
  edm::EDGetTokenT<TrackingParticleCollection> simToken_;
  edm::EDGetTokenT<edm::View<reco::Muon> > muonToken_;
  edm::EDGetTokenT<reco::MuonToTrackingParticleAssociator> muAssocToken_;
  edm::EDGetTokenT <std::vector< pat::PackedCandidate> > tokenPackedCandidate ;

  TrackingParticleSelector tpSelector_;
};

MuonAnalyser::MuonAnalyser(const edm::ParameterSet& pset)
{
  vtxToken_ = consumes<vector<Vertex> >(pset.getParameter<edm::InputTag>("primaryVertex"));
  simToken_ = consumes<TrackingParticleCollection>(pset.getParameter<InputTag>("simLabel"));
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

  genttree_ = fs->make<TTree>("gen", "gen");
  genttree_->Branch("genMuon", "TLorentzVector", &b_genMuon);
  genttree_->Branch("genMuon_isTightOptimized", &b_genMuon_isTightOptimized, "genMuon_isTightOptimized/O");
  genttree_->Branch("genMuon_isTight", &b_genMuon_isTight, "genMuon_isTight/O");
  genttree_->Branch("genMuon_isMedium", &b_genMuon_isMedium, "genMuon_isMedium/O");
  genttree_->Branch("genMuon_isLoose", &b_genMuon_isLoose, "genMuon_isLoose/O");
  genttree_->Branch("genMuon_noRecHitGEM", &b_genMuon_noRecHitGEM, "genMuon_noRecHitGEM/I");
  genttree_->Branch("genMuon_isME0Muon", &b_genMuon_isME0Muon, "genMuon_isME0Muon/O");
  genttree_->Branch("genMuon_isGEMMuon", &b_genMuon_isGEMMuon, "genMuon_isGEMMuon/O");
  genttree_->Branch("genMuon_isMuon", &b_genMuon_isMuon, "genMuon_isMuon/O");
  genttree_->Branch("genMuon_TrkIsolation03",&b_genMuon_TrkIso03,"genMuon_TrkIsolation03/F");
  genttree_->Branch("genMuon_TrkIsolation05",&b_genMuon_TrkIso05,"genMuon_TrkIsolation05/F");
  genttree_->Branch("genMuon_PFIsolation03",&b_genMuon_pfIso03,"genMuon_PFIsolation03/F");
  genttree_->Branch("genMuon_PFIsolation04",&b_genMuon_pfIso04,"genMuon_PFIsolation04/F");
  genttree_->Branch("genMuon_puppiIsoWithLep",&b_genMuon_puppiIsoWithLep,"genMuon_puppiIsoWithLep/F");
  genttree_->Branch("genMuon_puppiIsoWithoutLep",&b_genMuon_puppiIsoWithoutLep,"genMuon_puppiIsoWithoutLep/F");
  genttree_->Branch("genMuon_puppiIsoCombined",&b_genMuon_puppiIsoCombined,"genMuon_puppiIsoCombined/F");
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
  recottree_->Branch("recoMuon", "TLorentzVector", &b_recoMuon);  
  recottree_->Branch("recoMuon_pdgId", &b_recoMuon_pdgId, "recoMuon_pdgId/I");
  recottree_->Branch("recoMuon_signal", &b_recoMuon_signal, "recoMuon_signal/O");
  recottree_->Branch("recoMuon_isTightOptimized", &b_recoMuon_isTightOptimized, "recoMuon_isTightOptimized/O");
  recottree_->Branch("recoMuon_isTight", &b_recoMuon_isTight, "recoMuon_isTight/O");
  recottree_->Branch("recoMuon_isMedium", &b_recoMuon_isMedium, "recoMuon_isMedium/O");
  recottree_->Branch("recoMuon_isLoose", &b_recoMuon_isLoose, "recoMuon_isLoose/O");
  recottree_->Branch("recoMuon_isME0Muon", &b_recoMuon_isME0Muon, "recoMuon_isME0Muon/O");
  recottree_->Branch("recoMuon_isGEMMuon", &b_recoMuon_isGEMMuon, "recoMuon_isGEMMuon/O");
  recottree_->Branch("recoMuon_isMuon", &b_recoMuon_isMuon, "recoMuon_isMuon/O");
  recottree_->Branch("recoMuon_numberOfValidMuonGEMHits",&b_recoMuon_numberOfValidMuonGEMHits,"recoMuon_numberOfValidMuonGEMHits/I");
  recottree_->Branch("recoMuon_numberOfValidMuonME0Hits",&b_recoMuon_numberOfValidMuonME0Hits,"recoMuon_numberOfValidMuonME0Hits/I");

  recottree_->Branch("recoMuon_TrkIsolation03",&b_recoMuon_TrkIso03,"recoMuon_TrkIsolation03/F");
  recottree_->Branch("recoMuon_TrkIsolation05",&b_recoMuon_TrkIso05,"recoMuon_TrkIsolation05/F");
  recottree_->Branch("recoMuon_PFIsolation04",&b_recoMuon_PFIso04,"recoMuon_PFIsolation04/F");
  recottree_->Branch("recoMuon_PFIsolation03",&b_recoMuon_PFIso03,"recoMuon_PFIsolation03/F");
  recottree_->Branch("recoMuon_puppiIsoWithLep",&b_recoMuon_puppiIsoWithLep,"recoMuon_puppiIsoWithLep/F");
  recottree_->Branch("recoMuon_puppiIsoWithoutLep",&b_recoMuon_puppiIsoWithoutLep,"recoMuon_puppiIsoWithoutLep/F");
  recottree_->Branch("recoMuon_puppiIsoCombined",&b_recoMuon_puppiIsoCombined,"recoMuon_puppiIsoCombined/F");
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
  Handle<VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices); 
  if (vertices->empty()) { cout << "noPV" << endl; return; }
  auto pv0 = vertices->front();

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
  assoByHits->associateMuons(muonToSimColl, simToMuonColl, muonHandle, reco::InnerTk, simHandle);
  
  vector<const Muon*> signalMuons; signalMuons.clear();
  
  for (TrackingParticleCollection::size_type i=0; i<simHandle->size(); i++) {
    TrackingParticleRef simRef(simHandle, i);
    const TrackingParticle* simTP = simRef.get();
    if ( ! tpSelector_(*simTP) ) continue;

    b_genMuon = TLorentzVector(simRef->momentum().x(), simRef->momentum().y(), simRef->momentum().z(), simRef->energy() );
    b_genMuon_isTightOptimized = false;
    b_genMuon_isTight = false;
    b_genMuon_isMedium = false;
    b_genMuon_isLoose = false;
    b_genMuon_noRecHitGEM = -1;
    b_genMuon_isME0Muon = false;
    b_genMuon_isGEMMuon = false;
    b_genMuon_isMuon = false;
    b_genMuon_signal = false;
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

	b_genMuon_TrkIso03 = mu->isolationR03().sumPt/mu->pt();
	b_genMuon_TrkIso05 = mu->isolationR05().sumPt/mu->pt();
	b_genMuon_pfIso03 = (mu->pfIsolationR03().sumChargedHadronPt + TMath::Max(0.,mu->pfIsolationR03().sumNeutralHadronEt + mu->pfIsolationR03().sumPhotonEt - 0.5*mu->pfIsolationR03().sumPUPt))/mu->pt();
	b_genMuon_pfIso04 = (mu->pfIsolationR04().sumChargedHadronPt + TMath::Max(0.,mu->pfIsolationR04().sumNeutralHadronEt + mu->pfIsolationR04().sumPhotonEt - 0.5*mu->pfIsolationR04().sumPUPt))/mu->pt();

	puppiIso puppiIsoGen = getPuppiIso( mu, candidates);
	//cout << "pIso.combined "<< pIso.combined  <<endl;
  
    b_genMuon_puppiIsoWithLep    = puppiIsoGen.withLep;
    b_genMuon_puppiIsoWithoutLep = puppiIsoGen.withoutlep;
    b_genMuon_puppiIsoCombined   = puppiIsoGen.combined;
	
	b_genMuon_isTightOptimized = isTightMuonCustomOptimized(*mu, pv0);
	b_genMuon_isTight = isTightMuonCustom(*mu, pv0);
	b_genMuon_isMedium = muon::isMediumMuon(*mu);
	b_genMuon_isLoose = muon::isLooseMuon(*mu);
	b_genMuon_isME0Muon = mu->isME0Muon();
	b_genMuon_isGEMMuon = mu->isGEMMuon();
	b_genMuon_isMuon = mu->isMuon();

    b_genMuon_deltaXME0 = mu->dX(0,5, mu->ME0SegmentAndTrackArbitration); // position difference of track and segement
    b_genMuon_deltaYME0 = mu->dY(0,5, mu->ME0SegmentAndTrackArbitration);
    b_genMuon_deltaXErrME0 = mu->pullX(0,5, mu->ME0SegmentAndTrackArbitration); // delta X divided by segment error and propagation error
    b_genMuon_deltaYErrME0 = mu->pullY(0,5, mu->ME0SegmentAndTrackArbitration);
    b_genMuon_deltaDXDZME0 = mu->dDxDz(0,5, mu->ME0SegmentAndTrackArbitration); // slope difference of track and segment
    b_genMuon_deltaDYDZME0 = mu->dDyDz(0,5, mu->ME0SegmentAndTrackArbitration);
    b_genMuon_deltaDXDZErrME0 = mu->pullDxDz(0,5, mu->ME0SegmentAndTrackArbitration); //delta dXdZ divided by segment error and propagation error
    b_genMuon_deltaDYDZErrME0 = mu->pullDyDz(0,5, mu->ME0SegmentAndTrackArbitration);
	
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

    b_recoMuon_TrkIso03 = mu->isolationR03().sumPt/mu->pt();
    b_recoMuon_TrkIso05 = mu->isolationR05().sumPt/mu->pt();
    b_recoMuon_PFIso04 = (mu->pfIsolationR04().sumChargedHadronPt + TMath::Max(0.,mu->pfIsolationR04().sumNeutralHadronEt + mu->pfIsolationR04().sumPhotonEt - 0.5*mu->pfIsolationR04().sumPUPt))/mu->pt();
    b_recoMuon_PFIso03 = (mu->pfIsolationR03().sumChargedHadronPt + TMath::Max(0.,mu->pfIsolationR03().sumNeutralHadronEt + mu->pfIsolationR03().sumPhotonEt - 0.5*mu->pfIsolationR03().sumPUPt))/mu->pt();
    
	puppiIso puppiIsoRec = getPuppiIso( mu, candidates);
    
    b_recoMuon_puppiIsoWithLep    = puppiIsoRec.withLep;
    b_recoMuon_puppiIsoWithoutLep = puppiIsoRec.withoutlep;
    b_recoMuon_puppiIsoCombined   = puppiIsoRec.combined;

    b_recoMuon_isTightOptimized = isTightMuonCustomOptimized(*mu, pv0);
    b_recoMuon_isTight = isTightMuonCustom(*mu, pv0);
    b_recoMuon_isMedium = muon::isMediumMuon(*mu);
    b_recoMuon_isLoose = muon::isLooseMuon(*mu);
    b_recoMuon_isME0Muon = mu->isME0Muon();
    b_recoMuon_isGEMMuon = mu->isGEMMuon();
    b_recoMuon_isMuon = mu->isMuon();
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
  if ( !(mu.isGlobalMuon()) ) return false;
  if ( !(mu.isPFMuon() || mu.isTrackerMuon()) ) return false;
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
  if ( !(mu.globalTrack()->normalizedChi2() < 4.) ) return false; // < 10.
  if ( !(mu.globalTrack()->hitPattern().numberOfValidMuonHits() > 7) ) return false; // > 0
  if ( !(mu.numberOfMatchedStations() > 1) ) return false;
  if ( !(fabs(mu.muonBestTrack()->dxy(pv0.position())) < 0.05) ) return false; // < 0.2
  //if ( !(fabs(mu.muonBestTrack()->dz(pv0.position())) < 0.5) ) return false;
  if ( !(mu.innerTrack()->hitPattern().numberOfValidPixelHits() > 0) ) return false;
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

MuonAnalyser::puppiIso MuonAnalyser::getPuppiIso(const reco::Muon *mu, const vector< pat::PackedCandidate> *pcs) const
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
  const double dR_threshold = 0.4;
  
  double dR2_threshold = dR_threshold * dR_threshold;

  double val_PuppiWithLep    [3]= {0,0,0} ;
  double val_PuppiWithoutLep [3]= {0,0,0} ;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  // loop ever all the candidates, and accumulate PT deposit around the lepton.
  for( std::vector<pat::PackedCandidate>::const_iterator cand = pcs -> begin();
    cand != pcs->end();
    cand ++ )
  {
    // calc DR

    double d_eta = fabs( cand->eta() - mu->eta() ) ;
    double d_phi = fabs( cand->phi() - mu->phi() ) ; 
    d_phi = ( d_phi < M_PI ) ? d_phi : 2 * M_PI - d_phi ; 
    double dR2 = d_eta * d_eta  + d_phi * d_phi ;
    
    if( dR2 > dR2_threshold ) continue ;
    
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
      
      continue ;
    }
    
    // Check particleType dependent DR cut (remove overlapped candiadte)
    // The threshold values were taken from 'MuonPFIsolationSequence_cff.py'.
    if( pType == CH && dR2 < 0.0001*0.0001 ) continue ;
    if( pType == NH && dR2 < 0.01  *0.01   ) continue ;
    if( pType == PH && dR2 < 0.01  *0.01   ) continue ;

    // The candidate passed all the selection.
    // Now, add its PT to the variable with weight.

    val_PuppiWithLep   [ pType ] += cand -> pt() * cand -> puppiWeight() ;
    val_PuppiWithoutLep[ pType ] += cand -> pt() * cand -> puppiWeightNoLep();

   
  }// end of candidate LOOP.

  const double reliso_Puppi_withLep    = ( val_PuppiWithLep   [CH] + val_PuppiWithLep   [NH] + val_PuppiWithLep   [PH] ) / mu->pt() ;
  const double reliso_Puppi_withoutlep = ( val_PuppiWithoutLep[CH] + val_PuppiWithoutLep[NH] + val_PuppiWithoutLep[PH] ) / mu->pt() ;

  puppivalues.withLep    = reliso_Puppi_withLep;
  puppivalues.withoutlep = reliso_Puppi_withoutlep;
  puppivalues.combined   = _mix_fraction_ * reliso_Puppi_withLep + ( 1.0 - _mix_fraction_) * reliso_Puppi_withoutlep;
  
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

DEFINE_FWK_MODULE(MuonAnalyser);
