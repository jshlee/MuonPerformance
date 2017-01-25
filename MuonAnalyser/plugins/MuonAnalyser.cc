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

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/Associations/interface/MuonToTrackingParticleAssociator.h"
#include "SimMuon/MCTruth/interface/MuonToSimAssociatorByHits.h"
#include "SimTracker/Common/interface/TrackingParticleSelector.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"

#include <TTree.h>
#include <TFile.h>
#include <TLorentzVector.h>

class MuonAnalyser : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit MuonAnalyser(const edm::ParameterSet&);
  ~MuonAnalyser();

  bool isLooseMuonCustom(const reco::Muon& mu) const;
  bool isTightMuonCustom(const reco::Muon& mu, reco::Vertex pv0) const;
  int nGEMhit(const reco::Muon * mu) const;
  void treereset();

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

  TTree* genttree_;
  TTree* recottree_;
  TLorentzVector b_genMuon;
  bool b_genMuon_isTight, b_genMuon_isMedium, b_genMuon_isLoose, b_genMuon_isME0Muon, b_genMuon_isGEMMuon, b_genMuon_isMuon;
  int b_genMuon_noRecHitGEM;
  float b_genMuon_pfIso03; float b_genMuon_pfIso04;
  float b_genMuon_TrkIso05; float b_genMuon_TrkIso03;
  int b_genMuon_numberOfValidMuonGEMHits, b_genMuon_numberOfValidMuonME0Hits;
  
  TLorentzVector b_recoMuon;
  bool b_recoMuon_signal, b_recoMuon_isTight, b_recoMuon_isMedium, b_recoMuon_isLoose, b_recoMuon_isME0Muon, b_recoMuon_isGEMMuon;
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
  bool b_recoMuon_isMuon;
  int b_recoMuon_numberOfValidMuonGEMHits, b_recoMuon_numberOfValidMuonME0Hits;
  
  edm::EDGetTokenT<std::vector<reco::Vertex> > vtxToken_;
  edm::EDGetTokenT<TrackingParticleCollection> simToken_;
  edm::EDGetTokenT<edm::View<reco::Muon> > muonToken_;
  edm::EDGetTokenT<reco::MuonToTrackingParticleAssociator> muAssocToken_;

  TrackingParticleSelector tpSelector_;
};
using namespace std;
using namespace reco;
using namespace edm;
MuonAnalyser::MuonAnalyser(const edm::ParameterSet& pset)
{
  vtxToken_ = consumes<vector<Vertex> >(pset.getParameter<edm::InputTag>("primaryVertex"));
  simToken_ = consumes<TrackingParticleCollection>(pset.getParameter<InputTag>("simLabel"));
  muonToken_ = consumes<View<Muon> >(pset.getParameter<InputTag>("muonLabel"));
  muAssocToken_ = consumes<reco::MuonToTrackingParticleAssociator>(pset.getParameter<InputTag>("muAssocLabel"));

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
  genttree_->Branch("genMuon_numberOfValidMuonGEMHits",&b_genMuon_numberOfValidMuonGEMHits,"genMuon_numberOfValidMuonGEMHits/I");
  genttree_->Branch("genMuon_numberOfValidMuonME0Hits",&b_genMuon_numberOfValidMuonME0Hits,"genMuon_numberOfValidMuonME0Hits/I");

  recottree_ = fs->make<TTree>("reco", "reco");
  recottree_->Branch("recoMuon", "TLorentzVector", &b_recoMuon);  
  recottree_->Branch("recoMuon_pdgId", &b_recoMuon_pdgId, "recoMuon_pdgId/I");
  recottree_->Branch("recoMuon_signal", &b_recoMuon_signal, "recoMuon_signal/O");
  recottree_->Branch("recoMuon_isTight", &b_recoMuon_isTight, "recoMuon_isTight/O");
  recottree_->Branch("recoMuon_isMedium", &b_recoMuon_isMedium, "recoMuon_isMedium/O");
  recottree_->Branch("recoMuon_isLoose", &b_recoMuon_isLoose, "recoMuon_isLoose/O");
  recottree_->Branch("recoMuon_isME0Muon", &b_recoMuon_isME0Muon, "recoMuon_isME0Muon/O");
  recottree_->Branch("recoMuon_isGEMMuon", &b_recoMuon_isGEMMuon, "recoMuon_isGEMMuon/O");
  recottree_->Branch("recoMuon_noChamberMatch", &b_recoMuon_noChamberMatch, "recoMuon_noChamberMatch/I");
  recottree_->Branch("recoMuon_isMuon", &b_recoMuon_isMuon, "recoMuon_isMuon/O");
  recottree_->Branch("recoMuon_numberOfValidMuonGEMHits",&b_recoMuon_numberOfValidMuonGEMHits,"recoMuon_numberOfValidMuonGEMHits/I");
  recottree_->Branch("recoMuon_numberOfValidMuonME0Hits",&b_recoMuon_numberOfValidMuonME0Hits,"recoMuon_numberOfValidMuonME0Hits/I");

  recottree_->Branch("recoMuon_TrkIsolation03",&b_recoMuon_TrkIso03,"recoMuon_TrkIsolation03/F");
  recottree_->Branch("recoMuon_TrkIsolation05",&b_recoMuon_TrkIso05,"recoMuon_TrkIsolation05/F");
  recottree_->Branch("recoMuon_PFIsolation04",&b_recoMuon_PFIso04,"recoMuon_PFIsolation04/F");
  recottree_->Branch("recoMuon_PFIsolation03",&b_recoMuon_PFIso03,"recoMuon_PFIsolation03/F");
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
    b_genMuon_isTight = false;
    b_genMuon_isMedium = false;
    b_genMuon_isLoose = false;
    b_genMuon_noRecHitGEM = -1;
    b_genMuon_isME0Muon = false;
    b_genMuon_isGEMMuon = false;
    b_genMuon_isMuon = false;
    b_genMuon_numberOfValidMuonGEMHits = -1;
    b_genMuon_numberOfValidMuonME0Hits = -1;

    if ( simToMuonColl.find(simRef) != simToMuonColl.end() ) {
      vector<pair<RefToBase<Muon>, double> > MuRefV = simToMuonColl[simRef];      
      if ( !MuRefV.empty()) {
	const Muon* mu = MuRefV.begin()->first.get();
	signalMuons.push_back(mu);

        b_genMuon_TrkIso03 = mu->isolationR03().sumPt/mu->pt();
        b_genMuon_TrkIso05 = mu->isolationR05().sumPt/mu->pt();
        b_genMuon_pfIso03 = (mu->pfIsolationR03().sumChargedHadronPt + TMath::Max(0.,mu->pfIsolationR03().sumNeutralHadronEt + mu->pfIsolationR03().sumPhotonEt - 0.5*mu->pfIsolationR03().sumPUPt))/mu->pt();
        b_genMuon_pfIso04 = (mu->pfIsolationR04().sumChargedHadronPt + TMath::Max(0.,mu->pfIsolationR04().sumNeutralHadronEt + mu->pfIsolationR04().sumPhotonEt - 0.5*mu->pfIsolationR04().sumPUPt))/mu->pt();

	b_genMuon_isTight = isTightMuonCustom(*mu, pv0);
	b_genMuon_isMedium = muon::isMediumMuon(*mu);
	b_genMuon_isLoose = muon::isLooseMuon(*mu);
	b_genMuon_isME0Muon = mu->isME0Muon();
	b_genMuon_isGEMMuon = mu->isGEMMuon();
	b_genMuon_isMuon = mu->isMuon();
	const reco::Track* muonTrack = 0;  
	if ( mu->globalTrack().isNonnull() ) muonTrack = mu->globalTrack().get();
	else if ( mu->outerTrack().isNonnull()  ) muonTrack = mu->outerTrack().get();
	if (muonTrack){
	  b_genMuon_noRecHitGEM = nGEMhit(mu);
	  b_genMuon_numberOfValidMuonGEMHits = muonTrack->hitPattern().numberOfValidMuonGEMHits();
	  b_genMuon_numberOfValidMuonME0Hits = muonTrack->hitPattern().numberOfValidMuonME0Hits();
	}
	//if ( b_genMuon.Eta()>2.4 ){ cout << fabs(b_genMuon.Eta()) << "  " << mu->isME0Muon() << endl; } 
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
    
    //Loose
    b_recoMuon_global = mu->isGlobalMuon();
    b_recoMuon_pf = mu->isPFMuon();

    //Medium
    if ( mu->globalTrack().isNonnull() ){
      b_recoMuon_chi2 = mu->globalTrack()->normalizedChi2();
    }
    b_recoMuon_chi2pos = mu->combinedQuality().chi2LocalPosition;
    b_recoMuon_trkKink = mu->combinedQuality().trkKink;
    b_recoMuon_segcompati = muon::segmentCompatibility(*mu);

    //Tight
    b_recoMuon_nstations = mu->numberOfMatchedStations();
    if ( mu->globalTrack().isNonnull() ){
      b_recoMuon_nglobalhits = mu->globalTrack()->hitPattern().numberOfValidMuonHits();
    }
    if ( mu->muonBestTrack().isNonnull() ){
      b_recoMuon_trackdxy = fabs(mu->muonBestTrack()->dxy(pv0.position()));
      b_recoMuon_trackdz = fabs(mu->muonBestTrack()->dz(pv0.position()));
    }
    if ( mu->innerTrack().isNonnull() ){
      b_recoMuon_ninnerhits = mu->innerTrack()->hitPattern().numberOfValidPixelHits();
      b_recoMuon_trackerlayers = mu->innerTrack()->hitPattern().trackerLayersWithMeasurement();
    }

    recottree_->Fill();
  }

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

DEFINE_FWK_MODULE(MuonAnalyser);
