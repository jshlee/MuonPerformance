#include <memory>
// user include files
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

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/Associations/interface/MuonToTrackingParticleAssociator.h"
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
  bool b_genMuon_isTight, b_genMuon_isMedium, b_genMuon_isLoose, b_genMuon_isME0Muon, b_genMuon_isGEMMuon;
  int b_genMuon_noRecHitGEM;
  
  TLorentzVector b_recoMuon;
  bool b_recoMuon_signal, b_recoMuon_isTight, b_recoMuon_isMedium, b_recoMuon_isLoose, b_recoMuon_isME0Muon, b_recoMuon_isGEMMuon;
  int b_recoMuon_noChamberMatch;
  int b_recoMuon_noSegment, b_recoMuon_noSegmentDT, b_recoMuon_noSegmentCSC, b_recoMuon_noSegmentRPC, b_recoMuon_noSegmentGEM, b_recoMuon_noSegmentME0;
  int b_recoMuon_noRecHitGEM, b_recoMuon_noRecHitME0;

  bool b_recoMuon_global; bool b_recoMuon_pf;
  float b_recoMuon_chi2pos; float b_recoMuon_trkKink; float b_recoMuon_segcompati;
  float b_recoMuon_chi2; int b_recoMuon_nglobalhits; int b_recoMuon_nstations;
  float b_recoMuon_trackdxy; float b_recoMuon_trackdz;
  int b_recoMuon_ninnerhits; float b_recoMuon_trackerlayers;

  edm::EDGetTokenT<std::vector<reco::Vertex> > vtxToken_;
  //edm::EDGetTokenT<std::vector<PSimHit> > MuonGEMHitsToken_;
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
  //MuonGEMHitsToken_ = consumes<vector<PSimHit> >(pset.getParameter<edm::InputTag>("MuonGEMHits"));
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
  
  recottree_ = fs->make<TTree>("reco", "reco");
  recottree_->Branch("recoMuon", "TLorentzVector", &b_recoMuon);  
  recottree_->Branch("recoMuon_signal", &b_recoMuon_signal, "recoMuon_signal/O");
  recottree_->Branch("recoMuon_isTight", &b_recoMuon_isTight, "recoMuon_isTight/O");
  recottree_->Branch("recoMuon_isMedium", &b_recoMuon_isMedium, "recoMuon_isMedium/O");
  recottree_->Branch("recoMuon_isLoose", &b_recoMuon_isLoose, "recoMuon_isLoose/O");
  recottree_->Branch("recoMuon_isME0Muon", &b_recoMuon_isME0Muon, "recoMuon_isME0Muon/O");
  recottree_->Branch("recoMuon_isGEMMuon", &b_recoMuon_isGEMMuon, "recoMuon_isGEMMuon/O");
  recottree_->Branch("recoMuon_noChamberMatch", &b_recoMuon_noChamberMatch, "recoMuon_noChamberMatch/I");

  recottree_->Branch("recoMuon_noSegment", &b_recoMuon_noSegment, "recoMuon_noSegment/I");
  recottree_->Branch("recoMuon_noSegmentDT", &b_recoMuon_noSegmentDT, "recoMuon_noSegmentDT/I");
  recottree_->Branch("recoMuon_noSegmentCSC", &b_recoMuon_noSegmentCSC, "recoMuon_noSegmentCSC/I");
  recottree_->Branch("recoMuon_noSegmentRPC", &b_recoMuon_noSegmentRPC, "recoMuon_noSegmentRPC/I");
  recottree_->Branch("recoMuon_noSegmentGEM", &b_recoMuon_noSegmentGEM, "recoMuon_noSegmentGEM/I");
  recottree_->Branch("recoMuon_noSegmentME0", &b_recoMuon_noSegmentME0, "recoMuon_noSegmentME0/I");
  recottree_->Branch("recoMuon_noRecHitGEM", &b_recoMuon_noRecHitGEM, "recoMuon_noRecHitGEM/I");
  recottree_->Branch("recoMuon_noRecHitME0", &b_recoMuon_noRecHitME0, "recoMuon_noRecHitME0/I");

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

  // Handle<vector<PSimHit>> MuonGEMHits;
  // iEvent.getByToken(MuonGEMHitsToken_, MuonGEMHits);

  Handle<TrackingParticleCollection> simHandle;
  iEvent.getByToken(simToken_, simHandle);
  //const TrackingParticleCollection simColl = *(simHandle.product());

  Handle<View<Muon> > muonHandle;
  iEvent.getByToken(muonToken_, muonHandle);
  //View<Muon> muonColl = *(muonHandle.product());
  
  MuonToTrackingParticleAssociator const* assoByHits = nullptr;
  Handle<MuonToTrackingParticleAssociator> associatorBase;
  iEvent.getByToken(muAssocToken_, associatorBase);
  assoByHits = associatorBase.product();

  RefToBaseVector<Muon> Muons;
  for (size_t i = 0; i < muonHandle->size(); ++i) {
    Muons.push_back(muonHandle->refAt(i));
  }
  RefVector<TrackingParticleCollection> allTPs;
  for (size_t i = 0; i < simHandle->size(); ++i) {
    allTPs.push_back(TrackingParticleRef(simHandle,i));
  }
  
  reco::MuonToSimCollection muonToSimColl;
  reco::SimToMuonCollection simToMuonColl;  
  assoByHits->associateMuons(muonToSimColl, simToMuonColl, Muons, reco::GlobalTk, allTPs);
  
  vector<const Muon*> signalMuons;signalMuons.clear();
  
  for (TrackingParticleCollection::size_type i=0; i<simHandle->size(); i++) {
    TrackingParticleRef simRef(simHandle, i);
    const TrackingParticle* simTP = simRef.get();
    if ( ! tpSelector_(*simTP) ) continue;

    b_genMuon = TLorentzVector(simRef->momentum().x(), simRef->momentum().y(), simRef->momentum().z(), simRef->energy() );
    
    // GlobalPoint  simVtx(simRef->vertex().x(), simRef->vertex().y(), simRef->vertex().z());
    // GlobalVector simMom(simRef->momentum().x(), simRef->momentum().y(), simRef->momentum().z());
    // const double simDxy = -simVtx.x()*sin(simPhi)+simVtx.y()*cos(simPhi);
    // const double simDz  = simVtx.z() - (simVtx.x()*simMom.x()+simVtx.y()*simMom.y())*simMom.z()/simMom.perp2();    
    // const unsigned int nSimHits = simRef->numberOfHits();
    
    // Get sim-reco association for a simRef
    b_genMuon_isTight = false;
    b_genMuon_isMedium = false;
    b_genMuon_isLoose = false;
    b_genMuon_noRecHitGEM = 0;
    b_genMuon_isME0Muon = false;
    b_genMuon_isGEMMuon = false;
    
    vector<pair<RefToBase<Muon>, double> > MuRefV;
    if ( simToMuonColl.find(simRef) != simToMuonColl.end() ) {
      MuRefV = simToMuonColl[simRef];
      
      if ( !MuRefV.empty()) {
	const Muon* mu = MuRefV.begin()->first.get();
    b_genMuon_noRecHitGEM = nGEMhit(mu);
	signalMuons.push_back(mu);
	if ( isTightMuonCustom(*mu, pv0) ) b_genMuon_isTight = true;
	//if ( muon::isTightMuon(*mu, pv0) ) b_genMuon_isTight = true;
	if ( muon::isMediumMuon(*mu) )  b_genMuon_isMedium = true;
	if ( muon::isLooseMuon(*mu) )  b_genMuon_isLoose = true;	
	if ( mu->isME0Muon() )  b_genMuon_isME0Muon = true;	
	if ( mu->isGEMMuon() )  b_genMuon_isGEMMuon = true;	

	b_genMuon_noRecHitGEM = nGEMhit(mu);
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
    b_recoMuon_isTight = false;
    b_recoMuon_isMedium = false;
    b_recoMuon_isLoose = false;
    b_recoMuon_isME0Muon = false;
    b_recoMuon_isGEMMuon = false;

    for (auto signal : signalMuons){
      if (mu == signal){
	b_recoMuon_signal = true;
	break;
      }
    }
    if (!b_recoMuon_signal){
      //vector<pair<RefToBase<TrackingParticle>, double> > trkRefV;    
      if ( muonToSimColl.find(muRef) != muonToSimColl.end() ) {
	auto trkRefV = muonToSimColl[muRef];
	if ( !trkRefV.empty()) {
	  const TrackingParticle* trkParticle = trkRefV.begin()->first.get();
	  cout << "trkParticle " << trkParticle->pdgId()
	       << " pt = " << trkParticle->pt()
	       << " eta = " << trkParticle->eta()
	       << endl;
	}
      }
    }
    
    if ( isTightMuonCustom(*mu, pv0) ) b_recoMuon_isTight = true;
    //if ( muon::isTightMuon(*mu, pv0) ) b_recoMuon_isTight = true;
    if ( muon::isMediumMuon(*mu) )  b_recoMuon_isMedium = true;
    if ( muon::isLooseMuon(*mu) )  b_recoMuon_isLoose = true;
    if ( mu->isME0Muon() )  b_recoMuon_isME0Muon = true;
    if ( mu->isGEMMuon() )  b_recoMuon_isGEMMuon = true;

    const vector<MuonChamberMatch>& chambers = mu->matches();
    b_recoMuon_noChamberMatch = chambers.size();
    b_recoMuon_noSegment = 0;
    b_recoMuon_noSegmentDT = 0;
    b_recoMuon_noSegmentCSC = 0;
    b_recoMuon_noSegmentRPC = 0;
    b_recoMuon_noSegmentGEM = 0;
    b_recoMuon_noSegmentME0 = 0;

    b_recoMuon_noRecHitGEM = nGEMhit(mu);
    b_recoMuon_noRecHitME0 = 0;

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
	if (chamber->detector() == 5){ //me0
	  ++b_recoMuon_noSegmentME0;
	  auto me0Segment = (*(*segment).me0SegmentRef);
	  b_recoMuon_noRecHitME0 += me0Segment.nRecHits();
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
