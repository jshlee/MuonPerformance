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

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/Associations/interface/MuonToTrackingParticleAssociator.h"
#include "SimTracker/Common/interface/TrackingParticleSelector.h"

#include <TTree.h>
#include <TFile.h>
#include <TLorentzVector.h>

class MuonAnalyser : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit MuonAnalyser(const edm::ParameterSet&);
  ~MuonAnalyser();

  bool isLooseMuonCustom(const reco::Muon& mu) const;
  bool isTightMuonCustom(const reco::Muon& mu, reco::Vertex pv0) const;


private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

  TTree* genttree_;
  TTree* recottree_;
  TLorentzVector b_genMuon;
  bool b_genMuon_isTight, b_genMuon_isMedium, b_genMuon_isLoose;
  TLorentzVector b_recoMuon;
  bool b_recoMuon_signal, b_recoMuon_isTight, b_recoMuon_isMedium, b_recoMuon_isLoose;
  int b_recoMuon_noChamberMatch;
  int b_recoMuon_noSegment, b_recoMuon_noSegmentDT, b_recoMuon_noSegmentCSC, b_recoMuon_noSegmentRPC, b_recoMuon_noSegmentGEM, b_recoMuon_noSegmentME0;

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
  
  recottree_ = fs->make<TTree>("reco", "reco");
  recottree_->Branch("recoMuon", "TLorentzVector", &b_recoMuon);  
  recottree_->Branch("recoMuon_signal", &b_recoMuon_signal, "recoMuon_signal/O");
  recottree_->Branch("recoMuon_isTight", &b_recoMuon_isTight, "recoMuon_isTight/O");
  recottree_->Branch("recoMuon_isMedium", &b_recoMuon_isMedium, "recoMuon_isMedium/O");
  recottree_->Branch("recoMuon_isLoose", &b_recoMuon_isLoose, "recoMuon_isLoose/O");
  recottree_->Branch("recoMuon_noChamberMatch", &b_recoMuon_noChamberMatch, "recoMuon_noChamberMatch/I");

  recottree_->Branch("recoMuon_noSegment", &b_recoMuon_noSegment, "recoMuon_noSegment/I");
  recottree_->Branch("recoMuon_noSegmentDT", &b_recoMuon_noSegmentDT, "recoMuon_noSegmentDT/I");
  recottree_->Branch("recoMuon_noSegmentCSC", &b_recoMuon_noSegmentCSC, "recoMuon_noSegmentCSC/I");
  recottree_->Branch("recoMuon_noSegmentRPC", &b_recoMuon_noSegmentRPC, "recoMuon_noSegmentRPC/I");
  recottree_->Branch("recoMuon_noSegmentGEM", &b_recoMuon_noSegmentGEM, "recoMuon_noSegmentGEM/I");
  recottree_->Branch("recoMuon_noSegmentME0", &b_recoMuon_noSegmentME0, "recoMuon_noSegmentME0/I");

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
    
    vector<pair<RefToBase<Muon>, double> > MuRefV;
    if ( simToMuonColl.find(simRef) != simToMuonColl.end() ) {
      MuRefV = simToMuonColl[simRef];
      
      if ( !MuRefV.empty()) {
	const Muon* mu = MuRefV.begin()->first.get();
	signalMuons.push_back(mu);
	if ( muon::isTightMuon(*mu, pv0) ) b_genMuon_isTight = true;
	if ( muon::isMediumMuon(*mu) )  b_genMuon_isMedium = true;
	if ( muon::isLooseMuon(*mu) )  b_genMuon_isLoose = true;	
      }
    }
    genttree_->Fill();
  }

  for (size_t i = 0; i < muonHandle->size(); ++i) {
    const Muon* mu = muonHandle->refAt(i).get();    
    b_recoMuon = TLorentzVector(mu->momentum().x(), mu->momentum().y(), mu->momentum().z(), mu->energy() );
    b_recoMuon_signal = false;
    b_recoMuon_isTight = false;
    b_recoMuon_isMedium = false;
    b_recoMuon_isLoose = false;
    for (auto signal : signalMuons){
      if (mu == signal){
	b_recoMuon_signal = true;
	break;
      }
    }
    
    if ( muon::isTightMuon(*mu, pv0) ) b_recoMuon_isTight = true;
    if ( muon::isMediumMuon(*mu) )  b_recoMuon_isMedium = true;
    if ( muon::isLooseMuon(*mu) )  b_recoMuon_isLoose = true;

    const vector<MuonChamberMatch>& chambers = mu->matches();
    b_recoMuon_noChamberMatch = chambers.size();
    b_recoMuon_noSegment = 0;
    b_recoMuon_noSegmentDT = 0;
    b_recoMuon_noSegmentCSC = 0;
    b_recoMuon_noSegmentRPC = 0;
    b_recoMuon_noSegmentGEM = 0;
    b_recoMuon_noSegmentME0 = 0;

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
	  ++b_recoMuon_noSegmentGEM;
	}
      }
      
      for ( std::vector<reco::MuonSegmentMatch>::const_iterator segment = chamber->me0Matches.begin(); segment != chamber->me0Matches.end(); ++segment ){
	++b_recoMuon_noSegment;      
	if (chamber->detector() == 5){ //me0
	  ++b_recoMuon_noSegmentME0;
	}
      }
      
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

DEFINE_FWK_MODULE(MuonAnalyser);
