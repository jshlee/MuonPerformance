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
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/Associations/interface/MuonToTrackingParticleAssociator.h"
#include "SimTracker/Common/interface/TrackingParticleSelector.h"

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

class MuonAnalyser : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit MuonAnalyser(const edm::ParameterSet&);
  ~MuonAnalyser();

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

  TTree* ttree_;
  TLorentzVector b_genMuon;
  bool b_genMuon_isTight, b_genMuon_isMedium, b_genMuon_isLoose;
  TLorentzVector b_recoMuon;
  
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

  ttree_ = fs->make<TTree>();
  ttree_->Branch("genMuon", "TLorentzVector", &b_genMuon);
  ttree_->Branch("genMuon_isTight", &b_genMuon_isTight, "genMuon_isTight/O");
  ttree_->Branch("genMuon_isMedium", &b_genMuon_isMedium, "genMuon_isMedium/O");
  ttree_->Branch("genMuon_isLoose", &b_genMuon_isLoose, "genMuon_isLoose/O");
  ttree_->Branch("recoMuon", "TLorentzVector", &b_recoMuon);
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
  const TrackingParticleCollection simColl = *(simHandle.product());

  Handle<View<Muon> > muonHandle;
  iEvent.getByToken(muonToken_, muonHandle);
  View<Muon> muonColl = *(muonHandle.product());
  
  MuonToTrackingParticleAssociator const* assoByHits = nullptr;
  Handle<MuonToTrackingParticleAssociator> associatorBase;
  iEvent.getByToken(muAssocToken_, associatorBase);
  assoByHits = associatorBase.product();

  RefToBaseVector<Muon> Muons;
  for (size_t i = 0; i < muonHandle->size(); ++i) {
      Muons.push_back(muonHandle->refAt(i));
  }
  RefVector<TrackingParticleCollection> allTPs;
  for (size_t i = 0; i < simColl.size(); ++i) {
      allTPs.push_back(TrackingParticleRef(simHandle,i));
  }
  
  reco::MuonToSimCollection muonToSimColl;
  reco::SimToMuonCollection simToMuonColl;  
  assoByHits->associateMuons(muonToSimColl, simToMuonColl, Muons, reco::GlobalTk, allTPs);

  for(TrackingParticleCollection::size_type i=0; i<simColl.size(); i++) {
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
	if ( muon::isTightMuon(*mu, pv0) ) b_genMuon_isTight = true;
	if ( muon::isMediumMuon(*mu) )  b_genMuon_isMedium = true;
	if ( muon::isLooseMuon(*mu) )  b_genMuon_isLoose = true;
	
      }
    }
    ttree_->Fill();
  }
  
      // for (TrackingParticleCollection::size_type i=0; i<trackingParticles->size(); i++){
      //   TrackingParticleRef tpr(trackingParticles, i);
  //   TrackingParticle* tp=const_cast<TrackingParticle*>(tpr.get());
  //   TrackingParticle::Vector momentumTP;
  //   TrackingParticle::Point vertexTP;

  //   if (abs(simTP->pdgId()) != 13) continue;

  //   bool SignalMuon = false;
  //   if(simTP->status() != -99){
  //     //==== Pythia8 gen status : home.thep.lu.se/~torbjorn/pythia81html/ParticleProperties.html
  //     if ((*simTP->genParticle_begin())->numberOfMothers()>0)  {
  // 	if ((*simTP->genParticle_begin())->mother()->numberOfMothers()>0){
  // 	}
  //     }
  //     if ( ( (simTP->status()==1) && ( (*simTP->genParticle_begin())->numberOfMothers()==0 ) )  ||
  // 	   ( (simTP->status()==1) )      )    SignalMuon=true;
  //   } // END if(simTP->status() != -99)        
    
  //   if(SignalMuon){
      



  //   }
  // }

}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonAnalyser);


