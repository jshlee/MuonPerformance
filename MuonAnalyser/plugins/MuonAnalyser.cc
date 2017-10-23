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
#include "DataFormats/GeometryCommonDetAlgo/interface/GlobalErrorBase.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/ErrorFrameTransformer.h"
#include "RecoParticleFlow/PFProducer/interface/PFMuonAlgo.h"

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/Associations/interface/MuonToTrackingParticleAssociator.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimMuon/MCTruth/interface/MuonToSimAssociatorByHits.h"
#include "SimTracker/Common/interface/TrackingParticleSelector.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/GEMGeometry/interface/ME0Geometry.h"
#include "Geometry/GEMGeometry/interface/ME0Chamber.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMSuperChamber.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

#include <TTree.h>
#include <TFile.h>
#include <TLorentzVector.h>

using namespace std;
using namespace reco;
using namespace edm;

class MuonAnalyser : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
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
  bool isTightMod(const reco::VertexCollection* vertices, const SimVertex &simPVh, const reco::Muon *muon,
		  bool useIPxy, bool useIPz, bool debug);
  
  void collectTMVAvalues(const reco::Muon& mu, reco::Vertex pv0);
  int nGEMhit(const reco::Muon * mu) const;
  int nME0hit(const reco::Muon * mu) const;
  
  bool isME0MuonSelNew(reco::Muon muon, double dEtaCut, double dPhiCut, double dPhiBendCut);
  void setBranches(TTree *tree);
  void fillBranches(TTree *tree, TLorentzVector tlv, edm::RefToBase<reco::Muon> muref, bool isSignal, bool isFromLeadingVtx, int pdgId, int mpdgId);
  
private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

  TH1D* h_nevents;
  TH1D* h_vertex;
  
  TH1D* h_vertexBS;
  TH1D* h_vertex1D;
  TH1D* h_vertex1DBS;
  TH1D* h_vertex4D;
  TH1D* h_vertex4DBS;
  TH1D* h_gemStation;
  
  int b_nvertex;

  TTree* genttree_;
  TTree* recottree_;

  int b_pu_density, b_pu_numInteractions;
  
  TLorentzVector b_muon;
  bool b_muon_signal, b_muon_leadingVtx;
  int b_muon_pdgId;
  int b_muon_no;
 

  //Resolution #1/4
  float b_muon_DZresolution;
  
  float b_muon_Etaresolution;
  float b_muon_Phiresolution;
  float b_muon_pTresolution, b_muon_pTinvresolution;
   
  float b_muon_DXDYresolution;
  float b_muon_Charge;
  
  
  bool b_muon_isTightOptimized, b_muon_isTightCustom, b_muon_isTight, b_muon_isMedium, b_muon_isLoose;
  bool b_muon_isME0Muon, b_muon_isME0MuonLoose, b_muon_isME0MuonMedium, b_muon_isME0MuonTight;
  bool b_muon_isGEMMuon, b_muon_isGE11Muon, b_muon_isGE21Muon, b_muon_isRPCMuon, b_muon_isCaloMuon, b_muon_isTrackerMuon;
  bool b_muon_isGlobalMuon, b_muon_isStandAloneMuon, b_muon_isPFMuon;
  bool b_muon_isLooseMod;
  bool b_muon_isTightModNoIP, b_muon_isTightModIPxy, b_muon_isTightModIPz, b_muon_isTightModIPxyz;

  /*
  bool b_muon_tracker; bool b_muon_global; bool b_muon_pf;
  float b_muon_chi2pos; float b_muon_trkKink; float b_muon_segmentCompatibility;
  float b_muon_chi2; int b_muon_nglobalhits; int b_muon_nstations;
  float b_muon_trackdxy; float b_muon_trackdz;
  int b_muon_ninnerhits; float b_muon_trackerlayers; 
  float b_muon_innerquality; float b_muon_caloCompatibility;
  */
  float b_muon_poszPV0, b_muon_poszSimPV, b_muon_poszMuon;

  float b_muon_ME0segX, b_muon_ME0chamX;

  float b_muon_ME0deltaX, b_muon_ME0deltaY, b_muon_ME0deltaZ, b_muon_ME0deltaDXDY, b_muon_ME0deltaDXDZ, b_muon_ME0deltaDYDZ, b_muon_ME0pullX, b_muon_ME0pullY, b_muon_ME0dPhi, b_muon_ME0dEta, b_muon_ME0pullPhi, b_muon_ME0dPhiBend;
  int b_muon_ME0noRecHit;
  float b_muon_GE11deltaX, b_muon_GE11deltaY, b_muon_GE11deltaDXDZ, b_muon_GE11deltaDYDZ, b_muon_GE11pullX, b_muon_GE11pullY, b_muon_GE11dPhi, b_muon_GE11dEta, b_muon_GE11pullPhi;
  int b_muon_GE11noRecHit;
  float b_muon_GE21deltaX, b_muon_GE21deltaY, b_muon_GE21deltaDXDZ, b_muon_GE21deltaDYDZ, b_muon_GE21pullX, b_muon_GE21pullY, b_muon_GE21dPhi, b_muon_GE21dEta, b_muon_GE21pullPhi;
  int b_muon_GE21noRecHit;
  
  float b_muon_PFIso04; float b_muon_PFIso03;
  float b_muon_PFIso03ChargedHadronPt, b_muon_PFIso03NeutralHadronEt;
  float b_muon_PFIso03PhotonEt, b_muon_PFIso03PUPt;
  float b_muon_PFIso04ChargedHadronPt, b_muon_PFIso04NeutralHadronEt;
  float b_muon_PFIso04PhotonEt, b_muon_PFIso04PUPt;
  float b_muon_TrkIso05; float b_muon_TrkIso03;
  float b_muon_puppiIso, b_muon_puppiIsoNoLep;
  float b_muon_puppiIso_ChargedHadron, b_muon_puppiIso_NeutralHadron, b_muon_puppiIso_Photon;  
  float b_muon_puppiIsoNoLep_ChargedHadron, b_muon_puppiIsoNoLep_NeutralHadron, b_muon_puppiIsoNoLep_Photon;  
  bool b_muon_isMuon;
  int b_muon_numberOfValidMuonGEMHits, b_muon_numberOfValidMuonME0Hits;

  float b_muon_ipxySim, b_muon_ipzSim;
  float b_muon_tmva_bdt, b_muon_tmva_mlp;
  float b_muon_tmva_me0bdt;
  
  float b_muon_istracker, b_muon_isglobal, b_muon_ispf;
  float b_muon_chi2, b_muon_chi2pos, b_muon_trkKink, b_muon_segmentCompatibility, b_muon_nstations, b_muon_nglobalhits;
  float b_muon_trackdxy, b_muon_trackdz, b_muon_ninnerhits, b_muon_trackerlayers, b_muon_innerquality, b_muon_caloCompatibility;

  int b_muon_mpdgId;

  TMVA::Reader* bdt_;
  TMVA::Reader* mlp_;
  TMVA::Reader* me0bdt_;

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
  edm::EDGetTokenT<std::vector<PileupSummaryInfo>> putoken;
  edm::EDGetTokenT<edm::View<reco::Muon> > muonToken_;
  edm::EDGetTokenT<reco::MuonToTrackingParticleAssociator> muAssocToken_;
  edm::EDGetTokenT <std::vector< pat::PackedCandidate> > tokenPackedCandidate ;
  edm::EDGetTokenT<edm::ValueMap<float> > PUPPIIsolation_charged_hadrons_;
  edm::EDGetTokenT<edm::ValueMap<float> > PUPPIIsolation_neutral_hadrons_;
  edm::EDGetTokenT<edm::ValueMap<float> > PUPPIIsolation_photons_;
  edm::EDGetTokenT<edm::ValueMap<float> > PUPPINoLeptonsIsolation_charged_hadrons_;
  edm::EDGetTokenT<edm::ValueMap<float> > PUPPINoLeptonsIsolation_neutral_hadrons_;
  edm::EDGetTokenT<edm::ValueMap<float> > PUPPINoLeptonsIsolation_photons_;

  string tmvaWeight_;
  string tmvaWeightme0_;

  const GEMGeometry* gemGeo_;
  const ME0Geometry* me0Geo_;

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
  
  //tmvaWeight_   = (pset.getParameter<string>("tmvaWeightLabel"));
  tmvaWeight_ = edm::FileInPath(pset.getParameter<std::string>("tmvaWeightLabel")).fullPath();
  tmvaWeightme0_ = edm::FileInPath(pset.getParameter<std::string>("tmvaWeightLabelme0")).fullPath();
  simToken_ = consumes<TrackingParticleCollection>(pset.getParameter<InputTag>("simLabel"));
  simVertexToken_ = consumes<std::vector<SimVertex> >(pset.getParameter<edm::InputTag> ("simVertexCollection"));  
  putoken = consumes<std::vector<PileupSummaryInfo>>(edm::InputTag("addPileupInfo"));
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
					 tpset.getParameter<double>("ptMax"),
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

  // tmva booking
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

  me0bdt_ = new TMVA::Reader();
  me0bdt_->AddVariable("muon_ME0deltaX", &b_muon_ME0deltaX);  
  me0bdt_->AddVariable("muon_ME0deltaY", &b_muon_ME0deltaY);  
  me0bdt_->AddVariable("muon_ME0pullX", &b_muon_ME0pullX);  
  me0bdt_->AddVariable("muon_ME0pullY", &b_muon_ME0pullY);  
  //me0bdt_->AddVariable("muon_ME0dPhiBend", &b_muon_ME0dPhiBend);  
  me0bdt_->AddVariable("muon_ME0dPhi", &b_muon_ME0dPhi);  
  me0bdt_->AddVariable("muon_ME0dEta", &b_muon_ME0dEta);  
  //me0bdt_->AddVariable("muon_ME0pullPhi", &b_muon_ME0pullPhi);  
  me0bdt_->AddVariable("muon_ME0deltaDXDZ", &b_muon_ME0deltaDXDZ);  
  me0bdt_->AddVariable("muon_ME0deltaDYDZ", &b_muon_ME0deltaDYDZ);  
  //me0bdt_->AddVariable("muon_ME0noRecHit", &b_muon_ME0noRecHit);  
  me0bdt_->BookMVA("BDT", tmvaWeightme0_);
  
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  h_nevents = fs->make<TH1D>("nevents", "nevents", 1, 0, 1);
  h_vertex = fs->make<TH1D>("vertex reco vs sim", "vertex reco vs sim", 100, -1, 1);
  
  h_vertexBS   = fs->make<TH1D>("vertex reco with BS vs sim", "vertex reco with BS vs sim", 100, -1, 1);
  h_vertex1D   = fs->make<TH1D>("vertex reco 1D vs sim", "vertex reco 1D vs sim", 100, -1, 1);
  h_vertex1DBS = fs->make<TH1D>("vertex reco 1D with BS vs sim", "vertex reco 1D with BS vs sim", 100, -1, 1);
  h_vertex4D   = fs->make<TH1D>("vertex reco 4D vs sim", "vertex reco 4D vs sim", 100, -1, 1);
  h_vertex4DBS = fs->make<TH1D>("vertex reco 4D with BS vs sim", "vertex reco 4D with BS vs sim", 100, -1, 1);
  h_gemStation = fs->make<TH1D>("gemstation", "gemstation", 20, -10, 10);
  genttree_ = fs->make<TTree>("gen", "gen");
  setBranches(genttree_);
  recottree_ = fs->make<TTree>("reco", "reco");
  setBranches(recottree_);

}
MuonAnalyser::~MuonAnalyser(){}
void MuonAnalyser::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  h_nevents->Fill(0.5);

  edm::ESHandle<GEMGeometry> gemg;
  iSetup.get<MuonGeometryRecord>().get(gemg);
  gemGeo_ = &*gemg;

  edm::ESHandle<ME0Geometry> me0g;
  iSetup.get<MuonGeometryRecord>().get(me0g);
  me0Geo_ = &*me0g;

  Handle<VertexCollection> vertices;
  Handle<VertexCollection> vertices1D, vertices1DBS, vertices4D, vertices4DBS, verticesBS;
  
  iEvent.getByToken(vtxToken_, vertices); 
  iEvent.getByToken(vtx1DToken_,   vertices1D); 
  iEvent.getByToken(vtx1DBSToken_, vertices1DBS); 
  iEvent.getByToken(vtx4DToken_,   vertices4D); 
  iEvent.getByToken(vtx4DBSToken_, vertices4DBS); 
  iEvent.getByToken(vtxBSToken_,   verticesBS); 
  
  b_nvertex = vertices->size();
  
  vertexes_ = vertices.product();
  if (vertexes_->empty()) { cout << "no PV" << endl; return; }

  Handle<std::vector<SimVertex> > simVertexCollection;
  iEvent.getByToken(simVertexToken_, simVertexCollection);  
  simVertex_ = simVertexCollection->at(0);

  edm::Handle<std::vector <PileupSummaryInfo> > PupInfo;
  iEvent.getByToken(putoken, PupInfo);
  b_pu_density = 0; b_pu_numInteractions = 0;
  std::vector<PileupSummaryInfo>::const_iterator ipu;
  for (ipu = PupInfo->begin(); ipu != PupInfo->end(); ++ipu) {
    if ( ipu->getBunchCrossing() != 0 ) continue; // storing detailed PU info only for BX=0
    for (unsigned int i=0; i<ipu->getPU_zpositions().size(); ++i) {
      if ( abs((ipu->getPU_zpositions())[i] - simVertex_.position().z()) < 0.1 )
	++b_pu_density;
      ++b_pu_numInteractions;
    }
  }
  
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
  b_muon_no = 0;
  for (TrackingParticleCollection::size_type i=0; i<simHandle->size(); i++) {
    TrackingParticleRef simRef(simHandle, i);
    const TrackingParticle* simTP = simRef.get();
    //if (simTP->pt()<5||abs(simTP->eta())>2.4) continue;
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

    int mpdgId = 0;
    const reco::GenParticle& gen = *(simTP->genParticles()[0]);
    if ( gen.mother()) { mpdgId = gen.mother()->pdgId(); }

    fillBranches(genttree_, gentlv, muonRef, true, true, simTP->pdgId(), mpdgId);
  }

  
  // reco muon loop - for fake rate
  b_muon_no = 0;
  for (size_t i = 0; i < muonHandle->size(); ++i) {    
    edm::RefToBase<reco::Muon> muRef = muonHandle->refAt(i);
    const Muon* mu = muRef.get();
    //if (mu->pt()<5||abs(mu->eta())>2.4) continue;
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
    bool isFromLeadingVtx = false;
    if (simTP){
      isSignalMuon = abs(simTP->pdgId())==13;
      isFromLeadingVtx = !simTP->genParticles().empty() && (simTP->eventId().event() == 0) && (simTP->eventId().bunchCrossing() == 0);
    }

    int mpdgId = 0;
    if (simTP){
      if ( !simTP->genParticles().empty()) {
        const reco::GenParticle& gen = *(simTP->genParticles()[0]);
        if ( gen.mother() ) { mpdgId = gen.mother()->pdgId(); }//cout << mpdgId << endl; }
      }
    }

    fillBranches(recottree_, recotlv, muRef, isSignalMuon, isFromLeadingVtx, pdgId, mpdgId);
  }

}

void MuonAnalyser::fillBranches(TTree *tree, TLorentzVector tlv, edm::RefToBase<reco::Muon> muref, bool isSignal, bool isFromLeadingVtx, int pdgId, int mpdgId)
{
  b_muon = tlv;
  b_muon_signal = isSignal;
  b_muon_leadingVtx = isFromLeadingVtx;
  b_muon_pdgId = pdgId;
  b_muon_mpdgId = mpdgId;
  ++b_muon_no;
  reco::Vertex pv0 = vertexes_->at(0);

  //Resolution #2/4
  b_muon_DZresolution = -999;
  b_muon_Etaresolution = -999;
  b_muon_Phiresolution = -999;
  b_muon_pTresolution = -999; b_muon_pTinvresolution = -999;
  b_muon_DXDYresolution = -999;
  b_muon_Charge = -999;


  b_muon_isTightOptimized = -999; b_muon_isTightCustom = -999; b_muon_isTight = -999; b_muon_isMedium = -999; b_muon_isLoose = -999;
  b_muon_isME0Muon = -999; b_muon_isME0MuonLoose = -999; b_muon_isME0MuonMedium = -999; b_muon_isME0MuonTight = -999;
  b_muon_isGEMMuon = -999; b_muon_isGE11Muon = -999; b_muon_isGE21Muon = -999; b_muon_isRPCMuon = -999; b_muon_isCaloMuon = -999; b_muon_isTrackerMuon = -999;
  b_muon_isStandAloneMuon = -999;
  b_muon_isGlobalMuon = -999; b_muon_isStandAloneMuon = -999; b_muon_isPFMuon = -999;
  b_muon_isLooseMod = -999;
  b_muon_isTightModNoIP = -999; b_muon_isTightModIPxy = -999; b_muon_isTightModIPz = -999; b_muon_isTightModIPxyz = -999;

  b_muon_ME0segX = 100; b_muon_ME0chamX = 100;

  b_muon_ME0deltaX = 100; b_muon_ME0deltaY = 100; b_muon_ME0deltaZ = 100; b_muon_ME0deltaDXDZ = 100; b_muon_ME0deltaDYDZ = 100; b_muon_ME0deltaDXDY = 100;  b_muon_ME0noRecHit = 100; b_muon_ME0pullX = 100; b_muon_ME0pullY = 100; b_muon_ME0dPhi = 100; b_muon_ME0dEta = 100; b_muon_ME0pullPhi = 100; b_muon_ME0dPhiBend = 100;
  b_muon_GE11deltaX = 100; b_muon_GE11deltaY = 100; b_muon_GE11deltaDXDZ = 100; b_muon_GE11deltaDYDZ = 100; b_muon_GE11noRecHit = 100; b_muon_GE11pullX = 100; b_muon_GE11pullY = 100; b_muon_GE11dPhi = 100; b_muon_GE11dEta = 100; b_muon_GE11pullPhi = 100;
  b_muon_GE21deltaX = 100; b_muon_GE21deltaY = 100; b_muon_GE21deltaDXDZ = 100; b_muon_GE21deltaDYDZ = 100; b_muon_GE21noRecHit = 100; b_muon_GE21pullX = 100; b_muon_GE21pullY = 100; b_muon_GE21dPhi = 100; b_muon_GE21dEta = 100; b_muon_GE21pullPhi = 100;

  b_muon_istracker = -999;  b_muon_isglobal = -999;  b_muon_ispf = -999;
  b_muon_chi2pos = -999;  b_muon_trkKink = -999;  b_muon_segmentCompatibility = -999;
  b_muon_chi2 = -999;  b_muon_nglobalhits = -999;  b_muon_nstations = -999;
  b_muon_trackdxy = -999;  b_muon_trackdz = -999;
  b_muon_ninnerhits = -999;  b_muon_trackerlayers = -999;
  b_muon_innerquality = -999; b_muon_caloCompatibility = -999;
  b_muon_poszPV0 = -999; b_muon_poszSimPV = -999; b_muon_poszMuon = -999;
  b_muon_PFIso04 = -999;  b_muon_PFIso03 = -999;
  b_muon_PFIso03ChargedHadronPt = -999; b_muon_PFIso03NeutralHadronEt = -999;
  b_muon_PFIso03PhotonEt = -999; b_muon_PFIso03PUPt = -999;
  b_muon_PFIso04ChargedHadronPt = -999; b_muon_PFIso04NeutralHadronEt = -999;
  b_muon_PFIso04PhotonEt = -999; b_muon_PFIso04PUPt = -999;
  b_muon_TrkIso05 = -999;  b_muon_TrkIso03 = -999;
  b_muon_puppiIso = -999; b_muon_puppiIso_ChargedHadron = -999; b_muon_puppiIso_NeutralHadron = -999; b_muon_puppiIso_Photon = -999;
  b_muon_puppiIsoNoLep = -999; b_muon_puppiIsoNoLep_ChargedHadron = -999; b_muon_puppiIsoNoLep_NeutralHadron = -999; b_muon_puppiIsoNoLep_Photon = -999;  
  b_muon_isMuon = -999;
  b_muon_numberOfValidMuonGEMHits = -999; b_muon_numberOfValidMuonME0Hits = -999;

  b_muon_ipxySim = -999; b_muon_ipzSim = -999;
  b_muon_tmva_bdt = -999; b_muon_tmva_mlp = -999;
  b_muon_tmva_me0bdt = -999;
  
  const Muon* mu = muref.get();
  if (mu){

    //Resolution #3/4
    b_muon_Etaresolution = (b_muon.Eta()-mu->eta())/b_muon.Eta();
    b_muon_pTresolution = (b_muon.Pt()-mu->pt())/b_muon.Pt();
    b_muon_Phiresolution = (b_muon.Phi()-mu->phi())/b_muon.Phi();
    b_muon_Charge = mu->charge();
    
    //b_muon_DZ1resolution = (pv0.position().z()-mu->muonBestTrack()->vz())/pv0.position().z();
    b_muon_DZresolution = (simVertex_.position().z()-mu->muonBestTrack()->dz(pv0.position()))/simVertex_.position().z();
    //b_muon_DXDY1resolution = (pv0.position().xy()-mu->muonBestTrack()->vxy())/pv0.position().xy();
    //b_muon_DXDYresolution =(simVertex_.position().xy()-mu->muonBestTrack()->dxy(pv0.position()))/simVertex_.position().PxPy();


    
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
        
    b_muon_isTightOptimized = isTightMuonCustomOptimized(*mu, pv0);
    b_muon_isTightCustom = isTightMuonCustom(*mu, pv0);
    b_muon_isTight = muon::isTightMuon(*mu, pv0);
    b_muon_isMedium = muon::isMediumMuon(*mu);
    b_muon_isLoose = muon::isLooseMuon(*mu);

    b_muon_isME0Muon = mu->isME0Muon();

    double mom = mu->p();
    double dPhiCut_ = std::min(std::max(1.2/mom,1.2/100),0.056);
    double dPhiBendCut_ = std::min(std::max(0.2/mom,0.2/100),0.0096);
    b_muon_isME0MuonLoose = isME0MuonSelNew(*mu, 0.077, dPhiCut_, dPhiBendCut_);
    
    bool ipxy = false, ipz = false, validPxlHit = false, highPurity = false;
    if (mu->innerTrack().isNonnull()){
      ipxy = std::abs(mu->muonBestTrack()->dxy(pv0.position())) < 0.2;
      ipz = std::abs(mu->muonBestTrack()->dz((pv0.position()))) < 0.5;
      validPxlHit = mu->innerTrack()->hitPattern().numberOfValidPixelHits() > 0;
      highPurity = mu->innerTrack()->quality(reco::Track::highPurity);
    }
    // isMediumME0 - just loose with track requirements for now, this needs to be updated
    b_muon_isME0MuonMedium = isME0MuonSelNew(*mu, 0.077, dPhiCut_, dPhiBendCut_) && ipxy && validPxlHit && highPurity;

    // tighter cuts for tight ME0
    dPhiCut_ = std::min(std::max(1.2/mom,1.2/100),0.032);
    dPhiBendCut_ = std::min(std::max(0.2/mom,0.2/100),0.0041);
    b_muon_isME0MuonTight = isME0MuonSelNew(*mu, 0.048, dPhiCut_, dPhiBendCut_) && ipxy && ipz && validPxlHit && highPurity;
    
    b_muon_isGEMMuon = mu->isGEMMuon();
    b_muon_isRPCMuon = mu->isRPCMuon();
    b_muon_isCaloMuon = mu->isCaloMuon();
    b_muon_isTrackerMuon = mu->isTrackerMuon();
    b_muon_isMuon = mu->isMuon();
    b_muon_isGlobalMuon = mu->isGlobalMuon();
    b_muon_isStandAloneMuon = mu->isStandAloneMuon();
    b_muon_isPFMuon = mu->isPFMuon();
    
    b_muon_isLooseMod = isLooseMod(mu);
    b_muon_isTightModNoIP  = isTightMod(vertexes_, simVertex_, mu, false, false, false);
    b_muon_isTightModIPxy  = isTightMod(vertexes_, simVertex_, mu, true,  false, false);
    b_muon_isTightModIPz   = isTightMod(vertexes_, simVertex_, mu, false, true, false);
    b_muon_isTightModIPxyz = isTightMod(vertexes_, simVertex_, mu, true,  true, false);

    if (b_muon_isTightCustom != b_muon_isTightModIPxyz){
      if (abs(mu->eta()) > 1.4 && abs(mu->eta()) < 2.0) {
	cout << "b_muon_isTightCustom "<< b_muon_isTightCustom
	     << " b_muon_isTightModIPxyz "<< isTightMod(vertexes_, simVertex_, mu, true,  true, true)
	     << endl;
	cout << "isPFMuon "<< mu->isPFMuon()
	     << " isPFMuonMod "<< (isGlobalTightMuon(mu) || isTrackerTightMuon(mu) || isIsolatedMuon(mu))
	     << endl;
      }
    }
    
    float me0SegX = 100;
    float ge11SegX = 100;
    float ge21SegX = 100;

    ErrorFrameTransformer tran;

    for (auto chamber : mu->matches()){
      for (auto segment : chamber.me0Matches){
	if (chamber.detector() == 5){
	  auto me0Segment = (*segment.me0SegmentRef);
	  b_muon_ME0segX = segment.x;
	  b_muon_ME0chamX = chamber.x;
	  me0SegX = abs( chamber.x - segment.x );	  
	  if (me0SegX < abs(b_muon_ME0deltaX)){
	    b_muon_ME0deltaX    = ( chamber.x - segment.x );
	    b_muon_ME0deltaY    = ( chamber.y - segment.y );
	    b_muon_ME0pullX  = (chamber.x - segment.x) / std::sqrt(chamber.xErr + segment.xErr);
	    b_muon_ME0pullY  = (chamber.y - segment.y) / std::sqrt(chamber.yErr + segment.yErr);
	    //b_muon_ME0dPhi = atan(chamber.dXdZ) - atan(segment.dXdZ);	      
	    
       
       b_muon_ME0deltaDXDZ = ( chamber.dXdZ - segment.dXdZ );
	    b_muon_ME0deltaDYDZ = ( chamber.dYdZ - segment.dYdZ );
	    b_muon_ME0noRecHit  = me0Segment.nRecHits();

	    const ME0Chamber * me0Det = me0Geo_->chamber(me0Segment.me0DetId());

	    LocalPoint segLp(segment.x, segment.y, 0);
	    LocalPoint chmLp(chamber.x, chamber.y, 0);
        LocalVector chmLv(chamber.dXdZ, chamber.dYdZ ,1);
	    GlobalPoint segGp = me0Det->toGlobal(segLp);
	    GlobalPoint chmGp = me0Det->toGlobal(chmLp);

	    LocalError segLe(segment.xErr, segment.yErr, 0);
	    LocalError chmLe(chamber.xErr, chamber.yErr, 0);
	    GlobalError segGe = tran.transform(segLe, me0Det->surface());
	    GlobalError chmGe = tran.transform(chmLe, me0Det->surface());

        double segdPhi = me0Segment.deltaPhi();
        double chmdPhi = me0Det->computeDeltaPhi(chmLp, chmLv); 

        b_muon_ME0dPhiBend = segdPhi - chmdPhi;
	    b_muon_ME0pullPhi = ( chmGp.phi() - segGp.phi() ) / std::sqrt(chmGe.phierr(chmGp) + segGe.phierr(segGp) );
	    b_muon_ME0dPhi = deltaPhi(float(chmGp.phi()), float(segGp.phi()));
	    b_muon_ME0dEta = chmGp.eta()- segGp.eta();

	  }
	}
      }
      for (auto segment : chamber.gemMatches){
	if (chamber.detector() == 4){
	  auto gemSegment = (*segment.gemSegmentRef);
	  if (gemSegment.gemDetId().station() == 1){
	    b_muon_isGE11Muon = 1;
	    ge11SegX = abs( chamber.x - segment.x );
	    if (ge11SegX < abs(b_muon_GE11deltaX)){
	      b_muon_GE11deltaX    = ( chamber.x - segment.x );
	      b_muon_GE11deltaY    = ( chamber.y - segment.y );
	      b_muon_GE11pullX  = (chamber.x - segment.x) / std::sqrt(chamber.xErr + segment.xErr);
	      b_muon_GE11pullY  = (chamber.y - segment.y) / std::sqrt(chamber.yErr + segment.yErr);
	      //b_muon_GE11dPhi = atan(chamber.dXdZ) - atan(segment.dXdZ);	      
	      b_muon_GE11deltaDXDZ = ( chamber.dXdZ - segment.dXdZ );
	      b_muon_GE11deltaDYDZ = ( chamber.dYdZ - segment.dYdZ );
	      b_muon_GE11noRecHit  = gemSegment.nRecHits();

	      const GEMSuperChamber * ge11Det = gemGeo_->superChamber(gemSegment.gemDetId());

	      LocalPoint segLp(segment.x, segment.y, 0);
	      LocalPoint chmLp(chamber.x, chamber.y, 0);
	      GlobalPoint segGp = ge11Det->toGlobal(segLp);
	      GlobalPoint chmGp = ge11Det->toGlobal(chmLp);

	      LocalError segLe(segment.xErr, segment.yErr, 0);
	      LocalError chmLe(chamber.xErr, chamber.yErr, 0);
	      GlobalError segGe = tran.transform(segLe, ge11Det->surface());
	      GlobalError chmGe = tran.transform(chmLe, ge11Det->surface());

	      b_muon_GE11pullPhi = ( chmGp.phi() - segGp.phi() ) / std::sqrt(chmGe.phierr(chmGp) + segGe.phierr(segGp) );
	      b_muon_GE11dPhi = deltaPhi(float(chmGp.phi()), float(segGp.phi()));
	      b_muon_GE11dEta = chmGp.eta()- segGp.eta();	      
	    }
	  }
	  if (gemSegment.gemDetId().station() == 2){
	    b_muon_isGE21Muon = 1;
	    ge21SegX = abs( chamber.x - segment.x );	  
	    if (ge21SegX < abs(b_muon_GE21deltaX)){
	      b_muon_GE21deltaX    = ( chamber.x - segment.x );
	      b_muon_GE21deltaY    = ( chamber.y - segment.y );
	      b_muon_GE21pullX  = (chamber.x - segment.x) / std::sqrt(chamber.xErr + segment.xErr);
	      b_muon_GE21pullY  = (chamber.y - segment.y) / std::sqrt(chamber.yErr + segment.yErr);
	      //b_muon_GE21dPhi = atan(chamber.dXdZ) - atan(segment.dXdZ);	      
	      b_muon_GE21deltaDXDZ = ( chamber.dXdZ - segment.dXdZ );
	      b_muon_GE21deltaDYDZ = ( chamber.dYdZ - segment.dYdZ );
	      b_muon_GE21noRecHit  = gemSegment.nRecHits();
	     
	      const GEMSuperChamber * ge21Det = gemGeo_->superChamber(gemSegment.gemDetId());

	      LocalPoint segLp(segment.x, segment.y, 0);
	      LocalPoint chmLp(chamber.x, chamber.y, 0);
	      GlobalPoint segGp = ge21Det->toGlobal(segLp);
	      GlobalPoint chmGp = ge21Det->toGlobal(chmLp);

	      LocalError segLe(segment.xErr, segment.yErr, 0);
	      LocalError chmLe(chamber.xErr, chamber.yErr, 0);
	      GlobalError segGe = tran.transform(segLe, ge21Det->surface());
	      GlobalError chmGe = tran.transform(chmLe, ge21Det->surface());

	      b_muon_GE21pullPhi = ( chmGp.phi() - segGp.phi() ) / std::sqrt(chmGe.phierr(chmGp) + segGe.phierr(segGp) );
	      b_muon_GE21dPhi = deltaPhi(float(chmGp.phi()), float(segGp.phi()));
	      b_muon_GE21dEta = chmGp.eta()- segGp.eta();	      
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

      if (isSignal){
	for(auto i=muonTrack->recHitsBegin(); i!=muonTrack->recHitsEnd(); i++) {
	  DetId hitId = (*i)->geographicalId();
	  if (!(*i)->isValid() ) continue;
	  if (hitId.det()!=DetId::Muon) continue;
	
	  if (hitId.subdetId() == MuonSubdetId::GEM){
	  
	    if ((*i)->recHits().size()){
	      for (auto rh : (*i)->recHits()){
		GEMDetId gemRHId(rh->geographicalId());
		h_gemStation->Fill( gemRHId.region()*(gemRHId.station()*3 + gemRHId.layer()));
	      }
	    }
	    else {
	      GEMDetId gemRHId(hitId);
	      h_gemStation->Fill( gemRHId.region()*(gemRHId.station()*3 + gemRHId.layer()));
	    
	    }
	  }
	}
	//      cout << "(*i)->size()="<< (*i)->recHits().size()<< " det="<< hitId.det() << " subdet=" << hitId.subdetId() <<endl;
      }
      
    }
    
    GlobalPoint point(simVertex_.position().x(), simVertex_.position().y(), simVertex_.position().z());
    b_muon_ipxySim = abs(mu->muonBestTrack()->dxy(math::XYZPoint(point.x(),point.y(),point.z())));
    b_muon_ipzSim = abs(mu->muonBestTrack()->dz(math::XYZPoint(point.x(),point.y(),point.z())));

    collectTMVAvalues(*mu, pv0);
    b_muon_tmva_bdt = bdt_->EvaluateMVA("BDT");

    b_muon_tmva_me0bdt = me0bdt_->EvaluateMVA("BDT");
  }
  tree->Fill();
}


void MuonAnalyser::collectTMVAvalues(const reco::Muon& mu, reco::Vertex pv0) 
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

bool MuonAnalyser::isLooseMuonCustom(const reco::Muon& mu) const
{
  if ( !(mu.isPFMuon()) ) return false;
  if ( !(mu.isGlobalMuon() || mu.isTrackerMuon()) ) return false;
  return true;
}

bool MuonAnalyser::isTightMuonCustom(const reco::Muon& muon, reco::Vertex vtx) const
{
  if ( !muon.isPFMuon() || !muon.isGlobalMuon() ) return false;
  
  bool muID = muon::isGoodMuon(muon,muon::GlobalMuonPromptTight) && (muon.numberOfMatchedStations() > 1);
      
  bool hits = muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 &&
    muon.innerTrack()->hitPattern().numberOfValidPixelHits() > 0; 
  
  bool ip = fabs(muon.muonBestTrack()->dxy(vtx.position())) < 0.2;
  //&& fabs(muon.muonBestTrack()->dz(vtx.position())) < 0.5;
  
  return muID && hits && ip;
}

bool MuonAnalyser::isTightMuonCustomOptimized(const reco::Muon& muon, reco::Vertex vtx) const
{
  if ( !muon.isPFMuon() || !muon.isGlobalMuon() ) return false;
  if ( !(muon.globalTrack().isNonnull()) )  return false;
  if ( !(muon.muonBestTrack().isNonnull()) ) return false;
  if ( !(muon.innerTrack().isNonnull()) ) return false;
  if ( !(muon.globalTrack()->normalizedChi2() < 2.) ) return false; // < 10.
  if ( !(muon.globalTrack()->hitPattern().numberOfValidMuonHits() > 10) ) return false; // > 0
  if ( !(muon.numberOfMatchedStations() > 1) ) return false;
  if ( !(abs(muon.muonBestTrack()->dxy(vtx.position())) < 0.02) ) return false; // < 0.2
  //if ( !(abs(muon.muonBestTrack()->dz(vtx.position())) < 0.5) ) return false;
  if ( !(muon.innerTrack()->hitPattern().numberOfValidPixelHits() > 3) ) return false; // > 0
  if ( !(muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 7) ) return false; // > 5
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

bool MuonAnalyser::isME0MuonSelNew(reco::Muon muon, double dEtaCut, double dPhiCut, double dPhiBendCut)
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
            
	  const ME0Chamber * me0chamber = me0Geo_->chamber(chamber->id);
            
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

bool MuonAnalyser::isTightMod(const reco::VertexCollection* vertices, const SimVertex &simPVh, const reco::Muon *muon, bool useIPxy, bool useIPz, bool debug = false)
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
            
      if(vertices->size() !=0 && vtxCoord[0] > 0.5 && vtxCoord[0] < 1.5){//DY samples
                
	ipz = abs(muon->muonBestTrack()->dz((*vertices)[indexFinal].position())) < 0.5;
	if (debug)
	  std::cout<<"vx: "<<pointDY.x()<<" vy: "<<pointDY.y()<<" vz: "<<pointDY.z()<<" |Dz|: "<<muon->muonBestTrack()->dz((*vertices)[indexFinal].position())<<std::endl;
                
      }
      else if(vtxCoord[0] > 1.5 && vtxCoord[0] < 3.5){
                
	ipz = ipzSimBool;
	if (debug){
	  std::cout<<"vx: "<<point.x()<<" vy: "<<point.y()<<" vz: "<<point.z()<<" |Dz|: "<<ipz<<std::endl;
	  std::cout<<"ipzSimBool: "<<ipzSimBool <<std::endl;
	}
                
      }
            
    }
    else if(useIPz == false) ipz = true;
        
    bool validPxlHit = muon->innerTrack()->hitPattern().numberOfValidPixelHits() > 0;
    //bool validPxlHit = muon->innerTrack()->hitPattern().pixelLayersWithMeasurement(3,2) > 0;
    //bool validPxlHit = muon->innerTrack()->hitPattern().pixelLayersWithMeasurement(4,3) > 0;

    if (debug)
      std::cout<<trkLayMeas<<" "<<isGlb<<" "<<isPF<<" "<<chi2<<" "<<validHits<<" "<<matchedSt<<" "<<ipxy<<" "<<ipz<<" "<<validPxlHit<<std::endl;
        
    if(trkLayMeas && isGlb && isPF && chi2 && validHits && matchedSt && ipxy && ipz && validPxlHit) result = true;
        
  }
    
  return result;
}
void MuonAnalyser::setBranches(TTree *tree)
{
  tree->Branch("nvertex", &b_nvertex, "nvertex/I");
  tree->Branch("pu_density", &b_pu_density, "pu_density/I");
  tree->Branch("pu_numInteractions", &b_pu_numInteractions, "pu_numInteractions/I");
  tree->Branch("muon", "TLorentzVector", &b_muon);  
  tree->Branch("muon_no", &b_muon_no, "muon_no/I");
  tree->Branch("muon_pdgId", &b_muon_pdgId, "muon_pdgId/I");
  tree->Branch("muon_signal", &b_muon_signal, "muon_signal/O");
  tree->Branch("muon_leadingVtx", &b_muon_leadingVtx, "muon_leadingVtx/O");
 
  //Resolution #4/4
  tree->Branch("muon_pTresolution",&b_muon_pTresolution,"muon_pTresolution/F");
  tree->Branch("muon_Phiresolution",&b_muon_Phiresolution,"muon_Phiresolution/F");
  tree->Branch("muon_Etaresolution",&b_muon_Etaresolution,"muon_Etaresolution/F");
  tree->Branch("muon_DXDYresolution",&b_muon_DXDYresolution,"muon_DXDYresolution/F");
  tree->Branch("muon_DZresolution",&b_muon_DZresolution,"muon_DZresolution/F");
  tree->Branch("muon_Charge",&b_muon_Charge,"muon_Charge/F");
  
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
  tree->Branch("muon_isME0MuonLoose", &b_muon_isME0MuonLoose, "muon_isME0MuonLoose/O");
  tree->Branch("muon_isME0MuonMedium", &b_muon_isME0MuonMedium, "muon_isME0MuonMedium/O");
  tree->Branch("muon_isME0MuonTight", &b_muon_isME0MuonTight, "muon_isME0MuonTight/O");
  tree->Branch("muon_isGEMMuon", &b_muon_isGEMMuon, "muon_isGEMMuon/O");
  tree->Branch("muon_isGE11Muon", &b_muon_isGE11Muon, "muon_isGE11Muon/O");
  tree->Branch("muon_isGE21Muon", &b_muon_isGE21Muon, "muon_isGE21Muon/O");
  tree->Branch("muon_isRPCMuon", &b_muon_isRPCMuon, "muon_isRPCMuon/O");
  tree->Branch("muon_isCaloMuon", &b_muon_isCaloMuon, "muon_isCaloMuon/O");
  tree->Branch("muon_isLooseMod", &b_muon_isLooseMod, "muon_isLooseMod/O");
  tree->Branch("muon_isTightOptimized", &b_muon_isTightOptimized, "muon_isTightOptimized/O");
  tree->Branch("muon_isTightCustom", &b_muon_isTightCustom, "muon_isTightCustom/O");
  tree->Branch("muon_isTightModNoIP", &b_muon_isTightModNoIP, "muon_isTightModNoIP/O");
  tree->Branch("muon_isTightModIPxy", &b_muon_isTightModIPxy, "muon_isTightModIPxy/O");
  tree->Branch("muon_isTightModIPz", &b_muon_isTightModIPz, "muon_isTightModIPz/O");
  tree->Branch("muon_isTightModIPxyz", &b_muon_isTightModIPxyz, "muon_isTightModIPxyz/O");
  
  tree->Branch("muon_istracker", &b_muon_istracker, "muon_istracker/F");
  tree->Branch("muon_isglobal", &b_muon_isglobal, "muon_isglobal/F");
  tree->Branch("muon_ispf", &b_muon_ispf, "muon_ispf/F");
  tree->Branch("muon_normalizedChi2", &b_muon_chi2, "muon_normalizedChi2/F");
  tree->Branch("muon_chi2LocalPosition", &b_muon_chi2pos, "muon_chi2LocalPosition/F");
  tree->Branch("muon_trkKink", &b_muon_trkKink, "muon_trkKink/F");
  tree->Branch("muon_segmentCompatibility", &b_muon_segmentCompatibility, "muon_segmentCompatibility/F");
  tree->Branch("muon_numberOfValidMuonHits", &b_muon_nglobalhits, "muon_numberOfValidMuonHits/F");
  tree->Branch("muon_numberOfMatchedStations", &b_muon_nstations, "muon_numberOfMatchedStations/F");
  tree->Branch("muon_pv0pos_dxy", &b_muon_trackdxy, "muon_pv0pos_dxy/F");
  tree->Branch("muon_pv0pos_dz", &b_muon_trackdz, "muon_pv0pos_dz/F");
  tree->Branch("muon_numberOfValidPixelHits", &b_muon_ninnerhits, "muon_numberOfValidPixelHits/F");
  tree->Branch("muon_trackerLayersWithMeasurement", &b_muon_trackerlayers, "muon_trackerLayersWithMeasurement/F");
  tree->Branch("muon_innerquality", &b_muon_innerquality, "muon_innerquality/F");
  tree->Branch("muon_caloCompatibility", &b_muon_caloCompatibility, "muon_caloCompatibility/F");
  tree->Branch("muon_ipxySim", &b_muon_ipxySim, "muon_ipxySim/F");
  tree->Branch("muon_ipzSim", &b_muon_ipzSim, "muon_ipzSim/F");
  tree->Branch("muon_tmva_bdt", &b_muon_tmva_bdt, "muon_tmva_bdt/F");
  tree->Branch("muon_tmva_mlp", &b_muon_tmva_mlp, "muon_tmva_mlp/F");  
  tree->Branch("muon_tmva_me0bdt", &b_muon_tmva_me0bdt, "muon_tmva_me0bdt/F");

  tree->Branch("muon_ME0segX", &b_muon_ME0segX, "muon_ME0segX/F");  
  tree->Branch("muon_ME0chamX", &b_muon_ME0chamX, "muon_ME0chamX/F");  

  //Resolution DXYZ #5
  tree->Branch("muon_ME0deltaX", &b_muon_ME0deltaX, "muon_ME0deltaX/F");  
  tree->Branch("muon_ME0deltaY", &b_muon_ME0deltaY, "muon_ME0deltaY/F");  
  tree->Branch("muon_ME0pullX", &b_muon_ME0pullX, "muon_ME0pullX/F");  
  tree->Branch("muon_ME0pullY", &b_muon_ME0pullY, "muon_ME0pullY/F");  
  tree->Branch("muon_ME0dPhiBend", &b_muon_ME0dPhiBend, "muon_ME0dPhiBend/F");  
  tree->Branch("muon_ME0dPhi", &b_muon_ME0dPhi, "muon_ME0dPhi/F");  
  tree->Branch("muon_ME0dEta", &b_muon_ME0dEta, "muon_ME0dEta/F");  
  tree->Branch("muon_ME0pullPhi", &b_muon_ME0pullPhi, "muon_ME0pullPhi/F");  
  tree->Branch("muon_ME0deltaDXDZ", &b_muon_ME0deltaDXDZ, "muon_ME0deltaDXDZ/F");  
  tree->Branch("muon_ME0deltaDYDZ", &b_muon_ME0deltaDYDZ, "muon_ME0deltaDYDZ/F");  
  tree->Branch("muon_ME0noRecHit", &b_muon_ME0noRecHit, "muon_ME0noRecHit/I");  

  tree->Branch("muon_GE11deltaX", &b_muon_GE11deltaX, "muon_GE11deltaX/F");  
  tree->Branch("muon_GE11deltaY", &b_muon_GE11deltaY, "muon_GE11deltaY/F");  
  tree->Branch("muon_GE11pullX", &b_muon_GE11pullX, "muon_GE11pullX/F");  
  tree->Branch("muon_GE11pullY", &b_muon_GE11pullY, "muon_GE11pullY/F");  
  tree->Branch("muon_GE11dPhi", &b_muon_GE11dPhi, "muon_GE11dPhi/F");  
  tree->Branch("muon_GE11dEta", &b_muon_GE11dEta, "muon_GE11dEta/F");  
  tree->Branch("muon_GE11pullPhi", &b_muon_GE11pullPhi, "muon_GE11pullPhi/F");  
  tree->Branch("muon_GE11deltaDXDZ", &b_muon_GE11deltaDXDZ, "muon_GE11deltaDXDZ/F");  
  tree->Branch("muon_GE11deltaDYDZ", &b_muon_GE11deltaDYDZ, "muon_GE11deltaDYDZ/F");  
  tree->Branch("muon_GE11noRecHit", &b_muon_GE11noRecHit, "muon_GE11noRecHit/I");  
  
  tree->Branch("muon_GE21deltaX", &b_muon_GE21deltaX, "muon_GE21deltaX/F");  
  tree->Branch("muon_GE21deltaY", &b_muon_GE21deltaY, "muon_GE21deltaY/F");  
  tree->Branch("muon_GE21pullX", &b_muon_GE21pullX, "muon_GE21pullX/F");  
  tree->Branch("muon_GE21pullY", &b_muon_GE21pullY, "muon_GE21pullY/F");  
  tree->Branch("muon_GE21dPhi", &b_muon_GE21dPhi, "muon_GE21dPhi/F");  
  tree->Branch("muon_GE21dEta", &b_muon_GE21dEta, "muon_GE21dEta/F");  
  tree->Branch("muon_GE21pullPhi", &b_muon_GE21pullPhi, "muon_GE21pullPhi/F");  
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

  tree->Branch("muon_mpdgId", &b_muon_mpdgId, "muon_mpdgId/I");

}
DEFINE_FWK_MODULE(MuonAnalyser);
