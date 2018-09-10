// cd /cms/ldap_home/iawatson/scratch/GEM/CMSSW_10_1_5/src/ && eval `scramv1 runtime -sh` && eval `scramv1 runtime -sh` && scram b -j 10
// cd ../../.. && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval `scramv1 runtime -sh` && eval `scramv1 runtime -sh` && scram b -j 10
// system include files
#include <memory>
#include <cmath>
#include <iostream>
#include <sstream>
#include<map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartition.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"

#include "DataFormats/GEMDigi/interface/GEMGEBStatusDigiCollection.h"
#include "DataFormats/GEMDigi/interface/GEMVfatStatusDigiCollection.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Run.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TTree.h"

using namespace std;
using namespace edm;

class SliceTestAnalysis : public edm::EDAnalyzer {
public:
  explicit SliceTestAnalysis(const edm::ParameterSet&);
  ~SliceTestAnalysis();

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginJob() override;
  virtual void endJob() override;

  virtual void beginRun(Run const&, EventSetup const&) override;
  virtual void endRun(Run const&, EventSetup const&) override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<GEMRecHitCollection> gemRecHits_;
  edm::EDGetTokenT<GEMGEBStatusDigiCollection> gebStatusCol_;
  edm::EDGetTokenT<GEMVfatStatusDigiCollection> vfatStatusCol_;
  edm::EDGetTokenT<edm::View<reco::Muon> > muons_;
  edm::EDGetTokenT<reco::VertexCollection> vertexCollection_;
  edm::Service<TFileService> fs;

  MuonServiceProxy* theService_;
  edm::ESHandle<Propagator> propagator_;
  edm::ESHandle<TransientTrackBuilder> ttrackBuilder_;
  edm::ESHandle<MagneticField> bField_; 

  static const int MAXCHAMBERS = 36;
  static const int MAXLAYERS = 2;
  
  TH2D* h_firstStrip[MAXCHAMBERS][MAXLAYERS];
  TH2D* h_allStrips[MAXCHAMBERS][MAXLAYERS];
  TH2D* h_globalPosOnGem;
  TH1D* h_clusterSize, *h_totalStrips, *h_bxtotal;
  TH1D* h_inEta[MAXCHAMBERS][MAXLAYERS];
  TH1D* h_hitEta[MAXCHAMBERS][MAXLAYERS];
  TH1D* h_trkEta[MAXCHAMBERS][MAXLAYERS];

  TH1D* h_res_x, *h_res_y, *h_pull_x, *h_pull_y;

  TTree *t_hit;
  int b_run, b_lumi, b_event;
  int b_firstStrip, b_nStrips, b_chamber, b_layer, b_etaPartition, b_muonQuality, b_bx;
  float b_x, b_y, b_z;
  float b_mu_eta, b_mu_phi, b_mu_pt;
  float b_pull_x, b_pull_y, b_res_x, b_res_y;

  int nEvents, nMuonTotal, nGEMFiducialMuon, nGEMTrackWithMuon;
  int m_trkHitFound, m_notrkHitFound, m_noHitFound;
  
  int nTightMuons, nTightMuonsFid, nTightMuonsWithGEM;
  int b_nMuons, b_nMuonsWithGEMHit;
  int b_valid;

  int b_nGEMHits;

  // muon branches
  TH1D* h_mu_all[MAXCHAMBERS][MAXLAYERS];
  TH1D* h_mu_all_gem[MAXCHAMBERS][MAXLAYERS];
  int m_nhits, m_nvalidhits;
  int m_nbounds;
  int m_quality;
  float m_pt, m_eta, m_phi;
  vector<int> m_roll, m_chamber, m_layer; // hit info
  vector<float> m_resx, m_resy, m_pullx, m_pully;
  vector<int> m_in_roll, m_in_chamber, m_in_layer; // propagation bound info
  vector<float> m_in_globalPhi, m_in_globalEta, m_in_nearGemPhi, m_in_nearGemEta; // global info
  
  TTree *t_run;
  TTree *t_muon;
  TTree *t_event;
};

struct GEMMuonAssociation {
  int muonQuality;
  float mu_phi;
  float mu_eta;
  float mu_pt;

  float pull_x;
  float pull_y;
  float res_x;
  float res_y;

  int valid;
};

SliceTestAnalysis::SliceTestAnalysis(const edm::ParameterSet& iConfig) :
  nEvents(0),
  nMuonTotal(0),
  nGEMFiducialMuon(0),
  nGEMTrackWithMuon(0)
{ 
  nTightMuons = 0; nTightMuonsFid = 0; nTightMuonsWithGEM = 0;
  m_trkHitFound = 0; m_notrkHitFound = 0; m_noHitFound = 0;
  gemRecHits_ = consumes<GEMRecHitCollection>(iConfig.getParameter<edm::InputTag>("gemRecHits"));
  gebStatusCol_ = consumes<GEMGEBStatusDigiCollection>(iConfig.getParameter<edm::InputTag>("gebStatusCol"));
  vfatStatusCol_ = consumes<GEMVfatStatusDigiCollection>(iConfig.getParameter<edm::InputTag>("vfatStatusCol"));
  
  muons_ = consumes<View<reco::Muon> >(iConfig.getParameter<InputTag>("muons"));
  vertexCollection_ = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollection"));
  edm::ParameterSet serviceParameters = iConfig.getParameter<edm::ParameterSet>("ServiceParameters");
  theService_ = new MuonServiceProxy(serviceParameters);

  h_clusterSize=fs->make<TH1D>(Form("clusterSize"),"clusterSize",100,0,100);
  h_totalStrips=fs->make<TH1D>(Form("totalStrips"),"totalStrips",200,0,200);
  h_bxtotal=fs->make<TH1D>(Form("bx"),"bx",31,-15,15);

  h_globalPosOnGem = fs->make<TH2D>(Form("onGEM"), "onGEM", 100, -100, 100, 100, -100, 100);

  t_run = fs->make<TTree>("Run", "Run");
  t_run->Branch("run", &b_run, "run/I");
  
  t_event = fs->make<TTree>("Event", "Event");
  t_event->Branch("nMuons", &b_nMuons, "nMuons/I");
  t_event->Branch("nMuonsWithGEMHit", &b_nMuonsWithGEMHit, "nMuonsWithGEMHit/I");
  t_event->Branch("nGEMHits", &b_nGEMHits, "nGEMHits/I");
  t_event->Branch("run", &b_run, "run/I");
  t_event->Branch("lumi", &b_lumi, "lumi/I");

  h_res_x=fs->make<TH1D>(Form("res_x"),"res_x",100,-50,50);
  h_res_y=fs->make<TH1D>(Form("res_y"),"res_y",100,-50,50);
  h_pull_x=fs->make<TH1D>(Form("pull_x"),"pull_x",100,-50,50);
  h_pull_y=fs->make<TH1D>(Form("pull_y"),"pull_y",100,-50,50);

  t_muon = fs->make<TTree>("Muon", "Muon");
  for (int ichamber=0; ichamber<MAXCHAMBERS;++ichamber) {
    for (int ilayer=0; ilayer<MAXLAYERS;++ilayer) {
      h_mu_all[ichamber][ilayer] = fs->make<TH1D>(Form("ieta mu all ch %i lay %i",ichamber, ilayer+1),"mu_all",8,0.5,8.5);
      h_mu_all_gem[ichamber][ilayer] = fs->make<TH1D>(Form("ieta mu all gem ch %i lay %i",ichamber, ilayer+1),"mu_all_gem",8,0.5,8.5);
    }
  }

  t_muon->Branch("nhits", &m_nhits, "nhits/I")->SetTitle("n GEM hits associated to muon");
  t_muon->Branch("nvalidhits", &m_nvalidhits, "nvalidhits/I")->SetTitle("n GEM hits associated to muon, and muon can propagate to eta partition of hit");
  t_muon->Branch("nbounds", &m_nbounds, "nbounds/I")->SetTitle("times muon is in GEM eta partition bounds");
  t_muon->Branch("pt", &m_pt, "pt/F");
  t_muon->Branch("eta", &m_eta, "eta/F");
  t_muon->Branch("phi", &m_phi, "phi/F");
  t_muon->Branch("quality", &m_quality, "quality/I")->SetTitle("muon quality :: 0:noid 1:looseID 2:tightID");
  t_muon->Branch("roll", &m_roll);
  t_muon->Branch("chamber", &m_chamber);
  t_muon->Branch("layer", &m_layer);
  t_muon->Branch("resx", &m_resx);
  t_muon->Branch("resy", &m_resy);
  t_muon->Branch("pullx", &m_pullx);
  t_muon->Branch("pully", &m_pully);
  t_muon->Branch("in_roll", &m_in_roll);
  t_muon->Branch("in_chamber", &m_in_chamber);
  t_muon->Branch("in_layer", &m_in_layer);
  t_muon->Branch("in_globalPhi", &m_in_globalPhi);
  t_muon->Branch("in_globalEta", &m_in_globalEta);
  t_muon->Branch("in_nearGemPhi", &m_in_nearGemPhi);
  t_muon->Branch("in_nearGemEta", &m_in_nearGemEta);

  t_hit = fs->make<TTree>("Hit", "Hit");
  t_hit->Branch("run", &b_run, "run/I");
  t_hit->Branch("lumi", &b_lumi, "lumi/I");  

  t_hit->Branch("bx", &b_bx, "bx/I");
  t_hit->Branch("firstStrip", &b_firstStrip, "firstStrip/I");
  t_hit->Branch("nStrips", &b_nStrips, "nStrips/I");
  t_hit->Branch("chamber", &b_chamber, "chamber/I");
  t_hit->Branch("layer", &b_layer, "layer/I");
  t_hit->Branch("etaPartition", &b_etaPartition, "etaPartition/I");
  t_hit->Branch("muonQuality", &b_muonQuality, "muonQuality/I")->SetTitle("muonQuality -1:none 0:noid 1:looseID 2:tightID");
  t_hit->Branch("x", &b_x, "x/F");
  t_hit->Branch("y", &b_y, "y/F");
  t_hit->Branch("z", &b_z, "z/F");
  t_hit->Branch("pull_x", &b_pull_x, "pull_x/F");
  t_hit->Branch("pull_y", &b_pull_y, "pull_y/F");
  t_hit->Branch("res_x", &b_res_x, "res_x/F");
  t_hit->Branch("res_y", &b_res_y, "res_y/F");
  t_hit->Branch("valid", &b_valid, "valid/I");
  t_hit->Branch("mu_phi", &b_mu_phi, "mu_phi/F");
  t_hit->Branch("mu_eta", &b_mu_eta, "mu_eta/F");
  t_hit->Branch("mu_pt", &b_mu_pt, "mu_pt/F");

  for (int ichamber=0; ichamber<MAXCHAMBERS;++ichamber) {
    // for (int ichamber=27; ichamber<=30;++ichamber) {
    for (int ilayer=0; ilayer<MAXLAYERS;++ilayer) {
      h_firstStrip[ichamber][ilayer] = fs->make<TH2D>(Form("firstStrip ch %i lay %i",ichamber, ilayer+1),"firstStrip",384,1,385,8,0.5,8.5);
      h_firstStrip[ichamber][ilayer]->GetXaxis()->SetTitle("strip");
      h_firstStrip[ichamber][ilayer]->GetYaxis()->SetTitle("iEta");
      
      h_allStrips[ichamber][ilayer] = fs->make<TH2D>(Form("allStrips ch %i lay %i",ichamber, ilayer+1),"allStrips",384,1,385,8,0.5,8.5);
      h_allStrips[ichamber][ilayer]->GetXaxis()->SetTitle("strip");
      h_allStrips[ichamber][ilayer]->GetYaxis()->SetTitle("iEta");

      h_inEta[ichamber][ilayer] = fs->make<TH1D>(Form("inEta ch %i lay %i",ichamber, ilayer+1),"inEta",8,0.5,8.5);
      h_hitEta[ichamber][ilayer] = fs->make<TH1D>(Form("hitEta ch %i lay %i",ichamber, ilayer+1),"hitEta",8,0.5,8.5);
      h_trkEta[ichamber][ilayer] = fs->make<TH1D>(Form("trkEta ch %i lay %i",ichamber, ilayer+1),"trkEta",8,0.5,8.5);
    }
  }
}

SliceTestAnalysis::~SliceTestAnalysis()
{
  std::cout << "::: GEM Slice Test Results :::" << std::endl;
  std::cout << ": From " << nEvents << " events" << std::endl;
  std::cout << std::endl;
  std::cout << " # Muons " << nMuonTotal << std::endl;
  std::cout << " # FidMu " << nGEMFiducialMuon << std::endl;
  std::cout << " # GEMMu " << nGEMTrackWithMuon << std::endl;

  std::cout << " # trk hit found " << m_trkHitFound << std::endl;
  std::cout << " # hit     found " << m_notrkHitFound << std::endl;
  std::cout << " # hit not found " << m_noHitFound << std::endl;
  
  std::cout << " # tight Muons " << nTightMuons << std::endl;
  std::cout << " # tight FidMu " << nTightMuonsFid << std::endl;
  std::cout << " # tight GEMMu " << nTightMuonsWithGEM << std::endl;
  std::cout << std::endl;
}

void
SliceTestAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  nEvents++;

  b_nMuons = 0;
  b_nMuonsWithGEMHit = 0;
  b_nGEMHits = 0;
  
  b_run = iEvent.run();
  b_lumi = iEvent.luminosityBlock();

  edm::ESHandle<GEMGeometry> hGeom;
  iSetup.get<MuonGeometryRecord>().get(hGeom);
  const GEMGeometry* GEMGeometry_ = &*hGeom;

  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",ttrackBuilder_);
  // iSetup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAny",propagator_);
  // iSetup.get<IdealMagneticFieldRecord>().get(bField_); 
  theService_->update(iSetup);
  auto propagator = theService_->propagator("SteppingHelixPropagatorAny");
  
  edm::Handle<GEMRecHitCollection> gemRecHits;  
  iEvent.getByToken(gemRecHits_, gemRecHits);
  //std::cout << "gemRecHits->size() " << gemRecHits->size() <<std::endl;  
  
  edm::Handle<reco::VertexCollection> vertexCollection;
  iEvent.getByToken( vertexCollection_, vertexCollection );
  if(vertexCollection.isValid()) {
    vertexCollection->size();
    //    std::cout << "vertex->size() " << vertexCollection->size() <<std::endl;
  }

  // checking readout status
  edm::Handle<GEMGEBStatusDigiCollection> gebStatusCol;  
  iEvent.getByToken(gebStatusCol_, gebStatusCol);
  
  edm::Handle<GEMVfatStatusDigiCollection> vfatStatusCol;  
  iEvent.getByToken(vfatStatusCol_, vfatStatusCol);

  for (auto ch : GEMGeometry_->chambers()) {
    GEMDetId cId = ch->id();
    auto gebs = gebStatusCol->get(cId); 
    for (auto geb = gebs.first; geb != gebs.second; ++geb) {
      std::cout << "geb id " << cId <<std::endl;
      std::cout << "geb read no. vfats " << int(geb->getVwh())/3
		<< " InFu " << int(geb->getInFu())
		<<std::endl;
    }
    for (auto roll : ch->etaPartitions()) {
      GEMDetId rId = roll->id();
      auto vfats = vfatStatusCol->get(rId); 
      for (auto vfat = vfats.first; vfat != vfats.second; ++vfat) {
	std::cout << rId
		  << " vfat pos " << vfat->position()
		  << " quality " << int(vfat->quality())
		  << " flag " << int(vfat->flag())
		  <<std::endl;

      }
    }
  }
 
  Handle<View<reco::Muon> > muons;
  iEvent.getByToken(muons_, muons);

  std::map<GEMRecHit*,GEMMuonAssociation> muHits;

  for (size_t i = 0; i < muons->size(); ++i) {
    b_nMuons++;
    nMuonTotal++;

    m_nhits = 0;
    m_nvalidhits = 0;
    m_nbounds = 0;
    m_roll.clear(); m_chamber.clear(); m_layer.clear();
    m_resx.clear(); m_resy.clear(); m_pullx.clear(); m_pully.clear();
    m_in_roll.clear(); m_in_chamber.clear(); m_in_layer.clear();
    m_in_globalPhi.clear(); m_in_globalEta.clear(); m_in_nearGemPhi.clear(); m_in_nearGemEta.clear();
    
    edm::RefToBase<reco::Muon> muRef = muons->refAt(i);
    const reco::Muon* mu = muRef.get();

    if (mu->passed(reco::Muon::Selector::CutBasedIdTight)) {
      m_quality = 2;
      nTightMuons++;
    }
    else if (mu->passed(reco::Muon::Selector::CutBasedIdLoose))
      m_quality = 1;
    else
      m_quality = 0;

    if (m_quality != 2) continue;
    if (mu->pt() < 10) continue;
    
    const reco::Track* muonTrack = 0;  
    if ( mu->globalTrack().isNonnull() ) muonTrack = mu->globalTrack().get();
    else if ( mu->outerTrack().isNonnull()  ) muonTrack = mu->outerTrack().get();
    if (muonTrack) {

      std::set<double> detLists;

      reco::TransientTrack ttTrack = ttrackBuilder_->build(muonTrack);
      //bool onDet = false;

      for (auto chamber : GEMGeometry_->chambers()) {

	if (chamber->id().chamber() == 1) continue;
	if (mu->eta() * chamber->id().region() < 0 ) continue;

	TrajectoryStateOnSurface tsos = propagator->propagate(ttTrack.outermostMeasurementState(),
							      chamber->surface());
	//cout <<" " << chamber->id() <<" " << chamber->surface().position().z()<< endl;
	if (!tsos.isValid()) continue;
	GlobalPoint tsosGP = tsos.globalPosition();
	
	// find matching etaPartition
	for (auto etaPart : chamber->etaPartitions()) {
	  auto gemid = etaPart->id();
	  
	  const LocalPoint locPos = etaPart->toLocal(tsosGP);
	  const LocalPoint locPos2D(locPos.x(), locPos.y(), 0);
	  const BoundPlane& bps(etaPart->surface());

	  // checking if muon is within eta partition
	  if (bps.bounds().inside(locPos2D)) {
	    nGEMFiducialMuon++;
	    h_mu_all[gemid.chamber()][gemid.layer()-1]->Fill(gemid.roll());
	    h_globalPosOnGem->Fill(tsosGP.x(), tsosGP.y());
	    m_nbounds++;
	    //onDet = true;
	    h_mu_all[gemid.chamber()][gemid.layer()-1]->Fill(gemid.roll());
	    h_inEta[gemid.chamber()][gemid.layer()-1]->Fill(gemid.roll());
	    m_in_roll.push_back(gemid.roll());
	    m_in_chamber.push_back(gemid.chamber());
	    m_in_layer.push_back(gemid.layer());
	    m_in_globalPhi.push_back(tsosGP.phi());
	    m_in_globalEta.push_back(tsosGP.eta());

	    bool hitFound = false;
	    // finding gem hit within this etaPartition
	    for (auto hit = muonTrack->recHitsBegin(); hit != muonTrack->recHitsEnd(); hit++) {
	      if ((*hit)->geographicalId().det() == 2 && (*hit)->geographicalId().subdetId() == 4) {
		if ((*hit)->rawId() == gemid.rawId() ) {
		  hitFound = true;
		  m_trkHitFound++;
		  if (m_quality == 2 && mu->pt()>20) ++nGEMTrackWithMuon;
		  //GEMDetId gemid((*hit)->geographicalId());
		  //auto etaPart = GEMGeometry_->etaPartition(gemid);
		  h_mu_all_gem[gemid.chamber()][gemid.layer()-1]->Fill(gemid.roll());		  
		}
	      }
	    }

	    float gemEta = +99.0, gemPhi = +99.0;      
	    // if (!hitFound) {
	    // finding gem hit closest gem hit
	    auto recHitsRange = gemRecHits->get(gemid);
	    LocalPoint hitLocPos(-999,-999,-999);
	    cout << "hit gemid "<< gemid
		 <<" gemRecHits->size() " << gemRecHits->size()
		 << endl;

	    auto gemRecHit = recHitsRange.first;
	    for (auto hit = gemRecHit; hit != recHitsRange.second; ++hit) {
	      hitLocPos = hit->localPosition();
	      auto gemGlob = etaPart->toGlobal(hitLocPos);
	      cout << "hit loc "<< locPos
		   << endl;		
	      // pick closest hit
	      if (fabs(gemGlob.phi() - tsosGP.phi()) < gemPhi) {
		gemPhi = gemGlob.phi();
		gemEta = gemGlob.eta();
	      }
	    }
	    if (fabs(locPos.x() - hitLocPos.x()) < 20.0) {
	      m_notrkHitFound++;
	    }
	    else {
	      m_noHitFound++;
	    }
	      
	    cout << "track locPos "<< locPos
		 << " closest hit "<< hitLocPos
		 <<" hitfound " << hitFound
		 << gemid
		 << endl;

	    //   }
	    
	    m_in_nearGemPhi.push_back(gemPhi);
	    m_in_nearGemEta.push_back(gemEta);
	    
	  }
	}
      }
    
      
      // if (onDet) {
      // 	if (m_quality == 2) {
      // 	  ++nTightMuonsFid;
      // 	}
      // }
      
      // if (muonTrack->hitPattern().numberOfValidMuonGEMHits()) {
      // 	++b_nMuonsWithGEMHit;
      // 	//++nGEMTrackWithMuon;
	
      // 	if (m_quality == 2) {
      // 	  ++nTightMuonsWithGEM;
      // 	}
	
      // 	for (auto hit = muonTrack->recHitsBegin(); hit != muonTrack->recHitsEnd(); hit++) {
      // 	  if ( (*hit)->geographicalId().det() == 2 && (*hit)->geographicalId().subdetId() == 4) {
      // 	    GEMDetId gemid((*hit)->geographicalId());
      // 	    auto etaPart = GEMGeometry_->etaPartition(gemid);

      // 	    m_nhits++;

      // 	    TrajectoryStateOnSurface tsos = propagator->propagate(ttTrack.outermostMeasurementState(),etaPart->surface());
      // 	    if (!tsos.isValid()) continue;
      // 	    // GlobalPoint tsosGP = tsos.globalPosition();

      // 	    m_nvalidhits++;

      // 	    LocalPoint && tsos_localpos = tsos.localPosition();
      // 	    LocalError && tsos_localerr = tsos.localError().positionError();
      // 	    LocalPoint && dethit_localpos = (*hit)->localPosition();     
      // 	    LocalError && dethit_localerr = (*hit)->localPositionError();
      // 	    auto res_x = (dethit_localpos.x() - tsos_localpos.x());
      // 	    auto res_y = (dethit_localpos.y() - tsos_localpos.y()); 
      // 	    auto pull_x = (dethit_localpos.x() - tsos_localpos.x()) / 
      // 	      std::sqrt(dethit_localerr.xx() + tsos_localerr.xx());
      // 	    auto pull_y = (dethit_localpos.y() - tsos_localpos.y()) / 
      // 	      std::sqrt(dethit_localerr.yy() + tsos_localerr.yy());
	    
      // 	    h_res_x->Fill(res_x);
      // 	    h_res_y->Fill(res_y);
      // 	    h_pull_x->Fill(pull_x);
      // 	    h_pull_y->Fill(pull_y);

      // 	    m_roll.push_back(gemid.roll());
      // 	    m_chamber.push_back(gemid.chamber());
      // 	    m_layer.push_back(gemid.layer());

      // 	    m_resx.push_back(res_x);
      // 	    m_resy.push_back(res_y);
      // 	    m_pullx.push_back(pull_x);
      // 	    m_pully.push_back(pull_y);
	    
      // 	    int isvalid =  (*hit)->isValid();
      // 	    auto mup = tsos.globalMomentum();
      // 	    muHits[static_cast<GEMRecHit*>(*hit)] = {muonQuality: m_quality,
      // 						     mu_phi: mup.phi(), mu_eta: mup.eta(), mu_pt: mup.perp(),
      // 						     pull_x: pull_x, pull_y: pull_y,
      // 						     res_x: res_x, res_y: res_y, valid: isvalid};
      // 	  }
      // 	}
      //   }
      // }
    
      if (m_nhits > 0 or m_nbounds > 0) {
	m_pt = mu->pt();
	m_eta = mu->eta();
	m_phi = mu->phi();
	t_muon->Fill();
      }
    }
  } // Muon Loop

  int totalStrips = 0;

  for (auto ch : GEMGeometry_->chambers()) {
    for(auto roll : ch->etaPartitions()) {
      GEMDetId rId = roll->id();

      auto recHitsRange = gemRecHits->get(rId); 
      auto gemRecHit = recHitsRange.first;
      for (auto hit = gemRecHit; hit != recHitsRange.second; ++hit) {

	// cout << "hit gemid "<< rId
	//      <<" hit->firstClusterStrip() " << hit->firstClusterStrip()
	//      << endl;

	h_hitEta[rId.chamber()][rId.layer()-1]->Fill(rId.roll());
	h_firstStrip[rId.chamber()][rId.layer()-1]->Fill(hit->firstClusterStrip(), rId.roll());
	h_clusterSize->Fill(hit->clusterSize());
	h_bxtotal->Fill(hit->BunchX());
	for (int nstrip = hit->firstClusterStrip(); nstrip < hit->firstClusterStrip()+hit->clusterSize(); ++nstrip) {
	  totalStrips++;
	  h_allStrips[rId.chamber()][rId.layer()-1]->Fill(nstrip, rId.roll());
	}

	b_firstStrip = hit->firstClusterStrip();
	b_bx = hit->BunchX();
	b_nStrips = hit->clusterSize();
	b_chamber = rId.chamber();
	b_layer = rId.layer();
	b_etaPartition = rId.roll();
	b_muonQuality = -1;
	b_pull_x = -99;
	b_pull_y = -99;
	b_res_x = -99;
	b_res_y = -99;
	b_valid = -1;
	b_mu_phi = b_mu_eta = b_mu_pt = -99;
	for (const auto & kv : muHits) {
	  auto muHit = kv.first;
	  if (*hit == *muHit) {
	    auto assoc = kv.second;
	    b_muonQuality = assoc.muonQuality;
	    b_pull_x = assoc.pull_x;
	    b_pull_y = assoc.pull_y;
	    b_res_x = assoc.res_x;
	    b_res_y = assoc.res_y;
	    b_valid = assoc.valid;
	    b_mu_eta = assoc.mu_eta;
	    b_mu_phi = assoc.mu_phi;
	    b_mu_pt = assoc.mu_pt;
	  }
	}

	auto globalPosition = roll->toGlobal(hit->localPosition());
	b_x = globalPosition.x();
	b_y = globalPosition.y();
	b_z = globalPosition.z();

	t_hit->Fill();
	b_nGEMHits++;
      }
    }
  }
  h_totalStrips->Fill(totalStrips);

  t_event->Fill();
}

void SliceTestAnalysis::beginJob(){}
void SliceTestAnalysis::endJob(){}

void SliceTestAnalysis::beginRun(Run const& run, EventSetup const&){
  b_run = run.run();
  t_run->Fill();
}
void SliceTestAnalysis::endRun(Run const&, EventSetup const&){}

//define this as a plug-in
DEFINE_FWK_MODULE(SliceTestAnalysis);
