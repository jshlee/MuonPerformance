// cd /cms/ldap_home/iawatson/scratch/GEM/CMSSW_10_1_5/src/ && eval `scramv1 runtime -sh` && eval `scramv1 runtime -sh` && scram b -j 10
// cd ../../.. && eval `scramv1 runtime -sh` && eval `scramv1 runtime -sh` && scram b -j 10
// system include files
#include <memory>
#include <cmath>
#include <iostream>
#include <sstream>

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
  edm::EDGetTokenT<edm::View<reco::Muon> > muons_;
  edm::EDGetTokenT<reco::VertexCollection> vertexCollection_;
  edm::Service<TFileService> fs;

  MuonServiceProxy* theService_;
  edm::ESHandle<Propagator> propagator_;
  edm::ESHandle<TransientTrackBuilder> ttrackBuilder_;
  edm::ESHandle<MagneticField> bField_; 
  
  TH2D* h_firstStrip[36][2];
  TH2D* h_allStrips[36][2];
  TH2D* h_globalPosOnGem;
  TH1D* h_clusterSize, *h_totalStrips, *h_bxtotal;
  TH1D* h_inEta[36][2];
  TH1D* h_hitEta[36][2];
  TH1D* h_trkEta[36][2];

  TTree *t_hit;
  int b_run, b_lumi, b_event;
  int b_firstStrip, b_nStrips, b_chamber, b_layer, b_etaPartition, b_muonQuality;
  float b_x, b_y, b_z;

  int nEvents, nMuonTotal, nGEMFiducialMuon, nGEMTrackWithMuon;
  int b_nMuons, b_nMuonsWithGEMHit;

  int b_nGEMHits;

  TTree *t_run;
  TTree *t_event;
};

SliceTestAnalysis::SliceTestAnalysis(const edm::ParameterSet& iConfig) :
  nEvents(0),
  nMuonTotal(0),
  nGEMFiducialMuon(0),
  nGEMTrackWithMuon(0)
{ 
  gemRecHits_ = consumes<GEMRecHitCollection>(iConfig.getParameter<edm::InputTag>("gemRecHits"));
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

  t_hit = fs->make<TTree>("Hit", "Hit");
  t_hit->Branch("run", &b_run, "run/I");
  t_hit->Branch("lumi", &b_lumi, "lumi/I");
  t_hit->Branch("run", &b_run, "run/I");

  t_event = fs->make<TTree>("Event", "Event");
  t_event->Branch("nMuons", &b_nMuons, "nMuons/I");
  t_event->Branch("nMuonsWithGEMHit", &b_nMuonsWithGEMHit, "nMuonsWithGEMHit/I");
  t_event->Branch("nGEMHits", &b_nGEMHits, "nGEMHits/I");
  

  t_hit->Branch("firstStrip", &b_firstStrip, "firstStrip/I");
  t_hit->Branch("nStrips", &b_nStrips, "nStrips/I");
  t_hit->Branch("chamber", &b_chamber, "chamber/I");
  t_hit->Branch("layer", &b_layer, "layer/I");
  t_hit->Branch("etaPartition", &b_etaPartition, "etaPartition/I");
  t_hit->Branch("muonQuality", &b_muonQuality, "muonQuality/I")->SetTitle("muonQuality -1:none 0:low 1:high");
  t_hit->Branch("x", &b_x, "x/F");
  t_hit->Branch("y", &b_y, "y/F");
  t_hit->Branch("z", &b_z, "z/F");
  
  for (int ichamber=0; ichamber<36;++ichamber) {
  // for (int ichamber=27; ichamber<=30;++ichamber) {
    for (int ilayer=0; ilayer<2;++ilayer) {
      h_firstStrip[ichamber][ilayer] = fs->make<TH2D>(Form("firstStrip ch %i lay %i",ichamber, ilayer),"firstStrip",384,1,385,8,0.5,8.5);
      h_firstStrip[ichamber][ilayer]->GetXaxis()->SetTitle("strip");
      h_firstStrip[ichamber][ilayer]->GetYaxis()->SetTitle("iEta");
      
      h_allStrips[ichamber][ilayer] = fs->make<TH2D>(Form("allStrips ch %i lay %i",ichamber, ilayer),"allStrips",384,1,385,8,0.5,8.5);
      h_allStrips[ichamber][ilayer]->GetXaxis()->SetTitle("strip");
      h_allStrips[ichamber][ilayer]->GetYaxis()->SetTitle("iEta");

      h_inEta[ichamber][ilayer] = fs->make<TH1D>(Form("inEta ch %i lay %i",ichamber, ilayer),"inEta",8,0.5,8.5);
      h_hitEta[ichamber][ilayer] = fs->make<TH1D>(Form("hitEta ch %i lay %i",ichamber, ilayer),"hitEta",8,0.5,8.5);
      h_trkEta[ichamber][ilayer] = fs->make<TH1D>(Form("trkEta ch %i lay %i",ichamber, ilayer),"trkEta",8,0.5,8.5);
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

  edm::Handle<reco::VertexCollection> vertexCollection;
  iEvent.getByToken( vertexCollection_, vertexCollection );
  if(vertexCollection.isValid()) {
    vertexCollection->size();
    //    std::cout << "vertex->size() " << vertexCollection->size() <<std::endl;
  }

  Handle<View<reco::Muon> > muons;
  iEvent.getByToken(muons_, muons);
  //std::cout << "muons->size() " << muons->size() <<std::endl;

  std::vector<GEMRecHit*> goodMuHits;
  std::vector<GEMRecHit*> poorMuHits;

  for (size_t i = 0; i < muons->size(); ++i) {
    b_nMuons++;
    nMuonTotal++;
    
    edm::RefToBase<reco::Muon> muRef = muons->refAt(i);
    const reco::Muon* mu = muRef.get();

    // if (mu->isGEMMuon()) {
    //   std::cout << "isGEMMuon " <<std::endl;
    // }
    
    const reco::Track* muonTrack = 0;  
    if ( mu->globalTrack().isNonnull() ) muonTrack = mu->globalTrack().get();
    else if ( mu->outerTrack().isNonnull()  ) muonTrack = mu->outerTrack().get();
    if (muonTrack) {

      std::set<double> detLists;

      reco::TransientTrack ttTrack = ttrackBuilder_->build(muonTrack);
      bool onDet = false;

      // prepropagate to near the GEM region, to speedup per/etapart prop. w/o loss of generatlity
      TrajectoryStateOnSurface tsos_ = propagator->propagate(ttTrack.outermostMeasurementState(),
      							    GEMGeometry_->etaPartitions()[0]->surface());
      if (!tsos_.isValid()) continue;
      
      // GlobalPoint tsosGP = tsos.globalPosition();
      
      for (auto ch : GEMGeometry_->etaPartitions()) {
	TrajectoryStateOnSurface tsos = propagator->propagate(tsos_, // ttTrack.outermostMeasurementState(),
							      ch->surface());
	if (!tsos.isValid()) continue;
	
	GlobalPoint tsosGP = tsos.globalPosition();
	//if ( !detLists.insert( ch->surface().position().z() ).second ) continue;
	const LocalPoint pos = ch->toLocal(tsosGP);
	const LocalPoint pos2D(pos.x(), pos.y(), 0);
	const BoundPlane& bps(ch->surface());
	//cout << "tsos gp   "<< tsosGP << ch->id() <<endl;
	h_globalPosOnGem->Fill(tsosGP.x(), tsosGP.y());

	if (bps.bounds().inside(pos2D)) {
	  onDet = true;
	  auto gemid = ch->id();
	  h_inEta[gemid.chamber()][gemid.layer()]->Fill(gemid.roll());
	    // cout << " in chamber "<< ch->id() << " pos = "<<pos<< " R = "<<pos.mag() <<" inside "
	  //      <<  bps.bounds().inside(pos2D) <<endl;
	  
	  for (auto hit = muonTrack->recHitsBegin(); hit != muonTrack->recHitsEnd(); hit++) {
	    if ( (*hit)->geographicalId().det() == 2 && (*hit)->geographicalId().subdetId() == 4) {
	      if ((*hit)->rawId() == ch->id().rawId() ) {
		GEMDetId gemid((*hit)->geographicalId());
		auto etaPart = GEMGeometry_->etaPartition(gemid);
		// cout << "found it "<< gemid
		//      << " lp " << (*hit)->localPosition()
		//      << " gp " << etaPart->toGlobal((*hit)->localPosition())
		//      << endl;
	      }
	    }
	  }
	}
      }

      for (auto hit = muonTrack->recHitsBegin(); hit != muonTrack->recHitsEnd(); hit++) {
	if ( (*hit)->geographicalId().det() == 2 && (*hit)->geographicalId().subdetId() == 4) {
	  GEMDetId gemid((*hit)->geographicalId());
	  h_trkEta[gemid.chamber()][gemid.layer()]->Fill(gemid.roll());
	  poorMuHits.push_back(static_cast<GEMRecHit*>(*hit));
	}
      }

      if (onDet) ++nGEMFiducialMuon;
      
      if (muonTrack->hitPattern().numberOfValidMuonGEMHits()) {
	++b_nMuonsWithGEMHit;
	++nGEMTrackWithMuon;
	// std::cout << "numberOfValidMuonGEMHits->size() " << muonTrack->hitPattern().numberOfValidMuonGEMHits()
	// 	  << " recHitsSize " << muonTrack->recHitsSize()
	// 	  << " pt " << muonTrack->pt()
	// 	  <<std::endl;
	for (auto hit = muonTrack->recHitsBegin(); hit != muonTrack->recHitsEnd(); hit++) {
	  if ( (*hit)->geographicalId().det() == 2 && (*hit)->geographicalId().subdetId() == 4) {
	    //if ((*hit)->rawId() == ch->id().rawId() ) {
	    GEMDetId gemid((*hit)->geographicalId());
	    auto etaPart = GEMGeometry_->etaPartition(gemid);

	    TrajectoryStateOnSurface tsos = propagator->propagate(ttTrack.outermostMeasurementState(),etaPart->surface());
	    if (!tsos.isValid()) continue;	    
	    // GlobalPoint tsosGP = tsos.globalPosition();

	    // LocalPoint && tsos_localpos = tsos.localPosition();
	    // LocalError && tsos_localerr = tsos.localError().positionError();
	    // LocalPoint && dethit_localpos = (*hit)->localPosition();     
	    // LocalError && dethit_localerr = (*hit)->localPositionError();
	    // auto res_x = (dethit_localpos.x() - tsos_localpos.x());
	    // auto res_y = (dethit_localpos.y() - tsos_localpos.y()); 
	    // auto pull_x = (dethit_localpos.x() - tsos_localpos.x()) / 
	    //   std::sqrt(dethit_localerr.xx() + tsos_localerr.xx());
	    // auto pull_y = (dethit_localpos.y() - tsos_localpos.y()) / 
	    //   std::sqrt(dethit_localerr.yy() + tsos_localerr.yy());
	    
	    // cout << "gem hit "<< gemid<< endl;
	    // cout << " gp " << etaPart->toGlobal((*hit)->localPosition())<< endl;
	    // cout << " tsosGP "<< tsosGP << endl;
	    // cout << " res_x " << res_x
	    // 	 << " res_y " << res_y
	    // 	 << " pull_x " << pull_x
	    // 	 << " pull_y " << pull_y
	    //   << endl;
	    }
	  
	}
	// auto res = muonTrack->residuals();
	// for (unsigned int i = 0; i < muonTrack->recHitsSize(); ++i) {
	//   cout << " res x "<< res.residualX(i)
	//        << " res y "<< res.residualY(i)
	//        << " pull x "<< res.pullX(i)
	//        << " pull y "<< res.pullY(i)
	//        <<endl;
	// }
      }
    }
  }

  // if (poorMuHits.size() > 0) {
  //   std::cout << "POOR MU ";
  //   for (auto & hit : poorMuHits) std::cout << hit << " ";
  //   std::cout << std::endl;
  // }
  
  // if (gemRecHits->size()) {
  //   std::cout << "gemRecHits->size() " << gemRecHits->size() <<std::endl;
  // }
  int totalStrips = 0;

  for (auto ch : GEMGeometry_->chambers()) {
    for(auto roll : ch->etaPartitions()) {
      GEMDetId rId = roll->id();

      //std::cout << "rId " << rId <<std::endl;
      auto recHitsRange = gemRecHits->get(rId); 
      auto gemRecHit = recHitsRange.first;
      for (auto hit = gemRecHit; hit != recHitsRange.second; ++hit) {

	h_hitEta[rId.chamber()][rId.layer()]->Fill(rId.roll());
	h_firstStrip[rId.chamber()][rId.layer()-1]->Fill(hit->firstClusterStrip(), rId.roll());
	h_clusterSize->Fill(hit->clusterSize());
	h_bxtotal->Fill(hit->BunchX());
	for (int nstrip = hit->firstClusterStrip(); nstrip < hit->firstClusterStrip()+hit->clusterSize(); ++nstrip) {
	  totalStrips++;
	  h_allStrips[rId.chamber()][rId.layer()-1]->Fill(nstrip, rId.roll());
	}

	b_firstStrip = hit->firstClusterStrip();
	b_nStrips = hit->clusterSize();
	b_chamber = rId.chamber();
	b_layer = rId.layer();
	b_etaPartition = rId.roll();
	b_muonQuality = -1;
	for (auto muHit : poorMuHits) {
	  if (*hit == *muHit) {
	    b_muonQuality = 0;
	    // std::cout << "  POOR HIT" << std::endl;
	  }
	}

	for (auto muHit : goodMuHits) {
	  if (*hit == *muHit)
	    b_muonQuality = 1;
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
