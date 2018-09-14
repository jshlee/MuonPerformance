// cd /cms/ldap_home/iawatson/scratch/GEM/CMSSW_10_1_5/src/ && eval `scramv1 runtime -sh` && eval `scramv1 runtime -sh` && scram b -j 10
// cd ../../.. && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval `scramv1 runtime -sh` && eval `scramv1 runtime -sh` && scram b -j 10
// system include files
#include <memory>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>

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
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/MuonData/interface/MuonDigiCollection.h"
#include "DataFormats/GEMDigi/interface/GEMAMCStatusDigi.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Run.h"

#include "EventFilter/Utilities/interface/json.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TTree.h"

using namespace std;
using namespace edm;

class SliceTestEfficiencyAnalysis : public edm::EDAnalyzer {
public:
  explicit SliceTestEfficiencyAnalysis(const edm::ParameterSet&);
  ~SliceTestEfficiencyAnalysis();

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginJob() override;
  virtual void endJob() override;

  virtual void beginRun(Run const&, EventSetup const&) override;
  virtual void endRun(Run const&, EventSetup const&) override;

  const GEMEtaPartition* findEtaPartition(const GEMChamber*& chamber, GlobalPoint& tsosGP);

  // ----------member data ---------------------------
  edm::EDGetTokenT<GEMRecHitCollection> gemRecHits_;
  edm::EDGetTokenT<edm::View<reco::Muon> > muons_;
  edm::EDGetTokenT<reco::VertexCollection> vertexCollection_;
  //edm::EDGetTokenT<MuonDigiCollection<unsigned short,GEMAMCStatusDigi>> gemDigis_;
  edm::Service<TFileService> fs;

  double Latency_;

  MuonServiceProxy* theService_;
  edm::ESHandle<Propagator> propagator_;
  edm::ESHandle<TransientTrackBuilder> ttrackBuilder_;
  edm::ESHandle<MagneticField> bField_; 

  static const int MAXCHAMBERS = 36;
  static const int MAXLAYERS = 2;
  static const int MAXROLL = 8;

  int nGEMHitInMuontrack;

  TH1D* h_nGEMHitInMuontrack;
  TH2D* h_inMap;
  TH2D* h_hitMap;
  
  TH2D* h_hitLumiMap[MAXCHAMBERS][MAXLAYERS];
  TH2D* h_stripLumiMap[MAXCHAMBERS][MAXLAYERS];
  TH1D* h_inRoll[MAXCHAMBERS][MAXLAYERS];
  TH1D* h_inStrip[MAXCHAMBERS][MAXLAYERS];
  TH2D* h_inVfat[MAXCHAMBERS][MAXLAYERS];
  TH2D* h_inPos[MAXCHAMBERS][MAXLAYERS];

  TH1D* h_hitRoll[MAXCHAMBERS][MAXLAYERS];
  TH1D* h_hitStrip[MAXCHAMBERS][MAXLAYERS];
  TH2D* h_hitVfat[MAXCHAMBERS][MAXLAYERS];
  TH2D* h_hitPos[MAXCHAMBERS][MAXLAYERS];
  TH1D* h_hitNstrip[MAXCHAMBERS][MAXLAYERS][MAXROLL];
  TH1D* h_resX_etaPart[MAXCHAMBERS][MAXLAYERS][MAXROLL];
  TH1D* h_resY_etaPart[MAXCHAMBERS][MAXLAYERS][MAXROLL];
  TH1D* h_resPhi_etaPart[MAXCHAMBERS][MAXLAYERS][MAXROLL];
  TH1D* h_resX_xPart[MAXCHAMBERS][MAXLAYERS][3];
  TH1D* h_resY_xPart[MAXCHAMBERS][MAXLAYERS][3];
  TH1D* h_resPhi_xPart[MAXCHAMBERS][MAXLAYERS][3];

  TH1D* h_resX[MAXCHAMBERS][MAXLAYERS];
  TH1D* h_resY[MAXCHAMBERS][MAXLAYERS];
  TH1D* h_resPhi[MAXCHAMBERS][MAXLAYERS];
  TH2D* h_resXvsNstrip[MAXCHAMBERS][MAXLAYERS];
  TH1D* h_pullX[MAXCHAMBERS][MAXLAYERS];
  TH1D* h_pullY[MAXCHAMBERS][MAXLAYERS];
  TH2D* h_inPos_matched[MAXCHAMBERS][MAXLAYERS];
  TH2D* h_inPhiVsHitPhi[MAXCHAMBERS][MAXLAYERS];
  TH2D* h_inXVsHitX[MAXCHAMBERS][MAXLAYERS];
  TH2D* h_inStripVsHitStrip[MAXCHAMBERS][MAXLAYERS];
  TH2D* h_inPhiVsHitStrip[MAXCHAMBERS][MAXLAYERS];

  int b_run, b_lumi;
  int nEvents;
  int b_nMuons, b_nMuonsInGEMRegion, b_nGEMHits;
  int b_latency;
  
  TTree *t_event;
};

SliceTestEfficiencyAnalysis::SliceTestEfficiencyAnalysis(const edm::ParameterSet& iConfig) :
  nEvents(0)
{ 
  gemRecHits_ = consumes<GEMRecHitCollection>(iConfig.getParameter<edm::InputTag>("gemRecHits"));
  muons_ = consumes<View<reco::Muon> >(iConfig.getParameter<InputTag>("muons"));
  edm::ParameterSet serviceParameters = iConfig.getParameter<edm::ParameterSet>("ServiceParameters");
  //gemDigis_ = consumes<MuonDigiCollection<unsigned short,GEMAMCStatusDigi>>(iConfig.getParameter<edm::InputTag>("gemDigis"));
  theService_ = new MuonServiceProxy(serviceParameters);
  Latency_ = iConfig.getParameter<double>("latency");

  t_event = fs->make<TTree>("Event", "Event");
  t_event->Branch("nMuons", &b_nMuons, "nMuons/I");
  t_event->Branch("nMuonsInGEMRegion", &b_nMuonsInGEMRegion, "nMuonsInGEMRegion/I");
  t_event->Branch("nGEMHits", &b_nGEMHits, "nGEMHits/I");
  t_event->Branch("run", &b_run, "run/I");
  t_event->Branch("lumi", &b_lumi, "lumi/I");
  t_event->Branch("latency", &b_latency, "latency/I");

  h_nGEMHitInMuontrack = fs->make<TH1D>(Form("nGEMHitInMuontrack"),"nGEMHitInMuontrack",4,0,4);
  h_inMap = fs->make<TH2D>(Form("inMap"),"inMap",18,26.5,31,18,0,9);
  h_hitMap = fs->make<TH2D>(Form("hitMap"),"hitMap",18,26.5,31,18,0,9);
  for (int ichamber=27; ichamber<=30;++ichamber) {
    for (int ilayer=0; ilayer<MAXLAYERS;++ilayer) {
      h_hitLumiMap[ichamber][ilayer] = fs->make<TH2D>(Form("hitLumiMap ch %i lay %i",ichamber, ilayer+1),"hitLumiMap",700,0,700,32,1,9);
      h_stripLumiMap[ichamber][ilayer] = fs->make<TH2D>(Form("stripLumiMap ch %i lay %i",ichamber, ilayer+1),"stripLumiMap",700,0,700,400,0,400);
      h_inRoll[ichamber][ilayer] = fs->make<TH1D>(Form("inRoll ch %i lay %i",ichamber, ilayer+1),"inRoll",8,0.5,8.5);
      h_inStrip[ichamber][ilayer] = fs->make<TH1D>(Form("inStrip ch %i lay %i",ichamber, ilayer+1),"inStrip",384,0,384);
      h_inVfat[ichamber][ilayer] = fs->make<TH2D>(Form("inVfat ch %i lay %i",ichamber, ilayer+1),"inVfat",3,0.5,3.5,8,0.5,8.5);
      h_inPos[ichamber][ilayer] = fs->make<TH2D>(Form("inPos ch %i lay %i",ichamber, ilayer+1),"inPos",100,-70,120,100,-260,-110);

      h_hitRoll[ichamber][ilayer] = fs->make<TH1D>(Form("hitRoll ch %i lay %i",ichamber, ilayer+1),"hitRoll",8,0.5,8.5);
      h_hitStrip[ichamber][ilayer] = fs->make<TH1D>(Form("hitStrip ch %i lay %i",ichamber, ilayer+1),"hitStrip",384,0,384);
      h_hitVfat[ichamber][ilayer] = fs->make<TH2D>(Form("hitVfat ch %i lay %i",ichamber, ilayer+1),"hitVfat",3,0.5,3.5,8,0.5,8.5);
      h_hitPos[ichamber][ilayer] = fs->make<TH2D>(Form("hitPos ch %i lay %i",ichamber, ilayer+1),"hitPos",100,-70,120,100,-260,-110);
      for (int ieta=0; ieta<MAXROLL;++ieta) {
        h_hitNstrip[ichamber][ilayer][ieta] = fs->make<TH1D>(Form("hitNstrip ch %i lay %i ieta %i",ichamber, ilayer+1, ieta),"hitNstrip",30,0,30);
        h_resX_etaPart[ichamber][ilayer][ieta] = fs->make<TH1D>(Form("resX_etaPart ch %i lay %i ieta %i",ichamber, ilayer+1, ieta),"resX",500,-3,3);
        h_resY_etaPart[ichamber][ilayer][ieta] = fs->make<TH1D>(Form("resY_etaPart ch %i lay %i ieta %i",ichamber, ilayer+1, ieta),"resY",500,-15,15);
        h_resPhi_etaPart[ichamber][ilayer][ieta] = fs->make<TH1D>(Form("resPhi_etaPart ch %i lay %i ieta %i",ichamber, ilayer+1, ieta),"resPhi",500,-0.03,0.03);
      }
      for (int xrange=0; xrange<3;++xrange) {
        h_resX_xPart[ichamber][ilayer][xrange] = fs->make<TH1D>(Form("resX_xPart ch %i lay %i xrange %i",ichamber, ilayer+1, xrange),"resX",500,-3,3);
        h_resY_xPart[ichamber][ilayer][xrange] = fs->make<TH1D>(Form("resY_xPart ch %i lay %i xrange %i",ichamber, ilayer+1, xrange),"resY",500,-15,15);
        h_resPhi_xPart[ichamber][ilayer][xrange] = fs->make<TH1D>(Form("resPhi_xPart ch %i lay %i xrange %i",ichamber, ilayer+1, xrange),"resPhi",500,-0.03,0.03);
      }

      h_resX[ichamber][ilayer] = fs->make<TH1D>(Form("resX ch %i lay %i",ichamber, ilayer+1),"resX",500,-3,3);
      h_resY[ichamber][ilayer] = fs->make<TH1D>(Form("resY ch %i lay %i",ichamber, ilayer+1),"resY",500,-15,15);
      h_resPhi[ichamber][ilayer] = fs->make<TH1D>(Form("resPhi ch %i lay %i",ichamber, ilayer+1),"resPhi",500,-0.03,0.03);
      h_resXvsNstrip[ichamber][ilayer] = fs->make<TH2D>(Form("resXvsNstrip ch %i lay %i",ichamber, ilayer+1),"resXvsNstrip",250,-3,3,15,0,15);
      h_pullX[ichamber][ilayer] = fs->make<TH1D>(Form("pullX ch %i lay %i",ichamber, ilayer+1),"pullX",500,-3,3);
      h_pullY[ichamber][ilayer] = fs->make<TH1D>(Form("pullY ch %i lay %i",ichamber, ilayer+1),"pullY",500,-3,3);
      h_inPos_matched[ichamber][ilayer] = fs->make<TH2D>(Form("inPos_matched ch %i lay %i",ichamber, ilayer+1),"inPos",100,-70,120,100,-260,-110);
      h_inPhiVsHitPhi[ichamber][ilayer] = fs->make<TH2D>(Form("inPhiVsHitPhi ch %i lay %i",ichamber, ilayer+1),"inPos",100,-1.9,-1.1,100,-1.9,-1.1);
      h_inXVsHitX[ichamber][ilayer] = fs->make<TH2D>(Form("inXVsHitX ch %i lay %i",ichamber, ilayer+1),"inPos",100,-70,120,100,-70,120);
      h_inStripVsHitStrip[ichamber][ilayer] = fs->make<TH2D>(Form("inStripVsHitStrip ch %i lay %i",ichamber, ilayer+1),"inPos",400,0,400,400,0,400);
      h_inPhiVsHitStrip[ichamber][ilayer] = fs->make<TH2D>(Form("inPhiVsHitStrip ch %i lay %i",ichamber, ilayer+1),"inPos",100,-1.9,-1.1,400,0,400);
    }
  }
}

SliceTestEfficiencyAnalysis::~SliceTestEfficiencyAnalysis(){}

void
SliceTestEfficiencyAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  nEvents++;

  b_nMuons = 0;
  b_nMuonsInGEMRegion = 0;
  b_nGEMHits = 0;
  
  b_run = iEvent.run();
  b_lumi = iEvent.luminosityBlock();

  edm::ESHandle<GEMGeometry> hGeom;
  iSetup.get<MuonGeometryRecord>().get(hGeom);
  const GEMGeometry* GEMGeometry_ = &*hGeom;

  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",ttrackBuilder_);
  theService_->update(iSetup);
  auto propagator = theService_->propagator("SteppingHelixPropagatorAny");
  
  edm::Handle<GEMRecHitCollection> gemRecHits;  
  iEvent.getByToken(gemRecHits_, gemRecHits);
 
  Handle<View<reco::Muon> > muons;
  iEvent.getByToken(muons_, muons);

  /*
  edm::Handle<MuonDigiCollection<unsigned short,GEMAMCStatusDigi>> gemDigis;
  iEvent.getByToken(gemDigis_, gemDigis);

  b_latency = -1;
  for (auto g : *gemDigis) {
    for (auto a = g.second.first; a != g.second.second; ++a) {
      b_latency = a->Param1();
    }
  }
  //if (b_latency != Latency_) return;
  */

  if (b_run != 319347) return;
  if ( (b_lumi <50)
     || (b_lumi>210 && b_lumi<240)
     || (b_lumi>640) ) return;

  for (auto ch : GEMGeometry_->chambers()) {
    for(auto roll : ch->etaPartitions()) {
      auto rId = roll->id();
      auto recHitsRange = gemRecHits->get(rId);
      b_nGEMHits += recHitsRange.second - recHitsRange.first;
      for (auto hit = recHitsRange.first; hit != recHitsRange.second; ++hit) {
        int vfat = hit->firstClusterStrip()/128 +1;
        h_hitLumiMap[rId.chamber()][rId.layer()-1]->Fill(b_lumi,rId.roll()+vfat/4.);
        h_stripLumiMap[rId.chamber()][rId.layer()-1]->Fill(b_lumi,hit->firstClusterStrip());
      }
    }
  }
  //if (b_nGEMHits == 0) continue;

  for (size_t i = 0; i < muons->size(); ++i) {

    edm::RefToBase<reco::Muon> muRef = muons->refAt(i);
    const reco::Muon* mu = muRef.get();

    // muon id
    int muonId = 0;
    if (mu->passed(reco::Muon::Selector::CutBasedIdTight)) muonId = 2;
    else if (mu->passed(reco::Muon::Selector::CutBasedIdLoose)) muonId = 1;

    // tight and pt > 20 muon only
    if (muonId != 2) continue;
    if (mu->pt() < 20) continue;
    
    const reco::Track* muonTrack = 0;  
    if ( mu->globalTrack().isNonnull() ) muonTrack = mu->globalTrack().get();
    else if ( mu->outerTrack().isNonnull()  ) muonTrack = mu->outerTrack().get();
    if (!muonTrack) continue;
    b_nMuons++;

    nGEMHitInMuontrack = 0;
    for (auto hit = muonTrack->recHitsBegin(); hit != muonTrack->recHitsEnd(); hit++) {
      if ((*hit)->geographicalId().det() == 2 && (*hit)->geographicalId().subdetId() == 4) {
        nGEMHitInMuontrack++;
      }
    }
    h_nGEMHitInMuontrack->Fill(nGEMHitInMuontrack);

    reco::TransientTrack ttTrack = ttrackBuilder_->build(muonTrack);
    for (auto chamber : GEMGeometry_->chambers()) {
	  if (chamber->id().chamber() == 1) continue; // ignore chammber 1
	  if (mu->eta() * chamber->id().region() < 0 ) continue;

	  TrajectoryStateOnSurface tsos = propagator->propagate(ttTrack.outermostMeasurementState(),
	  					            chamber->surface());
	  if (!tsos.isValid()) continue;

	  GlobalPoint tsosGP = tsos.globalPosition();
      auto etaPart =  findEtaPartition(chamber, tsosGP);
      if (!etaPart) continue;

      auto gemid = etaPart->id();
      auto locPos = etaPart->toLocal(tsosGP);
      auto strip = (int) etaPart->strip(locPos);
      auto vfat = ((int) strip/128)+1;

	  h_inRoll[gemid.chamber()][gemid.layer()-1]->Fill(gemid.roll());
	  h_inStrip[gemid.chamber()][gemid.layer()-1]->Fill(strip);
	  h_inVfat[gemid.chamber()][gemid.layer()-1]->Fill(vfat, gemid.roll());
	  h_inPos[gemid.chamber()][gemid.layer()-1]->Fill(tsosGP.x(), tsosGP.y());
	  h_inMap->Fill(gemid.chamber()+gemid.layer()/2., gemid.roll());

      //Find hit
      float resX = 999;
      GEMRecHit closestHit;
	  auto recHitsRange = gemRecHits->get(gemid);
	  for (auto hit = recHitsRange.first; hit != recHitsRange.second; ++hit) {
	    LocalPoint hitLocPos = hit->localPosition();
        if ( fabs(hitLocPos.x() - locPos.x()) < fabs(resX) ) {
          resX = hitLocPos.x() - locPos.x();
          closestHit = (*hit);
        }
	  }
	  if (resX == 999) continue;
	  auto hitLocPos = closestHit.localPosition();
      auto hitGlobPos = etaPart->toGlobal(hitLocPos);
      auto hitStrip = closestHit.firstClusterStrip();
      auto resY = hitLocPos.y() - locPos.y();
      auto resPhi = hitGlobPos.phi() - tsosGP.phi();
	  h_inPhiVsHitPhi[gemid.chamber()][gemid.layer()-1]->Fill(tsosGP.phi(), hitGlobPos.phi());
	  h_inXVsHitX[gemid.chamber()][gemid.layer()-1]->Fill(tsosGP.x(), hitGlobPos.x());
	  h_inStripVsHitStrip[gemid.chamber()][gemid.layer()-1]->Fill(strip, hitStrip);
	  h_inPhiVsHitStrip[gemid.chamber()][gemid.layer()-1]->Fill(tsosGP.phi(), hitStrip);

	  if (resX > 5.0) continue;
      LocalError && locErr = tsos.localError().positionError();
      LocalError && hitLocErr = closestHit.localPositionError();
      auto pullX = resX / std::sqrt(hitLocErr.xx() + locErr.xx());
      auto pullY = resY / std::sqrt(hitLocErr.yy() + locErr.yy());

      //Filling histograms
      h_hitRoll[gemid.chamber()][gemid.layer()-1]->Fill(gemid.roll());
	  h_hitStrip[gemid.chamber()][gemid.layer()-1]->Fill(hitStrip);
	  h_hitVfat[gemid.chamber()][gemid.layer()-1]->Fill(((int)hitStrip/128)+1, gemid.roll());
	  h_hitPos[gemid.chamber()][gemid.layer()-1]->Fill(hitGlobPos.x(), hitGlobPos.y());
	  h_hitNstrip[gemid.chamber()][gemid.layer()-1][gemid.roll()-1]->Fill(closestHit.clusterSize());
	  h_hitMap->Fill(gemid.chamber()+gemid.layer()/2., gemid.roll());

      int x = 1;
      if (abs(locPos.x())<10) x = 0;
      else if (abs(locPos.x())>20) x = 2;
	  h_resX_etaPart[gemid.chamber()][gemid.layer()-1][gemid.roll()-1]->Fill(resX);
	  h_resY_etaPart[gemid.chamber()][gemid.layer()-1][gemid.roll()-1]->Fill(resY);
	  h_resPhi_etaPart[gemid.chamber()][gemid.layer()-1][gemid.roll()-1]->Fill(resPhi);
	  h_resX_xPart[gemid.chamber()][gemid.layer()-1][x]->Fill(resX);
	  h_resY_xPart[gemid.chamber()][gemid.layer()-1][x]->Fill(resY);
	  h_resPhi_xPart[gemid.chamber()][gemid.layer()-1][x]->Fill(resPhi);
	  h_resX[gemid.chamber()][gemid.layer()-1]->Fill(resX);
	  h_resY[gemid.chamber()][gemid.layer()-1]->Fill(resY);
	  h_resPhi[gemid.chamber()][gemid.layer()-1]->Fill(resPhi);
	  h_resXvsNstrip[gemid.chamber()][gemid.layer()-1]->Fill(resX,closestHit.clusterSize());
	  h_pullX[gemid.chamber()][gemid.layer()-1]->Fill(pullX);
	  h_pullY[gemid.chamber()][gemid.layer()-1]->Fill(pullY);
	  h_inPos_matched[gemid.chamber()][gemid.layer()-1]->Fill(tsosGP.x(), tsosGP.y());

    }
  }
  t_event->Fill();
}

const GEMEtaPartition* SliceTestEfficiencyAnalysis::findEtaPartition(const GEMChamber*& chamber, GlobalPoint& tsosGP){
  for (auto etaPart : chamber->etaPartitions()) {
    const LocalPoint locPos = etaPart->toLocal(tsosGP);
    const LocalPoint locPos2D(locPos.x(), locPos.y(), 0);
    const BoundPlane& bps(etaPart->surface());
    if (!bps.bounds().inside(locPos2D)) continue;
    return etaPart;
  }
  return nullptr;
}

void SliceTestEfficiencyAnalysis::beginJob(){}
void SliceTestEfficiencyAnalysis::endJob(){}

void SliceTestEfficiencyAnalysis::beginRun(Run const& run, EventSetup const&){}
void SliceTestEfficiencyAnalysis::endRun(Run const&, EventSetup const&){}

//define this as a plug-in
DEFINE_FWK_MODULE(SliceTestEfficiencyAnalysis);
