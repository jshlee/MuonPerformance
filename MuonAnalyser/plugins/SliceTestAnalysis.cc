// system include files
#include <memory>
#include <cmath>
#include <iostream>
#include <sstream>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

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

#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"

using namespace std;
using namespace edm;

class SliceTestAnalysis : public edm::EDAnalyzer {
public:
  explicit SliceTestAnalysis(const edm::ParameterSet&);
  ~SliceTestAnalysis(){};

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginJob() ;
  virtual void endJob() ;

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
  TH1D* h_clusterSize, *h_totalStrips, *h_bxtotal;
};

SliceTestAnalysis::SliceTestAnalysis(const edm::ParameterSet& iConfig)
{ 
  gemRecHits_ = consumes<GEMRecHitCollection>(iConfig.getParameter<edm::InputTag>("gemRecHits"));
  muons_ = consumes<View<reco::Muon> >(iConfig.getParameter<InputTag>("muons"));
  vertexCollection_ = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollection"));
  edm::ParameterSet serviceParameters = iConfig.getParameter<edm::ParameterSet>("ServiceParameters");
  theService_ = new MuonServiceProxy(serviceParameters);

  h_clusterSize=fs->make<TH1D>(Form("clusterSize"),"clusterSize",100,0,100);
  h_totalStrips=fs->make<TH1D>(Form("totalStrips"),"totalStrips",200,0,200);
  h_bxtotal=fs->make<TH1D>(Form("bx"),"bx",31,-15,15);
  
  //for (int ichamber=0; ichamber<36;++ichamber) {
  for (int ichamber=27; ichamber<=30;++ichamber) {
    for (int ilayer=0; ilayer<2;++ilayer) {
      h_firstStrip[ichamber][ilayer] = fs->make<TH2D>(Form("firstStrip ch %i lay %i",ichamber, ilayer),"firstStrip",384,1,385,8,0.5,8.5);
      h_firstStrip[ichamber][ilayer]->GetXaxis()->SetTitle("strip");
      h_firstStrip[ichamber][ilayer]->GetYaxis()->SetTitle("iEta");
      
      h_allStrips[ichamber][ilayer] = fs->make<TH2D>(Form("allStrips ch %i lay %i",ichamber, ilayer),"allStrips",384,1,385,8,0.5,8.5);
      h_allStrips[ichamber][ilayer]->GetXaxis()->SetTitle("strip");
      h_allStrips[ichamber][ilayer]->GetYaxis()->SetTitle("iEta");
    }
  }
}

void
SliceTestAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
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

  for (size_t i = 0; i < muons->size(); ++i) {    
    edm::RefToBase<reco::Muon> muRef = muons->refAt(i);
    const reco::Muon* mu = muRef.get();
    if (mu->isGEMMuon()) {
      std::cout << "isGEMMuon " <<std::endl;
    }
    const reco::Track* muonTrack = 0;  
    if ( mu->globalTrack().isNonnull() ) muonTrack = mu->globalTrack().get();
    else if ( mu->outerTrack().isNonnull()  ) muonTrack = mu->outerTrack().get();
    if (muonTrack) {

      std::set<double> detLists;

      reco::TransientTrack ttTrack = ttrackBuilder_->build(muonTrack);
      for (auto ch : GEMGeometry_->etaPartitions()) {
	//if ( !detLists.insert( ch->surface().position().z() ).second ) continue;
	
	TrajectoryStateOnSurface tsos = propagator->propagate(ttTrack.outermostMeasurementState(),ch->surface());
	if (!tsos.isValid()) continue;

	GlobalPoint tsosGP = tsos.globalPosition();
	const LocalPoint pos = ch->toLocal(tsosGP);
	const LocalPoint pos2D(pos.x(), pos.y(), 0);
	const BoundPlane& bps(ch->surface());
	//cout << "tsos gp   "<< tsosGP << ch->id() <<endl;

	if (bps.bounds().inside(pos2D)) {
	  cout << " in chamber "<< ch->id() << " pos = "<<pos<< " R = "<<pos.mag() <<" inside "
	       <<  bps.bounds().inside(pos2D) <<endl;
	  
	  for (auto hit = muonTrack->recHitsBegin(); hit != muonTrack->recHitsEnd(); hit++) {
	    if ( (*hit)->geographicalId().det() == 2 && (*hit)->geographicalId().subdetId() == 4) {
	      if ((*hit)->rawId() == ch->id().rawId() ) {
		GEMDetId gemid((*hit)->geographicalId());
		auto etaPart = GEMGeometry_->etaPartition(gemid);
		cout << "found it "<< gemid
		     << " lp " << (*hit)->localPosition()
		     << " gp " << etaPart->toGlobal((*hit)->localPosition())
		     << endl;
	      }
	    }
	  }
	}
      }
      
      if (muonTrack->hitPattern().numberOfValidMuonGEMHits()) {
	std::cout << "numberOfValidMuonGEMHits->size() " << muonTrack->hitPattern().numberOfValidMuonGEMHits()
		  << " recHitsSize " << muonTrack->recHitsSize()
		  << " pt " << muonTrack->pt()
		  <<std::endl;
	for (auto hit = muonTrack->recHitsBegin(); hit != muonTrack->recHitsEnd(); hit++) {
	  if ( (*hit)->geographicalId().det() == 2 && (*hit)->geographicalId().subdetId() == 4) {
	    //if ((*hit)->rawId() == ch->id().rawId() ) {
	    GEMDetId gemid((*hit)->geographicalId());
	    auto etaPart = GEMGeometry_->etaPartition(gemid);

	    TrajectoryStateOnSurface tsos = propagator->propagate(ttTrack.outermostMeasurementState(),etaPart->surface());
	    if (!tsos.isValid()) continue;	    
	    GlobalPoint tsosGP = tsos.globalPosition();

	    LocalPoint && tsos_localpos = tsos.localPosition();
	    LocalError && tsos_localerr = tsos.localError().positionError();
	    LocalPoint && dethit_localpos = (*hit)->localPosition();     
	    LocalError && dethit_localerr = (*hit)->localPositionError();
	    auto res_x = (dethit_localpos.x() - tsos_localpos.x());
	    auto res_y = (dethit_localpos.y() - tsos_localpos.y()); 
	    auto pull_x = (dethit_localpos.x() - tsos_localpos.x()) / 
	      std::sqrt(dethit_localerr.xx() + tsos_localerr.xx());
	    auto pull_y = (dethit_localpos.y() - tsos_localpos.y()) / 
	      std::sqrt(dethit_localerr.yy() + tsos_localerr.yy());
	    
	    cout << "gem hit "<< gemid<< endl;
	    cout << " gp " << etaPart->toGlobal((*hit)->localPosition())<< endl;
	    cout << " tsosGP "<< tsosGP << endl;
	    cout << " res_x " << res_x
		 << " res_y " << res_y
		 << " pull_x " << pull_x
		 << " pull_y " << pull_y
	      << endl;
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
      for ( auto hit = gemRecHit; hit != recHitsRange.second; ++hit ) {

	h_firstStrip[rId.chamber()][rId.layer()-1]->Fill(hit->firstClusterStrip(), rId.roll());
	h_clusterSize->Fill(hit->clusterSize());
	h_bxtotal->Fill(hit->BunchX());
	for (int nstrip = hit->firstClusterStrip(); nstrip < hit->firstClusterStrip()+hit->clusterSize(); ++nstrip) {
	  totalStrips++;
	  h_allStrips[rId.chamber()][rId.layer()-1]->Fill(nstrip, rId.roll());
	}
      }
    }
  }
  h_totalStrips->Fill(totalStrips);

  
}

void SliceTestAnalysis::beginJob(){}
void SliceTestAnalysis::endJob(){}

//define this as a plug-in
DEFINE_FWK_MODULE(SliceTestAnalysis);
