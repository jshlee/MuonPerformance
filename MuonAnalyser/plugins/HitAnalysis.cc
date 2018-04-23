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

#include "DataFormats/GEMDigi/interface/GEMDigiCollection.h"
#include "DataFormats/GEMDigi/interface/GEMPadDigiCollection.h"
#include "DataFormats/GEMDigi/interface/GEMCoPadDigiCollection.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartition.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"

using namespace std;
using namespace edm;

class HitAnalysis : public edm::EDAnalyzer {
public:
  explicit HitAnalysis(const edm::ParameterSet&);
  ~HitAnalysis();

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  edm::EDGetTokenT<GEMDigiCollection> gemDigiInput_;
  edm::EDGetTokenT<GEMPadDigiCollection> gemPadDigiInput_;
  edm::EDGetTokenT<GEMCoPadDigiCollection> gemCoPadDigiInput_;

  edm::Service<TFileService> fs;

  TH1D* h_nPadPerGEB[12];
  TH2D* h_nPadPerGEBPerBX[12];
  TH1D* h_nPadGEB[12];
  TH2D* h_nPadGEBPerBX[12];
  TH1D* h_nPadRangeGEB[12];
  TH2D* h_nPadRangeGEBPerBX[12];  
  TH1D* h_nCoPadPerGEB[12];
  TH2D* h_nCoPadPerGEBPerBX[12];
  TH2D* h_maxPadvsCoPad[12];
  
};
HitAnalysis::HitAnalysis(const edm::ParameterSet& iConfig)
{
  gemDigiInput_ = consumes<GEMDigiCollection>(iConfig.getParameter<edm::InputTag>("gemDigiInput"));
  gemPadDigiInput_ = consumes<GEMPadDigiCollection>(iConfig.getParameter<edm::InputTag>("gemPadDigiInput"));
  gemCoPadDigiInput_ = consumes<GEMCoPadDigiCollection>(iConfig.getParameter<edm::InputTag>("gemCoPadDigiInput"));

  for (int i=0; i<12;++i){
    h_nPadPerGEB[i] =fs->make<TH1D>(Form("nPadPerGEB M%i",i+1),"nPadPerGEB",100,0,100);
    h_nPadPerGEB[i]->GetXaxis()->SetTitle("no. Pads Per GEB");
    
    h_nPadGEB[i] =fs->make<TH1D>(Form("nPadGEB M%i",i+1),"nPadGEB",100,0,100);
    h_nPadGEB[i]->GetXaxis()->SetTitle("no. Pad Ranges per GEB");
    
    h_nPadRangeGEB[i] =fs->make<TH1D>(Form("nPadRangeGEB M%i",i+1),"nPadRangeGEB",100,0,100);
    h_nPadRangeGEB[i]->GetXaxis()->SetTitle("size of Pad Range");
    
    h_nCoPadPerGEB[i] =fs->make<TH1D>(Form("nCoPadPerGEB M%i",i+1),"nCoPadPerGEB",100,0,100);
    h_nCoPadPerGEB[i]->GetXaxis()->SetTitle("no. coPads Per GEB");
  
    h_nPadPerGEBPerBX[i] =fs->make<TH2D>(Form("nPadPerGEBPerBX M%i",i+1),"nPadPerGEBPerBX",50,-25,25,100,0,100);
    h_nPadPerGEBPerBX[i]->GetXaxis()->SetTitle("BX");
    h_nPadPerGEBPerBX[i]->GetYaxis()->SetTitle("no. Pads Per GEB");

    h_nPadGEBPerBX[i] =fs->make<TH2D>(Form("nPadGEBPerBX M%i",i+1),"nPadGEBPerBX",50,-25,25,100,0,100);
    h_nPadGEBPerBX[i]->GetXaxis()->SetTitle("BX");
    h_nPadGEBPerBX[i]->GetYaxis()->SetTitle("no. Pad Ranges per GEB");
    
    h_nPadRangeGEBPerBX[i] =fs->make<TH2D>(Form("nPadRangeGEBPerBX M%i",i+1),"nPadRangeGEBPerBX",50,-25,25,100,0,100);
    h_nPadRangeGEBPerBX[i]->GetXaxis()->SetTitle("BX");
    h_nPadRangeGEBPerBX[i]->GetYaxis()->SetTitle("size of Pad Range");

    h_nCoPadPerGEBPerBX[i] =fs->make<TH2D>(Form("nCoPadPerGEBPerBX M%i",i+1),"nCoPadPerGEBPerBX",50,-25,25,100,0,100);
    h_nCoPadPerGEBPerBX[i]->GetXaxis()->SetTitle("BX");
    h_nCoPadPerGEBPerBX[i]->GetYaxis()->SetTitle("no. coPads Per GEB");
    
    h_maxPadvsCoPad[i] =fs->make<TH2D>(Form("maxPadvsCoPad M%i",i+1),"maxPadvsCoPad",100,0,100,100,0,100);    
    h_maxPadvsCoPad[i]->GetXaxis()->SetTitle("max no. Pads");
    h_maxPadvsCoPad[i]->GetYaxis()->SetTitle("no. CoPads");
  }
}
HitAnalysis::~HitAnalysis(){}
void
HitAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<GEMDigiCollection> gem_digis;  
  edm::Handle<GEMPadDigiCollection> gempad_digis;
  edm::Handle<GEMCoPadDigiCollection> gemcopad_digis;
  iEvent.getByToken(gemDigiInput_, gem_digis);
  iEvent.getByToken(gemPadDigiInput_, gempad_digis);
  iEvent.getByToken(gemCoPadDigiInput_, gemcopad_digis);

  // GE21 18*2 chambers, 8 GEBs, bx or pads
  int nPadPerGEB[37][12][40] = {};
  int nPadGEB[37][12][384] = {};

  for(GEMPadDigiCollection::DigiRangeIterator cItr = gempad_digis->begin(); cItr != gempad_digis->end(); ++cItr){
    GEMDetId id = (*cItr).first; 
    auto range((*cItr).second);

    for (GEMPadDigiCollection::const_iterator digiIt = range.first; digiIt != range.second; digiIt++) {
      if (id.station() == 2){
	int m = (id.layer() - 1)*4 + (8-id.roll())/2 +1;
	// m1 = roll 7+8, m2 = roll 5+6, m3 = roll 3+4, m4 = roll 1+2
	// m5 = roll 7+8, m6 = roll 5+6, m7 = roll 3+4, m8 = roll 1+2
	// c1 = m1+m2
	// m9 = c1, m10 = c2 etc
	int c = (m-1)/2 +1;
	
	cout << "ge21 m "<< m
	     << " c " << c
	     << " id.layer() " << id.layer() 
	     << " id.roll() " << id.roll() 
	     << " pad digiIt->pad() " << digiIt->pad()
	     << " pad digiIt->bx() " << digiIt->bx()
	     << endl;	
	nPadPerGEB[18+id.region()*id.chamber()][m-1][15+digiIt->bx()]++;
	if (digiIt->bx() == 0) nPadGEB[18+id.region()*id.chamber()][m-1][digiIt->pad()]++;
	nPadPerGEB[18+id.region()*id.chamber()][8+c-1][15+digiIt->bx()]++;
	if (digiIt->bx() == 0) nPadGEB[18+id.region()*id.chamber()][8+c-1][digiIt->pad()]++;
	
      }
    }
  }
  for (int i = 0; i < 37; i++){
    for (int h = 0; h < 12; h++){
      
      if ( i == 18 ) continue;
      h_nPadPerGEB[h]->Fill(nPadPerGEB[i][h][15]);

      for (int j = -15; j < 25; j++){
	h_nPadPerGEBPerBX[h]->Fill(j,nPadPerGEB[i][h][15+j]);
      }
      
      int nContPadRanges = 0;
      int maxrange = 0;
      for (int k = 0; k < 384; k++){
	
	if (nPadGEB[i][h][k] && maxrange < 8 && k < 383){
	  maxrange++;
	}
	else if (maxrange){
	  h_nPadRangeGEB[h]->Fill(maxrange);
	  maxrange = 0;
	  nContPadRanges++;
	}	
      }
     
      h_nPadGEB[h]->Fill(nContPadRanges);

    }
  }

  // GE21 18*2 chambers, 8 GEBs, bx or pads
  int nCoPadPerGEB[37][12][40] = {};
  
  for(GEMCoPadDigiCollection::DigiRangeIterator cItr = gemcopad_digis->begin(); cItr != gemcopad_digis->end(); ++cItr){
    GEMDetId id = (*cItr).first; 
    auto range((*cItr).second);
    for (GEMCoPadDigiCollection::const_iterator digiIt = range.first; digiIt != range.second; digiIt++) {

      if (id.station() == 2){
	int m = (id.layer() - 1) + id.roll()/2;
	nCoPadPerGEB[18+id.region()*id.chamber()][m][15+digiIt->first().bx()]++;
      }
    }
  }

  for (int i = 0; i < 37; i++){
    for (int h = 0; h < 4; h++){
      
      if ( i == 18 ) continue;
      h_nCoPadPerGEB[h]->Fill(nCoPadPerGEB[i][h][15]);

      for (int j = -15; j < 25; j++){
	h_nCoPadPerGEBPerBX[h]->Fill(j,nCoPadPerGEB[i][h][15+j]);
      }
      
      if (nPadPerGEB[i][h][15] || nPadPerGEB[i][h+4][15]){
	h_maxPadvsCoPad[h]->Fill(std::max(nPadPerGEB[i][h][15], nPadPerGEB[i][h+4][15]) , nCoPadPerGEB[i][h][15] );
      }
    }
  }
  
}

void HitAnalysis::beginJob(){}
void HitAnalysis::endJob(){}

//define this as a plug-in
DEFINE_FWK_MODULE(HitAnalysis);
