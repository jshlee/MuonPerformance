#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"

class GEMSkim : public edm::EDFilter {
public:
  GEMSkim(const edm::ParameterSet&);
  ~GEMSkim() override;

private:
  void beginJob() override;
  bool filter( edm::Event &, edm::EventSetup const& ) override;
  void endJob() override;
  
  edm::EDGetTokenT<GEMRecHitCollection> gemRecHits_;
};

GEMSkim::GEMSkim(const edm::ParameterSet& iConfig)
{ 
  gemRecHits_ = consumes<GEMRecHitCollection>(iConfig.getParameter<edm::InputTag>("gemRecHits"));
}

bool
GEMSkim::filter(edm::Event & iEvent, edm::EventSetup const& iSetup)
{
  edm::Handle<GEMRecHitCollection> gemRecHits;  
  iEvent.getByToken(gemRecHits_, gemRecHits);

  if (gemRecHits->size()) return true;

  return false;
}
void GEMSkim::beginJob(){}
void GEMSkim::endJob(){}
GEMSkim::~GEMSkim(){}

//define this as a plug-in
DEFINE_FWK_MODULE(GEMSkim);
