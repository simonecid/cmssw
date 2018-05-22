// -*- C++ -*-
//
// Package:    l1tEtSumConvolutionCurves/MatchGenJetWithL1Objects
// Class:      MatchGenJetWithL1Objects
// 
/**\class MatchGenJetWithL1Objects MatchGenJetWithL1Objects.cc l1tEtSumConvolutionCurves/MatchGenJetWithL1Objects/plugins/MatchGenJetWithL1Objects.cc

 Description: Matches a gen jet with a L1T object and creates a tree

 Implementation:
    Receives tags of the various L1T objects
*/
//
// Original Author:  Simone Bologna
//         Created:  Wed, 12 Jul 2017 14:21:56 GMT
//
//


// system include files
#include <memory>
#include <signal.h>

// user include files
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/L1Trigger/interface/BXVector.h"

#include "DataFormats/L1Trigger/interface/EtSum.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "TTree.h"

#include <utility>
#include "L1TJetConvolutionCurves/MatchGenJetWithL1Objects/plugins/MatchingAlgorithms.h"


struct Particle {
  unsigned int id;
  float pt, eta, phi;
};

struct TriggerObject {
  unsigned int id;
  float pt, eta, phi;
  int hwQual;
};

class ComputeMHT : public edm::one::EDAnalyzer<edm::one::SharedResources> {
  public:
    explicit ComputeMHT(const edm::ParameterSet&);
    ~ComputeMHT();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    edm::EDGetTokenT< BXVector<l1t::EtSum> >* _l1tEtSumCollectionTag;

    TTree * _MHTTree;

    float _mht;
};

ComputeMHT::ComputeMHT(const edm::ParameterSet& iConfig)
{  
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  this -> _MHTTree = fs -> make<TTree>("mhtTree", "TTree with L1T MHT");
  this -> _MHTTree -> Branch("mht", &(this -> _mht), "mht/F");

  this -> _l1tEtSumCollectionTag = new edm::EDGetTokenT< BXVector<l1t::EtSum> >(consumes< BXVector<l1t::EtSum> > (iConfig.getParameter< edm::InputTag >("l1tEtSumCollectionTag")));
  
}

ComputeMHT::~ComputeMHT()
{
  if (this -> _l1tEtSumCollectionTag) delete this -> _l1tEtSumCollectionTag;
}


//
// member functions
//

// ------------ method called for each event  ------------
void
ComputeMHT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //Retrieving gen and l1t stuff

  edm::Handle < BXVector< l1t::EtSum > > l1tEtSumCollectionHandle;
  iEvent.getByToken(*(this -> _l1tEtSumCollectionTag), l1tEtSumCollectionHandle);

  for (auto l1tEtSumIterator = l1tEtSumCollectionHandle -> begin(0); l1tEtSumIterator != l1tEtSumCollectionHandle -> end(0); l1tEtSumIterator++ )
  {
    if (l1tEtSumIterator -> getType() == l1t::EtSum::EtSumType::kMissingHt)
      std::cout << "EtSum " << l1tEtSumIterator -> pt() << "\t" << "Type MHT: " << l1tEtSumIterator -> getType() << std::endl;
    if (l1tEtSumIterator -> getType() == l1t::EtSum::EtSumType::kMissingHtHF)
      std::cout << "EtSum " << l1tEtSumIterator -> pt() << "\t" << "Type MHTHF: " << l1tEtSumIterator -> getType() << std::endl;
  }

}

// ------------ method called once each job just before starting event loop  ------------
void ComputeMHT::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void ComputeMHT::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void ComputeMHT::fillDescriptions(edm::ConfigurationDescriptions& descriptions) 
{
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(ComputeMHT);
