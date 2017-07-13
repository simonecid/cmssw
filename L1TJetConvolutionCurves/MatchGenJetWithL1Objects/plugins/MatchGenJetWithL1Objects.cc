// -*- C++ -*-
//
// Package:    L1TJetConvolutionCurves/MatchGenJetWithL1Objects
// Class:      MatchGenJetWithL1Objects
// 
/**\class MatchGenJetWithL1Objects MatchGenJetWithL1Objects.cc L1TJetConvolutionCurves/MatchGenJetWithL1Objects/plugins/MatchGenJetWithL1Objects.cc

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

// user include files
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/JetReco/interface/GenJet.h"

#include "DataFormats/L1Trigger/interface/BXVector.h"

#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/L1Trigger/interface/Tau.h"


class MatchGenJetWithL1Objects : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:
    explicit MatchGenJetWithL1Objects(const edm::ParameterSet&);
    ~MatchGenJetWithL1Objects();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    void _getTokens(const edm::ParameterSet&);
    void _freeTokens();

    // Beautiful template function to perform matching between a gen jet and whatever L1T object
    template <class T> // <3
    void _matchGenJetWithL1Object
    (
      const edm::Event& iEvent,
      const edm::Handle< std::vector< reco::GenJet > > &,
      const edm::EDGetTokenT < BXVector < T > > & 
    );
  
    edm::EDGetTokenT< std::vector< reco::GenJet > > *_genJetCollectionTag;
    edm::EDGetTokenT< BXVector< l1t::Muon > > *_l1tMuonCollectionTag;
    edm::EDGetTokenT< BXVector< l1t::Jet > > *_l1tJetCollectionTag;
    //!TODO: It should splitted in two different objects in the Stage-2
    edm::EDGetTokenT< BXVector< l1t::EGamma > > *_l1tEGammaCollectionTag;
    //!TODO: Like this
    //edm::EDGetTokenT< BXVector< l1t::Electron > > _l1tElectronCollectionTag;
    //edm::EDGetTokenT< BXVector< l1t::Gamma > > _l1tGammaCollectionTag;
    edm::EDGetTokenT< BXVector< l1t::Tau > > *_l1tTauCollectionTag;
    
    // ----------member data ---------------------------
};

MatchGenJetWithL1Objects::MatchGenJetWithL1Objects(const edm::ParameterSet& iConfig)
{  
  this -> _getTokens(iConfig);
  usesResource("TFileService");
}

void MatchGenJetWithL1Objects::_getTokens(const edm::ParameterSet& iConfig)
{
  
  this -> _genJetCollectionTag = new edm::EDGetTokenT< std::vector< reco::GenJet > >(consumes< std::vector< reco::GenJet > > (iConfig.getParameter< edm::InputTag >("genJetCollectionTag")));
  // Taking the tag of the various L1T object collections
  // If a parameter is omitted that object will not be studied
  try
  {
    this -> _l1tMuonCollectionTag = new edm::EDGetTokenT< BXVector < l1t::Muon > >(consumes< BXVector< l1t::Muon > > (iConfig.getParameter< edm::InputTag >("l1tMuonCollectionTag")));
  } catch (std::exception const & ex) 
  {
    this -> _l1tMuonCollectionTag = NULL;
  }
  
  try
  {
    this -> _l1tJetCollectionTag = new edm::EDGetTokenT< BXVector < l1t::Jet > >(consumes< BXVector< l1t::Jet > > (iConfig.getParameter< edm::InputTag >("l1tJetCollectionTag")));
  } catch (std::exception const & ex) 
  {
    this -> _l1tJetCollectionTag = NULL;
  }
  
  try
  {
    this -> _l1tEGammaCollectionTag = new edm::EDGetTokenT< BXVector < l1t::EGamma > >(consumes< BXVector< l1t::EGamma > > (iConfig.getParameter< edm::InputTag >("l1tEGammaCollectionTag")));
  } catch (std::exception const & ex) 
  {
    this -> _l1tEGammaCollectionTag = NULL;
  }
  
  try
  {
    this -> _l1tTauCollectionTag = new edm::EDGetTokenT< BXVector < l1t::Tau > >(consumes< BXVector< l1t::Tau > > (iConfig.getParameter< edm::InputTag >("l1tTauCollectionTag")));
  } catch (std::exception const & ex) 
  {
    this -> _l1tTauCollectionTag = NULL;
  }
  
  return;
}

void MatchGenJetWithL1Objects::_freeTokens()
{
  if (this -> _genJetCollectionTag) delete this -> _genJetCollectionTag;
  if (this -> _l1tMuonCollectionTag) delete this -> _l1tMuonCollectionTag;
  if (this -> _l1tJetCollectionTag) delete this -> _l1tJetCollectionTag;
  if (this -> _l1tEGammaCollectionTag) delete this -> _l1tEGammaCollectionTag;
  if (this -> _l1tTauCollectionTag) delete this -> _l1tTauCollectionTag;
}

MatchGenJetWithL1Objects::~MatchGenJetWithL1Objects()
{
  this -> _freeTokens();
}


//
// member functions
//

// ------------ method called for each event  ------------
void
MatchGenJetWithL1Objects::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //Retrieving gen and l1t stuff
  edm::Handle < std::vector< reco::GenJet > > genJetCollectionHandle;
  iEvent.getByToken(*(this -> _genJetCollectionTag), genJetCollectionHandle);
  if (this -> _l1tMuonCollectionTag != NULL)
    this -> _matchGenJetWithL1Object<>(iEvent, genJetCollectionHandle, *(this -> _l1tMuonCollectionTag));
  if (this -> _l1tJetCollectionTag != NULL)
    this -> _matchGenJetWithL1Object<>(iEvent, genJetCollectionHandle, *(this -> _l1tJetCollectionTag));
  if (this -> _l1tEGammaCollectionTag != NULL)
    this -> _matchGenJetWithL1Object<>(iEvent, genJetCollectionHandle, *(this -> _l1tEGammaCollectionTag));
  if (this -> _l1tTauCollectionTag != NULL)
    this -> _matchGenJetWithL1Object<>(iEvent, genJetCollectionHandle, *(this -> _l1tTauCollectionTag));
}

template <class T> // <3
void MatchGenJetWithL1Objects::_matchGenJetWithL1Object
(
  const edm::Event& iEvent,
  const edm::Handle < std::vector< reco::GenJet > > & genJetCollectionHandle,
  const edm::EDGetTokenT < BXVector < T > > & l1tObjectCollectionTag 
)
{
  edm::Handle < BXVector< T > > l1tObjectCollectionHandle;
  iEvent.getByToken(l1tObjectCollectionTag, l1tObjectCollectionHandle);
}

// ------------ method called once each job just before starting event loop  ------------
void MatchGenJetWithL1Objects::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void MatchGenJetWithL1Objects::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void MatchGenJetWithL1Objects::fillDescriptions(edm::ConfigurationDescriptions& descriptions) 
{
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(MatchGenJetWithL1Objects);
