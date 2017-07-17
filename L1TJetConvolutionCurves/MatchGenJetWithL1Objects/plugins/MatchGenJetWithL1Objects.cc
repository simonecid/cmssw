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

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "TTree.h"

#include <utility>

struct PtEtaPhi {
  float pt, eta, phi;
};

class MatchGenJetWithL1Objects : public edm::one::EDAnalyzer<edm::one::SharedResources> {
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
    const std::vector < std::tuple < const T *, const reco::GenJet *, float > > 
    _matchGenJetWithL1Object
    (
      const edm::Handle< std::vector< reco::GenJet > > &,
      const edm::Handle < BXVector < T > > & 
    );
    
    // Another lovely template function to fill trees
    template <class T> // <3
    void
    _fillTreeWithMatchedPairs
    (
      TTree &,
      const std::vector < std::tuple <const T*, const reco::GenJet*, float > > &
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

    TTree * _l1tMuonGenJetTree;
    TTree * _l1tJetGenJetTree;
    TTree * _l1tEGammaGenJetTree;
    TTree * _l1tTauGenJetTree;
    TTree * _genJetTree;

    PtEtaPhi _genJetPtEtaPhi;
    PtEtaPhi _l1tObjectPtEtaPhi;
    float _deltaR2;
    unsigned int _eventNumber;
};

MatchGenJetWithL1Objects::MatchGenJetWithL1Objects(const edm::ParameterSet& iConfig):
_eventNumber(0)
{  
  this -> _getTokens(iConfig);
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  this -> _l1tMuonGenJetTree = fs -> make<TTree>("matchedL1TMuonGenJetTree", "TTree with generator-level jet / L1T Muon information");
  this -> _l1tJetGenJetTree = fs -> make<TTree>("matchedL1TJetGenJetTree", "TTree with generator-level jet / L1T Jet information");
  this -> _l1tEGammaGenJetTree = fs -> make<TTree>("matchedL1TEGammaGenJetTree", "TTree with generator-level jet / L1T EGamma information");
  this -> _l1tTauGenJetTree = fs -> make<TTree>("matchedL1TTauGenJetTree", "TTree with generator-level jet / L1T Tau information");
  this -> _genJetTree = fs -> make<TTree>("genJetTree", "TTree with generator-level jet information");

  this -> _l1tMuonGenJetTree -> Branch("eventNumber", & (this -> _eventNumber), "eventNumber/i");
  this -> _l1tMuonGenJetTree -> Branch("genJet_pt", &(this -> _genJetPtEtaPhi.pt), "genJet_pt/F");
  this -> _l1tMuonGenJetTree -> Branch("genJet_eta", &(this -> _genJetPtEtaPhi.eta), "genJet_eta/F");
  this -> _l1tMuonGenJetTree -> Branch("genJet_phi", &(this -> _genJetPtEtaPhi.phi), "genJet_phi/F");
  this -> _l1tMuonGenJetTree -> Branch("l1tMuon_pt", &(this -> _l1tObjectPtEtaPhi.pt), "l1tMuon_pt/F");
  this -> _l1tMuonGenJetTree -> Branch("l1tMuon_eta", &(this -> _l1tObjectPtEtaPhi.eta), "l1tMuon_eta/F");
  this -> _l1tMuonGenJetTree -> Branch("l1tMuon_phi", &(this -> _l1tObjectPtEtaPhi.phi), "l1tMuon_phi/F");
  this -> _l1tMuonGenJetTree -> Branch("deltaR2", &(this -> _deltaR2), "deltaR2/F");

  this -> _l1tJetGenJetTree -> Branch("eventNumber", & (this -> _eventNumber), "eventNumber/i");
  this -> _l1tJetGenJetTree -> Branch("genJet_pt", &(this -> _genJetPtEtaPhi.pt), "genJet_pt/F");
  this -> _l1tJetGenJetTree -> Branch("genJet_eta", &(this -> _genJetPtEtaPhi.eta), "genJet_eta/F");
  this -> _l1tJetGenJetTree -> Branch("genJet_phi", &(this -> _genJetPtEtaPhi.phi), "genJet_phi/F");
  this -> _l1tJetGenJetTree -> Branch("l1tJet_pt", &(this -> _l1tObjectPtEtaPhi.pt), "l1tJet_pt/F");
  this -> _l1tJetGenJetTree -> Branch("l1tJet_eta", &(this -> _l1tObjectPtEtaPhi.eta), "l1tJet_eta/F");
  this -> _l1tJetGenJetTree -> Branch("l1tJet_phi", &(this -> _l1tObjectPtEtaPhi.phi), "l1tJet_phi/F");
  this -> _l1tJetGenJetTree -> Branch("deltaR2", &(this -> _deltaR2), "deltaR2/F");

  this -> _l1tEGammaGenJetTree -> Branch("eventNumber", & (this -> _eventNumber), "eventNumber/i");
  this -> _l1tEGammaGenJetTree -> Branch("genJet_pt", &(this -> _genJetPtEtaPhi.pt), "genJet_pt/F");
  this -> _l1tEGammaGenJetTree -> Branch("genJet_eta", &(this -> _genJetPtEtaPhi.eta), "genJet_eta/F");
  this -> _l1tEGammaGenJetTree -> Branch("genJet_phi", &(this -> _genJetPtEtaPhi.phi), "genJet_phi/F");
  this -> _l1tEGammaGenJetTree -> Branch("l1tEGamma_pt", &(this -> _l1tObjectPtEtaPhi.pt), "l1tEGamma_pt/F");
  this -> _l1tEGammaGenJetTree -> Branch("l1tEGamma_eta", &(this -> _l1tObjectPtEtaPhi.eta), "l1tEGamma_eta/F");
  this -> _l1tEGammaGenJetTree -> Branch("l1tEGamma_phi", &(this -> _l1tObjectPtEtaPhi.phi), "l1tEGamma_phi/F");
  this -> _l1tEGammaGenJetTree -> Branch("deltaR2", &(this -> _deltaR2), "deltaR2/F");

  this -> _l1tTauGenJetTree -> Branch("eventNumber", & (this -> _eventNumber), "eventNumber/i");
  this -> _l1tTauGenJetTree -> Branch("genJet_pt", &(this -> _genJetPtEtaPhi.pt), "genJet_pt/F");
  this -> _l1tTauGenJetTree -> Branch("genJet_eta", &(this -> _genJetPtEtaPhi.eta), "genJet_eta/F");
  this -> _l1tTauGenJetTree -> Branch("genJet_phi", &(this -> _genJetPtEtaPhi.phi), "genJet_phi/F");
  this -> _l1tTauGenJetTree -> Branch("l1tTau_pt", &(this -> _l1tObjectPtEtaPhi.pt), "l1tTau_pt/F");
  this -> _l1tTauGenJetTree -> Branch("l1tTau_eta", &(this -> _l1tObjectPtEtaPhi.eta), "l1tTau_eta/F");
  this -> _l1tTauGenJetTree -> Branch("l1tTau_phi", &(this -> _l1tObjectPtEtaPhi.phi), "l1tTau_phi/F");
  this -> _l1tTauGenJetTree -> Branch("deltaR2", &(this -> _deltaR2), "deltaR2/F");

  //Used to detemine the prob that a jet will be misidentified binned in pt
  this -> _genJetTree -> Branch("eventNumber", & (this -> _eventNumber), "eventNumber/i");
  this -> _genJetTree -> Branch("genJet_pt", &(this -> _genJetPtEtaPhi.pt), "genJet_pt/F");
  this -> _genJetTree -> Branch("genJet_eta", &(this -> _genJetPtEtaPhi.eta), "genJet_eta/F");
  this -> _genJetTree -> Branch("genJet_phi", &(this -> _genJetPtEtaPhi.phi), "genJet_phi/F");

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
    std::cerr << ">>> l1tMuonCollectionTag not found, proceeding without computing the corresponding convolution." << std::endl;
    this -> _l1tMuonCollectionTag = NULL;
  }
  
  try
  {
    this -> _l1tJetCollectionTag = new edm::EDGetTokenT< BXVector < l1t::Jet > >(consumes< BXVector< l1t::Jet > > (iConfig.getParameter< edm::InputTag >("l1tJetCollectionTag")));
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> l1tJetCollectionTag not found, proceeding without computing the corresponding convolution." << std::endl;
    this -> _l1tJetCollectionTag = NULL;
  }
  
  try
  {
    this -> _l1tEGammaCollectionTag = new edm::EDGetTokenT< BXVector < l1t::EGamma > >(consumes< BXVector< l1t::EGamma > > (iConfig.getParameter< edm::InputTag >("l1tEGammaCollectionTag")));
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> l1tEGammaCollectionTag not found, proceeding without computing the corresponding convolution." << std::endl;
    this -> _l1tEGammaCollectionTag = NULL;
  }
  
  try
  {
    this -> _l1tTauCollectionTag = new edm::EDGetTokenT< BXVector < l1t::Tau > >(consumes< BXVector< l1t::Tau > > (iConfig.getParameter< edm::InputTag >("l1tTauCollectionTag")));
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> l1tTauCollectionTag not found, proceeding without computing the corresponding convolution." << std::endl;
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

  if (this -> _l1tMuonCollectionTag)
  {
    edm::Handle < BXVector< l1t::Muon > > l1tMuonCollectionHandle;
    iEvent.getByToken(*(this -> _l1tMuonCollectionTag), l1tMuonCollectionHandle);
    auto l1tMuonGenJetPairs = this -> _matchGenJetWithL1Object<>(genJetCollectionHandle, l1tMuonCollectionHandle);
    this -> _fillTreeWithMatchedPairs<>(*(this -> _l1tMuonGenJetTree), l1tMuonGenJetPairs);
  }
  if (this -> _l1tJetCollectionTag)
  {
    edm::Handle < BXVector< l1t::Jet > > l1tJetCollectionHandle;
    iEvent.getByToken(*(this -> _l1tJetCollectionTag), l1tJetCollectionHandle);
    auto l1tJetGenJetPairs = this -> _matchGenJetWithL1Object<>(genJetCollectionHandle, l1tJetCollectionHandle);
    this -> _fillTreeWithMatchedPairs<>(*(this -> _l1tJetGenJetTree), l1tJetGenJetPairs);
  }
  if (this -> _l1tEGammaCollectionTag)
  {
    edm::Handle < BXVector< l1t::EGamma > > l1tEGammaCollectionHandle;
    iEvent.getByToken(*(this -> _l1tEGammaCollectionTag), l1tEGammaCollectionHandle);
    auto l1tEGammaGenJetPairs = this -> _matchGenJetWithL1Object<>(genJetCollectionHandle, l1tEGammaCollectionHandle);
    this -> _fillTreeWithMatchedPairs<>(*(this -> _l1tEGammaGenJetTree), l1tEGammaGenJetPairs);
  }
  if (this -> _l1tTauCollectionTag)
  {
    edm::Handle < BXVector< l1t::Tau > > l1tTauCollectionHandle;
    iEvent.getByToken(*(this -> _l1tTauCollectionTag), l1tTauCollectionHandle);
    auto l1tTauGenJetPairs = this -> _matchGenJetWithL1Object<>(genJetCollectionHandle, l1tTauCollectionHandle);
    this -> _fillTreeWithMatchedPairs<>(*(this -> _l1tTauGenJetTree), l1tTauGenJetPairs);
  }

  for (const reco::GenJet & genJet : *genJetCollectionHandle)
  {
    this -> _genJetPtEtaPhi.pt = genJet.pt();
    this -> _genJetPtEtaPhi.eta = genJet.eta();
    this -> _genJetPtEtaPhi.phi = genJet.phi();
    this -> _genJetTree -> Fill();
  }

  this -> _eventNumber++;

}

// Another lovely template function to fill trees
template <class T> // <3
void
MatchGenJetWithL1Objects::_fillTreeWithMatchedPairs
(
  TTree & aTree,
  const std::vector < std::tuple <const T*, const reco::GenJet*, float > > & matchedL1TObjectJetPairs
)
{
  for (const auto & matchTuple : matchedL1TObjectJetPairs)
  {
    const T* l1tObject = std::get<0>(matchTuple);
    const reco::GenJet* genJet = std::get<1>(matchTuple);
    float deltaR2 = std::get<2>(matchTuple);
    this -> _genJetPtEtaPhi.pt = genJet -> pt();
    this -> _genJetPtEtaPhi.eta = genJet -> eta();
    this -> _genJetPtEtaPhi.phi = genJet -> phi();
    this -> _l1tObjectPtEtaPhi.pt = l1tObject -> pt();
    this -> _l1tObjectPtEtaPhi.eta = l1tObject -> eta();
    this -> _l1tObjectPtEtaPhi.phi = l1tObject -> phi();
    this -> _deltaR2 = deltaR2;
    aTree.Fill();
  }
}

template <class T> // <3
const std::vector <std::tuple < const T *, const reco::GenJet *, float > > 
MatchGenJetWithL1Objects::_matchGenJetWithL1Object
(
  const edm::Handle < std::vector< reco::GenJet > > & genJetCollectionHandle,
  const edm::Handle < BXVector < T > > & l1tObjectCollectionHandle 
)
{

  std::vector< std::tuple < const T *, const reco::GenJet * , float > > l1tObjectGenJetPairs;
  // for each object in the l1t collection we look for the closest jet in a wide range
  for (
    typename BXVector<T>::const_iterator bx0Iterator = l1tObjectCollectionHandle -> begin(0);
    bx0Iterator != l1tObjectCollectionHandle -> end(0);
    bx0Iterator++
  )
  {
    
    bool foundMatch = false;
    const T & l1tObject = *bx0Iterator;
    float dr2Min = 25; // i.e. dr = 5
    std::tuple<const T *, const reco::GenJet *, float > l1tObjectGenJetPair;
    
    for (const reco::GenJet & genJet : *genJetCollectionHandle)
    {
      float dr2 = reco::deltaR2(l1tObject, genJet);
      if (dr2 < dr2Min)
      {
        foundMatch = true;
        dr2Min = dr2;
        l1tObjectGenJetPair = std::make_tuple(&l1tObject, &genJet, dr2Min);
      }
    }
    
    //if we have found a compatible gen jet we push the object-jet pair in a vector which will be our result
    if (foundMatch)
      l1tObjectGenJetPairs.push_back(l1tObjectGenJetPair);
  }

  return l1tObjectGenJetPairs;
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
