// -*- C++ -*-
//
// Package: L1CaloTrigger
// Class: PrintMomentum
//
/**\class PrintMomentum PrintMomentum.cc L1Trigger/L1CaloTrigger/plugin/PrintMomentum.cc

Description: Dumps momentums

*/
//
// Original Simone Bologna
// Created: Mon Jul 02 2018
//


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/Phase2L1ParticleFlow/interface/PFCandidate.h"
#include "DataFormats/Phase2L1ParticleFlow/interface/PFCluster.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/L1Trigger/interface/L1Candidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"

#include "TH2F.h"

class PrintMomentum : public edm::EDAnalyzer {
   public:
      explicit PrintMomentum(const edm::ParameterSet&);
      ~PrintMomentum();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      /// Get the energy of a certain tower while correctly handling phi periodicity in case of overflow
      // <3 handy method to print stuff
      template <class Collection>
      void _print(const Collection & aCollection, const std::string & collectionName = "");

      edm::EDGetTokenT<BXVector<l1t::Jet>> *_phase1L1TJetFromPfCandidatesTag;
      edm::EDGetTokenT<std::vector<reco::GenJet>> *_genJetCollectionTag;
      edm::EDGetTokenT<BXVector<l1t::Jet>> *_phase1L1TJetFromPfClustersTag;

};

PrintMomentum::PrintMomentum(const edm::ParameterSet& iConfig)
{
  // Retrieving the pfCandidates and pfClusters if input tag has been provided
  try
  {
    this -> _genJetCollectionTag = new edm::EDGetTokenT< std::vector<reco::GenJet> >(consumes< std::vector<reco::GenJet> > (iConfig.getParameter< edm::InputTag >("genJetCollectionTag")));
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> genJetCollectionTag not found" << std::endl;
    this -> _genJetCollectionTag = NULL;
  }

  try
  {
    this -> _phase1L1TJetFromPfCandidatesTag = new edm::EDGetTokenT< BXVector<l1t::Jet> >(consumes< BXVector<l1t::Jet> > (iConfig.getParameter< edm::InputTag >("phase1L1TJetFromPfCandidatesTag")));
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> phase1L1TJetFromPfCandidatesTag not found" << std::endl;
    this -> _phase1L1TJetFromPfCandidatesTag = NULL;
  }
  
  try
  {
    this -> _phase1L1TJetFromPfClustersTag = new edm::EDGetTokenT< BXVector<l1t::Jet> >(consumes< BXVector<l1t::Jet> > (iConfig.getParameter< edm::InputTag >("phase1L1TJetFromPfClustersTag")));
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> phase1L1TJetFromPfClustersTag not found" << std::endl;
    this -> _phase1L1TJetFromPfClustersTag = NULL;
  }

}

PrintMomentum::~PrintMomentum()
{
  if (this -> _genJetCollectionTag) delete this -> _genJetCollectionTag;
  if (this -> _phase1L1TJetFromPfCandidatesTag) delete this -> _phase1L1TJetFromPfCandidatesTag;
  if (this -> _phase1L1TJetFromPfClustersTag) delete this -> _phase1L1TJetFromPfClustersTag;
}

void PrintMomentum::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // Retrieving the pfCandidates and pfClusters if input tag has been provided

  if (this -> _genJetCollectionTag) {
    edm::Handle < std::vector< reco::GenJet > > genJetCollectionHandle;
    iEvent.getByToken(*(this -> _genJetCollectionTag), genJetCollectionHandle);
    this -> _print<>(*genJetCollectionHandle, "Jets from Generator");
  }

  if (this -> _phase1L1TJetFromPfCandidatesTag) {
    edm::Handle < BXVector<l1t::Jet> > phase1L1TJetFromPfCandidatesHandle;
    iEvent.getByToken(*(this -> _phase1L1TJetFromPfCandidatesTag), phase1L1TJetFromPfCandidatesHandle);
    // dumping the data
    //std::cout << ">>>>>> DUMPING PFCANDIDATES <<<<<<" << std::endl;
    //for (auto pfCandidateIterator = pfCandidateCollectionHandle -> begin(); pfCandidateIterator != pfCandidateCollectionHandle -> end(); pfCandidateIterator++) 
    //{
    //  std::cout << pfCandidateIterator -> pt() << "\t" << pfCandidateIterator -> eta() << "\t" << pfCandidateIterator -> phi() << "\t" << std::endl;
    //}
    this -> _print<>(*phase1L1TJetFromPfCandidatesHandle, "Jets from Candidates");
  }

  if (this -> _phase1L1TJetFromPfClustersTag) {
    edm::Handle < BXVector<l1t::Jet> > phase1L1TJetFromPfClustersHandle;
    iEvent.getByToken(*(this -> _phase1L1TJetFromPfClustersTag), phase1L1TJetFromPfClustersHandle);
    // dumping the data
    //std::cout << ">>>>>> DUMPING PFCLUSTERS <<<<<<" << std::endl;
    //for (auto pfClusterIterator = pfClusterCollectionHandle -> begin(); pfClusterIterator != pfClusterCollectionHandle -> end(); pfClusterIterator++) 
    //{
    //  std::cout << pfClusterIterator -> pt() << "\t" << pfClusterIterator -> eta() << "\t" << pfClusterIterator -> phi() << "\t" << std::endl;
    //}
    this -> _print<>(*phase1L1TJetFromPfClustersHandle, "Jets from Clusters");
  }


  return;

}

template <class Collection>
void PrintMomentum::_print(const Collection & aCollection, const std::string & collectionName)
{
  std::cout << ">>> Dumping collection " << collectionName << std::endl;
  //Filling the calo grid with the primitives
  for (const auto & collectionItem : aCollection) 
  {
    std::cout << "(pt, eta, phi): (" << collectionItem.pt() << ", " << collectionItem.eta() << ", " << collectionItem.phi() << ")" << std::endl;
  }
  return;
}

void
PrintMomentum::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(PrintMomentum);