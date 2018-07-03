// -*- C++ -*-
//
// Package: L1CaloTrigger
// Class: L1TJetPhase1Producer
//
/**\class L1TJetPhase1Producer L1TJetPhase1Producer.cc L1Trigger/L1CaloTrigger/plugin/L1TJetPhase1Producer.cc

Description: Produces jets with sliding window algorithm using pfcluster and pfcandidates

*/
//
// Original Simone Bologna
// Created: Mon Jul 02 2018
//


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Phase2L1ParticleFlow/interface/PFCandidate.h"
#include "DataFormats/Phase2L1ParticleFlow/interface/PFCluster.h"
#include "DataFormats/L1Trigger/interface/L1Candidate.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"

class L1TJetPhase1Producer : public edm::EDProducer {
   public:
      explicit L1TJetPhase1Producer(const edm::ParameterSet&);
      ~L1TJetPhase1Producer();

   private:
      virtual void produce(edm::Event&, const edm::EventSetup&);
      edm::EDGetTokenT<std::vector<l1t::PFCandidate>> *_pfCandidateCollectionTag;
      edm::EDGetTokenT<std::vector<l1t::PFCluster>> *_pfClusterCollectionTag;

};

L1TJetPhase1Producer::L1TJetPhase1Producer(const edm::ParameterSet& iConfig)
{
  // Retrieving the pfCandidates and pfClusters if input tag has been provided
  try
  {
    this -> _pfCandidateCollectionTag = new edm::EDGetTokenT< std::vector<l1t::PFCandidate> >(consumes< std::vector<l1t::PFCandidate> > (iConfig.getParameter< edm::InputTag >("pfCandidateCollectionTag")));
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> pfCandidateCollectionTag not found" << std::endl;
    this -> _pfCandidateCollectionTag = NULL;
  }
  
  try
  {
    this -> _pfClusterCollectionTag = new edm::EDGetTokenT< std::vector<l1t::PFCluster> >(consumes< std::vector<l1t::PFCluster> > (iConfig.getParameter< edm::InputTag >("pfClusterCollectionTag")));
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> pfClusterCollectionTag not found" << std::endl;
    this -> _pfClusterCollectionTag = NULL;
  }
}

L1TJetPhase1Producer::~L1TJetPhase1Producer()
{
  if (this -> _pfCandidateCollectionTag) delete this -> _pfCandidateCollectionTag;
  if (this -> _pfClusterCollectionTag) delete this -> _pfClusterCollectionTag;
}

void L1TJetPhase1Producer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // Retrieving the pfCandidates and pfClusters if input tag has been provided

  if (this -> _pfCandidateCollectionTag) {
    edm::Handle < std::vector< l1t::PFCandidate > > pfCandidateCollectionHandle;
    iEvent.getByToken(*(this -> _pfCandidateCollectionTag), pfCandidateCollectionHandle);
    // dumping the data
    std::cout << ">>>>>> DUMPING PFCANDIDATES <<<<<<" << std::endl;
    for (auto pfCandidateIterator = pfCandidateCollectionHandle -> begin(); pfCandidateIterator != pfCandidateCollectionHandle -> end(); pfCandidateIterator++) 
    {
      std::cout << pfCandidateIterator -> pt() << "\t" << pfCandidateIterator -> eta() << "\t" << pfCandidateIterator -> phi() << "\t" << std::endl;
    }
  }

  if (this -> _pfClusterCollectionTag) {
    edm::Handle < std::vector< l1t::PFCluster > > pfClusterCollectionHandle;
    iEvent.getByToken(*(this -> _pfClusterCollectionTag), pfClusterCollectionHandle);
    // dumping the data
    std::cout << ">>>>>> DUMPING PFCLUSTERS <<<<<<" << std::endl;
    for (auto pfClusterIterator = pfClusterCollectionHandle -> begin(); pfClusterIterator != pfClusterCollectionHandle -> end(); pfClusterIterator++) 
    {
      std::cout << pfClusterIterator -> pt() << "\t" << pfClusterIterator -> eta() << "\t" << pfClusterIterator -> phi() << "\t" << std::endl;
    }
  }

  return;

}

DEFINE_FWK_MODULE(L1TJetPhase1Producer);