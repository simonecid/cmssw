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
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/Phase2L1ParticleFlow/interface/PFCandidate.h"
#include "DataFormats/Phase2L1ParticleFlow/interface/PFCluster.h"
#include "DataFormats/L1Trigger/interface/L1Candidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH2F.h"

#include <csignal>
#include <algorithm>

//UNCOMMENT TO CREATE DEBUG HISTO
//#define DEBUG

//class L1TJetPhase1Producer : public edm::EDProducer {
class L1TJetPhase1Producer : public edm::one::EDProducer<edm::one::SharedResources> {
   public:
      explicit L1TJetPhase1Producer(const edm::ParameterSet&);
      ~L1TJetPhase1Producer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void produce(edm::Event&, const edm::EventSetup&);
      /// Finds the seeds in the caloGrid, seeds are saved in a vector that contain the index in the TH2F of each seed
      std::vector<std::tuple<int, int>> _findSeeds(const TH2F & caloGrid, float seedThreshold);
      BXVector<l1t::Jet> _buildJetsFromSeeds(const TH2F & caloGrid, const std::vector<std::tuple<int, int>> & seeds);
      void _subtract9x9Pileup(const TH2F & caloGrid, l1t::Jet & l1tjet);
      /// Get the energy of a certain tower while correctly handling phi periodicity in case of overflow
      float _getTowerEnergy(const TH2F & caloGrid, int iEta, int iPhi);

      // <3 handy method to fill the calogrid with whatever type
      template <class TriggerPrimitive>
      void _fillCaloGrid(TH2F & caloGrid, const std::vector<TriggerPrimitive> & triggerPrimitives);

      edm::EDGetTokenT<std::vector<l1t::PFCandidate>> *_pfCandidateCollectionTag;
      edm::EDGetTokenT<std::vector<l1t::PFCluster>> *_pfClusterCollectionTag;
      TH2F* _caloGridPfCandidate;
      TH2F* _caloGridPfCluster;

      std::vector<double> _etaBinning;
      size_t _nBinsEta;
      unsigned int _nBinsPhi;
      double _phiLow;
      double _phiUp;
      unsigned int _jetIEtaSize;
      unsigned int _jetIPhiSize;
      double _seedPtThreshold;

};

L1TJetPhase1Producer::L1TJetPhase1Producer(const edm::ParameterSet& iConfig)
{

  usesResource("TFileService");
  edm::Service<TFileService> fs;

  this -> _etaBinning = iConfig.getParameter<std::vector<double> >("etaBinning");
  this -> _nBinsEta = this -> _etaBinning.size() - 1;
  this -> _nBinsPhi = iConfig.getParameter<unsigned int>("nBinsPhi");
  this -> _phiLow = iConfig.getParameter<double>("phiLow");
  this -> _phiUp = iConfig.getParameter<double>("phiUp");
  this -> _jetIEtaSize = iConfig.getParameter<unsigned int>("jetIEtaSize");
  this -> _jetIPhiSize = iConfig.getParameter<unsigned int>("jetIPhiSize");
  this -> _seedPtThreshold = iConfig.getParameter<double>("seedPtThreshold");

  // Retrieving the pfCandidates and pfClusters if input tag has been provided
  try
  {
    this -> _pfCandidateCollectionTag = new edm::EDGetTokenT< std::vector<l1t::PFCandidate> >(consumes< std::vector<l1t::PFCandidate> > (iConfig.getParameter< edm::InputTag >("pfCandidateCollectionTag")));
    #ifndef DEBUG
    this -> _caloGridPfCandidate = new TH2F("caloGridPfCandidate", "Calorimeter grid", this -> _nBinsEta, this-> _etaBinning.data(), this -> _nBinsPhi, this -> _phiLow, this -> _phiUp);
    #endif
    #ifdef DEBUG
    this -> _caloGridPfCandidate = fs -> make<TH2F>("caloGridPfCandidate", "Calorimeter grid", this -> _nBinsEta, this-> _etaBinning.data(), this -> _nBinsPhi, this -> _phiLow, this -> _phiUp);
    #endif
    this -> _caloGridPfCandidate -> GetXaxis() -> SetTitle("#eta");
    this -> _caloGridPfCandidate -> GetYaxis() -> SetTitle("#phi");
    produces<BXVector<l1t::Jet> >( "Phase1L1TJetFromPfCandidates" ).setBranchAlias("Phase1L1TJetFromPfCandidates");
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> pfCandidateCollectionTag not found" << std::endl;
    this -> _pfCandidateCollectionTag = NULL;
    this -> _caloGridPfCandidate = NULL;
  }
  
  try
  {
    this -> _pfClusterCollectionTag = new edm::EDGetTokenT< std::vector<l1t::PFCluster> >(consumes< std::vector<l1t::PFCluster> > (iConfig.getParameter< edm::InputTag >("pfClusterCollectionTag")));
    #ifndef DEBUG
    this -> _caloGridPfCluster = new TH2F("caloGridPfCluster", "Calorimeter grid", this -> _nBinsEta, this-> _etaBinning.data(), this -> _nBinsPhi, this -> _phiLow, this -> _phiUp);
    #endif
    #ifdef DEBUG
    this -> _caloGridPfCluster = fs -> make<TH2F>("caloGridPfCluster", "Calorimeter grid", this -> _nBinsEta, this-> _etaBinning.data(), this -> _nBinsPhi, this -> _phiLow, this -> _phiUp);
    #endif
    this -> _caloGridPfCluster -> GetXaxis() -> SetTitle("#eta");
    this -> _caloGridPfCluster -> GetYaxis() -> SetTitle("#phi");
    produces<BXVector<l1t::Jet> >( "Phase1L1TJetFromPfClusters" ).setBranchAlias("Phase1L1TJetFromPfClusters");
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> pfClusterCollectionTag not found" << std::endl;
    this -> _pfClusterCollectionTag = NULL;
    this -> _caloGridPfCluster = NULL;
  }

}

L1TJetPhase1Producer::~L1TJetPhase1Producer()
{
  if (this -> _pfCandidateCollectionTag) delete this -> _pfCandidateCollectionTag;
  if (this -> _pfClusterCollectionTag) delete this -> _pfClusterCollectionTag;
  #ifndef DEBUG
  if (this -> _caloGridPfCandidate) delete this -> _caloGridPfCandidate;
  if (this -> _caloGridPfCluster) delete this -> _caloGridPfCluster;
  #endif
}

float L1TJetPhase1Producer::_getTowerEnergy(const TH2F & caloGrid, int iEta, int iPhi)
{
  // We return the pt of a certain bin in the calo grid, taking account of the phi periodicity when overflowing (e.g. phi > phiSize), and returning 0 for the eta out of bounds

  int nBinsEta = caloGrid.GetNbinsX();
  int nBinsPhi = caloGrid.GetNbinsY();
  while (iPhi < 1) {
    iPhi += nBinsPhi;
  }
  while (iPhi > nBinsPhi) {
    iPhi -= nBinsPhi;
  }
  if (iEta < 1) {
    //std::cout << "Returning (pt, ieta, iphi): " << "(" << 0 << ", " << iEta << ", " << iPhi << ")" << std::endl;
    return 0;
  }
  if (iEta > nBinsEta) {
    //std::cout << "Returning (pt, ieta, iphi): " << "(" << 0 << ", " << iEta << ", " << iPhi << ")" << std::endl;
    return 0;
  }
    //std::cout << "Returning (pt, ieta, iphi): " << "(" << caloGrid.GetBinContent(iEta, iPhi) << ", " << iEta << ", " << iPhi << ")" << std::endl;
  return caloGrid.GetBinContent(iEta, iPhi);
}

void L1TJetPhase1Producer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // Retrieving the pfCandidates and pfClusters if input tag has been provided
  
#ifdef DEBUG
  this -> _caloGridPfCandidate -> Reset();
  this -> _caloGridPfCluster -> Reset();
  #include "L1Trigger/L1CaloTrigger/debug/CaloGridTestingCode.cdebug"
#endif

  if (this -> _pfCandidateCollectionTag) {
    edm::Handle < std::vector< l1t::PFCandidate > > pfCandidateCollectionHandle;
    iEvent.getByToken(*(this -> _pfCandidateCollectionTag), pfCandidateCollectionHandle);
    // dumping the data
    //std::cout << ">>>>>> DUMPING PFCANDIDATES <<<<<<" << std::endl;
    //for (auto pfCandidateIterator = pfCandidateCollectionHandle -> begin(); pfCandidateIterator != pfCandidateCollectionHandle -> end(); pfCandidateIterator++) 
    //{
    //  std::cout << pfCandidateIterator -> pt() << "\t" << pfCandidateIterator -> eta() << "\t" << pfCandidateIterator -> phi() << "\t" << std::endl;
    //}
    
    #ifndef DEBUG
    this -> _caloGridPfCandidate -> Reset();
    this -> _fillCaloGrid<>(*(this -> _caloGridPfCandidate), *pfCandidateCollectionHandle);
    #endif
    const auto seedsVector = this -> _findSeeds(*(this -> _caloGridPfCandidate), this -> _seedPtThreshold); // seedPtThreshold = 6
    auto l1jetVector = this -> _buildJetsFromSeeds(*(this -> _caloGridPfCandidate), seedsVector);
    
    std::unique_ptr< BXVector<l1t::Jet> > l1jetVectorPtr(new BXVector<l1t::Jet>(l1jetVector));
    iEvent.put(std::move(l1jetVectorPtr), "Phase1L1TJetFromPfCandidates");
  }

  if (this -> _pfClusterCollectionTag) {
    edm::Handle < std::vector< l1t::PFCluster > > pfClusterCollectionHandle;
    iEvent.getByToken(*(this -> _pfClusterCollectionTag), pfClusterCollectionHandle);
    // dumping the data
    //std::cout << ">>>>>> DUMPING PFCLUSTERS <<<<<<" << std::endl;
    //for (auto pfClusterIterator = pfClusterCollectionHandle -> begin(); pfClusterIterator != pfClusterCollectionHandle -> end(); pfClusterIterator++) 
    //{
    //  std::cout << pfClusterIterator -> pt() << "\t" << pfClusterIterator -> eta() << "\t" << pfClusterIterator -> phi() << "\t" << std::endl;
    //}
    
    #ifndef DEBUG
    this -> _caloGridPfCluster -> Reset();
    this -> _fillCaloGrid<>(*(this -> _caloGridPfCluster), *pfClusterCollectionHandle);
    #endif
    const auto seedsVector = this -> _findSeeds(*(this -> _caloGridPfCluster), this -> _seedPtThreshold); // seedPtThreshold = 6
    auto l1jetVector = this -> _buildJetsFromSeeds(*(this -> _caloGridPfCluster), seedsVector);
    
    std::unique_ptr< BXVector<l1t::Jet> > l1jetVectorPtr(new BXVector<l1t::Jet>(l1jetVector));
    iEvent.put(std::move(l1jetVectorPtr), "Phase1L1TJetFromPfClusters");
  }

  return;

}

void L1TJetPhase1Producer::_subtract9x9Pileup(const TH2F & caloGrid, l1t::Jet & l1tjet) {
  //For each jet we compute the 4 side bands
  //for (l1t::Jet & l1tjet: jetCollection) {
  //for (auto l1tjetIterator = jetCollection.begin(0); l1tjetIterator != jetCollection.end(0); l1tjetIterator++)
  // these variables host the total pt in each sideband and the total pileup contribution
  float topBandPt = 0;
  float leftBandPt = 0;
  float rightBandPt = 0;
  float bottomBandPt = 0;
  float pileUpEnergy;

  // hold the jet's x-y (and z, as I have to use it, even if 2D) location in the histo
  int xCenter, yCenter, zCenter;
  // Retrieving histo-coords for seed
  caloGrid.GetBinXYZ(caloGrid.FindFixBin(l1tjet.eta(), l1tjet.phi()), xCenter, yCenter, zCenter);
  
  // Computing pileup
  for (int x = -4; x <= 4; x++) {
    for (int y = 0; y < 3; y++) {
      // top band, I go up 5 squares to reach the bottom of the top band
      // +x scrolls from left to right, +y scrolls up
      topBandPt += this -> _getTowerEnergy(caloGrid, xCenter + x, yCenter + (5 + y));
      // left band, I go left 5 squares (-5) to reach the bottom of the top band
      // +x scrolls from bottom to top, +y scrolls left
      leftBandPt += this -> _getTowerEnergy(caloGrid, xCenter - (5 + y), yCenter + x);
      // right band, I go right 5 squares (+5) to reach the bottom of the top band
      // +x scrolls from bottom to top, +y scrolls right
      rightBandPt += this -> _getTowerEnergy(caloGrid, xCenter + (5 + y), yCenter + x);
      // right band, I go right 5 squares (+5) to reach the bottom of the top band
      // +x scrolls from bottom to top, +y scrolls right
      bottomBandPt += this -> _getTowerEnergy(caloGrid, xCenter + x, yCenter - (5 + y));
    }
  }
  // adding bands and removing the maximum band (equivalent to adding the three minimum bands)
  pileUpEnergy = topBandPt + leftBandPt + rightBandPt + bottomBandPt - std::max(topBandPt, std::max(leftBandPt, std::max(rightBandPt, bottomBandPt)));

  //preparing the new 4-momentum vector
  math::PtEtaPhiMLorentzVector ptVector;
  // removing pu contribution
  ptVector.SetPt(l1tjet.pt() - pileUpEnergy);
  ptVector.SetEta(l1tjet.eta());
  ptVector.SetPhi(l1tjet.phi());
  //updating the jet
  l1tjet.setP4(ptVector);
}

std::vector<std::tuple<int, int>> L1TJetPhase1Producer::_findSeeds(const TH2F & caloGrid, float seedThreshold) 
{
  int nBinsX = caloGrid.GetNbinsX();
  int nBinsY = caloGrid.GetNbinsY();

  std::vector<std::tuple<int, int>> seeds;

  int etaHalfSize = (int) this -> _jetIEtaSize/2;
  int phiHalfSize = (int) this -> _jetIPhiSize/2;

  // for each point of the grid check if it is a local maximum
  // to do so I take a point, and look if is greater than the points around it (in the 9x9 neighborhood)
  // to prevent mutual exclusion, I check greater or equal for points above and right to the one I am considering (including the top-left point)
  // to prevent mutual exclusion, I check greater for points below and left to the one I am considering (including the bottom-right point)

  for (int iPhi = 1; iPhi <= nBinsY; iPhi++)
  {
    for (int iEta = 1; iEta <= nBinsX; iEta++)
    {
      float centralPt = caloGrid.GetBinContent(iEta, iPhi);
      if (centralPt <= seedThreshold) continue;
      bool isLocalMaximum = true;

      // Scanning through the grid centered on the seed
      for (int etaIndex = -etaHalfSize; etaIndex <= etaHalfSize; etaIndex++)
      {
        for (int phiIndex = -phiHalfSize; phiIndex <= phiHalfSize; phiIndex++)
        {
          if ((etaIndex == 0) && (phiIndex == 0)) continue;
          if (phiIndex > 0) {
            // >=
            isLocalMaximum = ((isLocalMaximum) && (centralPt >= this -> _getTowerEnergy(caloGrid, iEta + etaIndex, iPhi + phiIndex)));
          } else if (phiIndex < 0) {
            // >
            isLocalMaximum = ((isLocalMaximum) && (centralPt > this -> _getTowerEnergy(caloGrid, iEta + etaIndex, iPhi + phiIndex)));
          } else if (etaIndex > 0) { //phiIndex = 0, left side is >, right is >=
            // >=
            isLocalMaximum = ((isLocalMaximum) && (centralPt >= this -> _getTowerEnergy(caloGrid, iEta + etaIndex, iPhi + phiIndex)));
          } else {
            // >
            isLocalMaximum = ((isLocalMaximum) && (centralPt > this -> _getTowerEnergy(caloGrid, iEta + etaIndex, iPhi + phiIndex)));
          }
        }
      }

      if (isLocalMaximum)
      {
        seeds.emplace_back(std::make_tuple(iEta, iPhi));
      }

    }
  }

  return seeds;
}

BXVector<l1t::Jet> L1TJetPhase1Producer::_buildJetsFromSeeds(const TH2F & caloGrid, const std::vector<std::tuple<int, int>> & seeds)
{

  // For each seed take a grid centered on the seed of the size specified by the user
  // Sum the pf in the grid, that will be the pt of the l1t jet. Eta and phi of the jet is taken from the seed.
  BXVector<l1t::Jet> jets;
  for (const auto& seed: seeds){
    int iEta = std::get<0>(seed);
    int iPhi = std::get<1>(seed);

    int etaHalfSize = (int) this -> _jetIEtaSize/2;
    int phiHalfSize = (int) this -> _jetIPhiSize/2;

    float ptSum = 0;
    // Scanning through the grid centered on the seed
    for (int etaIndex = -etaHalfSize; etaIndex <= etaHalfSize; etaIndex++)
    {
      for (int phiIndex = -phiHalfSize; phiIndex <= phiHalfSize; phiIndex++)
      {
        ptSum += this -> _getTowerEnergy(caloGrid, iEta + etaIndex, iPhi + phiIndex);
      }

    }

    // Creating a jet with eta phi centered on the seed and momentum equal to the sum of the pt of the components    
    //reco::Candidate::LorentzVector ptVector;
    math::PtEtaPhiMLorentzVector ptVector;
    ptVector.SetPt(ptSum);
    //ptVector.SetPtEtaPhiE(ptSum, caloGrid.GetXaxis() -> GetBinCenter(iEta), caloGrid.GetYaxis() -> GetBinCenter(iPhi), ptSum);
    ptVector.SetEta(caloGrid.GetXaxis() -> GetBinCenter(iEta));
    ptVector.SetPhi(caloGrid.GetYaxis() -> GetBinCenter(iPhi));
    l1t::Jet jet(ptVector);

    this -> _subtract9x9Pileup(caloGrid, jet);
    jets.push_back(0, jet);
  }

  return jets;
}

template <class TriggerPrimitive>
void L1TJetPhase1Producer::_fillCaloGrid(TH2F & caloGrid, const std::vector<TriggerPrimitive> & triggerPrimitives)
{
  //Filling the calo grid with the primitives
  for (auto primitiveIterator = triggerPrimitives.begin(); primitiveIterator != triggerPrimitives.end(); primitiveIterator++) 
  {
    caloGrid.Fill((float) primitiveIterator -> eta(), (float) primitiveIterator -> phi(), (float) primitiveIterator -> pt());
  }
  return;
}

void
L1TJetPhase1Producer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(L1TJetPhase1Producer);