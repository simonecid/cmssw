// -*- C++ -*-
//
// Package:    JetCalibration/ApplyCalibrationFactors
// Class:      ApplyCalibrationFactors
// 
/**\class ApplyCalibrationFactors ApplyCalibrationFactors.cc JetCalibration/ApplyCalibrationFactors/plugins/ApplyCalibrationFactors.cc

 Description: [one line class summary]

 Implementation:
    [Notes on implementation]
*/
//
// Original Author:  Simone Bologna
//         Created:  Wed, 19 Dec 2018 12:44:23 GMT
//
//


// system include files
#include <memory>

#include <csignal>

#include <algorithm>
#include <cmath>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/JetReco/interface/CaloJet.h"


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"


//
// class declaration
//

class ApplyCalibrationFactors : public edm::stream::EDProducer<> {
  public:
    explicit ApplyCalibrationFactors(const edm::ParameterSet&);
    ~ApplyCalibrationFactors();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    virtual void beginStream(edm::StreamID) override;
    virtual void produce(edm::Event&, const edm::EventSetup&) override;
    virtual void endStream() override;

    //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

    edm::EDGetTokenT < std::vector < reco::CaloJet > > _inputCollectionTag;
    std::vector<double> _absEtaBinning;
    size_t _nBinsEta;
    std::vector<edm::ParameterSet> _calibration;
    std::string _outputCollectionName;

    std::vector<std::vector<double>> _jetCalibrationFactorsBinnedInEta;
    std::vector<std::vector<double>> _jetCalibrationFactorsPtBins;

    // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
ApplyCalibrationFactors::ApplyCalibrationFactors(const edm::ParameterSet& iConfig):
  _inputCollectionTag(edm::EDGetTokenT< std::vector < reco::CaloJet > >(consumes< std::vector < reco::CaloJet > > (iConfig.getParameter< edm::InputTag >("inputCollectionTag"))))
{
  this -> _absEtaBinning = iConfig.getParameter<std::vector<double> >("absEtaBinning");
  this -> _nBinsEta = this -> _absEtaBinning.size() - 1;
  this -> _calibration = iConfig.getParameter<std::vector<edm::ParameterSet> >("calibration");
  this -> _outputCollectionName = iConfig.getParameter<std::string>("outputCollectionName");

  for (const auto & pset: this -> _calibration)
  {
    this -> _jetCalibrationFactorsPtBins.emplace_back(pset.getParameter<std::vector<double> >("l1tPtBins"));
    this -> _jetCalibrationFactorsBinnedInEta.emplace_back(pset.getParameter<std::vector<double> >("l1tCalibrationFactors"));
  }

  produces<std::vector<reco::CaloJet> >( this->_outputCollectionName ).setBranchAlias(this->_outputCollectionName);
}


ApplyCalibrationFactors::~ApplyCalibrationFactors()
{
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
ApplyCalibrationFactors::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle < std::vector < reco::CaloJet > > inputCollectionHandle;
  iEvent.getByToken(this -> _inputCollectionTag, inputCollectionHandle);

  std::unique_ptr < std::vector < reco::CaloJet > > calibratedCollectionPtr (new std::vector < reco::CaloJet > () );

  // for each candidate:
  // 1 get pt and eta
  // 2 run a dicotomic search on the eta vector to find the eta index
  // 3 fetch the corresponding calibration elements
  // 4 run a dicotomic search on the pt binning to find the pt index
  // 5 fetch the calibration factor
  // 6 update the candidate pt by applying the calibration factor
  // 7 store calibrated candidate in a new collection
  for (const auto & candidate: *inputCollectionHandle )
  {
    // 1 
    float pt = candidate.pt();
    float eta = candidate.eta();

    //2
    auto etaBin = upper_bound(this -> _absEtaBinning.begin(), this -> _absEtaBinning.end(), fabs(eta));
    int etaIndex = etaBin - this -> _absEtaBinning.begin() - 1;

    //3
    const std::vector<double> & l1tPtBins = this -> _jetCalibrationFactorsPtBins[etaIndex];
    const std::vector<double> & l1tCalibrationFactors = this -> _jetCalibrationFactorsBinnedInEta[etaIndex];

    //4
    auto ptBin = upper_bound(l1tPtBins.begin(), l1tPtBins.end(), pt);
    int ptBinIndex = ptBin - l1tPtBins.begin() - 1;

    //5
    const double & l1tCalibrationFactor = l1tCalibrationFactors[ptBinIndex];

    //6
    reco::Candidate::PolarLorentzVector candidateP4(candidate.polarP4());
    reco::CaloJet * newCandidate = candidate.clone();
    candidateP4.SetPt(candidateP4.pt() * l1tCalibrationFactor);
    newCandidate -> setP4(candidateP4);
    
    //7
    calibratedCollectionPtr -> emplace_back(*newCandidate);
    // clean up

    #ifdef DEBUG
    if (newCandidate -> pt() < 0) {
      std::cout << "PRE-CALIBRATION " << std::endl;
      std::cout << "\t Jet properties (pt, eta, phi, pile-up): " << candidate.pt() << "\t" << candidate.eta() << "\t" << candidate.phi() << "\t" << candidate.pileup() << std::endl;
      std::cout << "CALIBRATION " << std::endl;
      std::cout << "\t Using eta - pt - factor " << *etaBin << " - " << *ptBin << " - " << l1tCalibrationFactor << std::endl;
      std::cout << "POST-CALIBRATION " << std::endl;
      std::cout << "\t Jet properties (pt, eta, phi, pile-up): " << newCandidate -> pt() << "\t" << newCandidate -> eta() << "\t" << newCandidate -> phi() << "\t" << newCandidate -> pileup() << std::endl;
    }
    #endif

    delete newCandidate;
  }

  iEvent.put(std::move(calibratedCollectionPtr), this -> _outputCollectionName);

}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
ApplyCalibrationFactors::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
ApplyCalibrationFactors::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
ApplyCalibrationFactors::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
ApplyCalibrationFactors::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
ApplyCalibrationFactors::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
ApplyCalibrationFactors::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ApplyCalibrationFactors::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ApplyCalibrationFactors);
