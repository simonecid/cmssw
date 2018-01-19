
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include <tuple>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"

namespace MatchingAlgorithms{
  template <class TParticle, class TTrigger> // <3
  std::tuple<const TTrigger*, const TParticle*, float, int>
  matchParticleWithL1Object
  (
    const TParticle &particle,
    const edm::Handle<std::vector<TParticle> > &particleCollectionHandle,
    const edm::Handle<BXVector<TTrigger> > &l1tObjectCollectionHandle,
    float dr2Max,
    bool performCrossCheck
  );

  template <class TParticle, class TTrigger> // <3
  std::tuple<const TTrigger *, const TParticle *, float, int>
  matchL1ObjectWithParticle
  (
    const TTrigger &l1tObject,
    const edm::Handle<BXVector<TTrigger>> &l1tObjectCollectionHandle,
    const edm::Handle<std::vector<TParticle>> &particleCollectionHandle,
    float dr2Max,
    bool performCrossCheck
  );
}


template <class TParticle, class TTrigger> // <3
std::tuple<const TTrigger *,const TParticle *, float, int>
MatchingAlgorithms::matchParticleWithL1Object(
    const TParticle &particle,
    const edm::Handle<std::vector<TParticle>> &particleCollectionHandle,
    const edm::Handle<BXVector<TTrigger>> &l1tObjectCollectionHandle,
    float dr2Max,
    bool performCrossCheck)
{

  bool foundPreliminaryMatch = false;
  bool foundMatch = false;
  std::tuple<const TTrigger *,const TParticle *, float, int> l1tObjectGenJetPair;

  const TTrigger *bestMatch = NULL;
  std::get<1>(l1tObjectGenJetPair) = &particle;
  std::get<0>(l1tObjectGenJetPair) = bestMatch;

  float currentDR2Best = dr2Max;
  float currentDR2Min = 0;

  int matchQuality = 1;

  do
  {
    // do not delete those 2 lines
    foundPreliminaryMatch = false;
    foundMatch = false;

    // particle to object match
    for (
        typename BXVector<TTrigger>::const_iterator bx0Iterator = l1tObjectCollectionHandle->begin(0);
        bx0Iterator != l1tObjectCollectionHandle->end(0);
        bx0Iterator++)
    {
      float dr2 = reco::deltaR2(*bx0Iterator, particle);

      if ((dr2 > currentDR2Min) && (dr2 < currentDR2Best))
      {

        bestMatch = &(*bx0Iterator);
        std::get<0>(l1tObjectGenJetPair) = bestMatch;
        std::get<2>(l1tObjectGenJetPair) = dr2;
        std::get<3>(l1tObjectGenJetPair) = matchQuality;
        currentDR2Best = dr2;
        foundPreliminaryMatch = true;
      }
    }

    // if asked check the cross match if we have found a one-direction match
    if (foundPreliminaryMatch)
    {
      std::cout << "################## " << std::endl;
      std::cout << "Preliminary match: " << std::endl;
      std::cout << "Gen jet pt : " << particle.pt() << std::endl;
      std::cout << "l1tJet pt : " << bestMatch -> pt() << std::endl;
      std::cout << "deltaR2 : " << currentDR2Best << std::endl;
      if (performCrossCheck)
      // performing the cross check
      {
        // inverse matching
        const std::tuple<const TTrigger*, const TParticle*, float, int> l1tObjectGenJetPair2 =
            matchL1ObjectWithParticle<>(
                *bestMatch,
                l1tObjectCollectionHandle,
                particleCollectionHandle,
                dr2Max,
                false);

        // if inverse matching return the initial particle we are good
        if (std::get<1>(l1tObjectGenJetPair2) == &particle)
        {
          foundMatch = true;
          std::cout << "Preliminary match's cross-check at quality " << matchQuality << " was successful." << std::endl;
        }
        // if not we want to check if there is a secondary particle
        else
        {
          foundMatch = false;
          currentDR2Min = currentDR2Best;
          currentDR2Best = dr2Max;
          std::cout << "Preliminary match's cross-check at quality " << matchQuality << " failed, trying again." << std::endl;
          matchQuality++;
          std::get<0>(l1tObjectGenJetPair) = NULL;
        }
      }
      else
      {
        foundMatch = true;
      }
    }

  } while ((performCrossCheck) && (foundPreliminaryMatch) && (!foundMatch));

  return l1tObjectGenJetPair;
}

template <class TParticle, class TTrigger> // <3
std::tuple<const TTrigger *,const  TParticle *, float, int>
MatchingAlgorithms::matchL1ObjectWithParticle(
    const TTrigger &l1tObject,
    const edm::Handle<BXVector<TTrigger>> &l1tObjectCollectionHandle,
    const edm::Handle<std::vector<TParticle>> &particleCollectionHandle,
    float dr2Max,
    bool performCrossCheck)
{

  bool foundPreliminaryMatch = false;
  bool foundMatch = false;
  std::tuple<const TTrigger *, const TParticle *, float, int> l1tObjectGenJetPair;

  const TParticle *bestMatch = NULL;
  std::get<0>(l1tObjectGenJetPair) = &l1tObject;
  std::get<1>(l1tObjectGenJetPair) = bestMatch;

  float currentDR2Best = dr2Max;
  float currentDR2Min = 0;

  int matchQuality = 1;

  do
  {
    // do not delete those 2 lines
    foundPreliminaryMatch = false;
    foundMatch = false;

    // object to particle match
    for (
        typename std::vector<TParticle>::const_iterator bx0Iterator = particleCollectionHandle->begin();
        bx0Iterator != particleCollectionHandle->end();
        bx0Iterator++)
    {
      float dr2 = reco::deltaR2(*bx0Iterator, l1tObject);

      if ((dr2 > currentDR2Min) && (dr2 < currentDR2Best))
      {

        bestMatch = &(*bx0Iterator);
        std::get<1>(l1tObjectGenJetPair) = bestMatch;
        std::get<2>(l1tObjectGenJetPair) = dr2;
        std::get<3>(l1tObjectGenJetPair) = matchQuality;
        currentDR2Best = dr2;
        foundPreliminaryMatch = true;
      }
    }

    // if asked check the cross match if we have found a one-direction match
    if (foundPreliminaryMatch)
    {
      if (performCrossCheck)
      // performing the cross check
      {
        // inverse matching
        const std::tuple<const TTrigger*, const TParticle*, float, int> l1tObjectGenJetPair2 =
            matchParticleWithL1Object<>(
                *bestMatch,
                particleCollectionHandle,
                l1tObjectCollectionHandle,
                dr2Max,
                false);

        // if inverse matchin retgrun the initial object we are good
        if (std::get<0>(l1tObjectGenJetPair2) == &l1tObject)
        {
          foundMatch = true;
        }
        // if not we want to check if there is a secondary particle
        else
        {
          foundMatch = false;
          currentDR2Min = currentDR2Best;
          currentDR2Best = dr2Max;
          matchQuality++;
          std::get<1>(l1tObjectGenJetPair) = NULL;
        }
      }
      else
      {
        foundMatch = true;
      }
    }

  } while ((performCrossCheck) && (foundPreliminaryMatch) && (!foundMatch));

  return l1tObjectGenJetPair;
}