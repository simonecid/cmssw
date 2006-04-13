#ifndef HcalPedestals_h
#define HcalPedestals_h

/** 
\class HcalPedestals
\author Fedor Ratnikov (UMd)
POOL container to store Pedestal values 4xCapId
$Author: ratnikov
$Date: 2005/12/15 23:38:04 $
$Revision: 1.6 $
*/

#include <vector>
#include <algorithm>

#include "CondFormats/HcalObjects/interface/HcalPedestal.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"

// 
class HcalPedestals {
 public:
  HcalPedestals();
  ~HcalPedestals();
  /// get array of values for 4 capIds
  const HcalPedestal* getValues (HcalDetId fId) const;
  /// get value for given capId = 0..3
  float getValue (HcalDetId fId, int fCapId) const;
  /// get list of all available channels
  std::vector<HcalDetId> getAllChannels () const;
  /// check if data are sorted
  bool sorted () const {return mSorted;}
  /// fill values
  bool addValue (HcalDetId fId, const float fValues [4]);
  /// fill values
  bool addValue (HcalDetId fId, float fValue0, float fValue1, float fValue2, float fValue3);
  /// sort values by channelId  
  void sort ();
  // helper typedefs
  typedef HcalPedestal Item;
  typedef std::vector <Item> Container;
 private:
  Container mItems;
  bool mSorted;
};

#endif
