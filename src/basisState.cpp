#include <random>
#include <iostream>
#include "basisState.h"

basisState generateRandomState(int nOrbitals){
  basisState randomState(nOrbitals);
  double p;
  std::random_device rng;
  int newOccupation;
  double const normalizer=static_cast<double>(rng.max());
  for(int iOrbital=0;iOrbital<nOrbitals;++iOrbital){
    p=rng()/normalizer;
    if(p<0.5){
      newOccupation=1;
    }
    else{
      newOccupation=0;
    }
    randomState.setOccupation(iOrbital,newOccupation);
  }
  return randomState;
}

void printState(basisState const &a){
  int const d=a.getBasisSize();
  for(int i=0;i<d;++i)
    std::cout<<a[i];
}
