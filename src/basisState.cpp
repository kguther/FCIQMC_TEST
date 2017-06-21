#include <random>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "basisState.h"

basisState generateRandomState(int nOrbitals){
	//	Generate some random state with arbitrary number of electrons
  basisState randomState(nOrbitals);
  double p;
  std::random_device rng;
  int newOccupation;
  double const normalizer{static_cast<double>(rng.max())};
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

basisState generateRandomState(int nOrbitals, int nElectrons){
	// Generate some random state with fixed number of electrons
	basisState randomState(nOrbitals);
	int p;
	std::random_device rng;
	std::vector<int> buffer(nOrbitals);
	for(int i=0;i<nOrbitals;++i){
		buffer[i]=i;
		randomState.setOccupation(i,0);
	}
	int normalizer{rng.max()/nOrbitals};
	for(int iElectron=0;iElectron<nElectrons;++iElectron){
		p=rng()/normalizer;
		// Just in case...
		if(p>nOrbitals-iElectron) p=nOrbitals-iElectron;
		randomState.setOccupation(buffer[p],1);
		buffer.erase(buffer.begin()+p);
		normalizer=rng.max()/(nOrbitals-iElectron-1);
	}
	return randomState;
}

bool compareDeterminants(basisState const &a, basisState const &b){
  for(int i=0;i<a.getBasisSize();++i){
    if(a.getOccupation(i)>b.getOccupation(i))
      return false;
    if(a.getOccupation(i)<b.getOccupation(i))
      return true;
  }
  return false;
}

void printState(basisState const &a){
  int const d=a.getBasisSize();
  for(int i=0;i<d;++i)
    std::cout<<a[i];
}
