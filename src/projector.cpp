#include "projector.h"
#include <random>
#include <math.h>
#include <iostream>
#include <algorithm>

double abv(double x){
  if(x>0)
    return x;
  return -x;
}

int getWalkerNumber(std::vector<walker> const &sortedEnsemble, int pos){
  int numI{1};
  for(int j=pos-1;j>=0;--j){
    if(!(sortedEnsemble[j]==sortedEnsemble[pos]))
      break;
    ++numI;
  }
  for(unsigned int j=pos+1;j<sortedEnsemble.size();++j){
    if(!(sortedEnsemble[j]==sortedEnsemble[pos]))
      break;
    ++numI;
  }
  return numI;
}

projector::projector(hamiltonian const &HIn, std::vector<basisState> const &initialStates, parData const &parsIn):
  H(HIn),
  pars(parsIn),
  averagedShift(0.0)
{
  for(unsigned int i=0;i<initialStates.size();++i){
    for(int j=0;j<pars.getInitiatorThreshold()+5;++j){
      ensemble.push_back(walker(initialStates[i],1));
    }
  }
}

void projector::spawn(){
  basisState target, source;
  double p,pSpawn;
  std::random_device rng;
  double const timeStep=pars.getTimeStep();
  double const normalizer=static_cast<double>(rng.max());
  double randNum, matrixElement;
  int numSpawns;
  int walkerNumber;
  int const initiatorThreshold=pars.getInitiatorThreshold();
  for(unsigned int i=0;i<ensemble.size();++i){
    walkerNumber=std::count(ensemble.begin(),ensemble.end(),ensemble[i]);
    if(walkerNumber>initiatorThreshold){
      numSpawns=0;
      source=ensemble[i].getDeterminant();
      target=getRandomCoupledState(source,p);
      matrixElement=H.getMatrixElement(source,target);
      pSpawn=timeStep/p*abv(matrixElement);
      randNum=rng()/normalizer;
      numSpawns+=floor(pSpawn);
      pSpawn-=floor(pSpawn);
      if(randNum<pSpawn){
	++numSpawns;
      }
      int newSign;
      if(matrixElement>0){
	newSign=ensemble[i].getSign();
      }
      else{
	newSign=-ensemble[i].getSign();
      }
      for(int j=0;j<numSpawns;++j){
	newWalkers.push_back(walker(target,newSign));
      }
    }
  }
}

void projector::death(){
  double pDeath, matrixElement, randNum;
  basisState source;
  double const timeStep=pars.getTimeStep();
  double const S=pars.getS();
  std::random_device rng;
  double const normalizer{static_cast<double>(rng.max())};
  for(unsigned int i=0;i<ensemble.size();++i){
    source=ensemble[i].getDeterminant();
    matrixElement=H.getMatrixElement(source,source);
    pDeath=timeStep*(matrixElement-S);
    randNum=rng()/normalizer;
    if(randNum<abv(pDeath)){
      if(pDeath<0){
        newWalkers.push_back(walker(ensemble[i].getDeterminant(),-ensemble[i].getSign()));

      }
      else{
	ensemble.erase(ensemble.begin()+i);
	--i;

      }
    }
  }
}

void projector::annihilate(){
  //First, annihilate newly spawned/cloned walkers with eachother
  for(unsigned int i=0;i<newWalkers.size();++i){
    for(unsigned int j=i+1;j<newWalkers.size();++j){
      if(checkAnnihilation(newWalkers[i],newWalkers[j])){
	newWalkers.erase(newWalkers.begin()+j);
	//j>i -> newWalkers.begin()+i is still the same element
	newWalkers.erase(newWalkers.begin()+i);
	--i;
	break;
      }
    }
  }
  //Then, compare the remaining walkers to the ensemble
  for(unsigned int i=0;i<newWalkers.size();++i){
    for(unsigned int j=0;j<ensemble.size();++j){
      if(checkAnnihilation(newWalkers[i],ensemble[j])){
	ensemble.erase(ensemble.begin()+j);
	newWalkers.erase(newWalkers.begin()+i);
	--i;
	break;
      }
    }
  }
  //finally, merge the newWalkers into the ensemble
  //ensemble.insert(ensemble.end(),newWalkers.begin(),newWalkers.end());
  std::vector<walker>::iterator destination;
  //to keep the ensemble sorted, each new walker is inserted at the position corresponding to its binary value
  for(unsigned int m=0;m<newWalkers.size();++m){
    walker cWalker(newWalkers[m]);
    destination=std::lower_bound(ensemble.begin(),ensemble.end(),cWalker);
    ensemble.insert(destination,cWalker);
  }
  //and clear the newWalkers
  newWalkers.clear();
}

void projector::updateShift(){
  double S{pars.getS()};
  double const dt{pars.getTimeStep()*pars.getA()};
  //  std::cout<<ensemble.size()<<" vs "<<ensembleSizeBackup<<std::endl;
  S-=pars.getDamping()/dt*log(ensemble.size()/static_cast<double>(ensembleSizeBackup));
  pars.setS(S);
  ensembleSizeBackup=ensemble.size();
}

void projector::updateAveragedShift(unsigned int i){
  averagedShift = (averagedShift*i + pars.getS()*pars.getA())/(i+pars.getA());
}

void projector::fullProjection(unsigned int numSteps, unsigned int targetNumber){
  bool constN{false};
  bool equilibrated{false};
  unsigned int iEquil{500U};
  unsigned int iConstN{0U};
  int adaptionCounter{0};
  for(unsigned int i=0;i<numSteps;++i){
    prStep();
    if(ensemble.size()>=targetNumber && !constN){
      constN=true;
      iConstN=i;
      ensembleSizeBackup=ensemble.size();
    }
    if(constN){
      ++adaptionCounter;
      if(!equilibrated and (i-iConstN)>iEquil){
    	  equilibrated=true;
    	  iConstN = i;
	averagedShift = pars.getS();
      }
      if(adaptionCounter%pars.getA()==0){
	updateShift();
	std::cout<<"New shift: "<<pars.getS()<<std::endl;
	std::cout<<"Number of walkers: "<<ensemble.size()<<std::endl;
	if(equilibrated){
	  updateAveragedShift(i-iConstN);
	}
	adaptionCounter=0;
      }
    }
  }
  std::cout<<"Energy estimate: "<<averagedShift<<std::endl;
}
    
void printEnsemble(std::vector<walker> const &ensemble){
  std::cout<<"Printing ensemble\n";
  for(std::vector<walker>::const_iterator it=ensemble.begin();it!=ensemble.end();++it){
    printState(it->getDeterminant());
  }
}

