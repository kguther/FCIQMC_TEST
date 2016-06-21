#include "projector.h"
#include <random>
#include <iostream>

double abv(double x){
  if(x>0)
    return x;
  return -x;
}

projector::projector(hamiltonian const &HIn, basisState const &initialState, parData const &parsIn):
  H(HIn),
  pars(parsIn)
{
  ensemble.push_back(walker(initialState,1));
}

void projector::spawn(){
  basisState target, source;
  double p,pSpawn;
  std::random_device rng;
  double const timeStep=pars.getTimeStep();
  double const normalizer=static_cast<double>(rng.max());
  double randNum, matrixElement;
  int numSpawns;
  for(unsigned int i=0;i<ensemble.size();++i){
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

      std::cout<<"Spawned new walker with sign "<<newSign<<" at ";
      printState(target);
      std::cout<<std::endl;

    }
  }
}

void projector::death(){
  double pDeath, matrixElement, randNum;
  basisState source;
  double const timeStep=pars.getTimeStep();
  double const S=pars.getS();
  std::random_device rng;
  double const normalizer=static_cast<double>(rng.max());
  for(unsigned int i=0;i<ensemble.size();++i){
    source=ensemble[i].getDeterminant();
    matrixElement=H.getMatrixElement(source,source);
    pDeath=timeStep*(matrixElement-S);
    randNum=rng()/normalizer;
    if(randNum<abv(pDeath)){
      if(pDeath<0){
        newWalkers.push_back(walker(ensemble[i].getDeterminant(),-ensemble[i].getSign()));
	
	std::cout<<"Cloned walker\n";

      }
      else{
	ensemble.erase(ensemble.begin()+i);
	--i;

	std::cout<<"Killed walker\n";

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
	newWalkers.erase(newWalkers.begin()+i);
	
	std::cout<<"Annihilated walkers\n";
	
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
	
	std::cout<<"Annihilated walkers\n";
	
	--i;
	break;
      }
    }
  }
  //finally, merge the newWalkers into the ensemble
  ensemble.insert(ensemble.end(),newWalkers.begin(),newWalkers.end());
  //and clear the newWalkers
  newWalkers.clear();
}
    
void printEnsemble(std::vector<walker> const &ensemble){
  std::cout<<"Printing ensemble\n";
  for(std::vector<walker>::const_iterator it=ensemble.begin();it!=ensemble.end();++it){
    printState(it->getDeterminant());
  }
}
