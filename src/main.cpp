//FCIQMC frontend

#include "basisState.h"
#include "walker.h"
#include "parData.h"
#include "hamiltonian.h"
#include "projector.h"
#include <iostream>

void printState(basisState const &a);
hamiltonian generateHubbard(int dim, double U, double t);

int main(int argc, char **argv){
  int const L{6};
  int const sysSize{2*L};
  parData testData(0.5,0.005,3,.4,2);
  double U{6};
  double t{1};
  std::vector<int> occupations;
  basisState testState=generateRandomState(sysSize,L);
  std::cout<<"Initial ";
  printState(testState);
  std::cout<<std::endl;
  hamiltonian Hubbard=generateHubbard(sysSize,U,t);
  std::vector<basisState> testStateVec;
  testStateVec.push_back(testState);
  projector testSys(Hubbard,testStateVec,testData);
  /*
  double p;
  basisState testComp=getRandomCoupledState(testState,p);
  std::cout<<"Initial ";
  printState(testComp);
  std::cout<<std::endl;
  std::cout<<Hubbard.getMatrixElement(testComp,testState)<<std::endl;
  */
  testSys.fullProjection(4000U,300U);
  std::cout<<"Initial ";
  printState(testState);
  std::cout<<std::endl;
  
}

hamiltonian generateHubbard(int dim, double U, double t){
  hamiltonian H(dim);
  for(int i=0;i<dim-1;i+=2){
    //Hubbard interaction
    H.setMatrixElement(i,i+1,i,i+1,-U);
  }
  for(int i=0;i<dim-2;++i){
    //hopping, excluding the pbc term
    H.setMatrixElement(i,i+2,-t);
    H.setMatrixElement(i+2,i,-t);
  }
  for(int sigma=0;sigma<2;++sigma){
    //boundary term for pbc
    H.setMatrixElement(dim-1-sigma,1-sigma,-t);
    H.setMatrixElement(1-sigma,dim-1-sigma,-t);
  }
  return H;
}
