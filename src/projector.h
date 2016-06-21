#ifndef PROJECTOR_CLASS
#define PROJECTOR_CLASS

#include "walker.h"
#include "hamiltonian.h"
#include "parData.h"

//for testing purpose
void printEnsemble(std::vector<walker> const &ensemble);

class projector{
 public:
  projector(hamiltonian const &H, basisState const &initialState, parData const &pars);
  void prStep(){spawn();death();annihilate();}
  void spawn();
  void death();
  void annihilate();
  
  //for testing purpose
  void printDeterminants() const{printEnsemble(ensemble);}
 private:
  hamiltonian H;
  parData pars;
  std::vector<walker> ensemble;
  std::vector<walker> newWalkers;
  std::vector<basisState> initiators;
};



#endif
