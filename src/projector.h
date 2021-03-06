#ifndef PROJECTOR_CLASS
#define PROJECTOR_CLASS

#include "walker.h"
#include "hamiltonian.h"
#include "parData.h"

//for testing purpose
void printEnsemble(std::vector<walker> const &ensemble);

class projector{
 public:
  projector(hamiltonian const &H, std::vector<basisState> const &initialState, parData const &pars);
  void fullProjection(unsigned int numSteps, unsigned int targetNumber);
  void prStep(){spawn();death();annihilate();}
  void spawn();
  void death();
  void annihilate();
  void updateShift();
  void updateAveragedShift(unsigned int i);
  //for testing purpose
  void printDeterminants() const{printEnsemble(ensemble);}
 private:
  hamiltonian H;
  parData pars;
  std::vector<walker> ensemble;
  std::vector<walker> newWalkers;
  double averagedShift;
  int initiatorThreshold;
  int ensembleSizeBackup;
};



#endif
