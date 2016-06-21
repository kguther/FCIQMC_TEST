#ifndef HAMILTONIAN_CLASS
#define HAMILTONIAN_CLASS

#include "basisState.h"

class hamiltonian{
 public:
  hamiltonian(){}
  //dimension is the dimension of the single-particle hilbert space
  explicit hamiltonian(int dimension):d(dimension),oneBodyEntries(std::vector<double>(d*d,0.0)),twoBodyEntries(std::vector<double>(d*d*d*d,0.0)){}
  void setMatrixElement(int r, int s, double newEntry);
  void setMatrixElement(int p, int q, int r, int s, double newEntry);
  double getMatrixElement(basisState const &alpha, basisState const &beta);
  void printMatrix(int N);
 private:
  int d;
  //these are the coefficients of the second quantized hamiltonian
  std::vector<double> oneBodyEntries;
  std::vector<double> twoBodyEntries;
};

basisState getRandomCoupledState(basisState const &source, double &p);
int getFermiSign(basisState const &alpha, int annihilatorIndex, int creatorIndex);

#endif
