#ifndef BASIS_STATE_CLASS
#define BASIS_STATE_CLASS

#include <vector>

class basisState{
 public:
 basisState():nOrbitals_(0){}
 explicit basisState(int nOrbitals):nOrbitals_(nOrbitals),occupations_(std::vector<int>(nOrbitals)){}
 explicit basisState(std::vector<int> occupations):nOrbitals_(occupations.size()),occupations_(occupations){}
  void setOccupation(int iOrbital, int newOccupation){occupations_[iOrbital]=newOccupation;}
  void addParticle(int iOrbital){++(occupations_[iOrbital]);}
  void removeParticle(int iOrbital){--(occupations_[iOrbital]);}
  int getOccupation(int iOrbital) const{return occupations_[iOrbital];}
  const int& operator[](int iOrbital) const{return occupations_[iOrbital];}
  int getBasisSize() const{return nOrbitals_;}
 private:
  int nOrbitals_;
  std::vector<int> occupations_;
};

basisState generateRandomState(int nOrbitals);
inline bool operator==(basisState const &a, basisState const &b){
  if(a.getBasisSize()!=b.getBasisSize()){
    return false;
  }
  for(int i=0;i<a.getBasisSize();++i){
    if(a.getOccupation(i)!=b.getOccupation(i)){
      return false;
    }
  }
  return true;
}

void printState(basisState const &a);

#endif
