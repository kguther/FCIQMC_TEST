#ifndef WALKER_CLASS
#define WALKER_CLASS

#include "basisState.h"

class walker{
 public:
 walker(basisState determinant, int sign):determinant_(determinant),sign_(sign){}
  int getSign() const{return sign_;}
  basisState getDeterminant() const{return determinant_;}
 private:
  basisState determinant_;
  int sign_;
};

inline bool operator==(walker const &a, walker const &b){
  return (a.getDeterminant()==b.getDeterminant()) && (a.getSign()==b.getSign());
}

inline bool operator<(walker const &a, walker const &b){
  if(compareDeterminants(a.getDeterminant(),b.getDeterminant())){
    return true;
  }
  if(a.getDeterminant()==b.getDeterminant()){
    if(a.getSign()==b.getSign())
      return false;
    if(a.getSign()==-1)
      return true;
  }
  return false;
}

inline bool checkAnnihilation(walker const &a, walker const &b){
  return a.getDeterminant()==b.getDeterminant() && a.getSign()!=b.getSign();
}

#endif
