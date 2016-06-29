#ifndef POD_INPUT
#define POD_INPUT

class parData{
 public:
 parData(double S, double dt, int A, double xi):S_(S),dt_(dt),A_(A),xi_(xi){}
  double getTimeStep() const{return dt_;}
  double getS() const{return S_;}
  int getA() const{return A_;}
  double getDamping() const{return xi_;}
  void setS(double S){S_=S;}
 private:
  double S_;
  double dt_;
  int A_;
  double xi_;
};

#endif
