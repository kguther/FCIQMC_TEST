#ifndef POD_INPUT
#define POD_INPUT

class parData{
 public:
  //The damping xi and the timescale A determine how flexible the shift is:
  //Large A and small xi make the shift more rigid and make it harder to adjust the walker number -> shift has to be averaged
 parData(double S, double dt, int A, double xi, int I=0):S_(S),dt_(dt),A_(A),xi_(xi),I_(I){}
  double getTimeStep() const{return dt_;}
  double getS() const{return S_;}
  int getA() const{return A_;}
  double getDamping() const{return xi_;}
  double getInitiatorThreshold() const{return I_;}
  void setS(double S){S_=S;}
 private:
  double S_;
  double dt_;
  int A_;
  double xi_;
  int I_;
};

#endif
