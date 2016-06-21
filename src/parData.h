#ifndef POD_INPUT
#define POD_INPUT

class parData{
 public:
 parData(double S, double dt):S_(S),dt_(dt){}
  double getTimeStep() const{return dt_;}
  double getS() const{return S_;}
 private:
  double S_;
  double dt_;
};

#endif
