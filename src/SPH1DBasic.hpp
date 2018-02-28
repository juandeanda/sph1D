#ifndef SPH1DBasic_H
#define SPH1DBasic_H
#include<vector>
#include<array>
#include<iostream>
template<class T,class E>
class SPH1DBasic{
 private:
 double dt;
 public: 
 SPH1DBasic(double dt){
   this->dt = dt;
 }
 void SPH(T &particle,E &Eq){
     int index,it=0;
     double rho,dv,**q,sumVel,*Q,sumQ;
     T auxPar;
     q = new double *[particle.size()];
     Q = new double [particle.size()];
     for(int i=0;i<particle.size();i++){
       q[i] = new double [particle.size()];
     }
      do{
        Eq.ArvVonNeumann(q,particle);
        Eq.MometumEq(particle,q);
        Eq.Update(particle);
        Eq.rhoSPh(particle);
        Eq.CQ(particle,q,Q);
        Eq.CA(Q,particle);
        Eq.EoS(particle);
        it++;
      }while(it<30);
 }
};
#endif
