#ifndef SPH1DBasic_H
#define SPH1DBasic_H
#include<vector>
#include<array>
#include<iostream>
template<class T,class E>
class SPH1DBasic{
 private:
 double tf;
 public: 
 SPH1DBasic(double tf){
   this->tf = tf;
 }
 ~SPH1DBasic(){}
 void SPH(T &particle,E &Eq){
     int index,it=0;
     double rho,dv,**q,sumVel,*Q,sumQ,auxtf=0,dt;
     T auxPar;
     q = new double *[particle.size()];
     Q = new double [particle.size()];
     for(int i=0;i<particle.size();i++){
       q[i] = new double [particle.size()];
     }
      do{
        Eq.ArvVonNeumann(q,particle);
        dt= Eq.dtc(particle);
        if(dt<1e-3){
           dt = 0.15/30.0;
        }
        Eq.MometumEq(particle,q,dt);
        Eq.Update(particle,dt);
        Eq.rhoSPh(particle);
        Eq.CQ(particle,q,Q);
        Eq.CA(Q,particle,dt);
        Eq.EoS(particle);
        auxtf+=dt;
      }while(auxtf<tf);
      for(int i=0;i<particle.size();i++){
       delete[] q[i];
     }
      delete [] q;
      delete [] Q;
 }
};
#endif
