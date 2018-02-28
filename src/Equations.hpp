#ifndef Equations_H
#define Equations_H
#define gamma 1.4
#include<iostream>
template<class T, class K>
class Equations{
private:
 K kernel;
 double h,alpha,epsilon,dt;
public:
  Equations(double h,double alpha, double epsilon,double dt){
    this->h = h;
    this->alpha = alpha;
    this->epsilon = epsilon;
    this->dt = dt;
  }
  double divvk(double vk, double xi, double xk, double vi){
   return (vi-vk)*kernel.dwij(xi,xk,h);
 }
 void ArvVonNeumann(double **q,T &particle){
     int i,j,k;
    double vij,xij,rhom,c;
    for(i=0;i<particle.size();i++){
       for(k=0;k<particle.size();k++){
         vij=particle[i].vel[0]-particle[k].vel[0];
         xij=particle[i].pos[0]-particle[k].pos[0];
         rhom=0.5*(particle[i].rho+particle[k].rho);
         c = sqrt(gamma*particle[i].P/particle[i].rho);
       if((vij*xij)<0.0){
         q[i][k]= -1.0*alpha*h*c*(1/rhom)*((vij*xij)/(pow(xij,2)+(epsilon*pow(h,2))));
       }else{
        q[i][k]=0.0;
       }
     }
    
   }
 }
 void MometumEq(T &particle,double **q){
   int i,k,j;
   double sum,aux[particle.size()];
   for(i=0;i<particle.size();i++){
     sum=0.0;
     for(k=0;k<particle.size();k++){
        sum+=particle[k].mk*((particle[k].P/pow(particle[k].rho,2))+(particle[i].P/pow(particle[i].rho,2))+q[i][k])*kernel.dwij(particle[i].pos[0],particle[k].pos[0],h);
     }
     aux[i]=particle[i].vel[0]-(sum*dt);
   }
  for(i=0;i<particle.size();i++){
     particle[i].vel[0]=aux[i];
   }

}
void Update(T &particle){
  for(int i=0; i<particle.size();i++){
      particle[i].pos[0]=particle[i].pos[0]+particle[i].vel[0]*dt;
  }
}
void rhoSPh(T &particle){
   double sum;
   int k,i,j;
   for(i=0;i<particle.size();i++){
     sum=0.0;
     for(k=0;k<particle.size();k++){
        sum+= kernel.wij(particle[i].pos[0],particle[k].pos[0],h);
     }
     particle[i].rho = particle[i].mk*sum;
   }
 } 
 void CQ(T &particle, double **q,double *Q){
    int  i,k,j;
    double sum;
    for(i=0;i<particle.size();i++){
      sum=0.0;
       for(k=0;k<particle.size();k++){
         sum+= ((q[i][k])*(particle[i].vel[0]-particle[k].vel[0])*kernel.dwij(particle[i].pos[0],particle[k].pos[0],h));
       }
       Q[i]=(particle[i].mk/2.0)*sum;
    }
 }
 void CA(double *Q, T &particle){
    int i;
     for(i=0;i<particle.size();i++){
         particle[i].A[0]=particle[i].A[0]*(1.0+((dt*(gamma-1)*particle[i].rho*Q[i])/particle[i].P));
     }
 } 
 void EoS(T &particle){
   int i;
   for(i=0;i<particle.size();i++){
         particle[i].P= particle[i].A[0]*pow(particle[i].rho,gamma);
   }
 }
 double Interpolate1DRho(double x,T &particle){
    double sum=0.0;
     for(int k=0;k<particle.size();k++){
         sum+= particle[k].mk*kernel.wij(x,particle[k].pos[0],h);
     }
   return sum;
 }
 double Interpolate1DVel(double x,T &particle){
    double sum=0.0;
     for(int k=0;k<particle.size();k++){
         sum+= particle[k].mk*(particle[k].vel[0]/particle[k].rho)*kernel.wij(x,particle[k].pos[0],h);
     }
   return sum;
 }
 double Interpolate1DPre(double x,T &particle){
    double sum=0.0;
     for(int k=0;k<particle.size();k++){
         sum+= particle[k].mk*(particle[k].P/particle[k].rho)*kernel.wij(x,particle[k].pos[0],h);
     }
   return sum;
 }
};
#endif
