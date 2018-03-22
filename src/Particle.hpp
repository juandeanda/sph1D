#ifndef Particle_H
#define Particle_H
#include<vector>
using namespace std;
class Particle{
public:
vector<double> pos {0, 0, 0};
vector<double> vel {0,0,0};
vector<double> A {0,0,0};
double rho =0;
double P=0;
double mk=0;
double dv=0;
double c=0;

};
#endif
