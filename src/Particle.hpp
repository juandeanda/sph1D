#ifndef Particle_H
#define Particle_H
#include <CompactNSearch>
#include<array>
using namespace std;
class Particle{
public:
array<CompactNSearch::Real, 3> pos {0, 0, 0};
array<double,3> vel {0,0,0};
array<double,3> A {0,0,0};
double rho =0;
double P=0;
double mk=0;
double dv=0;

};
#endif
