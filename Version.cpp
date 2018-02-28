#include<iostream>
#include<string>
using namespace std;
int main(int argc,char **argv){
 string tipoDat =string(argv[1]);
 if(tipoDat=="--version"){
   cout<<"sph1D"<<" "<<"1.0"<<endl;
 }
return 0;
}
