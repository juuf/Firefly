//
//  main.cpp
//  test
//
//  Created by Julian on 26/04/2019.
//  Copyright © 2019 Julian. All rights reserved.
//

#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>
#include <fstream>
//#include <iterator>
//#include <ostream>
#include<numeric>
using namespace std;

// function declaration
double fireflySimulation(int f, double samplerate);
double median(vector<double> medi);

// global intialisation
//int in = 0;
//double interactions [1][11]; //interactions[0] = {1,2,3,4,5,6,7,8,9,10,11};
std::vector<double> E;
std::vector<double> C;
vector<double> domega(2,0);
int fire = 0;

int simlength = 8e4;
int sr =1e3;

// initial phase (temp2 = rand(2,1);)
vector<double> temp2={0.191519450378892, 0.622108771039832};
double pas_phi = temp2[0];
vector<double> phi_ext;

// initial frequency (temp = 2*2.^((rand(2,1)-0.5)*2);)
vector<double> temp={1.834587198568937,2.970523429156172};
double omega = temp[0];
//omega_ext= temp[1];
//int f = 0;
//vector<int> inter;
int delta_fire=1;
vector<int> df;
int ii;

#define pi 3.14159265358979323846


int main() {
//    local intialisation
    double omega_ext = temp[1];
    int f = 0;
    vector<double> omegas;
    vector<double> phase;
    
    phi_ext.push_back(temp2[1]);
    E.push_back(1);
    C.push_back(1);
    
    phase.push_back(pas_phi);
    omegas.push_back(omega);
    
    for (int i=2-1; i<simlength; i++) {
        ii=i;
        if (phi_ext[phi_ext.size()-1]<1) {
            phi_ext.push_back(phi_ext[phi_ext.size()-1]+omega_ext/sr); //tooth saw -> increase (phase=amplitude)
            if (phi_ext[phi_ext.size()-1]>=1){
                phi_ext[phi_ext.size()-1]=1;
            }
        }
        else{
        f=(f+1)%2;
//              std::cout << f;
            phi_ext.push_back(0);
        }
        phase.push_back(fireflySimulation(f,sr));
        omegas.push_back(omega);
    }
    
    
    
//    Output to file
 //   std::ofstream ff1("abosulte path/phi_ext.txt");
   //  for(vector<double>::const_iterator i = phi_ext.begin(); i != phi_ext.end(); ++i) {
    //     ff1 << *i << '\n';}
    // std::ofstream ff2("abosulte path//phi.txt");
     //for(vector<double>::const_iterator i = phase.begin(); i != phase.end(); ++i) {
        // ff2 << *i << '\n';}
     //std::ofstream ff3("abosulte path/omegas.txt");
     //for(vector<double>::const_iterator i = omegas.begin(); i != omegas.end(); ++i) {
        // ff3 << *i << '\n';}

//    std::cout<<median({1,2,3,4});
    
    return 0;
}





double fireflySimulation(int f, double sr) {
    // local variable declaration
    const double Alpha = 0.4;
    const double Beta = 0.4;
    const int FilterLength = 8;
    const double tref = 5;
    const int DoublingThreshold = 20;
    double phi;
    double phiref;
    double allzeros;
    vector<double> CE;
    vector<double> EE;
//    vector<int> res_;
    
    delta_fire++;
    if (pas_phi < 1){
        phi = pas_phi+omega/sr; //tooth saw -> increase (phase=amplitude)
        if (phi>=1) {phi=1;}
    }
    else{
        if (fire==1){
//            inter.push_back(ii);
            df.push_back(delta_fire);
            delta_fire = 0;
        }
        fire = (fire+1)%2; //fire only on every other peak
        phi = 0;
    }
    
    phiref = tref*omega*5e-4;
    
    if (phi_ext[phi_ext.size()-1] >= 1){    // if external node is at maximum
        if (f == 0){
//            pas_.in=pas_.in+1;
//            pas_.interactions(pas_.in,1:11) = [-1 i j 0 0 0 0 0 0 0 0];
        }
        if (f == 1 && fire == 1){ // node j pas_.fires and firing has already occured at this time.
//        pas_.in=pas_.in+1;
//        pas_.interactions(pas_.in,1:3) = [0 i j];
        }
    
        if (f == 1 && fire == 0){// if node j pas_.fires and firing has not already occured at this time.
//            phi_ext[phi_ext.size()-1]  = 1;

            if (phi > phiref){ // if the phase of node k is above phiref
                
                // PHASE COUPLING FUNCTION
                phi = phi+Alpha*(-sin(phi*2*pi))*abs(sin(phi*2*pi));
                if (phi<0){phi=0;} // trim phase to [0 1] range
//                else if (phi>=1){phi=1;}
                if (phi>=1){phi=1;}
        
                if (phi < phiref || phi > 1-phiref){
                    E.push_back(0);
                }
                else{
                    E.push_back((1-cos(2*pi*phi))/2); //else, error is calculated like this
                }

    
                if (E.size() < FilterLength){
                    C.push_back(median(E));
                }
                else{
                    for(int ii =0;ii<(FilterLength);ii++){
                        CE.push_back(E[E.size()-FilterLength+ii]);
                    }
                    C.push_back(median(CE));
                }

            // FREQUENCY COUPLING FUNCTION
            domega[0] = ((omega * pow(2,Beta*-sin(2*pi*phi)*C[C.size()-1]))+domega[0]*domega[1])/(domega[1]+1);
            domega[1] = domega[1]+1;
            if (domega[1] > DoublingThreshold){
                omega = 2*omega;
                domega = {0,0};
            }
//            pas_.in=pas_.in+1;
//            pas_.interactions(pas_.in,1:9) = [1 i j k 1 pas_.C(end) phi-pas_.phi Beta*sin(2*pi*phi) pas_.domega(1,1)];
        }
//        else{
////            pas_.in=pas_.in+1;
////            pas_.interactions(pas_.in,1:6) = [1 i j k 0 pas_.C(end)]; %col 5 == 0 means pas_.fire within tref
//            }
        }
    
        if (domega[1] > 0){ // if there has been freq adjustments
            omega = domega[0];    // Reachback pas_.firefly algorithm
            domega = {0,0};
        }

    
        if (E.size() < FilterLength){
            allzeros = std::accumulate(E.begin(), E.end(), 0.0);
        }
        else{
            for(int ii =0;ii<(FilterLength);ii++){
                EE.push_back(E[E.size()-FilterLength+ii]);
            }
            allzeros = std::accumulate(EE.begin(), EE.end(), 0.0);
        }
        if (allzeros == 0 && (df.size() > 8)){
//            std::sort(inter.begin(), inter.end());
//            std::reverse(inter.begin(), inter.end());
//            vector<int> res_(inter.size());
//            std::adjacent_difference(inter.begin(), inter.end()-(inter.size()-FilterLength), res_.begin()); // difference of vector elements
//            double tt = std::accumulate(res_.begin()+1,res_.end(),0); // sum of difference to take the mean afterwards
//            omega = -2000*(FilterLength-1)/tt;
            // better implementation
            vector<int> tmp_;
            for(int ii=0;ii<(FilterLength-1);ii++){
                tmp_.push_back(df[df.size()-FilterLength+1+ii]);
            }
            double ttt = std::accumulate(tmp_.begin(),tmp_.end(),0);
            omega = 2000*(FilterLength-1)/ttt;
        }
    }
    if (phi>1){phi=1;}
    
    pas_phi = phi;
    
    return phi;
}




double median(vector<double> medi)
{
    sort(medi.begin(), medi.end());     // sort values
    
    double tmedian;
    if (medi.size() % 2 == 0)           // even
        tmedian = (medi[medi.size() / 2 - 1] + medi[medi.size() / 2]) / 2;
    else                                // odd
        tmedian = medi[medi.size() / 2];
    
    return tmedian;
}
