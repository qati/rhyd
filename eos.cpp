#include "eos.hpp"
#include <math.h>
#include <iostream>
#include <cstdlib>

using namespace std;

EoS::EoS(dtype* kappa) : kappa(kappa)
{
    dtype Tm;
    v_kappa.resize(2001);
    for(uint i=0;i<2001;i++){
        Tm = (dtype)i+100;
        p = 0.21;
        for(T=100;T<=Tm;T+=dT){
            t  = T/200;
            I  = exp(-h1/t-h2/t/t)*(h0+f0*(tanh(f1*t+f2)+1)/(1+g1*t+g2*t*t));
            p += I*dT/T;
        }
        t  = Tm/200;
        I  = exp(-h1/t-h2/t/t)*(h0+f0*(tanh(f1*t+f2)+1)/(1+g1*t+g2*t*t));
        v_kappa.at(i) = 3+I/p;
    }
}

void EoS::set_mode(bool mode)
{
    use_QCD = mode;
    return;
}


void EoS::ck(const dtype& Tm)
{
    //cout<<*kappa<<endl;
    if (!use_QCD) return;
    if (Tm<=100){
        *kappa = v_kappa.at(0);
        return;
    }
    if (Tm>=2100){
        *kappa = v_kappa.at(2000);
        return;
    }
    modf((double)(Tm), &x1);
    y1 = v_kappa.at((uint)x1);
    modf(Tm+dT, &x2);
    y2 = v_kappa.at((uint)x2);
    m = (y2-y1)/(x2-x1);
    *kappa = (dtype)(m*(Tm-x1)+y1);
    if (*kappa!=*kappa){
        cerr<<"EoS::ck kappa nan! " << *kappa <<endl; 
        exit(1);
    }
    if (*kappa<0.1 || *kappa>20){
        cerr<<"EoS::ck kappa=" <<*kappa<<endl;
        exit(1);
    }
    return;
}
