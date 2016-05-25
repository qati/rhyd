#ifndef _EOS_HPP_
#define _EOS_HPP_


#include "matrix.hpp"

#define EOS_CONST_KAPPA 0
#define EOS_QCD_KAPPA 1

const dtype  h0 = 0.1396;
const dtype  h1 = -0.18;
const dtype  h2 = 0.035;
const dtype  f0 = 2.76;
const dtype  f1 = 6.79;
const dtype  f2 = -5.29;
const dtype  g1 = -0.47;
const dtype  g2 = 1.04;

const dtype dT  = 1;

class EoS{
private:
    dtype * kappa;
    dtype t, I, T, p;
    std::vector<dtype> v_kappa;
    uint index;
    double x1,x2,y1,y2, m;
    bool use_QCD;
public:
    EoS(dtype*);
    void ck(const dtype&);
    void set_mode(bool);
};





#endif
