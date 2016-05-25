#include "solverwa.hpp"
#include <cassert>
#include <iostream>
#include <cstdlib>
#include <cmath>


using namespace std;


void fw_calc(void * context, const uint& i, const uint& j)
{
    static_cast<solverwa*>(context)->analitic(i, j);
    static_cast<solverwa*>(context)->error(i, j);

    return;
}


solverwa::solverwa(io * FS) : solver(FS)
{
    fptr.push_back(&fw_calc);

    a.resize(p.row(), p.col());

    rd.resize(4);
    ad.resize(4);

    vector<vector<dtype>*> tmp;
    tmp.push_back(&ad);
    tmp.push_back(&rd);

    fs->add_hyd(&a);
    fs->add_vec(tmp);

    kappa = fs->get("kappa");
    iv.n0 = fs->get("n0");
    iv.p0 = fs->get("p0");

    X       = fs->get("X0");
    Y       = fs->get("Y0");
    V       = fs->get("V0");
    W       = fs->get("W0");


    p_s_corr = fs->get("p_s_corr");

}

solverwa::~solverwa()
{
}


void solverwa::analitic(const uint& i, const uint& j)
{
    cout<<"Calling virtual function testmod::analitic("<<i<<", "<<j<<")!"<<endl;
    assert(0);
    return;
}


void solverwa::error(const uint& i, const uint& j)
{
    for(ik=0;ik<4;ik++){
        rd.at(ik)     += fabs(p.at(ik, i, j) - a.at(ik, i, j));
        ad.at(ik)     += fabs(p.at(ik, i, j) - a.at(ik, i, j));
    }

    return;
}


void solverwa::start_step()
{
    cout<<"Calling virtual function testmod::start_step!"<<endl;
    assert(0);
}


void solverwa::end_step()
{
    cout<<"Calling virtual function testmod::end_step!"<<endl;
    assert(0);
}


void solverwa::init()
{
    uint j, k;

    start_step();

    p.at(0, i0, j0) = iv.n0;
    p.at(3, i0, j0) = iv.p0;
    for(uint i=1;i<=a.row();i++){
        for(j=1;j<=a.col();j++){
            analitic(i, j);
            for(k=0;k<4;k++) p.at(k, i, j) = a.at(k, i, j);
            calc_integral(i, j);
        }
    }
    calc_eps();

    nsolver.pq();

    fs->write(t);

    return;
}




/**
 * SOLPCONST MODULE.
 */
solpconst::solpconst(io * FS) : solverwa(FS)
{
    X   = V*t;
    Y   = W*t;
    t0  = t;
}


void solpconst::start_step()
{
    X = V*t;
    Y = W*t;

    for(ik=0;ik<4;ik++){
        ad.at(ik) = 0;
        rd.at(ik) = 0;
    }

    return;
}


void solpconst::end_step()
{
    for(ik=0;ik<4;ik++){
        rd.at(ik) /= a.row()*a.col();
        ad.at(ik) /= a.row()*a.col();
    }

    return;
}


void solpconst::analitic(const uint& i, const uint& j)
{
    x = st.x1 + (i-1)*st.dx;
    y = st.y1 + (j-1)*st.dy;

    tau0 = t0;
    tau  = sqrt(t*t-x*x-y*y);

    a.at(0, i, j) = iv.n0 * pow(tau0/tau, 2) * exp(-x*x/X/X/2-y*y/Y/Y/2);
    a.at(1, i, j) = x * V / X;
    a.at(2, i, j) = y * W / Y;
    a.at(3, i, j) = iv.p0 * pow(tau0/tau, 2+2/kappa);

    return;
}


/**
 * ACCSOL MODULE.
 */
accsol::accsol(io * FS) : solverwa(FS)
{
    t0 = t;
    if (FS->get("kappa")!=2.33333){
        cerr << "ACCSOL fatal error! Wrong kappa in opt file! Please update to 7/3!" << endl;
        assert(0);
    }
}


void accsol::start_step()
{
    for(ik=0;ik<4;ik++){
        ad.at(ik) = 0;
        rd.at(ik) = 0;
    }

    return;
}


void accsol::end_step()
{
    for(ik=0;ik<4;ik++){
        rd.at(ik) /= a.row()*a.col();
        ad.at(ik) /= a.row()*a.col();
    }

    return;
}


void accsol::analitic(const uint& i, const uint& j)
{
    x = st.x1 + (i-1)*st.dx;
    y = st.y1 + (j-1)*st.dy;

    tau  = sqrt(t*t-x*x-y*y);
    tau0  = sqrt(t0*t0-x*x-y*y);
    eta = atanh(sqrt(x*x+y*y)/t);

    /*if (isinf((double)tau) || isinf((double)eta) || isinf((double)tau0) || std::isnan((double)tau) || std::isnan((double)eta) || std::isnan((double)tau0)){
        cerr << "ACCSOL::ANALITIC inf/nan: (i, j)=(" << i << ", " << j << "); (x, y)=(" << x << ", " << y << "); tau="
            << tau<< "; tau0=" << tau0 << "; eta=" << eta << endl << endl;
        asert(0);
    }*/

    //cout << x << " " << y << " ; tau=" << tau <<" ; tau0=" << tau0 << "; eta= " <<eta<<endl;

    a.at(0, i, j) = iv.n0 * pow(tau0/tau, 3) * pow(cosh(eta/2), -1);
    a.at(1, i, j) = tanh(3*eta/2)*x/sqrt(x*x+y*y); //tanh/sqrt2
    a.at(2, i, j) = tanh(3*eta/2)*y/sqrt(x*x+y*y);
    a.at(3, i, j) = iv.p0 * pow(tau0/tau, 30/7) * pow(cosh(eta/2), -10/7);
    if (fabs(x)<DEF_NULL && fabs(y)<DEF_NULL) a.at(1,i,j)=a.at(2,i,j)=0;

    /*if (std::isnan((double)a.at(0, i, j)) || isinf((double)a.at(0, i, j))){
        cerr << "ACCSOL::ANALITIC inf/nan N: (i, j)=(" << i << ", " << j << "); (x, y)=(" << x << ", " << y << "); tau="
            << tau<< "; tau0=" << tau0 << "; eta=" << eta << ";; a.at(0, i, j) = " << a.at(0, i, j) << endl << endl;
        assert(0);
    }
    if (std::isnan((double)a.at(1, i, j)) || isinf((double)a.at(1, i, j))){
        cerr << "ACCSOL::ANALITIC inf/nan vx: (i, j)=(" << i << ", " << j << "); (x, y)=(" << x << ", " << y << "); tau="
            << tau<< "; tau0=" << tau0 << "; eta=" << eta << ";; a.at(1, i, j) = " << a.at(1, i, j) << endl << endl;
        assert(0);
    }
    if (std::isnan((double)a.at(2, i, j)) || isinf((double)a.at(2, i, j))){
        cerr << "ACCSOL::ANALITIC inf/nan vy: (i, j)=(" << i << ", " << j << "); (x, y)=(" << x << ", " << y << "); tau="
            << tau<< "; tau0=" << tau0 << "; eta=" << eta << ";; a.at(2, i, j) = " << a.at(2, i, j) << endl << endl;
        assert(0);
    }
    if (std::isnan((double)a.at(3, i, j)) || isinf((double)a.at(3, i, j))){
        cerr << "ACCSOL::ANALITIC inf/nan P: (i, j)=(" << i << ", " << j << "); (x, y)=(" << x << ", " << y << "); tau="
            << tau<< "; tau0=" << tau0 << "; eta=" << eta << ";; a.at(3, i, j) = " << a.at(3,  i, j) << endl << endl;
        assert(0);
    }*/
    return;
}
