#include "rhyd.hpp"
#include <cstdlib>
#include <iostream>
#include <stdexcept>


using namespace std;


rhyd::rhyd()
{
    p     = NULL;
    kappa = 0;

    eos = new EoS(&kappa);
}

rhyd::rhyd(const dtype& DT, const dtype& DX, const dtype& DY, hyd* h)
{
    set(DT, DX, DY, h);

    kappa = 0;

    eos = new EoS(&kappa);
}


void rhyd::set(const dtype& DT, const dtype& DX, const dtype& DY, hyd* h)
{
    p  = h;

    dt = DT;
    dx = DX;
    dy = DY;

    kappa = 0;

    r = p->row();
    c = p->col();

    q.resize(r, c);
    tp.resize(4);

    return;
}

void rhyd::set(const dtype& KAPPA)
{
    if (KAPPA <= 0) eos->set_mode(EOS_QCD_KAPPA);
    else eos->set_mode(EOS_CONST_KAPPA);
    kappa = KAPPA;
    return;
}


void rhyd::pq()
{
    uint i, j;
    for(i=1;i<=r;i++){
        for(j=1;j<=c;j++){
            T = p->at(3, i, j)/p->at(0, i, j);
            eos->ck(T);
            u0 = sqrt(1+pow(p->at(1, i, j),2)+pow(p->at(2, i, j), 2));
            q.at(0, i, j) = p->at(0, i, j) * u0;
            q.at(1, i, j) = (kappa+1)*p->at(3, i, j)*u0*u0-p->at(3, i, j);
            q.at(2, i, j) = (kappa+1)*p->at(3, i, j)*u0*p->at(1, i, j);
            q.at(3, i, j) = (kappa+1)*p->at(3, i, j)*u0*p->at(2, i, j);
        }
    }
    return;
}


void rhyd::qtp(const dtype& q0, const dtype& q1, const dtype& q2, const dtype& q3)
{
    try {
        if (q0!=q0) throw logic_error("0");

        t0  = pow(q2, 2) + pow(q3, 2);
        t1  = pow(1+kappa, 2)*pow(q1, 2) - 4*kappa*t0;

        if (t0<DEF_NULL || fabs(t1)<DEF_NULL){
            tp.at(0) = q0;
            tp.at(1) = 0;
            tp.at(2) = 0;
            tp.at(3) = q1/kappa;
            return;
        }
        if (t1<0) t1 = 0; //throw logic_error("1");
        t1  = sqrt(t1);
        t2  = -(1+kappa)*pow(q1, 2) + 2*kappa*t0;
        t3  = (t2 + q1*t1);
        t3 /= 2*(1+kappa)*t0*(pow(q1, 2)-t0);

        if (t3<0) t3 = 0; //throw logic_error("2");

        /** Calculating p from q! **/
        tp.at(0) = (1+kappa)*q0*((1-kappa)*q1+t1)/2/kappa * sqrt(t3);
        tp.at(1) = q2*sqrt(t3);
        tp.at(2) = sqrt(t3) * q3 ;
        tp.at(3) = ((1-kappa)*q1+t1)/2/kappa;

        if (tp.at(0)<DEF_NULL || tp.at(3)<DEF_NULL) T = 0;
        else T = tp.at(3)/tp.at(0);
        //if (T>100000){cout<< T<<" " <<tp.at(0)<<" " << tp.at(3)<<endl;}
        eos->ck(T);

    } catch(logic_error& e){
        string w(e.what());
        if (w=="0") non++;
        else if (w=="1"){
            cerr << endl << "======" << endl << "rhyd::qtp: t1 is less than 0 (" << t1 << ")!" << endl;
            exit(1);
        } else if (w=="2"){
            cerr << endl << "======" << endl << "rhyd::qtp: t3 is less than 0 (" << t3 << ")!" << endl;
            exit(1);
        } else {
            cerr << endl << "======" << endl << "rhyd::qtp: Unexpected error: " << w << "!"<<endl;
            exit(1);
        }
    }
    return;
}


void rhyd::qp(const uint& i, const uint& j)
{
    qtp(q.at(0, i, j), q.at(1, i, j), q.at(2, i, j), q.at(3, i, j));
    p->at(0, i, j) = tp.at(0);
    p->at(1, i, j) = tp.at(1);
    p->at(2, i, j) = tp.at(2);
    p->at(3, i, j) = tp.at(3);
    return;
}

void rhyd::flux_f(vector<dtype>& f, const dtype& q0, const dtype& q1, const dtype& q2, const dtype& q3)
{
    qtp(q0, q1, q2, q3);

    T = tp.at(3)/tp.at(0);
    eos->ck(T);

    u0 = sqrt(1+pow(tp.at(1), 2)+pow(tp.at(2), 2));

    f.at(0) = tp.at(0)*tp.at(1);
    f.at(1) = (1+kappa)*(tp.at(3)-tp.at(0))*tp.at(1)*u0;
    f.at(2) = (1+kappa)*(tp.at(3)-tp.at(0))*tp.at(1)*tp.at(1)+tp.at(3);
    f.at(3) = (1+kappa)*(tp.at(3)-tp.at(0))*tp.at(1)*tp.at(2);

    return;
}


dtype rhyd::get_T(){return T;}

void rhyd::flux_g(vector<dtype>& f, const dtype& q0, const dtype& q1, const dtype& q2, const dtype& q3)
{
    qtp(q0, q1, q2, q3);

    T = tp.at(3)/tp.at(0);
    eos->ck(T);

    u0 = sqrt(1+pow(tp.at(1), 2)+pow(tp.at(2), 2));

    f.at(0) = tp.at(0)*tp.at(2);
    f.at(1) = (1+kappa)*(tp.at(3)-tp.at(0))*tp.at(2)*u0;
    f.at(2) = (1+kappa)*(tp.at(3)-tp.at(0))*tp.at(2)*tp.at(1);
    f.at(3) = (1+kappa)*(tp.at(3)-tp.at(0))*tp.at(2)*tp.at(2)+tp.at(3);

    return;
}
