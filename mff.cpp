#include "mff.hpp"
#include <cstdlib>
#include <iostream>

using namespace std;


mff::mff()
{
}

mff::mff(const dtype& DT, const dtype& DX,const dtype& DY, hyd* P)
{
    mff::set(DT, DX, DY, P);
}


void mff::set(const dtype& DT, const dtype& DX,const dtype& DY, hyd* P)
{
    rhyd::set(DT, DX, DY, P);

    fi_x.resize(P->row(), P->col());
    fi_y.resize(P->row(), P->col());

    ql.resize(4);
    qm.resize(4);
    qr.resize(4);
    fl.resize(4);
    fm.resize(4);
    fr.resize(4);

    return;
}


void mff::boundary(const uint& i, const uint& j)
{
    for(uint k=0;k<4;k++){
        q.at(k, i, 0)   = 2*q.at(k, i, 1)-q.at(k, i, 2);
        q.at(k, i, c+1) = 2*q.at(k, i, c)-q.at(k, i, c-1);
        q.at(k, 0, j)   = 2*q.at(k, 1, j)-q.at(k, 2, j);
        q.at(k, r+1, j) = 2*q.at(k, r, j)-q.at(k, r-1, j);

        if (i==0)         q.at(k, 0, c+1)   = q.at(k, 1, c);
        if (j==0)         q.at(k, r+1, 0)   = q.at(k, r, 1);
        if (i==0 && j==0) q.at(k, 0, 0)     = q.at(k, 1, 1);
        if (j==c)         q.at(k, r+1, c+1) = q.at(k, r, c);
    }
    return;
}


void mff::add_callfunc(void *context, void(*func)(void*, const uint&, const uint&))
{
    call_context.push_back(context);
    call_func.push_back(func);

    return;
}


void mff::iflux(void (rhyd::*flux)(std::vector<dtype>&, const dtype&, const dtype&, const dtype&, const dtype&))
{
    uint i, j, k, l, ri, rj;
    dtype ds;
    dtype tmp;
    bool dir_x;

    if (flux==&mff::flux_f){
        ds          = dx;
        dir_x       = true;
    } else {
        ds          = dy;
        dir_x       = false;
    }
    for(i=0;i<=r;i++){
        for(j=0;j<=c;j++){
            boundary(i, j);
            boundary(i+1, j+1);

            if (dir_x){
                ri = i+1;
                rj = j;
            } else {
                ri = i;
                rj = j+1;
            }

            for(k=0;k<4;k++){
                ql.at(k) = q.at(k, i, j);
                qr.at(k) = q.at(k, ri, rj);
            }
            for(l=0;l<MUSTA_K;l++){
                (this->*flux)(fl, ql.at(0), ql.at(1), ql.at(2), ql.at(3));
                (this->*flux)(fr, qr.at(0), qr.at(1), qr.at(2), qr.at(3));

                for(k=0;k<4;k++){
                    qm.at(k) = (ql.at(k)+qr.at(k))/2 - (fr.at(k)-fl.at(k))*dt/ds/2;
                }
                (this->*flux)(fm, qm.at(0), qm.at(1), qm.at(2), qm.at(3));

                for(k=0;k<4;k++){
                    tmp = (fl.at(k)+2*fm.at(k)+fr.at(k) - (qr.at(k)-ql.at(k))*ds/dt )/4;
                    if (dir_x){
                        fi_x.at(k, i, j) = tmp;
                    } else {
                        fi_y.at(k, i, j) = tmp;
                    }
                }

                for(k=0;k<4 && (l+1)<MUSTA_K;k++){
                    tmp = (dir_x) ? fi_x.at(k, i, j) : fi_y.at(k, i, j);
                    ql.at(k) -= (tmp      - fl.at(k) ) * dt/ds;
                    qr.at(k) -= (fr.at(k) - tmp      ) * dt/ds;
                }
            }
        }
    }
    return;
}


void mff::init_max()
{
    vmax  = 0;
    kappamax = 0;
    return;
}


void mff::find_max(const uint& i, const uint& j)
{
    u0 = sqrt(1+pow(p->at(1, i, j),2)+pow(p->at(2, i, j),2));
    if (vmax<fabs(p->at(1, i, j)/u0)){
        vmax  = fabs(p->at(1, i, j)/u0);
        mdirx = true;
    }
    if (vmax<fabs(p->at(2, i, j)/u0)){
        vmax  = fabs(p->at(2, i, j)/u0);
        mdirx = false;
    }
    eos->ck(p->at(3, i, j)/p->at(0, i, j));
    if (kappamax<kappa){
        kappamax = kappa;
    }
    return;
}


void mff::update_dt()
{
    if (vmax==0){
        return;
    }
    dt = CFLN * ( (mdirx) ? dx : dy ) / (vmax+1/sqrt(kappamax));
    return;
}


dtype mff::navg(bool FX, const uint& k,  const uint& i, const uint& j)
{
    if (FX==FLUX_X)
        return (fi_x.at(k, i, j-1)+fi_x.at(k, i-1, j) + fi_x.at(k, i, j+1) + fi_x.at(k, i+1, j) ) /4;
    else
        return (fi_y.at(k, i, j-1)+fi_y.at(k, i-1, j) + fi_y.at(k, i, j+1) + fi_y.at(k, i+1, j) ) /4;
    cerr<<"mff::navg fatal error"<<endl;
    exit(1);
}


void mff::step()
{

    uint i, j, k;
    init_max();
    
    dtype tt;
    iflux(&mff::flux_f);
    for(i=1;i<=r;i++){
        for(j=1;j<=c;j++){
            for(k=0;k<4;k++){
                /*if (i!=1 & i!=r){
                    tt = navg(FLUX_X, k, i, j);
                    if (fabs(tt-fi_x.at(k, i, j))>15){
                        fi_x.at(k, i, j) = tt;
                    }
                    tt = navg(FLUX_X, k, i-1, j);
                    if (abs(tt-fi_x.at(k, i-1, j))>15){
                        fi_x.at(k, i-1, j) = tt;
                    }
                }*/
                q.at(k, i, j) -= (fi_x.at(k, i, j)-fi_x.at(k, i-1, j)) * dt/dx;
            }
        }
    }

    iflux(&mff::flux_g);
    non = 0;
    for(i=0;i<=r;i++){
        for(j=1;j<=c;j++){
            for(k=0;k<4;k++){
                /*if (j!=1 && j!=c && i!=0){
                    tt = navg(FLUX_Y, k, i, j);
                    if (abs(tt-fi_y.at(k, i, j))>15){
                        fi_y.at(k, i, j) = tt;
                    }
                    tt = navg(FLUX_Y, k, i, j-1);
                    if (abs(tt-fi_y.at(k, i, j-1))>15){
                        fi_y.at(k, i, j-1) = tt;
                    }
                }*/
                q.at(k, i, j) -= (fi_y.at(k, i, j)-fi_y.at(k, i, j-1)) *dt/dy;
            }

            /**
             * Call functions.
             */
            qp(i, j);
            find_max(i, j);
            if (!call_func.empty() || call_func.size()!=call_context.size()){
                for(k=0;k<call_func.size();k++){
                    call_func.at(k)(call_context.at(k), i, j);
                }
            }
        }

    }

    update_dt();
    return;
}


dtype mff::get_flux(bool FLUX, const uint& k, const uint& i, const uint& j)
{
    return ( (FLUX==FLUX_X) ? fi_x.at(k, i, j) : fi_y.at(k, i, j) );
}
