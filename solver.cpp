#include "solver.hpp"
#include <iostream>
#include <cassert>
#include <cstdlib>
#include <string>


using namespace std;


void fw_calc_integral(void* context, const uint& i, const uint& j)
{
    static_cast<solver*>(context)->calc_integral(i, j);
    static_cast<solver*>(context)->make_sol(i, j);
    return;
}


solver::solver(io * FS)
{
    st.x1 = FS->get("x1");
    st.x2 = FS->get("x2");
    st.y1 = FS->get("y1");
    st.y2 = FS->get("y2");
    st.dx = FS->get("dx");
    st.dy = FS->get("dy");

    p.resize( (fabs(st.x1)+fabs(st.x2))/st.dx+1, (fabs(st.y1)+fabs(st.y2))/st.dy+1 );

    nsolver.set(FS->get("dt"), FS->get("dx"), FS->get("dy"), &p);
    nsolver.set(FS->get("kappa"));

    LAST_V_ID = FS->get("LAST_V_ID");
    eps.resize(LAST_V_ID*3);
    eint.resize((LAST_V_ID+1)*3);

    ir = FS->get("int_rad");

    fs = FS;
    fs->set_sgrid(p.row(), p.col());
    fs->add_hyd(&p);
    fs->add_vec(vector<vector<dtype>*>(1, &eps));
    i0 = -st.x1/st.dx+1,
    j0 = -st.y1/st.dy+1;

    dt = fs->get("dt");
    t  = fs->get("t0");

    fptr.push_back(&fw_calc_integral);

    T0 = fs->get("T0");

    sol.resize(p.row(), p.col());

    cih.resize((p.row()+1)*(p.col()+1));


    hc   = 0;
    Tmax = 0;
    Tmin = 100000;
}

solver::~solver()
{
    fs = NULL;
}


void solver::gauss()
{
    dtype n0    = fs->get("n0"),
          p0    = fs->get("p0"),
          t0    = fs->get("t0");

    dtype R  = fs->get("R"),
          E2 = fs->get("E2"),
          E3 = fs->get("E3"),
          E4 = fs->get("E4");

    dtype p_s_corr = fs->get("p_s_corr");

    dtype s;

    uint i, j;
    uint  r = p.row(),
          c = p.col();

    dtype psi2 = fs->get("psi2"),
          psi3 = fs->get("psi3"),
          psi4 = fs->get("psi4");

    cout << "Initial condition!"
            << "psi2 = " << psi2 << endl
            << "psi3 = " << psi3 << endl
            << "psi4 = " << psi4 << endl
            << "p_s_corr = " << p_s_corr << endl
            << "====================================" <<endl;

    cout << "Integral radius: "<<ir<<endl;
    
    p.at(0, i0, j0) = n0;
    p.at(3, i0, j0) = p0;
    cout << R << " " << E2 <<" " <<E3 << " " <<E4<<endl;
    for(i=1;i<=r;i++){
        for(j=1;j<=c;j++){
            x = st.x1 + (i-1)*st.dx;
            y = st.y1 + (j-1)*st.dy;

            if (x==0 && y==0){
                s = 0;
            } else {
                fi = atan2(y, x);
                s = (x*x+y*y)  * pow(R, -2) * (1+E2*cos(2*(fi-psi2))+E4*cos(4*(fi-psi4))+E3*cos(3*(fi-psi3)));
            }
            
            
            if (s>-R*0.25 && s<R*0.25){
            
                p.at(0, i, j) = n0 * exp(-s/2);
                p.at(1, i, j) = 0;
                p.at(2, i, j) = 0;
                p.at(3, i, j) = p0*exp(-p_s_corr*s/2);
            } else {
                p.at(0, i, j) = 0;
                p.at(1, i, j) = 0;
                p.at(2, i, j) = 0;
                p.at(3, i, j) = 0;
            }

            calc_integral(i, j);
        }
    }
    calc_eps();
    
    cout<<"New integral radius: "<<ir <<endl;
    return;
}


void solver::init()
{
    gauss();

    fs->write(-1); //it's importatnt that t=-1, because in this way write will print the whole matrix

    //if gauss with 0 then we need to make te distribution smooth
    cout << "Calling make_smth.py for data file: " << fs->get_dfn()<<endl;
    string call = "python make_smth.py "+fs->get_dfn()+" "+to_string(t);
    system(call.c_str());
    fs->read_datfile_to_p(fs->get_dfn());
    cout <<" End of data smoothing"<<endl; 
    //writing the initial condition to file, but not the whole matrix
    fs->write(0, t);
    nsolver.pq();

    return;
}


void solver::calc_integral(const uint& i, const uint& j)
{
    uint k;

    x = st.x1 + (i-1)*st.dx;
    y = st.y1 + (j-1)*st.dy;

    if ((x*x+y*y)>ir*(st.x2*st.x2+st.y2*st.y2)) return;
    if (p.at(0,i,j)<DEF_NULL && p.at(3, i, j)<DEF_NULL) return;
    if (p.at(2,i,j)==0.0 && p.at(1, i,j)==0.0) return;

    if (x==0 && y==0){
        for(vector<dtype>::iterator it=eint.begin();it!=eint.end();++it){
            (*it) -= 1;
        }
    } else {
        fi = atan2(y, x);
        u0 = (pow(p.at(1, i, j), 2)+pow(p.at(2, i, j), 2));
        for(k=0;k<=LAST_V_ID;k++){
            eint.at(k)               += p.at(0, i, j) * cos(k*fi);
            eint.at(2*LAST_V_ID+2+k) += p.at(3, i, j) * cos(k*fi);
            //if ((p.at(1, i, j)!=0.0 && p.at(2, i,j)!=0.0) && (p.at(1, i, j)>DEF_NULL && p.at(2, i, j)>DEF_NULL)){
            if (p.at(0, i,j)>DEF_NULL && p.at(3,i,j)>DEF_NULL){
                eint.at(LAST_V_ID+1+k) += exp(-u0)    * cos(k*fi);}
        }
    }
    return;
}


void solver::null_integral()
{
    eint.clear();
    eint.resize((LAST_V_ID+1)*3);
    return;
}


void solver::calc_eps()
{
    for(uint i=0;i<LAST_V_ID;i++){
        if (eint.at(i+1)==eint.at(0) && eint.at(i+1)==0) eps.at(i) = 0;
        else eps.at(i) = eint.at(i+1)/eint.at(0);
    }
    for(uint i=LAST_V_ID;i<2*LAST_V_ID;i++){
        if (eint.at(i+2)==eint.at(LAST_V_ID+1) && eint.at(i+2)==0) eps.at(i) = 0;
        else eps.at(i) = eint.at(i+2)/eint.at(LAST_V_ID+1);
    }
    for(uint i=2*LAST_V_ID;i<3*LAST_V_ID;i++){
        if (eint.at(i+3)==eint.at(2*LAST_V_ID+2) && eint.at(i+3)==0) eps.at(i) = 0;
        else eps.at(i) = eint.at(i+3)/eint.at(2*LAST_V_ID+2);

    }

    return;
}


void solver::start_step()
{
    return;
}


void solver::end_step()
{
    return;
}

bool is_null(vector<dtype> v){
    int k=0;
    for(int i=0;i<4;i++){
        if (v.at(i)==0) k++;
    }
    return k==4;
}

void solver::make_sol(const uint& i, const uint& j)
{
    T = nsolver.get_T();//p.at(3, i, j)/p.at(0, i, j);
    if (Tmax<T) Tmax = T;
    if (Tmin>T && T>0) Tmin = T;

    return;
}


void solver::loop()
{
    uint N = 0;
    dtype tmax = t + fs->get("RTT");
    //dtype tws  = (dtype)fs->get("RTT")/fs->get("NOPTBW"),
          //ct   = 0;
    uint ws = fs->get("WS");
    //uint rc = p.row()*p.col();

    for(vector<void (*)(void*, const uint&, const uint&)>::iterator it=fptr.begin();it!=fptr.end();++it){
        nsolver.add_callfunc(this, *it);
    }

    for(uint i=0;t<tmax;++i){//&& hc<rc;i++){
        dt = nsolver.get_dt();
        t += dt;

        null_integral();
        start_step();
        nsolver.step();
        end_step();

        if (i%ws==0){
            fs->write_dt(dt);
            calc_eps();
            fs->write(t);
            fs->write_nos(N++);
            cout<<"param: kappa="<<fs->get("kappa")<<", p0="<<fs->get("p0")<<", pc="<<fs->get("p_s_corr")
                   <<"; "<<"time: "<<t<<"fm/c; hc: "<<hc<<"; Tmax="<<Tmax<<"; Tmin="<<Tmin<<"; nop: "<<nsolver.get_nop()
                                    <<"; non: "<<nsolver.get_non()<<endl;

        }
        /*Tmax = 0;
        ct += dt;
        if (ct>=tws){
            ct = 0;
            calc_eps();
            fs->write(t);
            fs->write_nos(N++);
        }*/
        /*for(uint k=1;k<=p.row();++k){
            for(uint l=1;l<=p.col();++l){
                T = p.at(3, k, l)/p.at(0, k, l);
                if (T<=T0 && !cih.at(k*p.col()+l)){
                    sol.at(0, k, l) = p.at(0, k, l);
                    sol.at(1, k, l) = p.at(1, k, l);
                    sol.at(2, k, l) = p.at(2, k, l);
                    sol.at(3, k, l) = p.at(3, k, l);
                    cih.at(k*p.col()+l) = 1;
                    hc++;
                }
            }
        }*/
    }
   /* int tt=0;
    for(uint i=1;i<=p.row();++i){
        for(uint j=1;j<=p.col();j++){
           if (sol.at(0, i, j)==0){
            cout<<" n=0 @ "<<i<<", "<<j<<endl;
            tt++;
           }
           if (sol.at(3, i, j)==0){
            cout<<" p=0 @ "<<i<<", "<<j<<endl;
           }
        }
    }
    cout<<endl<<"TTT: "<<tt<<endl;*/
    return;
}


hyd* solver::get_sol()
{
    return &sol;
}

dtype solver::get_t()
{
    return t;
}
