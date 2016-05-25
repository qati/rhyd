#include <memory>
#include <iostream>
#include "io.hpp"
#include <cassert>
#include <cstdlib>
#include <algorithm>
#include <cmath>


using namespace std;


io::io(vector<string>& files)
{
    open(files);
    dfn = files.at(0);
}

io::~io()
{
    for(vector<ofstream*>::iterator i=out.begin();i!=out.end();++i){
        (*i)->close();
    }
}


void io::set_sgrid(const uint& row, const uint& col)
{
    dtype x0 = get("x1"),
          y0 = get("y1"),
          dx = get("dx"),
          dy = get("dy");
    r = row;
    c = col;

    x.resize(r);
    y.resize(c);

    for(uint i=0;i<r;i++){
        x.at(i) = x0 + i*dx;
    }
    for(uint i=0;i<c;i++){
        y.at(i) = y0 + i*dy;
    }
    
    mod_r = 1;
    mod_c = 1;
    if (r>300){
        mod_r = 2;
    }
    if (c>300){
        mod_c = 2;
    }
    if (r>400){
        mod_r = 3;
    }
    if (c>400){
        mod_c = 4;
    }
    return;
}


void io::add_hyd(hyd *H1)
{
    h.push_back(H1);

    return;
}


void io::add_vec(vector<vector<dtype>*> dat)
{
    d.push_back(dat);

    return;
}


void io::add_m(vector<dtype>* MP, matrix* M)
{
    if (MP->size()!=M->row()){
        cerr<<"io::add_mvec: Wrong usage! MP.size="<<MP->size()<<"; M.row="<<M->row()<<endl;
    }
    m.push_back(M);
    mp.push_back(MP);

    return;
}


void io::open(vector<string>& files)
{
    ofstream * o;
    read(*files.begin());
    files.erase(files.begin(), ++files.begin());
    if (files.begin()->find("opt")!=string::npos){
        open(files);
        return;
    }
    for(vector<string>::const_iterator i=files.begin();i!=files.end();++i){
        try{
            o = new ofstream(i->c_str(), ios::out);
        } catch(exception& e){
            cerr<<"Can't open file: "<<*i<<"; Exception: "<<e.what()<<endl;
        }
        out.push_back(o);
    }
    return;
}


void io::read(string file)
{
    ifstream in(file.c_str(), ios::in);
    if (!in.is_open()){
        cerr << "File: " << file <<" not opened! Fatal error! Exit!" <<endl;
        exit(1);
    }
    string r2, s2, n2;
    while(!in.eof()){
        in>>r2;
        if (r2=="=" && r2!=s2){
            in>>n2;
            dat[s2] = (dtype)atof(n2.c_str());
        }
        s2 = r2;
    }
    in.close();
    return;
}


void io::read_datfile_to_p(string file)
{
    cout << "io::read_datfile_to_p start " << endl;
    ifstream in(file.c_str(), ios::in);
    if (!in.is_open()){
        cerr << "io::read_datfile_to_p: Can't open file: " << file <<"; Stopping.." <<endl;
        exit(1);
    }
    uint i, j;
    dtype rx, ry, rn, rvx, rvy, rp;
    dtype u0;
    for(i=1;i<=r;i++){
        for(j=1;j<=c;j++){
            in >> rx;
            //if (rx!="\n"){
                in >> ry >> rn >> rvx >> rvy >> rp;
               // cout<<rx <<" " <<ry <<" == " <<rn <<" "<<rvx<<" "<<rvy<<" "<<rp<<endl;
           // } else {
               // cout<<endl<<endl <<"newline:"<<rx<<endl;
            //}
            u0 = sqrt(1 + pow( h.at(0)->at(1, i, j), 2) + pow(h.at(0)->at(2, i, j),2));
           // if (u0 == 0 || abs(u0)<DEF_NULL){
             //   cerr << "io::read_datfile_to_p: Something wrong with speed! (i,j)=("<<i<<","<<j<<"); rvx="<<rvx<<"; rvy="<<rvy<<"; u0="<<u0<<endl;
               // exit(1);
            //}
            //rvx = rvx*u0;
            //rvy = rvy*u0;
            //if (abs(h.at(0)->at(1,i,j)-rvx)>0.1 || abs(h.at(0)->at(2,i,j)-rvy)>0.1){
             //   cerr << "Something wrong in reading smoth initial condition: speed doesn't match!" << endl;
              //  cerr << "(i,j)=("<<i<<","<<j<<"); rvx="<<rvx<<"; rvy="<<rvy<<"; u0="<<u0
               //      << "; p(1)="<<h.at(0)->at(1,i,j)<<"; p(2)="<<h.at(0)->at(2,i,j)<<endl;
               // exit(1);
            //}
            h.at(0)->at(0, i,j) = rn*u0;
            h.at(0)->at(3, i,j) = rp*u0;
            h.at(0)->at(1, i,j) = rvx*u0;
            h.at(0)->at(2, i,j) = rvy*u0;
            //if (h.at(0)->at(1,i,j)!=0){
              //  cout << i<< " "<< j<<endl;
            //}
        }
    } 
    in.close();
    return;
}

dtype io::get(const string& key)
{
    if (dat.find(key)==dat.end()){
        cerr<<"ERROR::io::get: Key '"<<key<<"' not found!"<<endl;
        assert(1);
    }
    return dat[key];
}


void io::write(const dtype& time, dtype t0)
{
    uint i, j;
    dtype u0;
    vector<hyd*>::iterator it;

    if (time==0){
        out.at(0)->seekp(0);
    }
    if (time==-1){
    for(i=1;i<=r;i++){
        for(j=1;j<=c;j++){
            if ((i%mod_r==0 && j%mod_c==0) || time==-1){
            (*out.at(0))<<x.at(i-1)<<" "<<y.at(j-1);
            for(it=h.begin();it!=h.end();++it){
                u0 = sqrt(1 + pow( (*it)->at(1, i, j), 2) + pow((*it)->at(2, i, j),2));
                (*out.at(0))<<" "<<(*it)->at(0, i, j)/u0<<" "<<(*it)->at(1, i, j)/u0<<" "<<(*it)->at(2, i, j)/u0<<" "<<(*it)->at(3, i, j)/u0;
            }
            (*out.at(0))<<endl;
        }}
        if (i%mod_r==0 || time==-1) (*out.at(0))<<endl;
    }
    (*out.at(0))<<endl;}
    
    if (time==-1) return;
    if (time!=0) t0 = time;
   
    vector<vector<dtype>*>::iterator itv;
    vector<dtype>::iterator itvv;
    for(i=0;i<d.size();i++){
        (*out.at(i+1))<<t0;
        for(itv=d.at(i).begin();itv!=d.at(i).end();++itv){
            for(itvv=(*itv)->begin();itvv!=(*itv)->end();++itvv){
                (*out.at(i+1))<<" "<<*itvv;
            }
        }
        (*out.at(i+1))<<endl;
    }

    return;
}

void io::write()
{
    uint i, j, fi, k;
    for(i=0;i<mp.size();++i){
        fi = d.size() + 1 + i;
        for(j=0;j<mp.at(i)->size();++j){
            (*out.at(fi))<<mp.at(i)->at(j);
            for(k=0;k<m.at(i)->col();++k){
                (*out.at(fi))<<" "<<m.at(i)->at(j, k);
            }
            (*out.at(fi))<<endl;
        }
    }
    return;
}

void io::write_nos(const uint& N)
{
    (**----out.end())<<N<<endl;
    return;
}

void io::write_dt(const dtype& dt)
{
    (*out.back())<<dt<<endl;
    return;
}

