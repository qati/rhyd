#include <iostream>
#include <algorithm>
#include <unistd.h>

#define SOLVER
//#define SOLPCONST
//#define ACCSOL

#include "solver.hpp"
#include "solverwa.hpp"
#include "mq.hpp"


#if defined SOLPCONST
typedef solpconst solt;
#elif defined ACCSOL
typedef accsol solt;
#else
typedef solver solt;
#endif

using namespace std;


int main(int argc, char ** argv)
{
    cout<<"Starting threads with following opt files:"<<endl;
    for(int i=1;i<argc;i++){
        cout<<"->"<<argv[i]<<endl;
    }

    vector< vector<string> > files;

    string fend;
    files.resize(argc>1 ? (argc-1) : 1);
    for(int i=0;i<(argc-1) || (argc<=1 && i<1);i++){
        files.at(i).push_back("opt/options");
        if (argc>1){
            string tmp(argv[i+1]);
            if (tmp.find("opt")==string::npos){
                cerr<<"Options files must be in opt dir!"<<endl;
                return 0;
            }
            fend = tmp.substr(tmp.find("__"));
            tmp.erase(tmp.find("__"));
            files.at(i).push_back(tmp);
        } else {
            fend = "";
        }
        files.at(i).push_back("data/dat"+fend);
        files.at(i).push_back("data/eps"+fend);
        files.at(i).push_back("data/vn"+fend);
        files.at(i).push_back("data/nfi"+fend);
        files.at(i).push_back("data/npt"+fend);
        #ifndef SOLVER
        files.at(i).push_back("data/error"+fend);
        #endif
        files.at(i).push_back("data/numofstep"+fend);
        files.at(i).push_back("data/dt"+fend);
    }
   

    vector<io*> pfstreams;
    vector<solt*> psolvers;
    vector<mq*> pmq;

    for(vector< vector<string> >::iterator it=files.begin();it!=files.end();++it){
        pfstreams.push_back(new io(*it));
        psolvers.push_back(new solt(pfstreams.back()));
        pmq.push_back(new mq(pfstreams.back()));
    }

    for(vector<solt*>::iterator it=psolvers.begin();it!=psolvers.end();++it){
        (*it)->init();
    }

    #pragma omp parallel for
    for(uint i=0;i<psolvers.size();++i){
        psolvers.at(i)->loop();
    }

    for(uint i=0;i<psolvers.size();++i){
        pmq.at(i)->add_sol(psolvers.at(i)->get_sol());
    }

    #pragma omp parallel for
    for(uint i=0;i<pmq.size();++i){
        pmq.at(i)->calc();
    }


    for(vector<mq*>::iterator it=pmq.begin();it!=pmq.end();++it){
        delete *it;
        *it = 0;
    }
    for(vector<solt*>::iterator it=psolvers.begin();it!=psolvers.end();++it){
        delete *it;
        *it = 0;
    }
    for(vector<io*>::iterator it=pfstreams.begin();it!=pfstreams.end();++it){
        delete *it;
        *it = 0;
    }

    return 0;
}
