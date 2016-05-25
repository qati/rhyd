#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <regex>
#include <sstream>
#include <cstdlib>

const int NC = 13; //1+3*4

using namespace std;


int number_of_lines(ifstream& f)
{  
    int l=0;
    string tmp;
    while(getline(f, tmp)){
        l++;
    }
    f.clear();
    f.seekg(0, ios::beg);
    return l;
}


string out_fs(string ifs)
{
   /* string r = regex_replace(ifs, regex("psi(3|4)\\ \\=\\ [[:digit:]]+(\\.[[:digit:]]+)?"), "");
    r.erase(r.begin(), r.begin()+2); 
    r = regex_replace(r, regex("\\\n"), "_");
    r = regex_replace(r, regex("\\ "), "");
     */ 
    /*string r = "", tmp;
    stringstream s(ifs);
    for(int j=0;getline(s, tmp); ++j){
        if (j>1){
            r += tmp + "\n";
        }   
    }*/
    string r(ifs);
    r.erase(r.begin()+r.find("_psi"), r.end());
    return r;
}


int main(int argc, char ** argv)
{
    cout << endl << "===================" << endl << "Running eps_agv" << endl;
    if (argc<=1){
        cerr << "Must pass files!" << endl;
    }
    
    cout << "Make average for the following files: " <<endl;
  //  for(int i=1;i<argc;i++){
    //    cout <<"->" << argv[i] << endl;    
   // }

	vector< string > files;

    string fend;
    files.resize(argc-2);
    for(int i=0;i<(argc-2);i++){
        string tmp(argv[i+1]);
        files.at(i) = "data/eps__"+tmp;
    }
    
    string selected_eps(argv[argc-1]);

    
    
    string out_file = "data/eps_avgPsi_Psi"+selected_eps+"__"+out_fs(string(argv[1]));    
    
    vector<ifstream> in(files.size());
    for(int i=0;i<files.size();++i){
        cout << "->" << files.at(i) << endl;
        in.at(i).open(files.at(i));
    }
    
    ofstream out(out_file); 
    
    int nol = number_of_lines(in.at(0));
   /*for(auto &i : in){
        if (nol!=number_of_lines(i)){
            cerr << "Eps file error! Not the same line number of 2 files!" << endl;
            return 1;
        }
    }*/
    
    string tmp;
    double s;
    for(int i=0;i<nol;++i){
        for(int j = 0; j < NC;++j){
            s = 0;
            for(auto& f : in){
                f >> tmp;
                s += (double)atof(tmp.c_str());
            }
            s /= in.size();
            out << s;
            if (j!=(NC-1)) out << " ";
        }
        out << endl;
    }
    
    for(auto &i: in){
      i.close();
    } 
    out.close();
    return 0;
}
