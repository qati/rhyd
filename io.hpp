#ifndef _IO_HPP_
#define _IO_HPP_


#include <fstream>
#include <string>
#include <map>
#include "matrix.hpp"


class io{
    private:
        /**
         * File streams.
         */
        std::vector<std::ofstream*> out;

        /**
         * Space grid and size.
         */
        std::vector<dtype> x, y;
        uint r, c;
    
        uint mod_r, mod_c;
        
        /**
         * Source variables.
         */
        std::vector<hyd*> h;
        std::vector< std::vector<std::vector<dtype>*> > d;
        std::vector<matrix*> m;
        std::vector< std::vector<dtype>*> mp;

        /**
         * Vector to file mode.
         */
        std::vector<uint> vtfm;
        uint lastvf;

        /**
         * Input.
         */
        std::map<std::string, dtype> dat;

        void read(std::string);
        
        std::string dfn;
        

    public:
        io(std::vector<std::string>&);
        ~io();


        void set_sgrid(const uint&, const uint&);
        void add_hyd(hyd*);
        void add_vec(std::vector<std::vector<dtype>*>);
        void add_m(std::vector<dtype>*, matrix*);

        void open(std::vector<std::string>&);

        dtype get(const std::string&);

        void write(const dtype&, dtype t0 = 0);
        void write();
        void write_nos(const uint&);
        void write_dt(const dtype&);

        void read_datfile_to_p(std::string);    
    
        std::string get_dfn(){return dfn;}
};


#endif // _IO_HPP_
