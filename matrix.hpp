#ifndef _MATRIX_HPP_
#define _MATRIX_HPP_


#include <vector>
#include <ostream>


typedef unsigned int uint;
typedef short int sint;
typedef double dtype;


#define HYD_VAR 4
#define NUM_OF_GHOST_CELL 2


const dtype DEF_NULL = 1e-14;


//first index row

class matrix{
    private:
        uint r,c;
        std::vector<dtype> v;
    public:
        matrix();
        matrix(uint, uint);
        matrix(const matrix&);

        void resize(uint, uint);

        const dtype & at(uint, uint) const;
        dtype & at(uint, uint);

        uint row()const{return r;}
        uint col()const{return c;}

        friend std::ostream& operator << (std::ostream&, const matrix&);
};


class hyd
{
    private:
        uint r, c;
        uint nghost;
        std::vector<dtype> v;
    public:
        hyd();
        hyd(uint, uint);
        hyd(const hyd&);

        void resize(uint, uint);

        const dtype & at(uint, uint, uint) const;
        dtype & at(uint, uint, uint);

        uint row()const{return r;}
        uint col()const{return c;}

        friend std::ostream& operator << (std::ostream&, const hyd&);
};

std::ostream& operator << (std::ostream& o, const std::vector<dtype>& vec);


inline dtype sgn(const dtype& v){
    return v<0 ? -1.0 : 1.0;
}

#endif // _MATRIX_HPP_
