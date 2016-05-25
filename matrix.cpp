#include "matrix.hpp"
#include <iostream>
#include <cassert>

using namespace std;

matrix::matrix() : r(0), c(0)
{
}

matrix::matrix(uint i, uint j) : r(i), c(j)
{
    v.resize(r*c);
}

matrix::matrix(const matrix& m)
{
    r = m.r;
    c = m.c;
    v.resize((r+NUM_OF_GHOST_CELL)*(c+NUM_OF_GHOST_CELL));
    for(uint i=0;i<(r+NUM_OF_GHOST_CELL);i++){
        for(uint j=0;j<(c+NUM_OF_GHOST_CELL);j++){
            at(i,j) = m.at(i,j);
        }
    }
}

void matrix::resize(uint i, uint j)
{
    r = i;
    c = j;
    v.resize((r+NUM_OF_GHOST_CELL)*(c+NUM_OF_GHOST_CELL));
}

const dtype & matrix::at(uint i, uint j) const
{
    if (i>=(r+NUM_OF_GHOST_CELL) || j>=(c+NUM_OF_GHOST_CELL)){
        cout<<"matrix::at const error!"<<endl<<"r:"<<r<<"; c:"<<c<<";"<<endl<<"    i:"<<i<<"j:"<<j<<endl;
        assert(i<r&&j<c);
    }
    return v.at(i*c+j);
}

dtype & matrix::at(uint i, uint j)
{
     if (i>=(r+NUM_OF_GHOST_CELL) || j>=(c+NUM_OF_GHOST_CELL)){
        cout<<"matrix::at error!"<<endl<<"r:"<<r<<"; c:"<<c<<";"<<endl<<"    i:"<<i<<"j:"<<j<<endl;
        assert(i<r&&j<c);
    }
    return v.at(i*c+j);
}



ostream& operator<<(ostream &o, const matrix& m)
{
    for(uint i=0;i<m.r;i++){
        for(uint j=0;j<m.c;j++){
            o<<m.at(i,j)<<" ";
        }
        o<<endl;
    }
    return o;
}



/**---------------------------------------------------HYD----------------------------------------**/

hyd::hyd() : r(0), c(0)
{
    nghost = NUM_OF_GHOST_CELL;
}

hyd::hyd(uint i, uint j) : r(i), c(j)
{
    nghost = NUM_OF_GHOST_CELL;
    v.resize(HYD_VAR*(r+nghost)*(c+nghost));
}

hyd::hyd(const hyd& h)
{
    r      = h.r;
    c      = h.c;
    nghost = NUM_OF_GHOST_CELL;
    v.resize(HYD_VAR*(r+nghost)*(c+nghost));

    for(int q=0;q<HYD_VAR;q++){
        for(uint i=0;i<(h.r+NUM_OF_GHOST_CELL);i++){
            for(uint j=0;j<(h.c+NUM_OF_GHOST_CELL);j++){
                v.at(j+i*(c+nghost)+q*(c+nghost)*(r+nghost)) = h.v.at(j+i*(c+nghost)+q*(c+nghost)*(r+nghost));
            }
        }
    }
}

void hyd::resize(uint row, uint col)
{
    r = row;
    c = col;
    v.resize(HYD_VAR*(r+nghost)*(c+nghost));
    return;
}

const dtype& hyd::at(uint q, uint i, uint j) const
{
    if (q>=HYD_VAR || i>=(r+nghost) || j>=(c+nghost)){
        cout<<"hyd::at const error!"<<endl<<"HYD_VAR: "<<HYD_VAR<<"; r:"<<r<<"; c:"<<c<<"; "<<"nghost:"<<nghost<<endl
            <<"    q:"<<q<<"; i:"<<i<<"; j:"<<j<<endl;
        assert(q<HYD_VAR&&i<r&&j<c);
    }
    return v.at(j+i*(c+nghost)+q*(r+nghost)*(c+nghost));
}

dtype& hyd::at(uint q, uint i, uint j)
{
    if (q>=HYD_VAR || i>=(r+nghost) || j>=(c+nghost)){
        cout<<"hyd::at error!"<<endl<<"HYD_VAR: "<<HYD_VAR<<"; r:"<<r<<"; c:"<<c<<"; "<<"nghost:"<<nghost<<endl
            <<"    q:"<<q<<"; i:"<<i<<"; j:"<<j<<endl;
        assert(q<HYD_VAR&&i<r&&j<c);
    }
    return v.at(j+i*(c+nghost)+q*(r+nghost)*(c+nghost));
}




ostream& operator<<(ostream &o, const hyd& h)
{
    for(uint q=0;q<HYD_VAR;q++){
        o<<"~~~~~~~~~~~~~~q="<<q<<"~~~~~~~~~~~~~~"<<endl;
        for(uint i=0;i<(h.r+NUM_OF_GHOST_CELL);i++){
            for(uint j=0;j<(h.c+NUM_OF_GHOST_CELL);j++){
                o<<h.at(q,i,j)<<" ";
            }
            o<<endl;
        }
        o<<endl;
    }
    return o;
}


std::ostream& operator << (std::ostream& o, const std::vector<dtype>& vec)
{
    uint i, n = vec.size();
    for(i=0;i<n;i++){
        o<<vec.at(i)<<"; ";
    }
    o<<std::endl;
    return o;
}

