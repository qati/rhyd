#ifndef _MQ_HPP_
#define _MQ_HPP_

#include "matrix.hpp"
#include "io.hpp"


class mq{
private:
    /**
     * Solution container.
     */
    hyd* p;
    dtype u0;

    /**
     * Temperature.
     */
    dtype T;

    /**
     * Particle mass.
     */
    dtype m;

    /**
     * IO module.
     */
    io * fs;

    /**
     * N integral at pt fi point.
     */
    dtype tN;
    matrix N;

    uint LAST_V_ID;

    /**
     * Vn parameters and pt parameters.
     */
    std::vector<dtype> pt;
    matrix vn, Npt;

    /**
     * Nfi functions and fi parameters.
     */
    std::vector<dtype> fi;
    matrix Nfi;


    /**
     * P space.
     */
    struct ps_t {dtype pt0, pt1, dpt, fi0, fi1, dfi;} ps;

    /**
     * Space.
     */
    dtype x, y;
    dtype x1, y1, dx, dy;

    /**
     * Calculate N(pt,fi).
     */
    void calcN(const dtype&, const dtype&);


    dtype pu, kappa;
    uint k, l;

    uint r, c;

public:
    mq(io*);

    /**
     * Add solution.
     */
    void add_sol(hyd*);

    /**
     * Measurable quantities.
     */
    void calc();
};

#endif // _MQ_HPP_
