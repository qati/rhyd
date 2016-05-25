#ifndef _RHYD_HPP_
#define _RHYD_HPP_

#include "matrix.hpp"
#include "eos.hpp"
#include <cmath>


class rhyd{
protected:
   /**
     * KAPPA.
     */
    dtype kappa;

    /**
     * @p: physical quantity: n, u1, u2, eps  & u0
     * @q: transformed quantity
     */
    hyd *p, q;
    dtype u0;


    /* Temperature.
     */
    dtype T;

    /**
     * EoS
     */
    EoS *eos;

    /**
     * Calculate fluxes from q1, q2, q3, q4, q5.
     */
    void flux_f(std::vector<dtype>&, const dtype&, const dtype&, const dtype&, const dtype&);
    void flux_g(std::vector<dtype>&, const dtype&, const dtype&, const dtype&, const dtype&);

    /**
     * Set variables.
     */
    void set(const dtype&, const dtype&,const dtype&, hyd*);

    /**
     * Makes qi to tp transformation.
     */
    void qtp(const dtype&, const dtype&, const dtype&, const dtype&);

    /**
     * Temp vars for qtp!
     **/
    std::vector<dtype> tp;
    dtype t0, t1, t2, t3;
    uint non;


     /**
     * Spacetime resolution.
     */
    dtype dt, dx, dy;

    /**
     * Number of rows and cols.
     */
    uint r, c;

    /**
     * Boundary condition.
     */
    virtual void boundary(const uint&, const uint&) = 0;

public:
    rhyd();
    rhyd(const dtype&, const dtype&,const dtype&, hyd*);

    /**
     * Set variables.
     */
    void set(const dtype&);

    /**
     * Makes p to q transformation.
     */
    void pq();

    /**
     * Makes q to p transformation.
     */
    void qp(const uint&, const uint&);

    /**
     * Time advance function.
     */
    virtual void step() = 0;

    /**
     * Get dt.
     * Get number of nans.
     * Get num of points.
     */
    dtype get_dt(){return dt;}
    uint get_non(){uint t=non; non=0; return t;}
    uint get_nop(){return r*c;}
    dtype get_T();
};


#endif // _RHYD_HPP_
