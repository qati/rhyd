#ifndef _MFF_HPP_
#define _MFF_HPP_


#include "rhyd.hpp"


#define MUSTA_K 24

const dtype CFLN = 0.3 ;

const bool FLUX_X = 1 ;
const bool FLUX_Y = 0 ;


class mff : public rhyd{
    /**
     * f_{i+1/2} flux.
     */
    hyd fi_x, fi_y;

    /**
     * TMP variables for MUSTA method.
     */
    std::vector<dtype> ql, qm, qr, fl, fm, fr;

    /**
     * Calculate fi with MUSTA FORCE K method.
     */
    void iflux(void (rhyd::*)(std::vector<dtype>&, const dtype&, const dtype&, const dtype&, const dtype&));

    /**
     * Function forwarder and obj.
     */
    std::vector<void*> call_context;
    std::vector<void (*)(void*, const uint&, const uint&)> call_func;

    /**
     * Boundary condition.
     */
    void boundary(const uint&, const uint&);

    /**
     * Find max of v.
     */
    dtype vmax, kappamax;
    bool mdirx;

    void init_max();
    void find_max(const uint&, const uint&);

    void update_dt();

    dtype navg(bool, const uint&, const uint&, const uint&);


    public:
        mff();
        mff(const dtype&, const dtype&,const dtype&, hyd*);

        /**
         * Set variables.
         */
        void set(const dtype& kappa){rhyd::set(kappa);}
        void set(const dtype&, const dtype&,const dtype&, hyd*);

        void add_callfunc(void*, void (*)(void*, const uint&, const uint&));

        void step();

        dtype get_flux(bool, const uint&, const uint&, const uint&);
};


#endif // _MFF_HPP_
