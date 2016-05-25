#ifndef _SOLVER_HPP_
#define _SOLVER_HPP_


#include "mff.hpp"
#include "io.hpp"


class solver{
protected:
    /**
     * Numeric solver.
     */
    mff nsolver;

    /**
     * Solution container.
     */
    hyd p;
    dtype u0;
    hyd sol;

    /**
     * IO module.
     */
    io * fs;

    /**
     * eps parameters.
     */
    uint LAST_V_ID;
    std::vector<dtype> eps;

    /**
     * Integrals to calculate eps parameters.
     */
    std::vector<dtype> eint;

    /**
     * Integral radius.
     */
    dtype ir;

    /**
     * x, y point.
     */
    dtype x, y, r, fi;

    /**
     * Space.
     */
    struct st_t {dtype x1, x2, y1, y2, dx, dy;} st;

    /**
     * P space.
     */
    struct ps_t {dtype pt0, pt1, dpt, fi0, fi1, dfi;} ps;

    /**
     * Index of grid point where n/eps is max & num of grid point.
     */
    uint i0, j0;

    /**
     * Time resolution.
     */
    dtype dt;

    /**
     * Time.
     */
    dtype t;

    /**
     * Hadronization temperature.
     */
    dtype T0, T, Tmax, Tmin;


    /**
     * Hadronized cells.
     */
    std::vector<bool> cih;
    uint hc;

    /**
      * Make sol.
      */
    void make_sol(const uint&, const uint&);

    /**
     * Calculate eps3 parameter.
     */
    void calc_integral(const uint&, const uint&);
    void null_integral();
    void calc_eps();

    /**
     * Start/end step.
     */
    virtual void start_step();
    virtual void end_step();

    /**
     * Functions.
     */
    std::vector<void (*)(void*, const uint&, const uint&)> fptr;

    /**
     * Initial distributions
     */
    void gauss();
    void woods_saxon();

    public:
        solver(io*);
        virtual ~solver();

        /**
         * Init function.
         * Sets up initial value.
         */
        virtual void init();

        /**
         * Loop function.
         */
        void loop();

        /**
         * Get solution.
         */
        hyd* get_sol();
        dtype get_t();

        /**
         * Forwarder functions.
         */
        friend void fw_calc_integral(void*, const uint&, const uint&);
};


#endif // _SOLVER_HPP_
