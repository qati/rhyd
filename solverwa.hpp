#ifndef _SOLVERWA_HPP_
#define _SOLVERWA_HPP_


#include "solver.hpp"


class solverwa : public solver{
    protected:
        /**
         * @a: analitic solution
         */
        hyd a;

        /**
         * Relative and absolute difference.
         */
        std::vector<dtype> rd, ad;

        /**
         * @XYVW: scale variables.
         */
        dtype X, Y, V, W;

        /**
         * Kappa.
         */
        dtype kappa;

        /**
         * Init values.
         */
        struct iv_t {dtype n0, p0;} iv;

        /**
         * Analitic solution.
         */
        virtual void analitic(const uint&, const uint&);

        /**
         * Calulate errors.
         */
        void error(const uint&, const uint&);

        /**
         * Start/End time step calc function.
         */
        virtual void start_step();
        virtual void end_step();

        /**
         * Loop index.
         */
        uint ik;

        /**
         * P correction.
         */
        dtype p_s_corr;

    public:
        solverwa(io*);
        virtual ~solverwa() = 0;

        /**
         * Initialization function.
         */
        void init();

         /**
         * Forwarder functions.
         */
        friend void fw_calc(void *, const uint&, const uint&);
};




/**
 * Exact 2+1D solution with const p.
 * @n   = n0*(tau/tau0)^3*exp(-s/2)
 * @u   = u0*(1, dX/X*x, dY/Y*y)
 * @p   = p0*(tau/tau0)^{3+3/kappa}
 * @XY  = X=dX*t, Y=dY*t
 * @tau = t/u0
 */
class solpconst : public solverwa{
    dtype t0, tau, tau0;

    void start_step();
    void end_step();
    void analitic(const uint&, const uint&);

    public:
        solpconst(io*);
};


class accsol : public solverwa{
    dtype tau, tau0, eta, t0;

    void start_step();
    void end_step();
    void analitic(const uint&, const uint&);

    public:
        accsol(io*);
};


#endif // _SOLVERWA_HPP_
