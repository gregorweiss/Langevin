#ifndef FORCE_H_
#define FORCE_H_

#include <cstdlib>
#include <ctime>
#include <math.h>
#include <vector>
#include <boost/foreach.hpp>

#include <boost/random/linear_congruential.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>

typedef boost::mt11213b generator_type;

#include "Particle.h"

/*****************************************************************************

Force is a scaffold class for all kind of forces. Each force can be 
constructed by inhereting from Force.

It has three overloaded definitions of operator() that are used for 
force computation. 
One of them is portected and virtual. It is the one implementing the respective 
force calculation.
The other two take either one particle_ptr for external and two for pair
forces and are publicly callable.

*******************************************************************************/

////////////////////////////////////////////////////////////////////////////
/////////////////// Scaffold for Potential and Force ///////////////////////
////////////////////////////////////////////////////////////////////////////

template<class RT>
class FuncOfparcl {
public:
    FuncOfparcl() :
            _myPi(3.141592654),
            _dir(3, 1.),
            _params(vector<double>()),
            _sign(1.),
            _zero(3, 0.) {}

    FuncOfparcl<RT>(const vector<double> &p) :
            _myPi(3.141592654),
            _dir(3, 1.),
            _params(p),
            _sign(1),
            _zero(3, 0.) {}

    FuncOfparcl(const vector<double> &p,
                const double s) :
            _myPi(3.141592654),
            _dir(3, 1.),
            _params(p),
            _sign(s),
            _zero(3, 0.) {}

    FuncOfparcl(const vector<double> &p,
                const double s,
                const parcl<double> &z) :
            _myPi(3.141592654),
            _dir(3, 1.),
            _params(p),
            _sign(s),
            _zero(z) {}

    FuncOfparcl(const parcl<double> &d,        // direction
                const vector<double> &p,    // parameters
                const double s,                        // sign
                const parcl<double> &z) :    // zero
            _myPi(3.141592654),
            _dir(d),
            _params(p),
            _sign(s),
            _zero(z) {}

    virtual ~FuncOfparcl() {};

    RT operator()(const parcl<double> &);

    RT operator()(const parcl<double> &, const parcl<double> &);

    RT operator()(const particle_ptr);

    RT operator()(const particle_ptr, const particle_ptr);

protected:

    virtual RT operator()(vector<parcl<double> > &vec);

    const parcl<double> &get_dir();

    double get_param(size_t n);

    double get_sign();

    const parcl<double> &get_zero();

    const double _myPi;

private:

    const parcl<double> _dir;
    const vector<double> _params;
    const double _sign;
    const parcl<double> _zero;

};

////////////////////////////////////////////////////////////////////////////
/////////////////// Scaffold for Potential and Force ///////////////////////
////////////////////////////////////////////////////////////////////////////

template<class RT>
RT FuncOfparcl<RT>::operator()(const particle_ptr p1) {
    vector<parcl<double> > &vec(*new vector<parcl<double> >(1, p1->last()));
    vec.push_back(parcl<double>(p1->last().size(), 0.)); // zero-particle

    const RT ret(this->operator()(vec));

    delete &vec;

    return ret;
}

template<class RT>
RT FuncOfparcl<RT>::operator()(const particle_ptr p1, const particle_ptr p2) {
    vector<parcl<double> > &vec(*new vector<parcl<double> >(1, p1->last()));
    vec.push_back(p2->last());

    const RT ret(this->operator()(vec));

    delete &vec;

    return ret;
}

template<class RT>
RT FuncOfparcl<RT>::operator()(const parcl<double> &pos1) {
    vector<parcl<double> > &vec(*new vector<parcl<double> >(1, pos1));
    vec.push_back(parcl<double>(pos1.size())); // zero-particle

    const RT ret(this->operator()(vec));

    delete &vec;

    return ret;
}

template<class RT>
RT FuncOfparcl<RT>::operator()(const parcl<double> &pos1, const parcl<double> &pos2) {
    vector<parcl<double> > &vec(*new vector<parcl<double> >(1, pos1));
    vec.push_back(pos2);

    const RT ret(this->operator()(vec));

    delete &vec;

    return ret;
}

template<class RT>
RT FuncOfparcl<RT>::operator()(vector<parcl<double> > &vec) { return RT(); }


template<class RT>
const parcl<double> &FuncOfparcl<RT>::get_dir() { return _dir; }

template<class RT>
double FuncOfparcl<RT>::get_param(const size_t n) { return _params.at(n); }

template<class RT>
double FuncOfparcl<RT>::get_sign() { return _sign; }

template<class RT>
const parcl<double> &FuncOfparcl<RT>::get_zero() { return _zero; }

////////////////////////////////////////////////////////////////////////////
///////////////////////// MOTHER class Force ///////////////////////////////
////////////////////////////////////////////////////////////////////////////

class Force : public FuncOfparcl<parcl<double> > {
public:
    Force() : FuncOfparcl<parcl<double> >() {}

    Force(const vector<double> &p) : FuncOfparcl<parcl<double> >(p) {}

    Force(const vector<double> &p,
          const double s) : FuncOfparcl<parcl<double> >(p, s) {}

    Force(const vector<double> &p,
          const double s,
          const parcl<double> z) : FuncOfparcl<parcl<double> >(p, s, z) {}

    Force(const parcl<double> &d,        // direction
          const vector<double> &p,    // parameters
          const double s,                        // sign
          const parcl<double> &z        // zero
    ) : FuncOfparcl<parcl<double> >(d, p, s, z) {}

    virtual ~Force() {};
};

//////////////////////////////////////////////////////////////////////////
/////////////////////       EXTERNAL FORCES     //////////////////////////
//////////////////////////////////////////////////////////////////////////

class Random : public Force {
public:
    Random(const vector<double> &p) : Force(p),
                                      stddev(get_param(0)),
                                      generator(static_cast<unsigned int>( std::time(0))),
                                      gauss_dist(0., stddev),
                                      gauss(generator, gauss_dist) {}

    ~Random() {}

protected:

    parcl<double> operator()(vector<parcl<double> > &vec);

private:

    const double stddev;
    generator_type generator;
    boost::normal_distribution<> gauss_dist;
    boost::variate_generator<generator_type &, boost::normal_distribution<> > gauss;
};

//////////////////// class DoubleWell : class Force //////////////////////

class DoubleWell : public Force {
public:
    DoubleWell(const parcl<double> &d,
               const vector<double> &p,
               const double s,
               const parcl<double> &z) : Force(d, p, s, z),
                                         a(get_zero() * 2),
                                         a_half(a * 0.5),
                                         b(get_param(0)),
                                         barrier(get_param(1)),
                                         bsqr(b * b),
                                         pre(4. * barrier / (bsqr * bsqr)) {};

    ~DoubleWell() {};

protected:

    parcl<double> operator()(vector<parcl<double> > &vec);

private:
    const parcl<double> a;
    const parcl<double> a_half;
    const double b;
    const double barrier;
    const double bsqr;
    const double pre;
};

//////////////////// class Bias : class Force //////////////////////

class Bias : public Force {
public:
    Bias(const parcl<double> &d,
         const vector<double> &p,
         const double s,
         const parcl<double> &z) : Force(d, p, s, z),
                                   _forceconstant(get_param(0)) {};

    ~Bias() {};

protected:

    parcl<double> operator()(vector<parcl<double> > &vec);

private:
    const double _forceconstant;
};

//////////////////// class GaussEnvelope : class Force //////////////////////

class GaussEnvelope : public Force {
public:
    GaussEnvelope(const parcl<double> &d,
                  const vector<double> &p,
                  const double s,
                  const parcl<double> &z) : Force(d, p, s, z),
                                            _m(3) {
        for (auto i = 0; i < (this->_m); ++i) {

            this->_c.push_back(2.);
            this->_c.push_back(2.);
            this->_c.push_back(2.);

            Eigen::VectorXd mu(2);
            Eigen::MatrixXd sigma(2, 2);

            mu(0) = 3.0;
            mu(1) = 8.;
            sigma(0, 0) = 2.;
            sigma(0, 1) = 0.;
            sigma(1, 0) = 0.;
            sigma(1, 1) = .85;
            Mvn gauss0(mu, sigma);
            _mvn.push_back(gauss0);

            mu(0) = 8.;
            mu(1) = 8.;
            sigma(0, 0) = 2.;
            sigma(0, 1) = 0.;
            sigma(1, 0) = 0.;
            sigma(1, 1) = .85;
            Mvn gauss1(mu, sigma);
            _mvn.push_back(gauss1);

            mu(0) = 5.75;
            mu(1) = 10.25;
            sigma(0, 0) = 2.;
            sigma(0, 1) = 0.;
            sigma(1, 0) = 0.;
            sigma(1, 1) = .85;
            Mvn gauss2(mu, sigma);
            _mvn.push_back(gauss2);
        }
    };

    ~GaussEnvelope() {};

protected:
    parcl<double> operator()(vector<parcl<double> > &vec);

private:
    const size_t _m;
    vector<double> _c;
    vector<Mvn> _mvn;
};

typedef boost::shared_ptr<Force> force_ptr;

const parcl<double>
update_pair(const particle_ptr &, const particle_ptr &, const vector<force_ptr> &, const double &);

const parcl<double>
update_extern(const particle_ptr &, const vector<force_ptr> &, const vector<double> &, const double &);

const parcl<double>
update_random(const particle_ptr &, const force_ptr &);

#endif /* FORCE_H_ */
