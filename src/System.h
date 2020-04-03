#ifndef SYSTEM_H_
#define SYSTEM_H_

#include <iostream>
#include <iomanip>

#include <cstdlib>
#include <vector>
#include <cmath>
#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>

#include <boost/random/linear_congruential.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>

typedef boost::mt11213b generator_type;

#include <omp.h>

#include "common.h"
#include "ArgParse.h"
#include "Particle.h"
#include "HIcouple.h"
#include "Force.h"

using std::vector;
using boost::shared_ptr;

////////////////////////////////////////////////////////////////////////////
/////////////////////////// MOTHER class System ////////////////////////////
////////////////////////////////////////////////////////////////////////////

/*
* class System implements langevin dynamics simulations. It can be used for
* setups of brownian dynamics systems with varying number of particles, particle 
* interaction as well as external forces and hydrodynamic interactions.
* class System implements the tools for system constuction and the integrated
* which implements the Ermak-McCammon algorithm.
*
* Definition of user flags at defined runtime:
* "-box"		default (double) 		10.0		// box length (for simplicity at the moment only cubic)
* "-dt"			default (double) 		1.e-3		// timestep for integration
* "-temp"		default (double) 		300.		// temperature
* "-prec"		default (unsigned)	9				// precision for writing member function print()
*	"-act"		default (double)		1.0			// parameter to turn on/off "actio" from the equations of motion
* "-react"	default (double)		1.0			// parameter to turn on/off "reactio" from the equations of motion
*/

class System {
public:
    System(ArgParse &input) :
            scale(),
            _dimension(2),
            _box(_dimension, input.flag<double>("-box", 100.0)),
            _freeze(_dimension, 1.),
            _kB(scale.kB),
            _timestep(input.flag<double>("-dt", 1.0e-3)),
            _temp(input.flag<double>("-temp", 300.)),
            _prec(input.flag<unsigned>("-prec", 9)),
            _actio(input.flag<double>("-act", 1.0)),
            _reactio(input.flag<double>("-react", 1.0)),
            lhs(vector<particle_ptr>()),
            cpl(new EmptyHI()),
            rhs_random(new Random(vector<double>(1, sqrt(2. * _timestep)))),
            rhs_pair(vector<vector<force_ptr> >(0)),
            rhs_extern(vector<vector<force_ptr> >(0)),
            state_extern(vector<vector<double> >(0)) {
        // Every system should have the opportunity for umbrella potential as used for
        // weighted histogram analysis method or spatial friction profiles e.g.
        input.flag<double>("-up", 2.0);
        input.flag<double>("-uk", 0.0);
    };

    virtual ~System() = default;

    // adds a particle
    void add_particle(const parcl<double> &,
                      const parcl<double> &,
                      const double);

    void add_particle(const parcl<double> &,
                      const parcl<double> &,
                      const double,
                      const double);

    // sets the hydrodynamics coupling
    template<class HItype>
    void set_coupling();

    // sets a pair force
    template<class Ftype>
    void add_pair_force(const size_t &n,
                        const parcl<double> d,
                        vector<double> &p,
                        const double s,
                        const parcl<double> z);

    // sets a external force
    template<class Ftype>
    void add_ext_force(const size_t &n,
                       const parcl<double> d,
                       vector<double> &p,
                       const double s,
                       const parcl<double> z);

    const parcl<double> operator()(const size_t, const size_t) const;

    parcl<double> operator()(const size_t, const size_t);

    // 'Herzstück' implementing the Ermak-McCammon algorithm
    void iterate();

    void copy(const size_t, boost::shared_ptr<vector<parcl<double> > >);

    void copy(const size_t, vector<parcl<double> > &);

    void reset();

    void reset(vector<parcl<double> > *);

    void print(unsigned, unsigned, unsigned);

protected:

    virtual bool init_input(ArgParse &input);

    const Scale scale;

    const size_t _dimension;
    const parcl<double> _box;
    const parcl<double> _freeze;

private:

    const double _kB;
    const double _timestep;
    const double _temp;
    const unsigned _prec;
    const double _actio;
    const double _reactio;

    vector<particle_ptr> lhs;
    couple_ptr cpl;

    force_ptr rhs_random;
    vector<vector<force_ptr> > rhs_pair;
    vector<vector<force_ptr> > rhs_extern;

    vector<vector<double> > state_extern;
};

////////////////////////////////////////////////////////////////////////////
///////////////////////// POCKET-LIGAND SYSTEMS ////////////////////////////
////////////////////////////////////////////////////////////////////////////

//////////////// class Single2D : class System /////////////////////

class Single2D : public System {
public:
    Single2D(ArgParse &input);

    ~Single2D() {};

protected:

    bool init_input(ArgParse &input);

private:

    bool started;
};

#endif /* SYSTEM_H_ */
