#ifndef PARTICLE_H_
#define PARTICLE_H_

#include <cstdlib>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>

#include "common.h"

using std::vector;
using boost::shared_ptr;
using boost::scoped_ptr;

///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// class Scale /////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

class Scale {
public:
    Scale() : tau(1.e-9),                                // ns = 1,
              lambda(1.e-9),                             // nm = 1,
              kBT(1.3806488e-23 * 300.),                 // kB T = 1,
              mu(kBT * pow(lambda, 2.) / pow(tau, 2.)),
              kB(1.3806488e-23 / kBT),
              T(1 / kB) {};

    ~Scale() {};

    double time(double);

    double length(double);

    double energy(double);

    double mass(double);

    double temperature(double);

    const double tau;     // time
    const double lambda;  // length
    const double kBT;     // energy
    const double mu;      // Mass
    const double kB;      // Boltzmann constant
    const double T;       // Temperature
};

///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// class Particle //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

class Particle {
public:

    Particle(const parcl<double> &init, const double radius, const double temp) :
            scale(),
            myPi(3.141592654),
            R(radius),
            T(temp / scale.T),
            eta(1.e-3 * 1.e-18 / scale.mu), // 1.e-3 is roughly the size of water viscosity at room temperature
            D(T / (6. * myPi * eta * R)),
            _freeze(3, 1.),
            coordinates(new vector<parcl<double> >(1, init)),
            prev_it(0),
            prev(coordinates->at(prev_it)) {};


    Particle(const parcl<double> &init, const parcl<double> &frz, const double radius, const double temp) :
            scale(),
            myPi(3.141592654),
            R(radius),
            T(temp / scale.T),
            eta(1.e-3 * 1.e-18 / scale.mu), // 1.e-3 is roughly the size of water viscosity at room temperature
            D(T / (6. * myPi * eta * R)),
            _freeze(frz),
            coordinates(new vector<parcl<double> >(1, init)),
            prev_it(0),
            prev(coordinates->at(prev_it)) {};

    ~Particle() {};

    parcl<double> operator()(size_t);

    parcl<double> last();

    void append(parcl<double>);

    void increment();

    size_t length();

    double diffusivity();

    double radius();

    double temperature();

    void copy(boost::shared_ptr<vector<parcl<double> > >);

    void copy(vector<parcl<double> > &);

    void clear();

private:

    Scale scale;
    const double myPi;
    const double R;
    const double T;
    const double eta;
    const double D;

    const parcl<double> _freeze;
    scoped_ptr<vector<parcl<double> > > coordinates;
    size_t prev_it;
    parcl<double> prev;
};

typedef boost::shared_ptr<Particle> particle_ptr;

#endif /* PARTICLE_H_ */
