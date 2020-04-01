
#include "Particle.h"

/////////////////////////////////////////////////////////////////////////
////////////////////////// class Scale //////////////////////////////////
/////////////////////////////////////////////////////////////////////////

double Scale::time(double var) { return var / tau; }

double Scale::length(double var) { return var / lambda; }

double Scale::energy(double var) { return var / kBT; }

double Scale::mass(double var) { return var / mu; }

double Scale::temperature(double var) { return var / T; }

/////////////////////////////////////////////////////////////////////////
////////////////////////// class Particle ///////////////////////////////
/////////////////////////////////////////////////////////////////////////

parcl<double> Particle::operator()(size_t t) { return this->coordinates->at(t); }

parcl<double> Particle::last() { return this->prev; }

void Particle::append(parcl<double> coord) { this->coordinates->push_back(coord.transpose(_freeze)); }

void Particle::increment() {
    this->prev_it++;
    this->prev = this->coordinates->at(prev_it);
}

size_t Particle::length() { return this->coordinates->size(); }

// Is turned of at the moment due to diffusivity tuning with "internal" temperature
//double Particle::mobility()
//{ return this->mu; }

double Particle::diffusivity() { return this->D; }

double Particle::radius() { return this->R; }

double Particle::temperature() { return this->T; }

void Particle::copy(boost::shared_ptr<vector<parcl<double> > > vec) {
    for (parcl<double> &val : *coordinates) {
        vec->push_back(val);
    }
}

void Particle::copy(vector<parcl<double> > &vec) {
    for (parcl<double> &val : *coordinates) {
        vec.push_back(val);
    }
}

void Particle::clear() {
    this->prev_it = 0;

    this->coordinates.reset(new vector<parcl<double> >(1, prev));
}

