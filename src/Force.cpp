
#include "Force.h"

//////////////////////////////////////////////////////////////////////////
///////////////////// 		EXTERNAL FORCES		//////////////////////////
//////////////////////////////////////////////////////////////////////////

/////////////////////// class Random : class Force ///////////////////////

parcl<double> Random::operator()(vector<parcl<double> > &vec) {
    parcl<double> ret(vec.at(0).size(), 0.);

    const vector<double>::iterator it_end(ret.end());
    vector<double>::iterator it(ret.begin());

    while (it != it_end) {
        *it = gauss();
        it++;
    }

    return ret;
}

//////////////////// class DoubleWell : class Force //////////////////////

parcl<double> DoubleWell::operator()(vector<parcl<double> > &vec) {
    parcl<double> dr(vec.at(0) - this->a_half);
    dr = dr.transpose(this->get_dir());

    double dist(dr * dr);
    dist = sqrt(dist);

    double ret(0.);
    if (dist != 0.0) {
        dr /= dr;
        ret = (-1.) * this->pre * (dist * dist * dist - (this->bsqr * dist));
    } else {
        dr = this->get_dir();
    }

    return dr * ret;
}

//////////////////// class Bias : class Force //////////////////////

parcl<double> Bias::operator()(vector<parcl<double> > &vec) {
    return this->get_dir() * _forceconstant;
}

//////////////////// class GaussEnvelope : class Force //////////////////////

parcl<double> GaussEnvelope::operator()(vector<parcl<double> > &vec) {

    Eigen::VectorXd x(2);
    x(0) = vec[0][0];
    x(1) = vec[0][1];

    double denominator(0.0);
    Eigen::VectorXd ret(2);
    for (auto i = 0; i < this->_c.size(); ++i) {
        denominator += _c[i] * _mvn[i].form(x);
        ret -= _c[i] * _mvn[i].gradient(x);
    }
    ret /= denominator;

    parcl<double> dr(2, 0.0);
    dr[0] = ret(0);
    dr[1] = ret(1);

    return dr;
}


//////////////////////////////////////////////////////////////////////
//////////////////////// Update Functions ////////////////////////////
//////////////////////////////////////////////////////////////////////

// Computes the positional shift dr from the pair forces including hydrodynamic coupling
const parcl<double> update_pair(const particle_ptr &p_i,
                                const particle_ptr &p_j,
                                const vector<force_ptr> &frcs,
                                const double &dt) {
    const size_t dimension(p_i->last().size());
    parcl<double> dr(dimension, 0.);

    for (force_ptr f : frcs) { dr += (f->operator()(p_i, p_j) * dt); }

    return dr;
}

// Computes the external forces on one particle including its mobility
const parcl<double> update_extern(const particle_ptr &p_i,
                                  const vector<force_ptr> &frcs,
                                  const vector<double> &states,
                                  const double &dt) {
    const size_t dimension(p_i->last().size());
    parcl<double> dr(dimension, 0.);

    for (size_t i = 0; i < frcs.size(); i++) { dr += (frcs.at(i)->operator()(p_i) * dt * states.at(i)); }

    return dr;
}

// Computes the random forces on one particle
const parcl<double> update_random(const particle_ptr &p_i,
                                  const force_ptr &f) {
    const size_t dimension(p_i->last().size());
    parcl<double> dr(dimension, 0.);

    dr += f->operator()(p_i);

    return dr;
}

