////////////////////////////////////////////////////////////////////////////
/////////////////////////// MOTHER class System ////////////////////////////
////////////////////////////////////////////////////////////////////////////

#include "System.h"

// This one is the core of System. It implements the
// ERMAK-McCAMMON ALGORITHM which is capable of intergrating brownian
// dynamics simulations WITH hydrodynamic interactions.
// First the hydrodynamic interaction are calculated 'cpl->gen( lhs )' 
// with the current particle positions 'lhs' the left-hand-side of 
// the Ermak-McCammon equation.
// Secondly the forces are calcluted within the FIRST LOOP and 
// SECOND LOOP where the force vector and random force vector are
// being filled.
// The THIRD LOOP calculates each summand of the Ermak-McCammon
// equation.
// So that the FORTH LOOP can calculate the new postitions.
void System::iterate() {
    size_t nop(lhs.size());

    cpl->gen(lhs);

    vector<parcl<double> > &random_dt(*new vector<parcl<double> >(nop, parcl<double>(_dimension, 0.)));

    vector<parcl<double> > &force_dt(*new vector<parcl<double> >(nop, parcl<double>(_dimension, 0.)));

    vector<vector<parcl<double> > > &forcemat_dt(
            *new vector<vector<parcl<double> > >(nop, vector<parcl<double> >(nop, parcl<double>(_dimension, 0.))));

    // FIRST LOOP
    // #pragma omp parallel for
    for (size_t i = 0; i < nop; i++) {
        random_dt.at(i) += update_random(lhs.at(i), rhs_random);

        force_dt.at(i) += update_extern(lhs.at(i), rhs_extern.at(i), state_extern.at(i), this->_timestep);

        for (size_t j = i + 1; j < nop; j++) {
            forcemat_dt.at(i).at(j) += update_pair(lhs.at(i), lhs.at(j), rhs_pair.at(i), this->_timestep);

            // Newton's Third Law "actio-reactio"
            // Variable _actio allows to turn off/neglect Newton's thrid law
            // But be careful, first of all due to PHYSICALTIY and
            // second of all because then the order of the particle in lhs
            // plays a crucial role ! Physical case for which this is approximatly
            // true is if the masses are very different.
            force_dt.at(i) += forcemat_dt.at(i).at(j) * _actio;
        }
    }

    // SECOND LOOP
    // #pragma omp parallel for
    for (size_t i = 0; i < nop; i++) {
        // Newton's Third Law "actio-reactio"
        // Variable _reactio allows to turn off/neglect Newton's thrid law
        // But be careful, first of all due to PHYSICALTIY and
        // second of all because then the order of the particle in lhs
        // plays a crucial role ! Physical case for which this is approximatly
        // true is if the masses are very different.
        for (int j = i - 1; j >= 0; j--) {
            force_dt.at(i) -= forcemat_dt.at(j).at(i) * _reactio;
        }
    }
    delete &forcemat_dt;

    vector<parcl<double> > &delta_dt(*new vector<parcl<double> >(nop, parcl<double>(_dimension, 0.)));

    // THIRD LOOP
    // #pragma omp parallel for
    for (size_t i = 0; i < nop; i++) {
        for (size_t j = 0; j <= i; j++) {
            double D_ij(cpl->operator()(i, j));
            delta_dt.at(i) += (force_dt.at(j) * (D_ij / (_kB * _temp)));

            double sigma_ij(cpl->sigma(i, j));
            delta_dt.at(i) += (random_dt.at(j) * sigma_ij);
        }

        for (size_t j = i+1; j < nop; j++) {
            parcl<double> D_ij(cpl->operator()(i, j));
            delta_dt.at(i) += (force_dt.at(j) * (D_ij / (_kB * _temp)));
        }
    }
    force_dt.clear();
    delete &force_dt;
    random_dt.clear();
    delete &random_dt;

    // FORTH LOOP
    // #pragma omp parallel for
    for (size_t i = 0; i < nop; i++) {
        parcl<double> r_t(lhs.at(i)->last());
        parcl<double> r_tdt(r_t + delta_dt.at(i));

        // periodic boundary condition
//		for ( size_t dim = 0; dim < _dimension; dim++ )
//		{	r_tdt[dim] -= floor( r_tdt[dim] / _box[dim]) * _box[dim]; }

        // reflective boundary condition
//		for ( size_t dim = 0; dim < _dimension; dim++ )
//		{	r_tdt[dim] -= floor( r_tdt[dim] / _box[dim] ) * ( 2.0 * ( r_tdt[dim] - _box[dim] ) ); }
//		{	r_tdt[dim] -= floor( r_tdt[dim] / _box[dim] ) * ( 2.0 * ( r_tdt[dim] - _box[dim] ) ); }

//		if ( delta_dt.at(i).at(2) > 0.2 )
//		{ std::cout << "To large timestep " << i << " " << delta_dt.at(i).at(2) << std::endl; }

        lhs.at(i)->append(r_tdt);
        lhs.at(i)->increment();
    }

    delete &delta_dt;
}

void System::copy(const size_t particle, boost::shared_ptr<vector<parcl<double> > > vec) {
    lhs.at(particle)->copy(vec);
}

void System::copy(const size_t particle, vector<parcl<double> > &vec) {
    lhs.at(particle)->copy(vec);
}

// reset() clears the system and frees the memory except for the last step
void System::reset() {
    BOOST_FOREACH(const particle_ptr &prtcl, lhs) {
                    prtcl->clear();
                }
}

// Does the same as reset() but inserts another coordinate as last step
void System::reset(vector<parcl<double> > *pos) {
    size_t n(0);
    BOOST_FOREACH(const particle_ptr &prtcl, lhs) {
                    prtcl->append(pos->at(n));
                    prtcl->increment();
                    prtcl->clear();

                    n++;
                }
}

// Prints the system's timeline.
void System::print(unsigned t,
                   unsigned freq,
                   unsigned nc) {
    for (unsigned dt = 0; dt < nc / freq; dt++) {

        const unsigned current_t((t - nc) + dt * freq);

        if ((current_t / 10) > 6) {

            std::cout << std::setprecision(_prec) << current_t * _timestep;

        } else {

            std::cout << current_t * _timestep;

        }

        BOOST_FOREACH (const particle_ptr &prtcl, lhs) {
                        for (size_t d = 0; d < _dimension; d++) {
                            std::cout << " " << prtcl->operator()(dt * freq).at(d);
                        }
                    }

        std::cout << std::endl;
    }

    BOOST_FOREACH (const particle_ptr &prtcl, lhs) {
                    prtcl->clear();
                }
}

// Adds a particle at th initial position 'pos' with hydrodynamic radius 'size'
void System::add_particle(const parcl<double> &pos,
                          const parcl<double> &frz,
                          const double size) {
    const particle_ptr prtcl(new Particle(pos, frz, size, _temp));
    lhs.push_back(prtcl);

    rhs_extern.push_back(vector<force_ptr>());
    rhs_pair.push_back(vector<force_ptr>());
    state_extern.push_back(vector<double>());
}

// Adds a particle at th initial position 'pos' with hydrodynamic radius 'size'
void System::add_particle(const parcl<double> &pos,
                          const parcl<double> &frz,
                          const double size,
                          const double temp) {
    const particle_ptr prtcl(new Particle(pos, frz, size, temp));
    lhs.push_back(prtcl);

    rhs_extern.push_back(vector<force_ptr>());
    rhs_pair.push_back(vector<force_ptr>());
    state_extern.push_back(vector<double>());
}

// Set the hydrodynamic interaction
template<class HItype>
void System::set_coupling() {
    cpl.reset(new HItype(lhs, _temp));
}

// Adds an external force of type 'Ftype' using the force parameters 'params' to 
// particle at postition 'n' in 'lhs' with defined 'sign'.
template<class Ftype>
void System::add_ext_force(const size_t &n,
                           const parcl<double> direct,
                           vector<double> &params,
                           const double sign,
                           const parcl<double> shift) {
    force_ptr frc(new Ftype(direct, params, sign, shift));
    rhs_extern.at(n).push_back(frc);
    state_extern.at(n).push_back(1.0);
}

// Adds a pair force of type 'Ftype' using the force parameters 'params' to 
// particle at postition 'n' in 'lhs' with defined 'sign'.
template<class Ftype>
void System::add_pair_force(const size_t &n,
                            const parcl<double> direct,
                            vector<double> &params,
                            const double sign,
                            const parcl<double> shift) {
    force_ptr frc(new Ftype(direct, params, sign, shift));
    rhs_pair.at(n).push_back(frc);
}

// Returns the coordinate of n'th particle at time t.
const parcl<double> System::operator()(const size_t n,
                                       const size_t t) const { return lhs.at(n)->operator()(t); }

parcl<double> System::operator()(const size_t n,
                                 const size_t t) { return lhs.at(n)->operator()(t); }


bool System::init_input(ArgParse &input) {
    return true;
}

////////////////////////////////////////////////////////////////////////////
///////////////////////// POCKET-LIGAND SYSTEMS ////////////////////////////
////////////////////////////////////////////////////////////////////////////

//////////////// class Single2D : class System /////////////////////

/*
* class Single2D sets up an child of class System. It contains a single ligand
* which is simulated in one dimension. It interacts with two types of potentials
* as external potential at four sites, such that the system is basically mirrored,
* basically representing a system with reflecting boundary in the middle. The
* external potentials are instances from class SSolvation : class Force and 
* class Niner : class Force.
*
*	Definition of user flags at runtime:
*	"-p1"		default	(double) 1.0 			// left droplet position
*	"-p2"		default	(double) 5.0 			// right (mirrored) droplet position
*	"-lig"	default (double) 3.0			// ligand position
*	"-a1" 	default (double) 0.4			// radius of (mirrored) droplet
*	"-a2" 	default (double) 0.4 			// radius of ligand
*	"-eps"	default (double) 1.0			// lennard jones strength of ligand
*	"-att"	default (double) 1.0			// parameter to tune solvation potential
*/

Single2D::Single2D(ArgParse &input) : System(input),
                                      started(init_input(input)) {
    // Direction:
    parcl<double> tmp(_dimension, 1.);

    // The standard Single2D system contains one particles.
    parcl<double> ligand(_dimension, input.flag<double>("-lig"));
    add_particle(ligand, tmp, input.flag<double>("-a1"));

    vector<double> params(0);
    add_ext_force< GaussEnvelope >(0, tmp, params,  1., 0.0);

    // Hydrodynamic interaction off
    set_coupling<NoHI>();
}

bool Single2D::init_input(ArgParse &input) {
    input.insert_flag<double>("-lig", "--ligand", 3.0, 1);    // ligand position
    input.insert_flag<double>("-a1", "--radius1", .4, 1);     // radius of particle

    return true;
}



