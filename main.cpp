#include <iostream>

#include "src/ArgParse.h"
#include "src/Particle.h"
#include "src/System.h"

int main(int argc, char *argv[]) {

    ArgParse &input(*new ArgParse(argc, argv));

    System &system(*new System(input));

    const double delta_time(input.flag<double>("-dt", 1e-2));
    const size_t nt(static_cast<size_t>( input.flag<double>("-t", 100.0) / delta_time ));
    const size_t freq(static_cast<size_t>( input.flag<double>("-of", 0.1) / delta_time ));
    const size_t pckg_size(1 * freq);

    input.insert_flag<double>("-lig", "--ligand", 3.0, 1);
    input.insert_flag<double>("-a1", "--radius1", .4, 1);

    parcl<double> tmp(2, 1.);

    // The standard Single2D system contains one particles.
    parcl<double> ligand(2, input.flag<double>("-lig"));
    system.add_particle(ligand, tmp, input.flag<double>("-a1"));

    vector<double> params(0);
    system.add_ext_force< GaussEnvelope >(0, tmp, params,  1., 0.0);

    // Hydrodynamic interaction off
    system.set_coupling<NoHI>();

    for (unsigned t = 0; t <= nt; t++) {
        system.iterate();

        if ((t > 0) && (t % pckg_size == 0)) {
            system.print(t, freq, pckg_size);
        }
    }

    delete &input;
    delete &system;

    return 0;
}
