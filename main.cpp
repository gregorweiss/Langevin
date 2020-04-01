#include <iostream>

#include "src/ArgParse.h"
#include "src/Particle.h"
#include "src/System.h"

int main(int argc, char *argv[]) {

    ArgParse &input(*new ArgParse(argc, argv));

    Single2D &system(*new Single2D(input));

    const double delta_time(input.flag<double>("-dt", 1e-2));
    const size_t nt(static_cast<size_t>( input.flag<double>("-t", 100.0) / delta_time ));
    const size_t freq(static_cast<size_t>( input.flag<double>("-of", 0.1) / delta_time ));
    const size_t pckg_size(1 * freq);

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
