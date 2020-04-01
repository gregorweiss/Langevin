#ifndef HICOUPLE_H_
#define HICOUPLE_H_

#include <stdlib.h>
#include <vector>
#include <math.h>
#include <boost/shared_ptr.hpp>

#include "common.h"
#include "Particle.h"

using std::vector;

class HIcouple {
public:
    HIcouple() :
            scale(),
            myPi(3.141592654),
            kB(scale.kB),
            temp(300. / scale.T),
            eta(1e-3 * scale.lambda * scale.tau / scale.kBT),
            kBToverPiEta(scale.kB * temp / (myPi * eta)),
            Dtensor(vector<vector<double> >()),
            Stensor(vector<vector<double> >()) {};

    HIcouple(vector<particle_ptr> prtcls, const double T) :
            scale(),
            myPi(3.141592654),
            kB(scale.kB),
            temp(T),
            eta(1e-3 * scale.lambda * scale.tau / scale.kBT),
            kBToverPiEta(kB * temp / (myPi * eta)),
            Dtensor(gen_symtensor(prtcls.size())),
            Stensor(gen_symtensor(prtcls.size())) {
        size_t nop(prtcls.size());
        for (size_t i = 0; i < nop; i++) {
            Dtensor[i][i] = prtcls.at(i)->diffusivity();
        }
    };

    virtual ~HIcouple() {};

    void gen(vector<particle_ptr> prtcls);

    double operator()(size_t i, size_t j);

    double sigma(size_t i, size_t j);

    double T();

protected:

    virtual double func(particle_ptr, particle_ptr) = 0;

    const Scale scale;

    const double myPi;

    const double kB;
    const double temp;
    const double eta;
    const double kBToverPiEta;

private:

    void gen_Stensor();

    vector<vector<double> > gen_symtensor(size_t);

    vector<vector<double> > Dtensor;
    vector<vector<double> > Stensor;
};

class EmptyHI : public HIcouple {
public:
    EmptyHI() : HIcouple() {}

    ~EmptyHI() {}

protected:
    double func(particle_ptr p_i, particle_ptr p_j);
};

class NoHI : public HIcouple {
public:
    NoHI(vector<particle_ptr> prtcls, const double T) : HIcouple(prtcls, T) {};

    ~NoHI() {};

protected:

    double func(particle_ptr, particle_ptr);
};

class RPY : public HIcouple {
public:
    RPY(vector<particle_ptr> &prtcls, const double T) : HIcouple(prtcls, T) {};

    ~RPY() {};

protected:

    double func(particle_ptr, particle_ptr);
};

typedef boost::shared_ptr<HIcouple> couple_ptr;

#endif /* HICOUPLE_H_ */
