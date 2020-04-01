
#include "HIcouple.h"

vector<vector<double>> HIcouple::gen_symtensor(size_t dim) {
    vector<vector<double>> ret(0);
    for (size_t j = 0; j < dim; j++) {
        ret.push_back(vector<double>(j + 1, 0.));
    }

    return ret;
}

void HIcouple::gen(vector<particle_ptr> prtcls) {
    for (size_t i = 0; i < prtcls.size(); i++) {
        for (size_t j = 0; j < i + 1; j++) {
            Dtensor[i][j] = kBToverPiEta * func(prtcls[i], prtcls[j]);
        }
    }

    gen_Stensor();
}


double HIcouple::operator()(size_t i, size_t j) {
    if (i >= j) {

        return Dtensor[i][j];

    } else if (i < j) {

        return Dtensor[j][i];

    }

    return 0.;
}

double HIcouple::sigma(size_t i, size_t j) {
    double ret(0.);
    if (i <= j) { ret = Stensor[i][j]; }

    return ret;
}

double HIcouple::T() { return temp; }

void HIcouple::gen_Stensor() {
    int nop(Dtensor.size());;
    for (int i = 0; i < nop; i++) {
        double sigma_ii = Dtensor[i][i];
        for (int k = 0; k < i - 1; k++) {
            sigma_ii -= Stensor[i][k] * Stensor[i][k];
        }
        Stensor[i][i] = sqrt(sigma_ii);

        for (int j = 0; j < i; j++) {
            double sigma_ij = Dtensor[i][j];
            for (int k = 0; k < j - 1; k++) {
                sigma_ij -= Stensor[i][k] * Stensor[j][k];
            }

            sigma_ij /= Stensor[j][j];

            Stensor[i][j] = sigma_ij;
        }
    }
}

double EmptyHI::func(particle_ptr p_i, particle_ptr p_j) {
    return 0.;
}

double NoHI::func(particle_ptr p_i, particle_ptr p_j) {
    if (p_i == p_j) {
        return p_i->diffusivity() / kBToverPiEta;

    } else {
        return 0.;

    }
}

double RPY::func(particle_ptr p_i, particle_ptr p_j) {
    double ai_plus_aj(p_i->radius() + p_j->radius());
    double ai_diff_aj(sqrt(pow(p_i->radius()
                               - p_j->radius(), 2.)));
    parcl<double> dr(p_i->last() - p_j->last());
    double dist(dr * dr);
    dist = sqrt(dist);

    double ret(0.);
    if (ai_plus_aj < dist) {

        double ai2(p_i->radius() * p_i->radius());
        double aj2(p_j->radius() * p_j->radius());

        double ratio(ai2 + aj2 / dist);

        double ratio_plus(1. + ratio / 3.);
        double ratio_minus(1. - ratio);

        ret = ratio_plus + ratio_minus * dist;
        ret /= 8. * myPi * dist;

    } else if ((ai_diff_aj < dist) && (dist <= ai_plus_aj)) {

        double first(16. * pow(dist, 3) * ai_plus_aj);
        first -= pow(pow(ai_diff_aj, 2) + 3. * pow(dist, 2), 2);
        first /= 32. * pow(dist, 3);

        double second(3 * pow(pow(ai_diff_aj, 2) - pow(dist, 2), 2));
        second /= 32. * pow(dist, 3);

        ret = first + second * dist;

    } else if (dist <= ai_diff_aj) {
        if (p_i->radius() >= p_j->radius()) {
            ret = 1. / (6. * myPi * p_i->radius());
        } else if (p_i->radius() < p_j->radius()) {
            ret = 1. / (6. * myPi * p_j->radius());
        }
    }

    return ret;
}


