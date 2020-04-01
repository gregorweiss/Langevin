#ifndef COMMON_H_
#define COMMON_H_

#include <cstdlib>

#include <cstdio>
#include <iostream>
#include <ostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstring>

#include <cmath>

#include <vector>
#include <iterator>
#include <boost/scoped_ptr.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/range/algorithm/remove_if.hpp>

using namespace std;

#include <Eigen/Dense>

class Mvn {
public:
    Mvn(const Eigen::VectorXd &mu,
        const Eigen::MatrixXd &s) : mean(mu),
                                    sigma(s),
                                    n(mu.rows()),
                                    sqrt2pi(std::sqrt(2 * M_PI)) {}

    ~Mvn() {}

    double form(const Eigen::VectorXd &x) {
        double quadform = (x - mean).transpose() *
                          sigma.inverse() * (x - mean);

        return exp(-0.5 * quadform);
    }

    double pdf(const Eigen::VectorXd &x) {
        double norm = std::pow(sqrt2pi, -0.5 * n) *
                      std::pow(sigma.determinant(), -0.5);

        return norm * form(x);
    }

    Eigen::VectorXd gradient(const Eigen::VectorXd &x) {
        return form(x) * sigma.inverse() * (x - mean);
    }

private:

    const Eigen::VectorXd mean;
    const Eigen::MatrixXd sigma;

    const double n;
    const double sqrt2pi;
};

////////////////////////////////////////////////////////////////////////////////
/////////////////////////// String Manipulation ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

const std::string trim(const std::string &pString);

const std::string trim2(const std::string &pString);

////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// Linear Algebra /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template<typename T>
const std::vector<std::vector<T> > transpose(const std::vector<std::vector<T> > &);

////////////////////////////////////////////////////////////////////////////////
//////////////////////////// Numerical Integration /////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template<typename T>
struct RulePtr {
    typedef const T (*type)(const typename std::vector<T>::iterator,
                            const typename std::vector<T>::iterator,
                            const T);
};

template<typename T>
const T Trapez(const typename std::vector<T>::iterator,
               const typename std::vector<T>::iterator,
               const T);

template<typename T>
const T Simpson(const typename std::vector<T>::iterator,
                const typename std::vector<T>::iterator,
                const T);

template<typename T>
const T integrate(typename RulePtr<T>::type,
                  std::vector<T> &,
                  const T);

template<typename T>
const T integrate(typename RulePtr<T>::type,
                  const typename std::vector<T>::iterator,
                  const typename std::vector<T>::iterator,
                  const T);


template<typename T>
std::vector<T> cumulative(typename RulePtr<T>::type,
                          const typename std::vector<T>::iterator,
                          const typename std::vector<T>::iterator,
                          const T);

template<typename T>
bool cumulative(typename RulePtr<T>::type,
                const typename std::vector<T>::iterator,
                const typename std::vector<T>::iterator,
                const T,
                std::vector<T> &res);


////////////////////////////////////////////////////////////////////////////////
//////////////////////////// template class parcl //////////////////////////////
////////////////////////////////////////////////////////////////////////////////


template<class PT>
class parcl {
public:

    parcl() :
            _q(new typename std::vector<PT>(3, 0.)) {};

    parcl(const typename std::vector<PT> &vec) :
            _q(new typename std::vector<PT>(vec)) {};

    parcl(const size_t S) :
            _q(new typename std::vector<PT>(S, 0.)) {};

    parcl(const size_t S, const PT V) :
            _q(new typename std::vector<PT>(S, V)) {};

    parcl(const parcl<PT> &other) :
            _q(new typename std::vector<PT>(other.size(), 0.)) { *this = other; };

    ~parcl() {};


    const PT &operator[](const size_t) const;

    PT &operator[](const size_t);

    const PT &at(const size_t) const;

    PT &at(const size_t);

    size_t size() const;

    size_t size();

    typename std::vector<PT>::iterator begin();

    typename std::vector<PT>::const_iterator begin() const;

    typename std::vector<PT>::iterator end();

    typename std::vector<PT>::const_iterator end() const;

    parcl<PT> &operator=(const parcl<PT> &);

    const parcl<PT> &operator+=(const parcl<PT> &) const;

    parcl<PT> &operator+=(const parcl<PT> &);

    parcl<PT> &operator-=(const parcl<PT> &);

    parcl<PT> &operator*=(const PT &);

    parcl<PT> &operator/=(const PT &);

    parcl<PT> &operator/=(const parcl<PT> &);

    const parcl<PT> operator+(const parcl<PT> &) const;

    const parcl<PT> operator-(const parcl<PT> &) const;

    const parcl<PT> operator*(const PT &) const;

    bool operator==(const parcl<PT> &);

    bool operator!=(const parcl<PT> &);

    const PT operator*(const parcl<PT> &) const; // Scalar Product
    const parcl<PT> transpose(const parcl<PT> &) const; // Matrix Product
    const parcl<PT> transpose(const typename std::vector<PT> &) const; // Matrix Product

    const parcl<PT> operator/(const parcl<PT> &) const; // Normalization

private:

    boost::scoped_ptr<typename std::vector<PT> > _q;

};

////////////////////////////////////////////////////////////////////////////////
//////////////////////////// operations on parcls //////////////////////////////
////////////////////////////////////////////////////////////////////////////////

std::ostream &operator<<(std::ostream &, const parcl<double> &);

template<typename T>
const T parcl_dist(const parcl<T> &,
                   const parcl<T> &);

template<typename T>
const T parcl_dist(const parcl<T> &,
                   const parcl<T> &,
                   parcl<T> *);

template<typename T>
int sgn(T val);

template<typename T>
void crossproduct(T vec1, T vec2, T &res);

template<typename T, typename P>
void crossproduct(T vec1, T vec2, P &res);

template<typename T, typename P>
void subtract(T vec1, T vec2, P res);

template<typename T>
void geocenter(std::vector<T> vecs, T &res);

template<typename T>
T scalarproduct(parcl<T> &vec1, parcl<T> &vec2);

template<typename T>
T norm(parcl<T> &vec);

/////////////////////////////////////////////////////////////////////////
//////////////////// template class parcl ///////////////////////////////
/////////////////////////////////////////////////////////////////////////

template<class PT>
const PT &parcl<PT>::operator[](const size_t n) const { return _q->operator[](n); }

template<class PT>
PT &parcl<PT>::operator[](const size_t n) { return _q->operator[](n); }

template<class PT>
const PT &parcl<PT>::at(const size_t n) const { return _q->at(n); }

template<class PT>
PT &parcl<PT>::at(const size_t n) { return _q->at(n); }

template<class PT>
size_t parcl<PT>::size() const { return this->_q->size(); }

template<class PT>
size_t parcl<PT>::size() { return this->_q->size(); }

template<class PT>
typename vector<PT>::iterator parcl<PT>::begin() { return _q->begin(); };

template<class PT>
typename vector<PT>::const_iterator parcl<PT>::begin() const { return _q->begin(); };

template<class PT>
typename vector<PT>::iterator parcl<PT>::end() { return _q->end(); };

template<class PT>
typename vector<PT>::const_iterator parcl<PT>::end() const { return _q->end(); };

// Overloading assignement operator
// Take const-reference of rhs & Return non-const reference of lhs
// Note: no need for check to self-assignment, because no memory is deallocated.
template<class PT>
parcl<PT> &parcl<PT>::operator=(const parcl<PT> &right) {
    if (right.size() != this->size()) { return *this; }

    const typename vector<PT>::iterator lhs_end(this->_q->end());
    const typename vector<PT>::const_iterator rhs_end(right.end());

    typename vector<PT>::iterator lhs(this->_q->begin());
    typename vector<PT>::const_iterator rhs(right.begin());

    while ((lhs != lhs_end) && (rhs != rhs_end)) {
        *lhs = *rhs;

        lhs++;
        rhs++;
    }

    return *this;
}


template<class PT>
const parcl<PT> &parcl<PT>::operator+=(const parcl<PT> &right) const {
    if (right.size() != this->size()) { return *this; }

    const typename vector<PT>::iterator lhs_end(this->_q->end());
    const typename vector<PT>::const_iterator rhs_end(right.end());

    typename vector<PT>::iterator lhs(this->_q->begin());
    typename vector<PT>::const_iterator rhs(right.begin());

    while ((lhs != lhs_end) && (rhs != rhs_end)) {
        *lhs += *rhs;

        lhs++;
        rhs++;
    }

    return *this;
}


template<class PT>
parcl<PT> &parcl<PT>::operator+=(const parcl<PT> &right) {
    if (right.size() != this->size()) { return *this; }

    const typename vector<PT>::iterator lhs_end(this->_q->end());
    const typename vector<PT>::const_iterator rhs_end(right.end());

    typename vector<PT>::iterator lhs(this->_q->begin());
    typename vector<PT>::const_iterator rhs(right.begin());

    while ((lhs != lhs_end) && (rhs != rhs_end)) {
        *lhs += *rhs;

        lhs++;
        rhs++;
    }

    return *this;
}

template<class PT>
parcl<PT> &parcl<PT>::operator-=(const parcl<PT> &right) { return (*this += (right * static_cast<PT>(-1.))); }


template<class PT>
parcl<PT> &parcl<PT>::operator*=(const PT &right) {
    const typename vector<PT>::iterator lhs_end(this->_q->end());
    typename vector<PT>::iterator lhs(this->_q->begin());

    while (lhs != lhs_end) {
        *lhs *= right;

        lhs++;
    }

    return *this;
}

template<class PT>
parcl<PT> &parcl<PT>::operator/=(const PT &right) {
    const typename vector<PT>::iterator lhs_end(this->_q->end());
    typename vector<PT>::iterator lhs(this->_q->begin());

    while (lhs != lhs_end) {
        *lhs /= right;

        lhs++;
    }

    return *this;
}

/*
* brief: this is also for normalization
*/
template<class PT>
parcl<PT> &parcl<PT>::operator/=(const parcl<PT> &right) {
    const typename vector<PT>::const_iterator rhs_end(right.end());
    typename vector<PT>::const_iterator rhs(right.begin());

    PT norm(static_cast<PT>(0.));
    while ((rhs != rhs_end)) {
        norm += pow(*rhs, 2.);

        rhs++;
    }

    return *this /= sqrt(norm);
}

template<class PT>
const parcl<PT> parcl<PT>::operator+(const parcl<PT> &right) const { return parcl<PT>(*this) += right; }

template<class PT>
const parcl<PT> parcl<PT>::operator-(const parcl<PT> &right) const { return parcl<PT>(*this) -= right; }

template<class PT>
const parcl<PT> parcl<PT>::operator*(const PT &right) const { return parcl<PT>(*this) *= right; }


// Scalar Product
template<class PT>
const PT parcl<PT>::operator*(const parcl<PT> &right) const {
    if (right.size() != this->size()) { return 0.; }

    const typename vector<PT>::iterator lhs_end(this->_q->end());
    const typename vector<PT>::const_iterator rhs_end(right.end());

    typename vector<PT>::iterator lhs(this->_q->begin());
    typename vector<PT>::const_iterator rhs(right.begin());

    PT ret(static_cast<PT>(0.));
    while ((lhs != lhs_end) && (rhs != rhs_end)) {
        ret += (*lhs) * (*rhs);

        lhs++;
        rhs++;
    }

    return ret;
}

template<class PT>
bool parcl<PT>::operator==(const parcl<PT> &right) {
    if (right.size() != this->size()) { return false; }

    const typename vector<PT>::iterator lhs_end(this->_q->end());
    const typename vector<PT>::const_iterator rhs_end(right.end());

    typename vector<PT>::iterator lhs(this->_q->begin());
    typename vector<PT>::const_iterator rhs(right.begin());

    int count(0);
    while ((lhs != lhs_end) && (rhs != rhs_end) && (*lhs == *rhs)) {
        count++;
        lhs++;
        rhs++;
    }

    return (count == 3);
}

template<class PT>
bool parcl<PT>::operator!=(const parcl<PT> &right) { return !(*this == right); }


template<class PT>
const parcl<PT> parcl<PT>::transpose(const parcl<PT> &right) const {
    if (right.size() != this->size()) { return *this; }

    const typename vector<PT>::iterator lhs_end(this->_q->end());
    const typename vector<PT>::const_iterator rhs_end(right.end());

    typename vector<PT>::iterator lhs(this->_q->begin());
    typename vector<PT>::const_iterator rhs(right.begin());

    parcl<PT> ret(this->size());

    size_t d(0);
    while ((lhs != lhs_end) && (rhs != rhs_end)) {
        ret[d] = ((*lhs) * (*rhs));

        lhs++;
        rhs++;
        d++;
    }

    return parcl<PT>(ret);
}

template<class PT>
const parcl<PT> parcl<PT>::transpose(const vector<PT> &right) const { return parcl<PT>(right).transpose(*this); }

template<class PT>
const parcl<PT> parcl<PT>::operator/(const parcl<PT> &right) const { return parcl<PT>(*this) /= right; }

/////////////////////////////////////////////////////////////////////////

template<typename T>
const T parcl_dist(const parcl<T> &vec1, const parcl<T> &vec2) {
    parcl<T> ret(vec1 - vec2);

    T dist = ret * ret;
    dist = sqrt(dist);

    return dist;
}

template<typename T>
const T parcl_dist(const parcl<T> &vec1, const parcl<T> &vec2, parcl<double> *retvec) {
    *retvec = vec1 - vec2;

    T dist = (*retvec) * (*retvec);
    dist = sqrt(dist);

    *retvec /= *retvec;

    return dist;
}

template<typename T>
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

template<typename T>
void crossproduct(T vec1, T vec2, T &res) {
    res[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
    res[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
    res[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
}

template<typename T, typename P>
void crossproduct(T vec1, T vec2, P &res) {
    res[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
    res[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
    res[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
}

template<typename T, typename P>
void subtract(T vec1, T vec2, P res) {
    res[0] = vec1[0] - vec2[0];
    res[1] = vec1[1] - vec2[1];
    res[2] = vec1[2] - vec2[2];
}

template<typename T>
void geocenter(std::vector<T> vecs, T &res) {

    size_t number(vecs.size());

    for (size_t n = 0; n < number; n++) { res += vecs[n]; }

    res /= number;
}

template<typename T>
T scalarproduct(parcl<T> &vec1, parcl<T> &vec2) {

    T val = 0;

    val += vec1[0] * vec2[0];
    val += vec1[1] * vec2[1];
    val += vec1[2] * vec2[2];

    return val;
}


template<typename T>
T norm(parcl<T> &vec) {

    T abs = 0;

    for (size_t dim = 0; dim < 3; dim++) {
        abs += vec[dim] * vec[dim];
    }

    return sqrt(abs);
}

#endif /* COMMON_H_ */
