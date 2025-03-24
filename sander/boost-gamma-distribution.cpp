#include <boost/random/gamma_distribution.hpp>

/// non-module alias for amrand (defined in resamplekin.f90)
extern "C" void amrand_c_alias(double *out);

namespace {

/***
 * Wraps default amber random engine
 *
 * Implements interface of boost Uniform Random Number Generator
 * https://www.boost.org/doc/libs/1_66_0/doc/html/boost_random/reference.html#boost_random.reference.concepts.uniform_random_number_generator
 * */
class AmrandGenerator {
public:
    typedef double result_type;

    result_type operator()() {
        double x;
        amrand_c_alias(&x);
        return x;
    }

    result_type min() const { return 0.0; }

    result_type max() const { return 1.0; }
};

}
extern "C" {
/**
 * Returns random value distributed as gamma_distribution of order \a rank
 * Uses default Amber random generator (amrand)
 * */
double gamdev(int rank) {
    typedef boost::random::gamma_distribution<double>::param_type param_type;
    static boost::random::gamma_distribution<double> gamma_distribution;
    static AmrandGenerator urng;
    return gamma_distribution(urng, param_type(rank, 1));
}
};
