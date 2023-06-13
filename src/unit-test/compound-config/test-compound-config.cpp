#define BOOST_TEST_MODULE TestCompoundConfig
#include <boost/test/included/unit_test.hpp>

// number of testing cycles to run
int TESTS = 2000;
// the seed for the entropy source
uint SEED = 42;
// changes the max random value to the U_LONG32 random value
#define RAND_MAX = ULONG_LONG_MAX;

// tests the lookup functions when reading in from file
BOOST_AUTO_TEST_CASE(testLookups)
{
    // creates the entropy source
    srand(SEED);
    // tracks the test cases
    std::vector<long long> tests;
    for (int test = 0; test < TESTS; test++)
    {
        tests.push_back(rand());
    }
}

// tests the ability to set correctly
BOOST_AUTO_TEST_CASE(testDynamics)
{

}

// tests the ability to read out correctly from sets
BOOST_AUTO_TEST_CASE(testLookupDynamics)
{

}