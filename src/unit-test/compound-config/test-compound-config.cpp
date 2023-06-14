#define BOOST_TEST_MODULE TestCompoundConfig
#include <boost/test/included/unit_test.hpp>
#include <compound-config/compound-config.hpp>

// number of testing cycles to run
int TESTS = 2000;
// the seed for the entropy source
uint SEED = 42;
// changes the max random value to the U_LONG32 random value
#undef RAND_MAX
#define RAND_MAX = ULONG_LONG_MAX

// static YAML file names we want to load in for the test
std::vector<std::string> FILES = {
    "example.yaml" // include source here if pulled from a repo somewhere
};


// we are only testing things in config
namespace config {
// tests the lookup functions when reading in from file
BOOST_AUTO_TEST_CASE(testLookups)
{
    for (auto FILE:FILES) 
    {
        // reads the YAML file into CompoundConfig
        CompoundConfig cConfig = CompoundConfig(FILE, "yaml");
        // reads in the YAML file independently of CompoundConfig
        YAML::Node testRef = YAML::LoadFile(FILE);

        // gets the root CompoundConfigNode
        config::CompoundConfigNode root = cConfig.getRoot();

        // goes through all keys and compares the values.
        // figure this out. you may need a helper file.
        for (auto node: testRef)
        {
            // unpacks the key + expected value
            std::string key, expected;
            key = node.first.as<std::string>();
            expected = node.second.as<std::string>();

            // unpacks actual value
            std::string actual;
            root.lookupValue(key, actual);

            // compares to the config
            BOOST_CHECK_EQUAL(actual, expected);
        }
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
} // namespace config