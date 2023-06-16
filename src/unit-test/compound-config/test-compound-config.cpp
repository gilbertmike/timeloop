#define BOOST_TEST_MODULE TestCompoundConfig

#include <iostream>
#include <boost/test/included/unit_test.hpp>
#include <compound-config/compound-config.hpp>

// makes sure for a certain type CCN agrees with YNode
template <typename T>
bool testScalarLookup(config::CompoundConfigNode& CNode, YAML::Node& YNode, std::string& key)
{
    T expectedScalar = YNode[key].as<T>();
    T actualScalar;
    CNode.lookupValue(key, actualScalar);

    return expectedScalar == actualScalar;
}

bool testSequenceLookup(config::CompoundConfigNode& CNode, YAML::Node& YNode, std::string& key)
{
    std::vector<std::string> expectedSeq = YNode[key].as<std::vector<std::string>>();
    std::vector<std::string> actualSeq;
    CNode.lookupArrayValue(key, actualSeq);

    return expectedSeq == actualSeq;
}

// number of testing cycles to run
int TESTS = 2000;
// the seed for the entropy source
uint SEED = 42;
// changes the max random value to the U_LONG32 random value
#undef RAND_MAX
#define RAND_MAX = ULONG_LONG_MAX

// the location of the test files
std::string TEST_LOC = "../unit-test/compound-config/tests/";

// static YAML file names we want to load in for the test
std::vector<std::string> FILES = {
    "example.yaml" // include source here if pulled from a repo somewhere
};


// we are only testing things in config
namespace config {
// tests the lookup functions when reading in from file
BOOST_AUTO_TEST_CASE(testStaticLookups)
{
    // goes through all testing files
    for (std::string FILE:FILES) 
    {
        // reads the YAML file into CompoundConfig
        CompoundConfig cConfig = CompoundConfig(TEST_LOC + FILE, "yaml");
        // reads in the YAML file independently of CompoundConfig
        YAML::Node testRef = YAML::LoadFile(TEST_LOC + FILE);

        // gets the root CompoundConfigNode
        config::CompoundConfigNode root = cConfig.getRoot();

        // goes through all keys and compares the values.
        // figure this out. you may need a helper file.
        for (auto node: testRef)
        {
            // extracts the key
            std::string key = node.first.as<std::string>();

            // tests all lookups with generic functions based on YAML type
            switch(node.Type())
            {
                // null should pull out the same thing as scalar
                case YAML::NodeType::Null:
                // tests all possible scalar output values
                case YAML::NodeType::Scalar:
                    // tests precision values
                    BOOST_CHECK(testScalarLookup<double>(root, node, key));
                    BOOST_CHECK(testScalarLookup<bool>(root, node, key));
                    BOOST_CHECK(testScalarLookup<int>(root, node, key));
                    BOOST_CHECK(testScalarLookup<unsigned int>(root, node, key));
                    // implicitly tests the lookupValueLongOnly
                    BOOST_CHECK(testScalarLookup<long long>(root, node, key));
                    BOOST_CHECK(testScalarLookup<unsigned long long>(root, node, key));
                    // tests floating point values
                    BOOST_CHECK(testScalarLookup<double>(root, node, key));
                    BOOST_CHECK(testScalarLookup<float>(root, node, key));
                    // tests strings
                    // TODO:: This doesn't compile figure it out later
                    // BOOST_CHECK(testScalarLookup<const char *>(root, node, key));
                    BOOST_CHECK(testScalarLookup<std::string>(root, node, key));
                    break;
                case YAML::NodeType::Sequence:
                    break;
                case YAML::NodeType::Map:
                    break;
                case YAML::NodeType::Undefined:
                    break;
                default:
                    throw std::runtime_error(
                        std::string("Wrong YAML FORMAT: ") + FILE
                    );
                    break;
            }
        }
    }
}

// tests the ability to set correctly
BOOST_AUTO_TEST_CASE(testSetters)
{
    std::cout << "not yet implemented" << std::endl;
}

// tests the ability to read out correctly from sets
BOOST_AUTO_TEST_CASE(testDynamicLookups)
{
    std:: cout << "not yet implemented" << std::endl;
}
} // namespace config