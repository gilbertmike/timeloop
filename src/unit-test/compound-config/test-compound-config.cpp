#define BOOST_TEST_MODULE TestCompoundConfig

#include <iostream>
#include <boost/test/included/unit_test.hpp>
#include <compound-config/compound-config.hpp>

// number of testing cycles to run
int TESTS = 2000;
// the seed for the entropy source
uint SEED = 42;
// changes the max random value to the U_LONG32 random value
#undef RAND_MAX
#define RAND_MAX = ULONG_LONG_MAX

// the location of the test files
std::string TEST_LOC = "./src/unit-test/compound-config/tests/";

// static YAML file names we want to load in for the test
std::map<std::string, std::vector<std::string>> FILES = {
    // https://github.com/Accelergy-Project/timeloop-accelergy-exercises
    {
        "accelergy-project/2020.ispass/timeloop/01/", {
            "1level.arch.yaml",
            "conv1d-1level.map.yaml",
            "conv1d.prob.yaml",
            "timeloop-model.ART_summary.yaml",
            "timeloop-model.ART.yaml",
            "timeloop-model.ERT_summary.yaml",
            "timeloop-model.ERT.yaml",
            "timeloop-model.flattened_architecture.yaml",
        }
    }
};

// makes sure for a certain type CCN agrees with YNode. Defaults to false if conversion is wrong.
template <typename T>
bool testScalarLookup(config::CompoundConfigNode& CNode, YAML::Node& YNode, const std::string& key)
{
    try {
        T expectedScalar = YNode[key].as<T>();
        T actualScalar;
        CNode.lookupValue(key, actualScalar);

        return expectedScalar == actualScalar;
    } catch(const YAML::TypedBadConversion<T>& e) {
        return false;
    }
}

// makes sure sequences agree in CCN and YNode
bool testSequenceLookup(config::CompoundConfigNode& CNode, YAML::Node& YNode, const std::string& key)
{
    std::vector<std::string> expectedSeq = YNode[key].as<std::vector<std::string>>();
    std::vector<std::string> actualSeq;
    CNode.lookupArrayValue(key, actualSeq);

    return expectedSeq == actualSeq;
}

// forward declaration
bool testMapLookup(config::CompoundConfigNode& CNode, YAML::Node&YNode);
// accesses the key for a map because C++ dislikes things not being owned by other things
bool testMapLookup(config::CompoundConfigNode& CNode, YAML::Node&YNode, const std::string& key)
{
    auto nextCNode = CNode.lookup(key);
    auto nextYNode = YNode[key];
    return testMapLookup(nextCNode, nextYNode);
}

// tests the CCN lookup functions provided a given root node. Treats all input nodes as maps.
bool testMapLookup(config::CompoundConfigNode& CNode, YAML::Node&YNode)
{
    // defines return value namespace
    bool ret = true;

    // goes through all keys and compares the values.
    for (auto nodeIterVal: YNode)
    {
        std::cout << "looped" << nodeIterVal.second.Type() << std::endl;
        // whether or not comparison at this node has passed
        bool nodePass = false;
        // extracts the key
        const std::string key = nodeIterVal.first.as<std::string>();

        // tests all lookups with generic functions based on the YAMLType of the YAML map value
        switch(nodeIterVal.second.Type())
        {
            // null should pull out the same thing as scalar
            case YAML::NodeType::Null:
                break;
            // tests all possible scalar output values
            case YAML::NodeType::Scalar:
                // tests precision values
                std::cout << nodePass << std::endl;
                BOOST_CHECK(nodePass = nodePass || testScalarLookup<double>(CNode, YNode, key));
                BOOST_CHECK(nodePass = nodePass || testScalarLookup<bool>(CNode, YNode, key));
                BOOST_CHECK(nodePass = nodePass || testScalarLookup<int>(CNode, YNode, key));
                BOOST_CHECK(nodePass = nodePass || testScalarLookup<unsigned int>(CNode, YNode, key));
                // implicitly tests the lookupValueLongOnly
                BOOST_CHECK(nodePass = nodePass || testScalarLookup<long long>(CNode, YNode, key));
                BOOST_CHECK(nodePass = nodePass || testScalarLookup<unsigned long long>(CNode, YNode, key));
                // tests floating point values
                BOOST_CHECK(nodePass = nodePass || testScalarLookup<double>(CNode, YNode, key));
                BOOST_CHECK(nodePass = nodePass || testScalarLookup<float>(CNode, YNode, key));
                // tests strings
                // TODO:: This doesn't compile figure it out later
                // BOOST_CHECK(testScalarLookup<const char *>(root, node, key));
                BOOST_CHECK(nodePass = nodePass || testScalarLookup<std::string>(CNode, YNode, key));
                break;
            case YAML::NodeType::Sequence:
                BOOST_CHECK(nodePass = testSequenceLookup(CNode, YNode, key));
                break;
            case YAML::NodeType::Map:
                BOOST_CHECK(nodePass = testMapLookup(CNode, YNode, key));
                break;
            case YAML::NodeType::Undefined:
                throw std::runtime_error(
                    std::string("None of our files should contain undefineds")
                );
                break;
            default:
                throw std::runtime_error(
                    std::string("You should not be here in testMapLookup")
                );
                break;
        }
        ret = ret && nodePass;
    }

    return ret;
}

// we are only testing things in config
namespace config {
// tests the lookup functions when reading in from file
BOOST_AUTO_TEST_CASE(testStaticLookups)
{
    // marker for test
    std::cout << "Beginning Static Lookups Test:\n---" << std::endl;
    // goes through all testing dirs
    for (auto FILEPATH:FILES) 
    {
        // calculates DIR relative location and extracts file's name
        std::string DIR = TEST_LOC + FILEPATH.first;
        std::vector<std::string> FILENAMES = FILEPATH.second;

        for (std::string FILE:FILENAMES)
        {
            // calculates filepath
            std::string FILEPATH = DIR + FILE;
            // debug printing info
            std::cout << "Now testing: " + FILEPATH << std::endl;
            // reads the YAML file into CompoundConfig and gets root
            CompoundConfig cConfig = CompoundConfig({FILEPATH});
            CompoundConfigNode root = cConfig.getRoot();
            // reads in the YAML file independently of CompoundConfig to serve as test reference
            YAML::Node ref = YAML::LoadFile(FILEPATH);

            // tests the entire file
            testMapLookup(root, ref);
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