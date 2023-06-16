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
        // predeclares values
        T expectedScalar, actualScalar;
        // value resolution
        expectedScalar = YNode[key].as<T>();
        actualScalar = CNode.lookupValue(key, expectedScalar);

        // equality check
        return expectedScalar == actualScalar;
    } catch(const YAML::TypedBadConversion<T>& e) {
        // defaults to false on bad conversion
        return false;
    }
}
// foreward declaration
bool nodeEq(config::CompoundConfigNode CNode, YAML::Node YNode, 
                    const std::string& key, YAML::NodeType::value TYPE);
// forward declaration
bool testMapLookup(config::CompoundConfigNode& CNode, YAML::Node&YNode);
// makes sure sequences agree in CCN and YNode
bool testSequenceLookup(config::CompoundConfigNode& CNode, YAML::Node& YNode, const std::string& key)
{
    bool equal = true;
    // grabs next values so they're owned somewhere
    auto childCNode = CNode.lookup(key);
    auto childYNode = YNode[key];

    // goes through all elements in the sequence
    for (int i = 0; (std::size_t) i < childYNode.size(); i++)
    {
        // unpacks element
        auto nextCNode = childCNode[i];
        auto nextYNode = childYNode[i];
        // only works because values are always associated by maps
        equal = equal && testMapLookup(nextCNode, nextYNode);
    }

    return equal;
}
// fetches child values as C++ doesn't like temporary values
bool testMapLookup(config::CompoundConfigNode& CNode, YAML::Node&YNode, const std::string& key)
{
    auto childCNode = CNode.lookup(key);
    auto childYNode = YNode[key];
    return testMapLookup(childCNode, childYNode);
}
// tests the CCN lookup functions provided a given root node. Treats all input nodes as maps.
bool testMapLookup(config::CompoundConfigNode& CNode, YAML::Node&YNode)
{
    // defines return value namespace
    bool equal = true;

    // goes through all keys and compares the values.
    for (auto nodeMapPair: YNode)
    {
        // extracts the key
        const std::string key = nodeMapPair.first.as<std::string>();
        // gets nodeEq result
        bool res = nodeEq(CNode, YNode, key, nodeMapPair.second.Type());
        // tests all lookups for this node
        equal = equal && res;
    }

    return equal;
}
// ensures that CNode and YNode are equal
bool nodeEq(config::CompoundConfigNode CNode, YAML::Node YNode, 
                    const std::string& key, YAML::NodeType::value TYPE)
{
    // namespace of if a node is correct
    bool nodePass = false;

    // determines what check to do based off child node type
    switch(TYPE)
    {
        // null should pull out the same thing as scalar
        case YAML::NodeType::Null:
            break;
        // tests all possible scalar output values
        case YAML::NodeType::Scalar:
            nodePass = nodePass ||
                        // tests precision values
                       testScalarLookup<double>(CNode, YNode, key) ||
                       testScalarLookup<bool>(CNode, YNode, key) ||
                       testScalarLookup<int>(CNode, YNode, key) ||
                       testScalarLookup<unsigned int>(CNode, YNode, key) ||
                       testScalarLookup<long long>(CNode, YNode, key) ||
                       testScalarLookup<unsigned long long>(CNode, YNode, key) ||
                       // tests floating points
                       testScalarLookup<double>(CNode, YNode, key) ||
                       testScalarLookup<float>(CNode, YNode, key) ||
                       // tests strings
                       // TODO:: This doesn't compile figure it out later
                       // BOOST_CHECK(testScalarLookup<const char *>(root, node, key));
                       testScalarLookup<std::string>(CNode, YNode, key);
            break;
        case YAML::NodeType::Sequence:
            nodePass = testSequenceLookup(CNode, YNode, key);
            break;
        case YAML::NodeType::Map:
            nodePass = testMapLookup(CNode, YNode, key);
            break;
        case YAML::NodeType::Undefined:
            break;
        default:
            break;
    }
    return nodePass;
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
            BOOST_CHECK(testMapLookup(root, ref));
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