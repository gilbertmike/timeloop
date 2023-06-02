// represents YAML maps
#include <map>
// for type-safe unions of YAML types
#include <variant>
// for YAML arrays
#include <vector>
// for std errors from variants
#include <stdexcept>
#include <memory>
#include <iostream>

namespace structured_config {
class CCRet; // forward definition

// Literal value types possible in a YAML file
using YAMLLiteral = std::variant<
    std::string, int, float
>;
// YAML vector representation
using YAMLVector = std::vector<std::unique_ptr<CCRet>>;
// YAML map representation
using YAMLMap = std::map<std::string, std::unique_ptr<CCRet>>;
// YAML value types possible
using YAMLType = std::variant<
    YAMLLiteral,
    YAMLVector,
    YAMLMap
>;

class CCRet
{
    private:
        YAMLType data_ = nullptr;

    public:
        /** VALUE RESOLUTION **/
        // unpacks a literal CCRet
        inline YAMLLiteral& GetValue()
        { return std::get<YAMLLiteral>(data_); } 
        // resolves a map CCRet
        inline CCRet& At(const std::string& key)
        { return *std::get<YAMLMap>(data_).at(key); }
        // resolves a list CCRet
        inline CCRet& At(YAMLVector::size_type index)
        { return *std::get<YAMLVector>(data_).at(index); }

        /** CONTAINER FUNCTIONS **/
        // We could internally cast everything to unsigned long long to cut down
        // on code size. 
        bool lookupValue(const char *name, bool &value) const;
        bool lookupValue(const char *name, int &value) const;
        bool lookupValue(const char *name, unsigned int &value) const;
        bool lookupValueLongOnly(const char *name, long long &value) const; // Only for values with an L like 123L
        bool lookupValueLongOnly(const char *name, unsigned long long &value) const; // Only for values with an L like 123L
        bool lookupValue(const char *name, long long &value) const;
        bool lookupValue(const char *name, unsigned long long &value) const;
        bool lookupValue(const char *name, double &value) const;
        bool lookupValue(const char *name, float &value) const;
        bool lookupValue(const char *name, const char *&value) const;
        // this is the only one we cant cast to unsigned long long
        bool lookupValue(const char *name, std::string &value) const;

        // Michael added this, no clue what it does for now
        template<typename... ArgsT>
        void EmplaceBack(ArgsT&&... args)
        {
            std::get<YAMLVector>(data_).push_back(
                std::make_unique<CCRet>(std::forward<ArgsT>(args)...)
            );
        }

        size_t Size() const
        {
            return std::visit(
                [] (auto&& data)
                {
                    using DataT = std::decay_t<decltype(data)>;
                    if constexpr (std::is_same_v<DataT, YAMLLiteral>)
                    {
                        return (size_t)1;
                    }
                    else
                    {
                        return data.size();
                    }
                },
                data_
            );
        }

        // Michael added this, no clue what it does for now
        template<typename T>
        static CCRet Literal(const T& val)
        {
            auto ret = CCRet();
            ret.data_ = val;
            return ret;
        }

        static CCRet Vector()
        {
            auto ret = CCRet();
            ret.data_ = YAMLVector();
            return ret;
        }

        static CCRet Map()
        {
            auto ret = CCRet();
            ret.data_ = YAMLMap();
            return ret;
        }

    private:
        CCRet() : data_() {}
};
}