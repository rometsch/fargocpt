/**
 * @file config.cpp
 */

#include <algorithm>
#include <string>
#include <type_traits>
#include <typeinfo>

#include "LowTasks.h"
#include "config.h"

#define UNITS_DEFAULT_DOMAIN units::domains::astronomy
#include "units/units.hpp"

#include "yaml-cpp/yaml.h"

namespace config {

Config cfg;

static YAML::Node lowercased_node(const YAML::Node &node) {
    YAML::Node lnode;
    for (auto it=node.begin(); it!=node.end(); ++it) {
        std::string key = it->first.as<std::string>();
        std::string lkey = lowercase(key);
        auto &subnode = it->second;
        if (subnode.IsMap()) {
            lnode[lkey] = lowercased_node(subnode);
        } else {
            lnode[lkey] = YAML::Node(subnode);
        }
    }
    return lnode;
}

static bool has_unit(const std::string val){
    auto q = units::measurement_from_string(val);
    const bool ret =  (q.units().unit_type_count() > 0);
    return ret;
}

Config::Config(const char *filename) { 
    load_file(filename); }

Config::Config(const YAML::Node &n)
{
    m_root = std::make_shared<YAML::Node>(lowercased_node(n));
};

Config Config::get_subconfig(const char *key)
{
    const auto &n = *m_root;
    return Config((n[lowercase(key)]));
}


std::vector<Config> Config::get_planet_config()
{
    const auto &n = *m_root;
    std::vector<Config> planets;
    std::cout << "creating planet configs" << std::endl;
    if (contains("planets")) {
        auto & node = n["planets"];
        for (auto &planet_node : n["planets"]) {
            planets.emplace_back(planet_node);
        }
        std::cout << "finished planet configs" << std::endl;
    }
    std::cout << "returning planet configs" << std::endl;
    return planets;
}

void Config::print() {
    std::cout << *m_root << std::endl;
}

void Config::load_file(const char *filename)
{
	std::ifstream infile(filename);
    YAML::Node node = YAML::Load(infile);
    m_root = std::make_shared<YAML::Node>(lowercased_node(node));
}

bool Config::get_flag(const char *key)
{
    const YAML::Node & root = *m_root;
    const std::string lkey = lowercase(key);
    bool ret = root[lkey].as<bool>();
    return ret;
}

bool Config::get_flag(const char *key, const bool default_value)
{
    if (contains(key)) {
	return get_flag(key);
    } else {
	return default_value;
    }
}

char Config::get_first_letter_lowercase(const char *key, const char *default_value)
{
    std::string value = get<std::string>(key, default_value);
    if (value.length() == 0) {
        value = std::string(default_value);
    }
    return tolower(value[0]);
}

char Config::get_first_letter_lowercase(const char *key)
{
    std::string value = get<std::string>(key);
    return tolower(value[0]);
}

bool Config::contains(const char *key)
{
    const YAML::Node & root = *m_root;
    if (root[lowercase(key)]) {
        return true;
    } else {
        return false;
    }
}

bool Config::contains(const std::string &key)
{
    const YAML::Node & root = *m_root;
    if (root[lowercase(key)]) {
        return true;
    } else {
        return false;
    }
}

template <typename F> 
typename std::enable_if<!std::is_arithmetic<F>::value,F>::type 
parse_units(const std::string &val) {
    const F rv = val;
    return rv;
}

template <typename F> 
typename std::enable_if<std::is_arithmetic<F>::value,F>::type 
parse_units(const std::string &val) {
    auto q = units::measurement_from_string(val);
    const F rv = (F) q.convert_to_base().value();
    if (has_unit(val)) {
    std::cout << "have value string '" << val <<"'" << std::endl;
    std::cout << "measurement is " << units::to_string(q) << std::endl;
    std::cout << "Parsing " << typeid(F).name() << " " << val << " as number" << std::endl;
    std::cout << "Result " << rv << " from " << units::to_string(q) << std::endl;
    }
    return rv;
}


template <typename F> F Config::get(const char *key)
{
    const auto &root = *m_root;
    std::string lkey = lowercase(key);
    const std::string val = root[lkey].as<std::string>();
    const F rv = parse_units<F>(val);
    return rv;
}

template <typename F> F Config::get(const char *key, const F &default_value)
{
    const std::string lkey = lowercase(key);
    F ret;
    if (contains(lkey)) {
	    ret = get<F>(key);
    } else {
	    ret = default_value;
    }
    return ret;
}

template double Config::get(const char *key);
template int Config::get(const char *key);
template unsigned int Config::get(const char *key);
template std::string Config::get(const char *key);

template double Config::get(const char *key, const double &d);
template int Config::get(const char *key, const int &d);
template unsigned int Config::get(const char *key, const unsigned int &d);
template std::string Config::get(const char *key, const std::string &d);

}