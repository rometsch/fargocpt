/**
 * @file config.cpp
 */

#include <algorithm>
#include <string>
#include <type_traits>

#include "LowTasks.h"
#include "config.h"

#include "units.h"
#include "yaml-cpp/yaml.h"

namespace config {

Config cfg;

static YAML::Node lowercased_node(const YAML::Node &node) {
    YAML::Node lnode;
    for (auto it=node.begin(); it!=node.end(); ++it) {
        std::string key = it->first.as<std::string>();
        std::string lkey = lowercase(key);
        if (it->second.IsMap()) {
            lnode[lkey] = lowercased_node(it->second);
        } else {
            lnode[lkey] = it->second;
        }
    }
    return lnode;
}

Config::Config(const char *filename) { 
    load_file(filename); }

Config::Config(const YAML::Node &n)
{
    m_root = std::make_shared<YAML::Node>(YAML::Clone(lowercased_node(n)));
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
    if (contains("planets")) {
        auto & node = n["planets"];
        for (auto &planet_node : n["planets"]) {
            planets.emplace_back(YAML::Clone(planet_node));
        }
    }
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
    bool rv;
    const YAML::Node & root = *m_root;
    const std::string lkey = lowercase(key);
    try {
        rv = root[lkey].as<bool>();
    } catch (YAML::TypedBadConversion<bool> const &) {
        die("Conversion from yaml failed for key '%s'\n", key);
    }
    return rv;
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

template <typename T> T Config::get(const char *key)
{
    if (!contains(key)) {
        die("Required parameter '%s' missing!\n", key);
    }
    T rv;
    const auto &root = *m_root;
    std::string lkey = lowercase(key);
    try {
        rv = root[lkey].as<T>();
    } catch (YAML::TypedBadConversion<T> const &) {
        die("Conversion from yaml failed for key '%s'\n", key);
    }
    return rv;
}

template <typename T> T Config::get(const char *key, const units::precise_unit& unit)
{
    if (!contains(key)) {
        die("Required parameter '%s' missing!\n", key);
    }
    T rv;
    const auto &root = *m_root;
    std::string lkey = lowercase(key);
    const std::string val = root[lkey].as<std::string>();
    if (units::has_unit(val)) {
        rv = units::parse_units<T>(val, unit);
    } else {
        try {
            rv = root[lkey].as<T>();
        } catch (YAML::TypedBadConversion<T> const &) {
            die("Conversion from yaml failed for key '%s'\n", key);
        }
    }
    return rv;
}


template <typename T> T Config::get(const char *key, const T &default_value)
{
    T ret;
    const std::string lkey = lowercase(key);
    if (contains(key)) {
        ret = get<T>(key);
    } else {
	    ret = default_value;
    }
    return ret;
}

template <typename T> T stoT(const std::string &val);

template <> int stoT(const std::string &val) {
    return std::stoi(val);
}

template <> unsigned int stoT(const std::string &val) {
    return std::stoul(val);
}

template <> double stoT(const std::string &val) {
    return std::stod(val);
}

template <typename T> T Config::get(const char *key, 
                            const T &default_value, 
                            const units::precise_unit& unit) 
{    
    T ret;
    const std::string lkey = lowercase(key);
    std::string val;
    const auto &root = *m_root;
    if (contains(key)) {
        val = root[lkey].as<std::string>();
        if (units::has_unit(val)) {
            ret = units::parse_units<T>(val, unit);
        } else {
            ret = stoT<T>(val);
        }
    } else {
	    ret = default_value;
    }
    return ret;
};

template <typename T> T Config::get(const char *key, 
                            const std::string &default_value, 
                            const units::precise_unit& unit) 
{
    T rv;
    const std::string lkey = lowercase(key);
    std::string val;
    if (contains(lkey)) {
        const auto &root = *m_root;
        val = root[lkey].as<std::string>();        
    } else {
	    val = default_value;
    }
    if (units::has_unit(val)) {
        rv = units::parse_units<T>(val, unit);
    } else {
        rv = stoT<T>(val);
    }
    return rv;
};

template double Config::get(const char *key);
template int Config::get(const char *key);
template unsigned int Config::get(const char *key);
template std::string Config::get(const char *key);

template double Config::get(const char *key, const double &d);
template int Config::get(const char *key, const int &d);
template unsigned int Config::get(const char *key, const unsigned int &d);
template std::string Config::get(const char *key, const std::string &d);

template double Config::get(const char *key, const double &d, const units::precise_unit& unit);
template int Config::get(const char *key, const int &d, const units::precise_unit& unit);
template unsigned int Config::get(const char *key, const unsigned int &d, const units::precise_unit& unit);

template double Config::get(const char *key, const std::string &d, const units::precise_unit& unit);
template int Config::get(const char *key, const std::string &d, const units::precise_unit& unit);
template unsigned int Config::get(const char *key, const std::string &d, const units::precise_unit& unit);

template double Config::get(const char *key, const units::precise_unit& unit);
template int Config::get(const char *key, const units::precise_unit& unit);
template unsigned int Config::get(const char *key, const units::precise_unit& unit);



}