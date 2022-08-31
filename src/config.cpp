/**
 * @file config.cpp
 */

#include <algorithm>
#include <string>
#include <type_traits>
#include <cmath>

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

Config::Config(const std::string &filename) { 
    load_file(filename); }

Config::Config(const YAML::Node &n)
{
    m_root = std::make_shared<YAML::Node>(YAML::Clone(lowercased_node(n)));
    m_default = std::make_shared<YAML::Node>();
};

Config Config::get_subconfig(const std::string &key)
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

void Config::print_default() {
    std::cout << *m_default << std::endl;
}

void Config::load_file(const std::string &filename)
{
	std::ifstream infile(filename);
    YAML::Node node = YAML::Load(infile);
    m_root = std::make_shared<YAML::Node>(lowercased_node(node));
    m_default = std::make_shared<YAML::Node>();
}

template <typename T> bool isnan(const T &x) {
    return std::isnan(x);
}

template <> bool isnan(const std::string &x) {
    if (lowercase(x) == "nan") {
        return true;
    } else {
        return false;
    }
}


static bool string_decide(const std::string &val) {
    bool rv;
    const std::string l = lowercase(val);
    if (l == "yes" || l == "true") {
        rv = true;
    } else if (l == "no" || l == "false") {
        rv = false;
    } else {
        die("Invalid argument for a flag : '%s'!\nValid choices are: yes/no, true/false", val.c_str());
    }
    return rv;
}

bool Config::get_flag(const std::string &key)
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

bool Config::get_flag(const std::string &key, const bool default_value)
{
    (*m_default)[key] = default_value;
    if (contains(key)) {
	return get_flag(key);
    } else {
	return default_value;
    }
}

bool Config::get_flag(const std::string &key, const std::string &default_value)
{
    (*m_default)[key] = default_value;
    if (contains(key)) {
	return get_flag(key);
    } else {
	return string_decide(default_value);
    }
}

bool Config::get_flag(const std::string &key, const char *default_value)
{
    (*m_default)[key] = default_value;
    return get_flag(key, std::string(default_value));
}

char Config::get_first_letter_lowercase(const std::string &key, const std::string &default_value)
{
    (*m_default)[key] = default_value;
    std::string value = get<std::string>(key, default_value);
    if (value.length() == 0) {
        value = std::string(default_value);
    }
    return tolower(value[0]);
}

char Config::get_first_letter_lowercase(const std::string &key)
{
    std::string value = get<std::string>(key);
    return tolower(value[0]);
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

template <typename T> T Config::get(const std::string &key)
{
    if (!contains(key)) {
        die("Required parameter '%s' missing!\n", key.c_str());
    }
    T rv;
    const auto &root = *m_root;
    std::string lkey = lowercase(key);
    try {
        rv = root[lkey].as<T>();
    } catch (YAML::TypedBadConversion<T> const &) {
        die("Conversion from yaml failed for key '%s'\n", key.c_str());
    }
    if (isnan(rv)) {
        die("Value for key '%s' was converted to nan!\n", key.c_str());
    }
    return rv;
}


template <typename T> T Config::get(const std::string &key, const units::precise_unit& unit)
{
    if (!contains(key)) {
        die("Required parameter '%s' missing!\n", key.c_str());
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
            die("Conversion from yaml failed for key '%s'\n", key.c_str());
        }
    }
    if (isnan(rv)) {
        die("Value for key '%s' was converted to nan!\n", key.c_str());
    }
    return rv;
}


template <typename T> T Config::get(const std::string &key, const T &default_value)
{
    (*m_default)[key] = default_value;
    T rv;
    const std::string lkey = lowercase(key);
    if (contains(key)) {
        rv = get<T>(key);
    } else {
	    rv = default_value;
    }
    if (isnan(rv)) {
        die("Value for key '%s' was converted to nan!\n", key.c_str());
    }
    return rv;
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

template <typename T> T Config::get(const std::string &key, 
                            const T &default_value, 
                            const units::precise_unit& unit) 
{    
    (*m_default)[key] = default_value;
    T rv;
    const std::string lkey = lowercase(key);
    std::string val;
    const auto &root = *m_root;
    if (contains(key)) {
        val = root[lkey].as<std::string>();
        if (units::has_unit(val)) {
            rv = units::parse_units<T>(val, unit);
        } else {
            rv = stoT<T>(val);
        }
    } else {
	    rv = default_value;
    }
    if (isnan(rv)) {
        die("Value for key '%s' was converted to nan!\n", key.c_str());
    }
    return rv;
};

template <typename T> T Config::get(const std::string &key, 
                            const std::string &default_value, 
                            const units::precise_unit& unit) 
{
    (*m_default)[key] = default_value;
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
    if (isnan(rv)) {
        die("Value for key '%s' was converted to nan!\n", key.c_str());
    }
    return rv;
};

template double Config::get(const std::string &key);
template int Config::get(const std::string &key);
template unsigned int Config::get(const std::string &key);
template std::string Config::get(const std::string &key);

template double Config::get(const std::string &key, const double &d);
template int Config::get(const std::string &key, const int &d);
template unsigned int Config::get(const std::string &key, const unsigned int &d);
template std::string Config::get(const std::string &key, const std::string &d);

template double Config::get(const std::string &key, const double &d, const units::precise_unit& unit);
template int Config::get(const std::string &key, const int &d, const units::precise_unit& unit);
template unsigned int Config::get(const std::string &key, const unsigned int &d, const units::precise_unit& unit);

template double Config::get(const std::string &key, const std::string &d, const units::precise_unit& unit);
template int Config::get(const std::string &key, const std::string &d, const units::precise_unit& unit);
template unsigned int Config::get(const std::string &key, const std::string &d, const units::precise_unit& unit);

template double Config::get(const std::string &key, const units::precise_unit& unit);
template int Config::get(const std::string &key, const units::precise_unit& unit);
template unsigned int Config::get(const std::string &key, const units::precise_unit& unit);



}