/**
 * @file config.cpp
 */

#include <algorithm>
#include <string>
#include <type_traits>
#include <cmath>
#include <numeric>
#include <typeinfo>

#include "LowTasks.h"
#include "config.h"

#include "units.h"

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

Config::Config() {
}

Config::Config(const std::string &filename) { 
    load_file(filename); 
}

Config::Config(const YAML::Node &n)
{
    m_root = std::make_shared<YAML::Node>(YAML::Clone(lowercased_node(n)));
    m_default = std::make_shared<YAML::Node>();
}

Config Config::get_subconfig(const std::string &key)
{
    m_visited_keys.insert(lowercase(key));
    const auto &n = *m_root;
    return Config((n[lowercase(key)]));
}


std::vector<Config> Config::get_nbody_config()
{
    m_visited_keys.insert("nbody");
    const auto &n = *m_root;
    std::vector<Config> bodies;
    if (contains("nbody")) {
        for (auto &new_node : n["nbody"]) {
            bodies.emplace_back(YAML::Clone(new_node));
        }
    }
    return bodies;
}

void Config::print() {
    std::cout << *m_root << std::endl;
}

void Config::print_default() {
    std::cout << *m_default << std::endl;
}

void Config::add_to_default(const std::string key, const std::string type, const bool unit_support) {
    YAML::Node & defs = *m_default;
    if (!defs[key]) {
        defs[key] = YAML::Node();
    }
    if (!defs[key]["type"]) {
        defs[key]["type"] = type;
    }
    if (!defs[key]["unitsupport"]) {
        defs[key]["unitsupport"] = unit_support ? "yes" : "no";
    }
}

template <typename T> void Config::add_to_default(const std::string key, const std::string type, const bool unit_support, const T default_value) {
    YAML::Node & defs = *m_default;
    if (!defs[key]) {
        defs[key] = YAML::Node();
    }
    if (!defs[key]["default"]) {
        defs[key]["default"] = default_value;
    }
    if (!defs[key]["type"]) {
        defs[key]["type"] = type;
    }
    if (!defs[key]["unitsupport"]) {
        defs[key]["unitsupport"] = unit_support ? "yes" : "no";
    }
}

std::string Config::visited_keys() {
    std::string rv;
    for (auto &key : m_visited_keys) {
        rv += key + ", ";
    }
    return rv;
}

std::string Config::unknown_keys() {
    std::string rv = std::accumulate(std::next(m_unknown_keys.begin()), m_unknown_keys.end(), *m_unknown_keys.begin(), 
        [](std::string a, std::string b) {
            return a + ", " + b; // append ", " between words
        }
    );
    return rv;
}

bool Config::are_all_keys_visited() {
    const auto &n = *m_root;
    for (auto it=n.begin(); it!=n.end(); ++it) {
        std::string key = it->first.as<std::string>();
        if (m_visited_keys.find(key) == m_visited_keys.end()) {
            m_unknown_keys.insert(key);
        }
    }
    return m_unknown_keys.size() == 0;
}

void Config::exit_on_unknown_key() {
    if (!are_all_keys_visited()) {
        die("Unknown key(s) found in config file: '%s'\nMaybe there is a typo?\n", unknown_keys().c_str());
    }
}

void Config::write_default(const std::string &filename) {
    const auto & default_values = *m_default;
    const std::string outstr = YAML::Dump(default_values);
    std::ofstream out(filename);
    out << outstr;
    out.close();
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
    bool rv = false;
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
        die("Conversion from yaml failed for key '%s'\n", key.c_str());
    }
    return rv;
}

bool Config::get_flag(const std::string &key, const bool default_value)
{
    m_visited_keys.insert(lowercase(key));
    add_to_default(key, typeid(bool).name(), false, default_value);
    if (contains(key)) {
	return get_flag(key);
    } else {
	return default_value;
    }
}

bool Config::get_flag(const std::string &key, const std::string &default_value)
{
    m_visited_keys.insert(lowercase(key));
    add_to_default(key, typeid(bool).name(), false, default_value);
    if (contains(key)) {
	return get_flag(key);
    } else {
	return string_decide(default_value);
    }
}

bool Config::get_flag(const std::string &key, const char *default_value)
{
    m_visited_keys.insert(lowercase(key));
    add_to_default(key, typeid(bool).name(), false, default_value);
    return get_flag(key, std::string(default_value));
}

char Config::get_first_letter_lowercase(const std::string &key, const std::string &default_value)
{
    m_visited_keys.insert(lowercase(key));
    add_to_default(key, typeid(std::string).name(), false, default_value);
    std::string value = get<std::string>(key, default_value);
    if (value.length() == 0) {
        value = std::string(default_value);
    }
    return (char)tolower(value[0]);
}

char Config::get_first_letter_lowercase(const std::string &key)
{
    m_visited_keys.insert(lowercase(key));
    add_to_default(key, typeid(std::string).name(), false);
    std::string value = get<std::string>(key);
    return (char)tolower(value[0]);
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
    m_visited_keys.insert(lowercase(key));
    add_to_default(key, typeid(T).name(), false);
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
    m_visited_keys.insert(lowercase(key));
    add_to_default(key, typeid(T).name(), true);
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
    m_visited_keys.insert(lowercase(key));
    add_to_default(key, typeid(T).name(), false, default_value);
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
    m_visited_keys.insert(lowercase(key));
    add_to_default(key, typeid(T).name(), true, default_value);
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
}

template <typename T> T Config::get(const std::string &key, 
                            const std::string &default_value, 
                            const units::precise_unit& unit) 
{
    m_visited_keys.insert(lowercase(key));
    add_to_default(key, typeid(T).name(), true, default_value);
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
}

std::string Config::get_lowercase(const std::string &key, const std::string &default_value) {
    m_visited_keys.insert(lowercase(key));
    add_to_default(key, typeid(std::string).name(), false, default_value);
    std::string rv;
    const std::string lkey = lowercase(key);
    if (contains(key)) {
        rv = get<std::string>(key);
    } else {
	    rv = default_value;
    }
    return lowercase(rv);
}

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

template void Config::add_to_default(const std::string key, const std::string type, const bool unit_support, const int default_value);
template void Config::add_to_default(const std::string key, const std::string type, const bool unit_support, const unsigned int default_value);
template void Config::add_to_default(const std::string key, const std::string type, const bool unit_support, const bool default_value);
template void Config::add_to_default(const std::string key, const std::string type, const bool unit_support, const std::string default_value);
template void Config::add_to_default(const std::string key, const std::string type, const bool unit_support, const char* default_value);

}
