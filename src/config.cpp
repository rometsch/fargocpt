/**
 * @file config.cpp
 */

#include <algorithm>
#include <string>

#include "LowTasks.h"
#include "config.h"

#include "json/json.hpp"

using json = nlohmann::json;

template <typename T> T json_caster(const json &j)
{
    T ret;
    if (j.is_string()) {
	const std::string val = j;
	std::stringstream ss(val);
	ss >> ret;
	std::string rest;
	ss >> rest;
	if (rest.length() > 0) {
	    die("Could not fully parse value '%s' to %s", val.c_str(),
		typeid(ret).name());
	}
    } else {
	ret = j.get<T>();
    }
    return ret;
}

Config::Config(const char *filename) { load_file(filename); }

Config::Config(const json &j)
{
    m_j = std::make_shared<json>(j);
    insert_lowercase_keys();
};

Config Config::get_subconfig(const char *key)
{
    const auto &j = *m_j;
    return Config(json(j[lowercase(key)]));
}

std::vector<Config> Config::get_planet_config()
{
    const auto &j = *m_j;
    std::vector<Config> planets;
    if (contains("planets")) {
	for (auto &j_planet : j["planets"]) {
	    planets.push_back(Config(j_planet));
	}
    }
    return planets;
}

void Config::load_file(const char *filename)
{
    try {
	m_j.reset(new json());
	std::ifstream infile(filename);
	infile >> *m_j;
	insert_lowercase_keys();
    } catch (const nlohmann::detail::parse_error &ex) {
	std::cerr << "Could not parse json file." << std::endl;
    }
}

void Config::insert_lowercase_keys()
{
    auto &j = *m_j;
    for (auto &el : m_j->items()) {
	j[lowercase(el.key())] = el.value();
    }
}

static bool string_decide(const std::string &des)
{
    bool ret = false;
    std::string s = std::string(des);
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    if (s == "yes" || s == "y" || s == "true" || s == "1" || s == "on") {
	ret = true;
    } else if (s == "no" || s == "n" || s == "false" || s == "0 " ||
	       s == "off") {
	ret = false;
    } else {
	die("Can not parse '%s' to a boolean value!", des.c_str());
    }
    return ret;
}

bool Config::get_flag(const char *key)
{
    bool ret = false;
    const auto &j = *m_j;
    json val;
    auto &el = j[lowercase(key)];
    if (el.contains("value")) {
	val = el["value"];
    } else {
	val = el;
    }
    if (val.is_string()) {
	ret = string_decide(val);
    } else if (val.is_boolean()) {
	ret = val;
    }
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

bool Config::contains(const char *key)
{
    const bool ret = m_j->contains(lowercase(key));
    return ret;
}

bool Config::contains(const std::string &key)
{
    const bool ret = m_j->contains(lowercase(key));
    return ret;
}

template <typename T> T Config::get(const char *key)
{
    const auto &j = *m_j;
    json val;
    auto &el = j[lowercase(key)];
    if (el.contains("value")) {
	val = el["value"];
    } else {
	val = el;
    }
    const T ret = json_caster<T>(val);
    return ret;
}

// template <typename T> bool Config::contains(const T &key)
// {
//     const bool ret = m_j.contains(lowercase(key));
//     return ret;
// }
template <typename T> T Config::get(const char *key, const T &default_value)
{
    const std::string lkey = lowercase(key);
    if (contains(lkey)) {
	T ret = get<T>(key);
	return ret;
    } else {
	return default_value;
    }
}

// template bool Config::contains(const char* key);
// template bool Config::contains(const std::string& key);

template double Config::get(const char *key);
template int Config::get(const char *key);
template unsigned int Config::get(const char *key);
template std::string Config::get(const char *key);

template double Config::get(const char *key, const double &d);
template int Config::get(const char *key, const int &d);
template unsigned int Config::get(const char *key, const unsigned int &d);
template std::string Config::get(const char *key, const std::string &d);
