/**
 * @file config.cpp
 */

#include <algorithm>
#include <string> 

#include "config.h"
#include "LowTasks.h"

Config::Config(const char *filename) { load_file(filename); }

void Config::load_file(const char *filename)
{
    try {
	std::ifstream infile(filename);
	infile >> m_j;
    insert_lowercase_keys();
    } catch (const nlohmann::detail::parse_error &ex) {
	std::cerr << "Could not parse json file." << std::endl;
    }
}

void Config::insert_lowercase_keys() {
	for (auto &el : m_j.items()) {
	    m_j[lowercase(el.key())] = el.value();
	}
}

AutoTypeInterface Config::get(const char *key)
{
    std::cout << "Getting config value " << key << std::endl;
    AutoTypeInterface wrapper = AutoTypeInterface(m_j[lowercase(key)]["value"]);
    return wrapper;
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
    auto val = m_j[lowercase(key)]["value"];
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