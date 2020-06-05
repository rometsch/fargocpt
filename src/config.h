#ifndef CONFIG_H
#define CONFIG_H

#include "parameters.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "LowTasks.h"
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
        die("Could not fully parse value '%s' to %s", val.c_str(), typeid(ret).name());
    }
    } else {
	ret = j.get<T>();
    }
    return ret;
}

class AutoTypeInterface
{
  private:
    json m_j;

  public:
    AutoTypeInterface(const json &j) : m_j(j) {}

    template <typename T> operator T() { return json_caster<T>(m_j); }
};

class Config
{
  private:
    json m_j;

  public:
    Config(){};
    Config(const char *filename);
    Config(const json &j) : m_j(j) { insert_lowercase_keys(); };

    ~Config(){};

    void load_file(const char *filename);
    void insert_lowercase_keys();

    AutoTypeInterface get(const char *key);

    template <typename T> T get(const char *key, const T &default_value);

    std::string get(const char *key, const char *default_value);

    bool get_flag(const char *key);
    bool get_flag(const char *key, const bool default_value);

    template <typename T> bool contains(const T &key);

    Config get_subconfig(const char *key)
    {
	return Config(json(m_j[lowercase(key)]));
    }
};

template <typename T> bool Config::contains(const T &key)
{
    const bool ret = m_j.contains(lowercase(key));
    return ret;
}

template <typename T> T Config::get(const char *key, const T &default_value)
{
    const std::string lkey = lowercase(key);
    if (contains(lkey)) {
	T ret = get(key);
	return ret;
    } else {
	return default_value;
    }
}

#endif // CONFIG_H
