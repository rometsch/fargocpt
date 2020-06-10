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

    template <typename T> T get(const char *key);

    template <typename T> T get(const char *key, const T &default_value);

    std::string get(const char *key, const char *default_value);

    bool get_flag(const char *key);
    bool get_flag(const char *key, const bool default_value);

    bool contains(const char* key);
    bool contains(const std::string& key);

    Config get_subconfig(const char *key)
    {
	return Config(json(m_j[lowercase(key)]));
    }
};

#endif // CONFIG_H
