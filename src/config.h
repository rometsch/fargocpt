#ifndef CONFIG_H
#define CONFIG_H

#include "parameters.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "LowTasks.h"

#include "nlohmann/json_fwd.hpp"

using json = nlohmann::json;

class Config
{
  private:
    std::shared_ptr<json> m_j;

  public:
    Config(){};
    Config(const char *filename);
    Config(const json &j);
    ~Config(){};

    void load_file(const char *filename);
    void insert_lowercase_keys();

    template <typename T> T get(const char *key);

    template <typename T> T get(const char *key, const T &default_value);

    std::string get(const char *key, const char *default_value);

    bool get_flag(const char *key);
    bool get_flag(const char *key, const bool default_value);

    bool contains(const char *key);
    bool contains(const std::string &key);

    Config get_subconfig(const char *key);
};

#endif // CONFIG_H
