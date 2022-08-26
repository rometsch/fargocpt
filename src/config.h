#ifndef CONFIG_H
#define CONFIG_H

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "LowTasks.h"

#include "yaml-cpp/yaml.h"
namespace config {


class Config
{
  private:
    std::shared_ptr<YAML::Node> m_root;

  public:
    Config(){};
    Config(const char *filename);
    Config(const YAML::Node &n);
    ~Config(){};

    void load_file(const char *filename);

    template <typename T> T get(const char *key);

    template <typename T> T get(const char *key, const T &default_value);

    template <typename T> T get(const char *key, 
                                const T &default_value, 
                                const units::precise_unit& unit);

    template <typename T> T get(const char *key, 
                                const std::string &default_value, 
                                const units::precise_unit& unit);

    std::string get(const char *key, const char *default_value);

    bool get_flag(const char *key);
    bool get_flag(const char *key, const bool default_value);

    char get_first_letter_lowercase(const char *key);
    char get_first_letter_lowercase(const char *key, const char *default_value);

    bool contains(const char *key);
    bool contains(const std::string &key);

    Config get_subconfig(const char *key);
    std::vector<Config> get_planet_config();

    void print();
};

extern Config cfg;

}
#endif // CONFIG_H
