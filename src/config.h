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
    Config(const std::string& filename);
    Config(const YAML::Node &n);
    ~Config(){};

    void load_file(const std::string& filename);

    template <typename T> T get(const std::string &key);

    template <typename T> T get(const std::string &key, const T &default_value);

    template <typename T> T get(const std::string &key, const units::precise_unit& unit);

    template <typename T> T get(const std::string &key, 
                                const T &default_value, 
                                const units::precise_unit& unit);

    template <typename T> T get(const std::string &key, 
                                const std::string &default_value, 
                                const units::precise_unit& unit);

    bool get_flag(const std::string &key);
    bool get_flag(const std::string &key, const bool default_value);

    char get_first_letter_lowercase(const std::string &key);
    char get_first_letter_lowercase(const std::string &key, const std::string &default_value);

    bool contains(const std::string &key);

    Config get_subconfig(const std::string &key);
    std::vector<Config> get_planet_config();

    void print();
};

extern Config cfg;

}
#endif // CONFIG_H
