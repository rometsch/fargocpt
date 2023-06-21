#pragma once

#include <cstdint>
#include <string>

namespace options
{
void usage(int argc, char **argv);
void parse(int argc, char **argv);

extern bool memory_usage;
extern std::string parameter_file;
extern bool disable;
extern int max_iteration_number;
} // namespace options
