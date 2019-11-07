#pragma once

#include <cstdint>

namespace options
{
void usage(int argc, char **argv);
void parse(int argc, char **argv);

extern bool memory_usage;
extern char *parameter_file;
extern bool disable;
} // namespace options
