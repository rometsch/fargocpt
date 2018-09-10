#ifndef _CONFIG_H_
#define _CONFIG_H_

namespace config {

int read_config_from_file(const char *filename);
void print_parsed_config();
int key_exists(const char *key);
bool value_as_bool(const char *key);
int value_as_int(const char *key);
unsigned int value_as_unsigned_int(const char *key);
double value_as_double(const char *key);
const char *value_as_string(const char *key);
bool value_as_bool_default(const char *key, bool defvalue);
int value_as_int_default(const char *key, int defvalue);
unsigned int value_as_unsigned_int_default(const char *key, unsigned int defvalue);
double value_as_double_default(const char *key, double defvalue);
const char *value_as_string_default(const char *key, const char *defvalue);

}

#endif
