/**
 * @file config.cpp
 * @author Tobias Mueller <Tobias_Mueller@twam.info>
 */

/***************************************************************************
 *   Copyright (C) 2009 by Tobias Mueller                                  *
 *   Tobias_Mueller@twam.info                                              *
 *                                                                         *
 *   based on code from git (http://git-scm.com/):                         *
 *     Linus Torvalds, 2005                                                *
 *     Johannes Schindelin, 2005                                           *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "config.h"
#include "LowTasks.h"
#include "logging.h"
#include <ctype.h>
#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

namespace config
{

#define MAXNAME (256)

/// current config file pointer
static FILE *file_ptr;
/// filename of config file
static const char *file_name;
/// current parsed line numbers of config file
static int file_linenumber;
/// end of config file reached?
static int file_EOF;

typedef struct {
    char *key;
    char *value;
} t_entry;

typedef struct {
    t_entry *data;
    unsigned int size;
} t_list;

static t_list list = {NULL, 0};

static void add_to_config_list(const char *key, const char *value);
static inline int iskeychar(int c);
static int get_next_char(void);
static char *parse_value(void);
static int get_value(char *key, unsigned int len);
static int parse_file();
static void die_bad_config(const char *key);
static int parse_unit_factor(const char *end, unsigned long *val);
static const char *value(const char *key);
static int parse_long(const char *value, long *ret);
static int parse_unsigned_long(const char *value, unsigned long *ret);
static int parse_double(const char *value, double *ret);
/*static unsigned long as_unsigned_long(const char *key, const char *value);*/
static int as_int(const char *key, const char *value);
static unsigned int as_unsigned_int(const char *key, const char *value);
static double as_double(const char *key, const char *value);
static int as_bool_or_int(const char *key, const char *value, int *is_bool);
static bool as_bool(const char *key, const char *value);

/**
	adds a key/value pair to the config list

	\param key key name of config entry
	\param value value of config entry
*/
static void add_to_config_list(const char *key, const char *value)
{
    list.data = (t_entry *)realloc(list.data, (++list.size) * sizeof(t_entry));

	const unsigned int key_length = strlen(key) + 1;
	const unsigned int val_length = strlen(value) + 1;
    list.data[list.size - 1].key =
	(char *)std::malloc((key_length) * sizeof(char));
    list.data[list.size - 1].value =
	(char *)std::malloc((val_length) * sizeof(char));

	if(list.data[list.size - 1].key == nullptr){
		die("add_to_config_list malloc failed %s\n", key);
	}
	if(list.data[list.size - 1].value == nullptr){
		die("add_to_config_list malloc failed %s\n", value);
	}

	strncpy(list.data[list.size - 1].key, key, key_length);
	strncpy(list.data[list.size - 1].value, value, val_length);
}

void free_config_list()
{

    for (unsigned int i = 0; i < list.size; ++i) {
	free(list.data[i].key);
	free(list.data[i].value);
    }
    free(list.data);
}

/**
	check if a char is a valid key character

	\param c char to check
	\returns 0 if no key char, otherwise != 0
*/
static inline int iskeychar(int c) { return isalnum(c) || c == '-'; }

/**
	read next char from config file

	\returns next char from config file
*/
static int get_next_char(void)
{
    int c;
    FILE *f;

    c = '\n';
    if ((f = file_ptr) != NULL) {
	c = fgetc(f);
	if (c == '\r') {
	    /* DOS like systems */
	    c = fgetc(f);
	    if (c != '\n') {
		ungetc(c, f);
		c = '\r';
	    }
	}
	if (c == '\n')
	    file_linenumber++;
	if (c == EOF) {
	    file_EOF = 1;
	    c = '\n';
	}
    }

    return c;
}

/**
	read the the value from the config file

	\returns value (as string)
*/
static char *parse_value(void)
{
    static char value[1024];
    int quote = 0, comment = 0, space = 0;
    unsigned int len = 0;

    for (;;) {
	int c = get_next_char();
	if (len >= sizeof(value) - 1)
	    return NULL;
	if (c == '\n') {
	    if (quote)
		return NULL;
	    value[len] = 0;
	    return value;
	}
	if (comment)
	    continue;
	if (isspace(c) && !quote) {
	    if (len)
		space++;
	    continue;
	}
	if (!quote) {
	    if (c == ';' || c == '#') {
		comment = 1;
		continue;
	    }
	}
	for (; space; space--)
	    value[len++] = ' ';
	if (c == '\\') {
	    c = get_next_char();
	    switch (c) {
	    case '\n':
		continue;
	    case 't':
		c = '\t';
		break;
	    case 'b':
		c = '\b';
		break;
	    case 'n':
		c = '\n';
		break;
	    /* Some characters escape as themselves */
	    case '\\':
	    case '"':
		break;
	    /* Reject unknown escape sequences */
	    default:
		return NULL;
	    }
	    value[len++] = c;
	    continue;
	}
	if (c == '"') {
	    quote = 1 - quote;
	    continue;
	}
	value[len++] = c;
    }
}

/**
	read the rest of the key name and value from the the config file

	\param key key correspondig to value
	\param len number of chars already read into key
	\returns 0 if successful
*/
static int get_value(char *key, unsigned int len)
{
    int c;
    char *value;

    /* Get the full key */
    for (;;) {
	c = get_next_char();
	if (file_EOF)
	    break;
	if (!iskeychar(c))
	    break;
	key[len++] = tolower(c);
	if (len >= MAXNAME)
	    return -1;
    }
    key[len] = 0;
    while (c == ' ' || c == '\t')
	c = get_next_char();

    value = NULL;

    if (c != '\n') {
	/* uncomment this line and comment the ungetc line if you want to use
	 * "parameter = value" style in your config file */
	/*if (c != '=')
		return -1;*/
	ungetc(c, file_ptr);
	value = parse_value();
	if (!value)
	    return -1;
    }

    if (value != NULL)
	add_to_config_list(key, value);

    return 0;
}

/**
	parse the config file
*/
static int parse_file()
{
    int comment = 0;
    int baselen = 0;
    static char var[MAXNAME];

    /* U+FEFF Byte Order Mark in UTF8 */
    static const unsigned char *utf8_bom = (unsigned char *)"\xef\xbb\xbf";
    const unsigned char *bomptr = utf8_bom;

    for (;;) {
	int c = get_next_char();
	if (bomptr && *bomptr) {
	    /* We are at the file beginning; skip UTF8-encoded BOM
	     * if present. Sane editors won't put this in on their
	     * own, but e.g. Windows Notepad will do it happily. */
	    if ((unsigned char)c == *bomptr) {
		bomptr++;
		continue;
	    } else {
		/* Do not tolerate partial BOM. */
		if (bomptr != utf8_bom)
		    break;
		/* No BOM at file beginning. Cool. */
		bomptr = NULL;
	    }
	}
	if (c == '\n') {
	    if (file_EOF)
		return 0;
	    comment = 0;
	    continue;
	}
	if (comment || isspace(c))
	    continue;
	if (c == '#' || c == ';') {
	    comment = 1;
	    continue;
	}
	if (!isalpha(c))
	    break;
	var[baselen] = tolower(c);
	if (get_value(var, baselen + 1) < 0)
	    break;
    }

    die("bad config file line %d in %s", file_linenumber, file_name);

    return 0;
}

/**
	read config file

	\param filename config file to read
	\returns 0 if successful
*/
int read_config_from_file(const char *filename)
{
    int ret;
    FILE *f = fopen(filename, "r");

    ret = -1;
    if (f) {
	file_ptr = f;
	file_name = filename;
	file_linenumber = 1;
	file_EOF = 0;
	ret = parse_file();
	fclose(f);
	f = NULL;
	file_name = NULL;
    }

    if (f != NULL) {
	fclose(f);
    }
    return ret;
}

/**
	print error

	\param key key on which error occured
*/
static void die_bad_config(const char *key)
{
    if (file_name)
	die("bad config value for '%s' in %s", key, file_name);
    die("bad config value for '%s'", key);
}

/**
	Get a unit factor like (k)ilo, (M)ega, (G)iga

	\param end pointer to red unit factor
	\param val value to multiply unit factor with
	\returns 1 on success, 0 on error
*/
static int parse_unit_factor(const char *end, unsigned long *val)
{
    if (!*end)
	return 1;
    else if (!strcasecmp(end, "k")) {
	*val *= 1000;
	return 1;
    } else if (!strcasecmp(end, "M")) {
	*val *= 1000 * 1000;
	return 1;
    } else if (!strcasecmp(end, "G")) {
	*val *= 1000 * 1000 * 1000;
	return 1;
    }
    return 0;
}

/**
	Get a value for a given key

	\param key key to search for
	\returns value corresponding to key
*/
static const char *value(const char *key)
{
    char keylower[MAXNAME];
    unsigned int i;

    for (i = 0; i <= strlen(key); ++i) {
	keylower[i] = tolower(key[i]);
    }

    for (i = 0; i < list.size; ++i) {
	if (strcmp(keylower, list.data[i].key) == 0)
	    return list.data[i].value;
    }

    die("bad config file: no setting for '%s' found", key);

    return NULL;
}

/**
	Parse a value as long

	\param value value to parse
	\param ret return value
	\returns 1 on success, 0 on error
*/
static int parse_long(const char *value, long *ret)
{
    if (value && *value) {
	char *end;
	long val = strtol(value, &end, 0);
	unsigned long factor = 1;
	if (!parse_unit_factor(end, &factor))
	    return 0;
	*ret = val * factor;
	return 1;
    }
    return 0;
}

/**
	Parse a value as unsigned long

	\param value value to parse
	\param ret return value
	\returns 1 on success, 0 on error
*/
static int parse_unsigned_long(const char *value, unsigned long *ret)
{
    if (value && *value) {
	char *end;
	unsigned long val = strtoul(value, &end, 0);
	if (!parse_unit_factor(end, &val))
	    return 0;
	*ret = val;
	return 1;
    }
    return 0;
}

/**
	Parse a value as double

	\param value value to parse
	\param ret return value
	\returns 1 on success, 0 on error
*/
static int parse_double(const char *value, double *ret)
{
    if (value && *value) {
	char *end;
	double val = strtod(value, &end);
	unsigned long factor = 1;
	if (!parse_unit_factor(end, &factor))
	    return 0;
	*ret = val * factor;
	return 1;
    }
    return 0;
}

/**
	Get a value as unsigned long

	\param key key
	\param value value to read
	\returns value
*/
/*static unsigned long as_unsigned_long(const char *key, const char *value)
{
	unsigned long ret;
	if (!parse_unsigned_long(value, &ret))
		die_bad_config(key);
	return ret;
}*/

/**
	Get a value as int

	\param key key
	\param value value to read
	\returns value
*/
static int as_int(const char *key, const char *value)
{
    long ret = 0;
    if (!parse_long(value, &ret))
	die_bad_config(key);
    return ret;
}

/**
	Get a value as unsigned int

	\param key key
	\param value value to read
	\returns value
*/
static unsigned int as_unsigned_int(const char *key, const char *value)
{
    unsigned long ret = 0;
    if (!parse_unsigned_long(value, &ret))
	die_bad_config(key);
    return ret;
}

/**
	Get a value as double

	\param key key
	\param value value to read
	\returns value
*/
static double as_double(const char *key, const char *value)
{
    double ret = 0;
    if (!strcmp(value, "NaN")) {
	logging::print_master(
	    (std::string(LOG_WARNING) + "Warning: '" + std::string(key) +
	     "' set to the unphysical 1e300 instead of NaN. NaN is not supported for parameters!\n")
		.c_str());
	return 1e300;
    }
    if (!parse_double(value, &ret))
	die_bad_config(key);
    return ret;
}

/**
	Get a value as bool or int

	\param key key
	\param value value to read
	\param is_bool 1 is value is a boolean, otherwise 0
	\returns value
*/
static int as_bool_or_int(const char *key, const char *value, int *is_bool)
{
    *is_bool = 1;
    if (!value)
	return 1;
    if (!*value)
	return 0;
    if (!strcasecmp(value, "true") || !strcasecmp(value, "yes") ||
	!strcasecmp(value, "on"))
	return 1;
    if (!strcasecmp(value, "false") || !strcasecmp(value, "no") ||
	!strcasecmp(value, "off"))
	return 0;
    *is_bool = 0;
    return as_int(key, value);
}

/**
	Get a value as bool

	\param key key
	\param value value to read
	\returns value
*/
static bool as_bool(const char *key, const char *value)
{
    int discard;
    return !!as_bool_or_int(key, value, &discard);
}

/**
	Get a value as bool to a corresponding key

	\param key key
	\returns value
*/
bool value_as_bool(const char *key) { return as_bool(key, value(key)); }

/**
	Get a value as int to a corresponding key

	\param key key
	\returns value
*/
int value_as_int(const char *key) { return as_int(key, value(key)); }

/**
	Get a value as unsigned int to a corresponding key

	\param key key
	\returns value
*/
unsigned int value_as_unsigned_int(const char *key)
{
    return as_unsigned_int(key, value(key));
}

/**
	Get a value as double to a corresponding key

	\param key key
	\returns value
*/
double value_as_double(const char *key) { return as_double(key, value(key)); }

/**
	Get a value as string to a corresponding key

	\param key key
	\returns value
*/
const char *value_as_string(const char *key) { return value(key); }

/**
	Get a value as bool to a corresponding key if available, else set to
   default

	\param key key
	\param defvalue default value
	\returns value
*/
bool value_as_bool_default(const char *key, bool defvalue)
{
    if (key_exists(key)) {
	return value_as_bool(key);
    } else {
	return defvalue;
    }
}

/**
	Get a value as int to a corresponding key if available, else set to
   default

	\param key key
	\param defvalue default value
	\returns value
*/
int value_as_int_default(const char *key, int defvalue)
{
    if (key_exists(key)) {
	return value_as_int(key);
    } else {
	return defvalue;
    }
}

/**
	Get a value as unsigned int to a corresponding key if available, else
   set to default

	\param key key
	\param defvalue default value
	\returns value
*/
unsigned int value_as_unsigned_int_default(const char *key,
					   unsigned int defvalue)
{
    if (key_exists(key)) {
	return value_as_unsigned_int(key);
    } else {
	return defvalue;
    }
}

/**
	Get a value as double to a corresponding key  if available, else set to
   default

	\param key key
	\param defvalue default value
	\returns value
*/
double value_as_double_default(const char *key, double defvalue)
{
    if (key_exists(key)) {
	return value_as_double(key);
    } else {
	return defvalue;
    }
}

/**
	Get a value as string to a corresponding key  if available, else set to
   default

	\param key key
	\param defvalue default value
	\returns value
*/
const char *value_as_string_default(const char *key, const char *defvalue)
{
    if (key_exists(key)) {
	return value_as_string(key);
    } else {
	return defvalue;
    }
}

/**
	Get a value as t_boundary_condition to a corresponding key  if
   available, else set to default

	\param key key
	\param defvalue default value
	\returns t_boundary_condition
*/
parameters::t_damping_type
value_as_boudary_damping_default(const char *key, const char *defvalue)
{
    char *string_key;
    if (key_exists(key)) {
	string_key = (char *)value_as_string(key);
    } else {
	string_key = (char *)defvalue;
    }

    parameters::t_damping_type boundary_condition;
    switch (tolower(*string_key)) {
    case 'n':
	boundary_condition = parameters::t_damping_type::damping_none;
	break;
    case 'i':
	boundary_condition = parameters::t_damping_type::damping_initial;
	break;
    case 'y': // for legacy compatibility
	boundary_condition = parameters::t_damping_type::damping_initial;
	break;
    case 'm':
	boundary_condition = parameters::t_damping_type::damping_mean;
	break;
    case 'z':
	boundary_condition = parameters::t_damping_type::damping_zero;
	break;
    case 'v':
	boundary_condition = parameters::t_damping_type::damping_visc;
	break;
    default:
	boundary_condition = parameters::t_damping_type::damping_none;
    }
    return boundary_condition;
}

/**
	check if a key exists

	\param key key to check for
	\returns 1 if exists, otherwise 0
*/
int key_exists(const char *key)
{
    char keylower[MAXNAME];
    unsigned int i;

    for (i = 0; i <= strlen(key); ++i) {
	keylower[i] = tolower(key[i]);
    }

    for (i = 0; i < list.size; ++i) {
	if (strcmp(keylower, list.data[i].key) == 0)
	    return 1;
    }

    return 0;
}

/**
	prints the parsed config file. mainly for debugging purposes
*/
void print_parsed_config()
{
    unsigned int i;

    for (i = 0; i < list.size; ++i) {
	printf("[%u]: %s = %s\n", i, list.data[i].key, list.data[i].value);
    }
}

} // namespace config
