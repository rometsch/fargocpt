#ifndef OPTIONS_H
#define OPTIONS_H

namespace options {

void usage(int argc, char** argv);
void parse(int argc, char** argv);

extern bool restart;
extern unsigned int restart_from;
extern bool memory_usage;
extern char *parameter_file;
extern bool disable;
}

#endif // OPTIONS_H
