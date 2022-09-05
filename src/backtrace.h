#ifndef BACKTRACE_H
#define BACKTRACE_H

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#ifndef __USE_GNU
#define __USE_GNU
#endif

#define BACKTRACE_MAXDEPTH 50

void PrintTrace(bool = true);

#endif // BACKTRACE_H
