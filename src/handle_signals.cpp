#include <csignal>
#include <functional>
#include "backtrace.h"
#include "handle_signals.h"
#include "global.h"
#include "logging.h"
#include "LowTasks.h"

static void handleSIGUSR1(__attribute__((unused)) int signum) { PRINT_SIG_INFO = 1; }

static void handleSIGUSR2(__attribute__((unused)) int signum) { PrintTrace(); }

static void handleSIGTERM(__attribute__((unused)) int signum) { SIGTERM_RECEIVED = 1; }


static void registerSIGUSR1() {
    struct sigaction sa;

    sa.sa_handler = handleSIGUSR1;
    sa.sa_flags = 0;
    sigemptyset(&sa.sa_mask);

    if (sigaction(SIGUSR1, &sa, NULL) == -1) {
        perror("sigaction SIGUSR1");
        exit(1);
    }
}

static void registerSIGUSR2() {
    struct sigaction sa;

    sa.sa_handler = handleSIGUSR2;
    sa.sa_flags = 0;
    sigemptyset(&sa.sa_mask);

    if (sigaction(SIGUSR2, &sa, NULL) == -1) {
        perror("sigaction SIGUSR2");
        exit(1);
    }
}

static void registerSIGTERM() {
    struct sigaction sa;

    sa.sa_handler = handleSIGTERM;
    sa.sa_flags = 0;
    sigemptyset(&sa.sa_mask);

    if (sigaction(SIGTERM, &sa, NULL) == -1) {
        perror("sigaction SIGTERM");
        exit(1);
    }
}


void register_signal_handlers()
{
    registerSIGUSR1();
    registerSIGUSR2();
	SIGTERM_RECEIVED = 0;
    registerSIGTERM();
}
