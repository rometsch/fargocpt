#pragma once

#include <cstdint>

namespace start_mode
{
enum StartMode { mode_start, mode_restart, mode_auto, mode_debug, mode_none };
extern StartMode mode;
extern std::int32_t restart_from;
extern std::int32_t restart_debug;

void configure_start_mode();
} // namespace start_mode
