#pragma once

#include <cstdint>

namespace start_mode
{
enum StartMode { mode_start, mode_restart, mode_auto, mode_none };
extern StartMode mode;
extern std::int32_t restart_from;

void configure_start_mode();
std::int32_t get_latest_output_num();
} // namespace start_mode
