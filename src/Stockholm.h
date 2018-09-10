#ifndef _STOCKHOLM_H_
#define _STOCKHOLM_H_

#include "types.h"
#include "data.h"

Force ComputeForceStockholm(t_data &data, double x, double y, double rsmoothing, double mass);
Pair MassInOut(t_data &data, double a);
void UpdateLogStockholm(t_data &data, int outputnb, double time);

#endif
