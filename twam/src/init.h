#ifndef _INIT_H_
#define _INIT_H_

#include "data.h"
#include "types.h"

void init_radialarrays(void);
void resize_radialarrays(unsigned int size);

void init_physics(t_data &data);
void init_shakura_sunyaev(t_data &data);

void init_gas_density(t_data &data);
void init_gas_energy(t_data &data);
void init_gas_velocities(t_data &data);

#endif
