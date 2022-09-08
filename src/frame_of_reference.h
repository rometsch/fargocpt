#ifndef FRAME_OF_REFERENCE_H
#define FRAME_OF_REFERENCE_H

#include "data.h"
#include "types.h"

namespace frame_of_reference{ 

extern double planet_corot_ref_old_x;
extern double planet_corot_ref_old_y;

extern Pair IndirectTerm;
extern Pair IndirectTermDisk;
extern Pair IndirectTermPlanets;

void init_corotation(t_data &data);
void handle_corotation(t_data &data, const double dt);

}

#endif // FRAME_OF_REFERENCE_H
