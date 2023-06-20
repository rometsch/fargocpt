#pragma once

#include "data.h"
#include "types.h"

namespace refframe{ 

extern double planet_corot_ref_old_x;
extern double planet_corot_ref_old_y;
extern double OmegaFrame;
extern double FrameAngle;

extern Pair IndirectTerm;
extern Pair IndirectTermDisk;
extern Pair IndirectTermPlanets;

void ComputeIndirectTermFully();

void ComputeIndirectTermNbody(t_data &data, const double current_time, const double dt);
void ComputeIndirectTermDisk(t_data &data);
void ComputeIndirectTermNbodyEuler(t_data &data);


void init_corotation(t_data &data);
void handle_corotation(t_data &data, const double dt);

}
