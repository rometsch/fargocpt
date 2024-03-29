#pragma once

#include "data.h"

// Computations of physical quantities from other physical quantities.

namespace compute {
    
    void midplane_density(t_data &data, const double current_time);
    void kappa_eff(t_data &data);
    void toomreQ(t_data &data);
}
