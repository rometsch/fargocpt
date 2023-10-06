#ifndef CONSTANTS_H
#define CONSTANTS_H

/// gravitational constant in code units
// const double G = 1.0;

/// Universal Gas Constant in code units
// const double R = 1.0;

// /// central mass in code units
// const double M = 1.0;

// overlap cells (Zeus-like overlap kernel. 2:transport; 2: source, 1:viscous
// stress)
//const unsigned int CPUOVERLAP = 5;
// overlap cells (Zeus-like overlap kernel. With leapfrog: 2:transport; 2x2: source, 1x2:viscous
// stress)
const unsigned int CPUOVERLAP = 7; // I don't need 8 cells to get bitwise identical result with fargo.par
const unsigned int GHOSTCELLS_A = 2;
const unsigned int GHOSTCELLS_B = 1;

namespace constants
{

class t_constant
{
  private:
    char *m_symbol;
    double m_code_value;
    double m_cgs_value;
    char *m_cgs_unit_symbol;

  public:
    t_constant();
    ~t_constant();

    // setter
    void set_symbol(const char *);
    void set_code_value(double);
    void set_cgs_value(double);
    void set_cgs_unit_symbol(const char *);

    // getter
    const char *get_symbol(void) const;
    double get_code_value(void) const;
    double get_cgs_value(void) const;
    const char *get_cgs_unit_symbol(void) const;

    // operator
    inline operator const double &() const { return m_code_value; }
    /* inline operator double&() { return m_cgs_value; } */
};

/// gravitational constant
extern t_constant _G;
extern double global_cgs_G;

/// Boltzmann constant
extern t_constant k_B;
/// molecular mass
extern t_constant m_u;
/// Planck constant
extern t_constant h;
/// speed of light
extern t_constant c;
/// specific gas constant
extern t_constant _R;
/// Stefan-Boltzmann constant
extern t_constant sigma;

/// electron volt
extern t_constant eV;
/// electron mass
extern t_constant m_e;
/// Hydrogen atom mass
extern t_constant m_H;

// _G and _R are always defined with _ and are linked to G and R with compilied
// without NDEBUG otherswise to const double 1
#ifndef NDEBUG
extern t_constant &G;
extern t_constant &R;
#else
const double R = 1.0;
const double G = 1.0;
#endif

void initialize_constants();
void calculate_constants_in_code_units();
void print_constants();

} // namespace constants

#endif // CONSTANTS_H
