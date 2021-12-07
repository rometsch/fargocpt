#ifndef POLARGRID_H
#define POLARGRID_H

#include "units.h"
#include <stddef.h>
#include <stdio.h>

class t_data;

class t_polargrid
{
  private:
    units::t_unit *m_unit;

    /// scalar or vector quantity (scalar uses grid from 0 to Nrad-1, vector
    /// from 0 to Nrad)
    bool m_scalar;
    /// Name of the PolarGrid (can be "dens", "vrad", "vtheta" or "label").
    char *m_name;
    /// write min/max values on 1D output
    bool m_write_max_max_1D;
    /// write 1D?
    bool m_write_1D;
    /// write 2D?
    bool m_write_2D;
    /// calculate grid even if m_write1D/m_write2D is false when write is called
    bool m_calculate_on_write;
    /// callback function to be called before write operations
    void (*m_do_before_write)(t_data &, unsigned int, bool);
    /// set all entries to zero after write operation
    bool m_clear_after_write;
    /// for mass flow we want to integrate azimuthally instead of averaging
    bool m_integrate_1D_for_write;

  public:
    t_polargrid();
    ~t_polargrid();

    // setter
    void set_name(const char *name);
    void set_unit(units::t_unit &unit);
    void set_size(ptrdiff_t size_radial, ptrdiff_t size_azimuthal);
    void set_scalar(bool value);
    void set_vector(bool value);
    inline void set_write_1D(bool value) { m_write_1D = value; }
    inline void set_write_2D(bool value) { m_write_2D = value; }
    inline void set_write(bool value, bool write1D)
    {
	set_write_1D(value && write1D);
	set_write_2D(value);
    }
    inline void set_do_before_write(void (*value)(t_data &, unsigned int, bool))
    {
	m_do_before_write = value;
    }
    inline void set_clear_after_write(bool value)
    {
	m_clear_after_write = value;
    }
    inline void set_integrate_azimuthally_for_1D_write(bool value)
    {
	m_integrate_1D_for_write = value;
    }

    // getter
    inline ptrdiff_t get_max_radial() const
    {
	return m_scalar ? Nrad - 1 : Nrad;
    }
    inline ptrdiff_t get_max_azimuthal() const { return Nsec - 1; }
    inline ptrdiff_t get_size_radial() const
    {
	return m_scalar ? Nrad : Nrad + 1;
    }
    inline ptrdiff_t get_size_azimuthal() const { return Nsec; }
    inline char *get_name() const { return m_name; }
    inline units::t_unit *get_unit() const { return m_unit; }

    inline bool get_write_1D() const { return m_write_1D; }
    inline bool get_write_2D() const { return m_write_2D; }
    inline bool get_write() const { return m_write_1D || m_write_2D; }
    inline bool get_clear_after_write() const { return m_clear_after_write; }
    inline bool get_integrate_azimuthally_for_1D_write() const
    {
	return m_integrate_1D_for_write;
    }

    inline bool is_scalar() const { return m_scalar; }
    inline bool is_vector() const { return !m_scalar; }

    void clear();

    void write_polargrid(unsigned int number, t_data &data, bool debug);
    // 2D read/write
    void write2D(const unsigned int number, const bool debug) const;
    void read2D(const char *filename);
    void read2D(unsigned int number, bool debug);

    // 1D read/write
    void write1D(unsigned int number) const;
    void read1D(const char *filename, bool skip_min_max);
    void read1D(unsigned int number, bool skip_min_max);

    unsigned int bytes_needed_1D() const;
    unsigned int bytes_needed_2D() const;

    unsigned int get_memory_usage(ptrdiff_t size_radial,
				  ptrdiff_t size_azimuthal) const;

    inline unsigned int cell(ptrdiff_t nRadial, ptrdiff_t nAzimuthal) const
    {
	return nAzimuthal + (nRadial * Nsec);
    }

    inline double &operator()(ptrdiff_t nRadial, ptrdiff_t nAzimuthal)
    {
	return Field[cell(nRadial, nAzimuthal)];
    }
    inline double operator()(ptrdiff_t nRadial, ptrdiff_t nAzimuthal) const
    {
	return Field[cell(nRadial, nAzimuthal)];
    }

    double get_max() const;
    t_polargrid &operator*=(double);
    t_polargrid &operator/=(double);

    /// Radial size of the grid, in number of zones
    ptrdiff_t Nrad;
    /// Azimuthal size of the grid, in number of zones
    ptrdiff_t Nsec;
    /// Pointer to the array of Nrad*Nsec reals (e.g., density, etc.)
    double *Field;
    /// Should be 1 if grid is already initialised (with data), else 0
    int Initialised;
    /// 1 if grid should write out each nth timeStep
    int WriteOut;
};

// for compatibility
typedef t_polargrid PolarGrid;
#endif // POLARGRID_H
