#ifndef RADIALGRID_H
#define RADIALGRID_H

#include "units.h"
#include <stddef.h>
#include <stdlib.h>
#include <string>

class t_data;
class t_polargrid;

class t_radialgrid
{
  private:
    units::t_unit *m_unit;

    /// scalar or vector quantity (scalar uses grid from 0 to Nrad-1, vector
    /// from 0 to Nrad)
    bool m_scalar;
    /// Name of the PolarGrid (can be "dens", "vrad", "vtheta" or "label").
    char *m_name;
    /// write 1D?
    bool m_write_1D;
    /// calculate grid even if m_write1D/m_write2D is false when write is called
    bool m_calculate_on_write;
    /// callback function to be called before write operations
    void (*m_do_before_write)(t_data &, unsigned int, bool);
    /// set all entries to zero after write operations
    bool m_clear_after_write;
    /// Radial size of the grid, in number of zones
    unsigned int m_size_radial;
    /// pointer to actual grid data
    double *m_data;

  public:
    t_radialgrid();
    ~t_radialgrid();

    // setter
    void set_name(const char *name);
    void set_unit(units::t_unit &unit);
    void set_size(ptrdiff_t size_radial);
    void set_scalar(bool value);
    void set_vector(bool value);
    inline void set_write_1D(bool value) { m_write_1D = value; }
    inline void set_write(bool value) { set_write_1D(value); }
    inline void set_clear_after_write(bool value)
    {
	m_clear_after_write = value;
    }
    inline void set_do_before_write(void (*value)(t_data &, unsigned int, bool))
    {
	m_do_before_write = value;
	m_calculate_on_write = true;
    }

    // getter
    inline ptrdiff_t get_max_radial() const
    {
	return m_scalar ? m_size_radial - 1 : m_size_radial;
    }
    inline ptrdiff_t get_size_radial() const
    {
	return m_scalar ? m_size_radial : m_size_radial + 1;
    }
    inline char *get_name() const { return m_name; }
    inline units::t_unit *get_unit() const { return m_unit; }
    inline bool get_write_1D() const { return m_write_1D; }

    inline bool is_scalar() const { return m_scalar; }
    inline bool is_vector() const { return !m_scalar; }

    void clear();

	void write_radialgrid(unsigned int number, t_data &data);
    void write(std::string filename, unsigned int number, t_data &data,
	       bool one_file, bool force_write);
    // 1D read/write
    void write1D(unsigned int timestep) const;
    void write1D(std::string filename, bool one_file) const;
    void read1D(const char *filename);
    void read1D(unsigned int number);

    unsigned int bytes_needed_1D();

    unsigned int get_memory_usage(ptrdiff_t size_radial);

    inline unsigned int cell(ptrdiff_t n_radial) const { return n_radial; }

    inline double &operator()(ptrdiff_t n_radial) { return m_data[n_radial]; }
    inline double operator()(ptrdiff_t n_radial) const
    {
	return m_data[n_radial];
    }

    t_radialgrid &operator=(const t_polargrid &);

    /*
		    t_polargrid& operator*=(double);
		    t_polargrid& operator/=(double);
    */
};

#endif // RADIALGRID_H
