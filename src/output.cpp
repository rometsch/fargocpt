#include "output.h"
#include "Force.h"
#include "LowTasks.h"
#include "Pframeforce.h"
#include "SideEuler.h"
#include "Theo.h"
#include "constants.h"
#include "global.h"
#include "logging.h"
#include "options.h"
#include "parameters.h"
#include "particles/particles.h"
#include "quantities.h"
#include "start_mode.h"
#include "stress.h"
#include "util.h"
#include "viscosity/viscosity.h"
#include "frame_of_reference.h"
#include "simulation.h"
#include "fld.h"

#include <dirent.h>

#include "unistd.h" // for access()
#include <cfloat>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <sys/stat.h>
#include <sys/statvfs.h>

namespace output
{

std::string snapshot_dir = "";
std::string last_snapshot_dir = "";
std::string outdir = "";

const std::map<const std::string, const int> quantities_file_column_v2_4 = {
    {"time step", 0},
    {"analysis time step", 1},
    {"physical time", 2},
    {"mass", 3},
    {"radius", 4},
    {"angular momentum", 5},
    {"total energy", 6},
    {"internal energy", 7},
    {"kinematic energy", 8},
    {"potential energy", 9},
    {"radial kinetic energy", 10},
    {"azimuthal kinetic energy", 11},
    {"eccentricity", 12},
    {"periastron", 13},
    {"viscous dissipation", 14},
    {"luminosity", 15},
    {"pdivv", 16},
	{"inner boundary mass inflow", 17},
	{"inner boundary mass outflow", 18},
	{"outer boundary mass inflow", 19},
	{"outer boundary mass outflow", 20},
	{"wave damping inner mass creation", 21},
	{"wave damping inner mass removal", 22},
	{"wave damping outer mass creation", 23},
	{"wave damping outer mass removal", 24},
	{"density floor mass creation", 25},
	{"aspect ratio", 26},
	{"indirect term nbody x", 27},
	{"indirect term nbody y", 28},
	{"indirect term disk x", 29},
	{"indirect term disk y", 30}
	};
static const auto quantities_file_column = quantities_file_column_v2_4;

const std::map<const std::string, const std::string> quantities_file_variables =
	{{"physical time", "time"},
	{"mass", "mass"},
	{"radius", "length"},
	{"angular momentum", "angular_momentum"},
	{"total energy", "energy"},
	{"internal energy", "energy"},
	{"kinematic energy", "energy"},
	{"potential energy", "energy"},
	{"viscous dissipation", "power"},
	{"luminosity", "power"},
	{"pdivv", "pressure per time"},
	{"radial kinetic energy", "energy"},
	{"azimuthal kinetic energy", "energy"},
	{"inner boundary mass inflow", "mass"},
	{"inner boundary mass outflow", "mass"},
	{"outer boundary mass inflow", "mass"},
	{"outer boundary mass outflow", "mass"},
	{"wave damping mass creation", "mass"},
	{"wave damping mass removal", "mass"},
	{"density floor mass creation", "mass"},
	{"wave damping inner mass creation", "mass"},
	{"wave damping inner mass removal", "mass"},
	{"wave damping outer mass creation", "mass"},
	{"wave damping outer mass removal", "mass"},
	{"time step", "1"},
	{"analysis time step", "1"},
	{"omega frame", "frequency"},
	{"lost mass", "mass"},
	{"frame angle", "frequency"},
	{"eccentricity", "1"},
	{"periastron", "1"},
	{"aspect ratio", "1"},
	{"indirect term nbody x", "acceleration"},
	{"indirect term nbody y", "acceleration"},
	{"indirect term disk x", "acceleration"},
	{"indirect term disk y", "acceleration"}};


void check_free_space(t_data &data)
{
    DIR *directory_pointer;
    struct statvfs fiData;

    // check if output directory exists
    if ((directory_pointer = opendir(outdir.c_str())) == nullptr) {
	logging::print_master(LOG_ERROR "Output directory %s doesn't exist!\n",
			      outdir.c_str());
	die("Not output directory!");
	return; // needed so that compuler understands directory_pointer !=
		// nullptr
    }

    closedir(directory_pointer);

    unsigned long int space_needed = 0l;
    unsigned int number_of_files = 0l;

    // go thru all polargrids in output_polargrids and check if they are going
    // to be put out
    for (unsigned int i = 0; i < t_data::N_POLARGRID_TYPES; ++i) {
	if (data[(t_data::t_polargrid_type)i].get_write_1D()) {
	    space_needed += data[(t_data::t_polargrid_type)i].bytes_needed_1D();
	    number_of_files += 1;
	}

	if (data[(t_data::t_polargrid_type)i].get_write_2D()) {
	    space_needed += data[(t_data::t_polargrid_type)i].bytes_needed_2D();
	    number_of_files += 1;
	}
    }

    space_needed *= parameters::NTOT / parameters::NINTERM;
    number_of_files *= parameters::NTOT / parameters::NINTERM;

    logging::print_master(LOG_INFO "Output information:\n");
    logging::print_master(LOG_INFO "   Output directory: %s\n", outdir.c_str());
    logging::print_master(LOG_INFO "    Number of files: %u\n",
			  number_of_files);
    logging::print_master(LOG_INFO "  Total output size: %.2f GB\n",
			  (double)space_needed / 1024.0 / 1024.0 / 1024.0);

    if ((statvfs(outdir.c_str(), &fiData)) < 0) {
	logging::print_master(
	    LOG_WARNING
	    "Couldn't stat filesystem. You have to check for enough free space manually!\n");
    } else {
	// free space of device, more precisely number of space available to
	// non-priv processes (like us)
	unsigned long int free_space = fiData.f_bavail * fiData.f_frsize;

	logging::print_master(LOG_INFO "    Space Available: %.2f GB\n",
			      (double)free_space / 1024.0 / 1024.0 / 1024.0);

	if (space_needed > free_space) {
	    logging::print_master(
		LOG_WARNING
		"There is not enough space for all outputs! The program will fail at same point!\n");
	}
    }
}

static void register_output(const std::string &snapshot_id)
{
    if (CPU_Master) {
	const std::string filename = outdir + "snapshots/list.txt";
	std::ofstream output_list(filename, std::ios_base::app);
	output_list << snapshot_id << std::endl;
	output_list.close();
    }
}

static void copy_parameters_to_snapshot_dir()
{
    if (CPU_Master) {
	const std::string src_file = options::parameter_file;
	const std::string dst_file = snapshot_dir + "/config.yml";
	std::filesystem::copy_file(src_file, dst_file);
    }
}

void write_output_version()
{
    if (CPU_Master) {
	const std::string filename = outdir + "/fargocpt_output_v1_2";
	std::ofstream versionfile(filename, std::ios_base::app);
	versionfile.close();
    }
}

void cleanup_autosave()
{
    const auto s = last_snapshot_dir;
    const auto l = std::string("autosave").length();
    if (s.length() > l) {
	const auto ts = s.substr(s.length() - l);
	if (ts.compare("autosave") == 0) {
	    std::filesystem::remove_all(s);
	}
    }
}

void write_full_output(t_data &data, const std::string &snapshot_id,
		       const bool register_snapshot)
{

    snapshot_dir = outdir + "snapshots/" + snapshot_id;
    delete_directory_if_exists(snapshot_dir);
    ensure_directory_exists(snapshot_dir);
    MPI_Barrier(MPI_COMM_WORLD);

    // Enable output of Qplus / Qminus for bitwise exact restarting.
	if (parameters::bitwise_exact_restarting && !parameters::Locally_Isothermal) {
    if (!data[t_data::QPLUS].get_write()) {
	data[t_data::QPLUS].set_write(true, false);
    }
    if (!data[t_data::QMINUS].get_write()) {
	data[t_data::QMINUS].set_write(true, false);
    }
	}

    if (parameters::variableGamma) {
	if (!data[t_data::GAMMAEFF].get_write()) {
	    data[t_data::GAMMAEFF].set_write(true, false);
	    data[t_data::MU].set_write(true, false);
	    data[t_data::GAMMA1].set_write(true, false);
	}
    }

	if (parameters::radiative_diffusion_test_module) {
		fld2::Ka.set_name("Ka");
		fld2::Ka.write2D();
		fld2::Kb.set_name("Kb");
		fld2::Kb.write2D();
		fld2::A.set_name("A");
		fld2::A.write2D();
		fld2::B.set_name("B");
		fld2::B.write2D();
		fld2::C.set_name("C");
		fld2::C.write2D();
		fld2::D.set_name("D");
		fld2::D.write2D();
		fld2::E.set_name("E");
		fld2::E.write2D();
		fld2::Trad.set_name("Trad");
		fld2::Trad.write2D();
		fld2::Erad.set_name("Erad");
		fld2::Erad.write2D();
		fld2::Xold.set_name("Xold");
		fld2::Xold.write2D();
	} else {
		fld::Ka.set_name("Ka");
		fld::Ka.write2D();
		fld::Kb.set_name("Kb");
		fld::Kb.write2D();
		fld::A.set_name("A");
		fld::A.write2D();
		fld::B.set_name("B");
		fld::B.write2D();
		fld::C.set_name("C");
		fld::C.write2D();
		fld::D.set_name("D");
		fld::D.write2D();
		fld::E.set_name("E");
		fld::E.write2D();
		fld::Trad.set_name("Trad");
		fld::Trad.write2D();
		fld::Erad.set_name("Erad");
		fld::Erad.write2D();
		fld::Xold.set_name("Xold");
		fld::Xold.write2D();
	}

    // write polar grids
    output::write_grids(data, sim::N_snapshot, sim::N_hydro_iter, sim::PhysicalTime);
    // write planet data
    data.get_planetary_system().write_planets(0);
	data.get_massflow_tracker().write_to_file();

    // write misc stuff (important for resuming)
    output::write_misc();
    // write particles
    if (parameters::integrate_particles) {
	particles::write();
    }

    if (register_snapshot) {
    // write time info for coarse output
    output::write_snapshot_time();
	register_output(snapshot_id);
    }

    copy_parameters_to_snapshot_dir();

    MPI_Barrier(MPI_COMM_WORLD);
}

void write_grids(t_data &data, int index, int iter, double phystime)
{
    logging::print_master(
	LOG_INFO "Writing output %s, Timestep Number %d, Physical Time %f.\n",
	snapshot_dir.c_str(), index, iter, phystime);

    // go thru all grids and write them
    for (unsigned int i = 0; i < t_data::N_POLARGRID_TYPES; ++i) {
	data[(t_data::t_polargrid_type)i].write_polargrid(data);
    }

    // go thru all grids and write them
    for (unsigned int i = 0; i < t_data::N_RADIALGRID_TYPES; ++i) {
	data[(t_data::t_radialgrid_type)i].write_radialgrid(index, data);
    }
}

/**

*/
void write_quantities(t_data &data, bool force_update)
{
    FILE *fd = 0;
    std::string filename = outdir + "monitor/Quantities.dat";
    auto fd_filename = filename.c_str();
    static bool fd_created = false;

    if (CPU_Master) {

	// check if file exists and we restarted
	if ((start_mode::mode == start_mode::mode_restart) && !(fd_created)) {
	    if (access(fd_filename, W_OK) != -1) {
		fd_created = true;
	    }
	}

	// open logfile
	if (!fd_created) {
	    fd = fopen(fd_filename, "w");
	} else {
	    fd = fopen(fd_filename, "a");
	}
	if (fd == NULL) {
	    logging::print_master(
		LOG_ERROR "Can't write 'Quantities.dat' file. Aborting.\n");
	    PersonalExit(1);
	}

	if (!fd_created) {
	    // print header
	    fprintf(fd, "#FargoCPT quantities file\n");
	    fprintf(fd, "#version: 2.4\n");
	    fprintf(fd, "%s",
		    text_file_variable_description(quantities_file_column,
						   quantities_file_variables)
			.c_str());
	    fd_created = true;
	}
    }

    double quantities_limit_radius;
    if (parameters::quantities_radius_limit < 0.0) {
	auto &primary = data.get_planetary_system().get_planet(0);
	// distance to primary is distance to secondary for the primary
	quantities_limit_radius = primary.get_distance_to_primary() *
				  primary.get_dimensionless_roche_radius();
    } else {
	quantities_limit_radius = parameters::quantities_radius_limit;
    }

    const auto disk_quantities = reduce_disk_quantities(
	data, sim::N_snapshot, force_update, quantities_limit_radius);
    const double disk_eccentricity = disk_quantities[0];
    const double disk_periastron = disk_quantities[1];

    // computate absolute deviation from start values (this has to be done on
    // all nodes!)

    const double totalMass =
	quantities::gas_total_mass(data, quantities_limit_radius);

    if (totalMass <= 0.0) { // If roche lobe is smaller than RMIN
	quantities::gas_total_mass(data, RMAX);
    }

    const double diskRadius = quantities::gas_disk_radius(data, totalMass);
    const double totalAngularMomentum =
	quantities::gas_angular_momentum(data, quantities_limit_radius);
    const double internalEnergy =
	quantities::gas_internal_energy(data, quantities_limit_radius);
    const double qplus =
	quantities::gas_viscous_dissipation(data, quantities_limit_radius);
    const double qminus =
	quantities::gas_luminosity(data, quantities_limit_radius);
    const double kinematicEnergy =
	quantities::gas_kinematic_energy(data, quantities_limit_radius);
    const double radialKinematicEnergy =
	quantities::gas_radial_kinematic_energy(data, quantities_limit_radius);
    const double azimuthalKinematicEnergy =
	quantities::gas_azimuthal_kinematic_energy(data,
						   quantities_limit_radius);

    if (!parameters::body_force_from_potential) {
	CalculateNbodyPotential(data, sim::PhysicalTime);
    }

	const double gravitationalEnergy =
	-quantities::gas_reduce_mass_average(data, data[t_data::POTENTIAL], quantities_limit_radius);

    const double totalEnergy =
	internalEnergy + kinematicEnergy + gravitationalEnergy;

	if(!(parameters::heating_star_enabled || parameters::self_gravity)){
	quantities::compute_aspectratio(data, sim::N_snapshot, force_update);
	}
	const double scale_height =
	quantities::gas_reduce_mass_average(data, data[t_data::ASPECTRATIO], quantities_limit_radius);

    double pdivv_total = 0.0;
    double InnerBoundaryInflow = 0.0;
    double InnerBoundaryOutflow = 0.0;
    double OuterBoundaryInflow = 0.0;
    double OuterBoundaryOutflow = 0.0;
	double InnerWaveDampingMassCreation = 0.0;
	double OuterWaveDampingMassCreation = 0.0;
	double InnerWaveDampingMassRemoval = 0.0;
	double OuterWaveDampingMassRemoval = 0.0;
    double FloorMassCreation = 0.0;

    MPI_Reduce(&data.pdivv_total, &pdivv_total, 1, MPI_DOUBLE, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    MPI_Reduce(&MassDelta.InnerBoundaryInflow, &InnerBoundaryInflow, 1, MPI_DOUBLE, MPI_SUM,
	       0, MPI_COMM_WORLD);
    MPI_Reduce(&MassDelta.InnerBoundaryOutflow, &InnerBoundaryOutflow, 1, MPI_DOUBLE, MPI_SUM,
	       0, MPI_COMM_WORLD);
    MPI_Reduce(&MassDelta.OuterBoundaryInflow, &OuterBoundaryInflow, 1, MPI_DOUBLE, MPI_SUM,
	       0, MPI_COMM_WORLD);
    MPI_Reduce(&MassDelta.OuterBoundaryOutflow, &OuterBoundaryOutflow, 1, MPI_DOUBLE, MPI_SUM,
	       0, MPI_COMM_WORLD);
	MPI_Reduce(&MassDelta.InnerWaveDampingMassCreation, &InnerWaveDampingMassCreation, 1,
	       MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&MassDelta.OuterWaveDampingMassCreation, &OuterWaveDampingMassCreation, 1,
		   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&MassDelta.InnerWaveDampingMassRemoval, &InnerWaveDampingMassRemoval, 1,
	       MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&MassDelta.OuterWaveDampingMassRemoval, &OuterWaveDampingMassRemoval, 1,
		   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&MassDelta.FloorMassCreation, &FloorMassCreation, 1, MPI_DOUBLE, MPI_SUM,
	       0, MPI_COMM_WORLD);

    if (CPU_Master) {
	// print to logfile
	fprintf(
	    fd,
		"%u\t%u\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\n",
	    sim::N_snapshot, sim::N_monitor, sim::PhysicalTime, totalMass, diskRadius,
	    totalAngularMomentum, totalEnergy, internalEnergy, kinematicEnergy,
	    gravitationalEnergy, radialKinematicEnergy,
	    azimuthalKinematicEnergy, disk_eccentricity, disk_periastron, qplus,
	    qminus, pdivv_total, 
		InnerBoundaryInflow, InnerBoundaryOutflow, 
		OuterBoundaryInflow, OuterBoundaryOutflow, 
		InnerWaveDampingMassCreation, InnerWaveDampingMassRemoval,
		OuterWaveDampingMassCreation, OuterWaveDampingMassRemoval,
	    FloorMassCreation, 
		scale_height,
		refframe::IndirectTermPlanets.x,
		refframe::IndirectTermPlanets.y,
		refframe::IndirectTermDisk.x,
		refframe::IndirectTermDisk.y
		);

	// close file
	fclose(fd);
    }
    // set mass delta to 0
    MassDelta.reset();
}

/**
	log misc. data
*/

void write_misc()
{
    if (!CPU_Master) {
	return;
    }

    std::ofstream wf;
    const std::string filename = snapshot_dir + "/misc.bin";

    wf = std::ofstream(filename, std::ios::out | std::ios::binary);

    if (!wf.is_open()) {
	logging::print_master(
	    LOG_ERROR "Can't write '%s' file in \"write_misc\". Aborting.\n",
	    filename.c_str());
	PersonalExit(1);
    }

	misc_entry misc;
	memset(&misc, 0, sizeof(misc_entry));

	
    misc.timestep = sim::N_snapshot;
    misc.nTimeStep = sim::N_monitor;
    misc.PhysicalTime = sim::PhysicalTime;
    misc.OmegaFrame = refframe::OmegaFrame;
    misc.FrameAngle = refframe::FrameAngle;
    misc.last_dt = sim::last_dt;
    misc.N_iter = sim::N_hydro_iter;

    wf.write((char *)&misc, sizeof(misc));

    wf.close();
}

std::string get_version(std::string filename)
{
    std::string version;
    std::ifstream infile(filename);

    if (infile.fail()) {
	logging::print_master(LOG_ERROR "Error: File %s cannot be opened!\n",
			      filename.c_str());
	PersonalExit(1);
    }

    std::string line_start;
    while (infile >> line_start) {
	// is it the version line
	if (line_start == "#version:") {
	    infile >> version;
	    return version;
	} else if (line_start.substr(0, 1) == "#") {
	    // were still in the header, jump to next line
	    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	} else {
	    // not a comment line anymore, so nothing was found
	    break;
	}
    }
    // nothing found up until here, assume dafault version 1.
    return "1";
}

static std::string unit_descriptor(double value, std::string unit)
{
    // produce a string containing the pair of value and unit as
    // a string such as '1.7823468234...e16 g'
    // i.e. the number with format #.16e
    std::stringstream us;
    us.precision(16);
    us << std::scientific << value << " " << unit;
    return us.str();
}

std::string text_file_variable_description(
    const std::map<const std::string, const int> &variables,
    const std::map<const std::string, const std::string> &units)
{
    // construct a header string describing each variable in
    // its own line including the column and its unit. e.g.
    // #variable: 1 | PhysicalTime | s

    std::map<int, std::string> vars_by_column;
    for (auto const &ent : variables) {
	std::string name = ent.first;
	int column = ent.second;
	vars_by_column[column] = name;
    }

    // build a map with all code units
    // need to do this runtime since units can be set in
    // parameter file
    // could also do this somewhere in initialization...
    std::map<std::string, std::string> unit_descriptors = {
	{"mass", units::mass.get_cgs_factor_symbol()},
	{"mass accretion rate",
	 units::mass_accretion_rate.get_cgs_factor_symbol()},
	{"angular_momentum", units::angular_momentum.get_cgs_factor_symbol()},
	{"time", units::time.get_cgs_factor_symbol()},
	{"energy", units::energy.get_cgs_factor_symbol()},
	{"frequency", unit_descriptor(1.0 / units::time, "1/s")},
	{"1", "1"},
	{"length", units::length.get_cgs_factor_symbol()},
	{"velocity", units::velocity.get_cgs_factor_symbol()},
	{"acceleration", units::acceleration.get_cgs_factor_symbol()},
	{"power", units::power.get_cgs_factor_symbol()},
	{"specific power", unit_descriptor(units::power.get_cgs_factor() /
					       units::length.get_cgs_factor() /
					       units::length.get_cgs_factor(),
					   "erg cm2/s/g")},
	{"potential", units::potential.get_cgs_factor_symbol()},
	{"pressure per time", unit_descriptor(units::pressure.get_cgs_factor() /
						  units::time.get_cgs_factor(),
					      "dyn/cm/s")},
	{"torque", unit_descriptor(units::torque.get_cgs_factor(),
				   units::torque.get_cgs_symbol())}};

    std::string var_descriptor;
    for (auto const &ent : vars_by_column) {
	std::string column = std::to_string(ent.first);
	std::string name = ent.second;
	std::string unit = units.at(name);
	var_descriptor += "#variable: " + column + " | " + name + " | " +
			  unit_descriptors[unit] + "\n";
    }
    return var_descriptor;
}

int load_misc()
{
    std::string filename;
    filename = snapshot_dir + "/misc.bin";

    std::ifstream rf(filename, std::ios::in | std::ios::binary);

    if (!rf.is_open()) {
	logging::print_master(
	    LOG_ERROR "Can't read '%s' file in \"get_misc\". Aborting.\n",
	    filename.c_str());
	PersonalExit(1);
    }

    misc_entry misc;

    rf.read((char *)&misc, sizeof(misc_entry));

    sim::N_snapshot = misc.timestep;
    sim::N_monitor = misc.nTimeStep;
    sim::PhysicalTime = misc.PhysicalTime;
    refframe::OmegaFrame = misc.OmegaFrame;
    refframe::FrameAngle = misc.FrameAngle;
    sim::last_dt = misc.last_dt;
    sim::N_hydro_iter = misc.N_iter;

    rf.close();
    return misc.timestep;
}

void write_torques(t_data &data, bool force_update)
{
    const int ns = data[t_data::SIGMA].Nsec;
    const double *cell_center_x = CellCenterX->Field;
    const double *cell_center_y = CellCenterY->Field;

    // do everything for all planets/stars
    for (unsigned int n_planet = 0;
	 n_planet < data.get_planetary_system().get_number_of_planets();
	 ++n_planet) {
	t_planet &planet = data.get_planetary_system().get_planet(n_planet);

	const double x = planet.get_x();
	const double y = planet.get_y();
	const double mass = planet.get_mass();

	// calculate smoothing length only once if not dependend on radius

	for (unsigned int n_radial = Zero_or_active; n_radial < Max_or_active;
	     ++n_radial) {

	    data[t_data::TORQUE_1D](n_radial) = 0.0;
	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal <= data[t_data::SIGMA].get_max_azimuthal();
		 ++n_azimuthal) {
		// calculate smoothing length if dependend on radius
		// i.e. for thickness smoothing with scale height at cell
		// location

		const double smooth =
		    compute_smoothing(data, n_radial, n_azimuthal);
		const int cell_id = n_azimuthal + n_radial * ns;
		const double xc = cell_center_x[cell_id];
		const double yc = cell_center_y[cell_id];
		const double cellmass =
		    Surf[n_radial] * data[t_data::SIGMA](n_radial, n_azimuthal);
		const double dx = xc - x;
		const double dy = yc - y;
		const double dist_sm_2 =
		    std::pow(dx, 2) + std::pow(dy, 2) + std::pow(smooth, 2);
		const double dist_sm_3 = dist_sm_2 * std::sqrt(dist_sm_2);
		const double inv_dist_sm_3 = 1.0 / dist_sm_3;

		const double Fx =
		    constants::G * cellmass * dx * inv_dist_sm_3 * mass;
		const double Fy =
		    constants::G * cellmass * dy * inv_dist_sm_3 * mass;

		const double Torque = x * Fy - y * Fx;

		data[t_data::TORQUE](n_radial, n_azimuthal) = Torque;
		data[t_data::TORQUE_1D](n_radial) += Torque;
	    }
	}

	const std::string name =  "torque_planet_1D_" + std::to_string(n_planet);
	data[t_data::TORQUE_1D].set_name(name.c_str());

	if (force_update == false) {
	    data[t_data::TORQUE_1D].write1D();
	}
    }
}

void write_1D_info(t_data &data)
{
	if (!CPU_Master) {
		return;
	}

    for (int i = 0; i < t_data::N_POLARGRID_TYPES; ++i) {
	if (data[t_data::t_polargrid_type(i)].get_write_1D()) {

	    int Nr = GlobalNRadial;

	    if (data[t_data::t_polargrid_type(i)].is_vector()) {
			Nr += 1;
		}

		const std::string name = std::string(data[t_data::t_polargrid_type(i)].get_name());
	    const std::string filename_info = outdir + name + "1D.info";

	    std::ofstream info_ofs(filename_info);
		info_ofs.precision(std::numeric_limits<double>::max_digits10);
	    info_ofs << "# version 0.1" << std::endl;
	    info_ofs
		<< "# " << data[t_data::t_polargrid_type(i)].get_name()
		<< " 1d radial, in first line alternating: radii | quantity | minimum quantity | maximum quantity"
		<< std::endl;
	    info_ofs << "# values at time in timestepCoarse.dat" << std::endl;
	    info_ofs << "Nr = " << Nr << std::endl;
	    std::string unit;
	    if (data[t_data::t_polargrid_type(i)].get_unit() != NULL) {
		unit = std::string(data[t_data::t_polargrid_type(i)]
				       .get_unit()
				       ->get_cgs_symbol());
	    } else {
		unit = "1";
	    }
	    info_ofs << "unit = " << unit << std::endl;
	    if (data[t_data::t_polargrid_type(i)].get_unit() != NULL) {
		info_ofs << "code_units_to_cgs_factor = "
			 << data[t_data::t_polargrid_type(i)].get_unit()->get_cgs_factor()
			 << std::endl;
	    } else {
		info_ofs << "code_units_to_cgs_factor = " << 1.0 << std::endl;
	    }
	    info_ofs << "bigendian = " << is_big_endian() << std::endl;
	    info_ofs.close();

	}
    }
}

/**
	Calculates eccentricity, semi major axis and periastron averaged with
   the radial cells
*/
std::vector<double> reduce_disk_quantities(t_data &data, unsigned int timestep,
					   bool force_update,
					   const double quantitiy_radius)
{
    double local_eccentricity = 0.0;
    double disk_eccentricity = 0.0;
	double local_mass = 0.0;
	double global_mass = 0.0;
    // double semi_major_axis = 0.0;
    // double local_semi_major_axis = 0.0;
    double periastron = 0.0;
    double local_periastron = 0.0;

    // calculate eccentricity, semi_major_axis and periastron grid
	quantities::calculate_disk_ecc_peri(data, timestep, force_update);

    // Loop thru all cells excluding GHOSTCELLS & CPUOVERLAP cells (otherwise
    // they would be included twice!)
	const unsigned int Nphi = data[t_data::SIGMA].get_size_azimuthal();
	#pragma omp parallel for collapse(2) reduction(+ : local_eccentricity, local_periastron, local_mass)
	for (unsigned int nr = radial_first_active; nr < radial_active_size; ++nr) {
	for (unsigned int naz = 0; naz < Nphi; ++naz) {
		if (Rmed[nr] <= quantitiy_radius) {
		// eccentricity and semi major axis weighted with cellmass
		const double cell_mass = data[t_data::SIGMA](nr, naz) * Surf[nr];
		local_mass += cell_mass;
		local_eccentricity +=
			data[t_data::ECCENTRICITY](nr, naz) *
			cell_mass;

		// local_semi_major_axis +=
		// data[t_data::SEMI_MAJOR_AXIS](n_radial, n_azimuthal) *
		// local_mass;
		local_periastron +=
			data[t_data::PERIASTRON](nr, naz) *
			cell_mass;
	    }
	}
    }

    // synchronize threads
    MPI_Reduce(&local_eccentricity, &disk_eccentricity, 1, MPI_DOUBLE, MPI_SUM,
	       0, MPI_COMM_WORLD);
    // MPI_Allreduce(&local_semi_major_axis, &semi_major_axis, 1, MPI_DOUBLE,
    // MPI_SUM, MPI_COMM_WORLD);
    MPI_Reduce(&local_periastron, &periastron, 1, MPI_DOUBLE, MPI_SUM, 0,
	       MPI_COMM_WORLD);

	MPI_Reduce(&local_mass, &global_mass, 1, MPI_DOUBLE, MPI_SUM, 0,
		  MPI_COMM_WORLD);

	if (global_mass > 0.0) {
	disk_eccentricity /= global_mass;
	periastron /= global_mass;
    } else {
	disk_eccentricity = 0.0;
	periastron = 0.0;
    }

    std::vector<double> rv = {disk_eccentricity, periastron};

    return rv;
}

void write_lightcurves(t_data &data, unsigned int timestep, bool force_update)
{
    // calculate luminosity
    quantities::calculate_radial_luminosity(data, timestep, force_update);

    // calculate dissipation
    quantities::calculate_radial_dissipation(data, timestep, force_update);

    double *luminosity_values =
	new double[parameters::lightcurves_radii.size()];
    double *dissipation_values =
	new double[parameters::lightcurves_radii.size()];

    for (unsigned int i = 0; i < parameters::lightcurves_radii.size(); ++i) {
	luminosity_values[i] = 0.0;
	dissipation_values[i] = 0.0;
    }

    // all cpu except the lowest one expects data
    if (CPU_Rank > 0) {
	MPI_Recv(luminosity_values, parameters::lightcurves_radii.size(),
		 MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, NULL);
	MPI_Recv(dissipation_values, parameters::lightcurves_radii.size(),
		 MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, NULL);
    }

	/// TODO: openMP parallel
    unsigned int current_lightcurves_bin = 0;
    for (unsigned int n_radial = radial_first_active;
	 n_radial < radial_active_size; ++n_radial) {
	while ((current_lightcurves_bin <
		parameters::lightcurves_radii.size() - 1) &&
	       (parameters::lightcurves_radii[current_lightcurves_bin] <
		Rmed[n_radial])) {
	    current_lightcurves_bin++;
	}
	luminosity_values[current_lightcurves_bin] +=
	    data[t_data::LUMINOSITY_1D](n_radial);
	dissipation_values[current_lightcurves_bin] +=
	    data[t_data::DISSIPATION_1D](n_radial);
    }

    // all cpu except the hightest one sends data
    if (CPU_Rank < CPU_Highest) {
	MPI_Send(luminosity_values, parameters::lightcurves_radii.size(),
		 MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD);
	MPI_Send(dissipation_values, parameters::lightcurves_radii.size(),
		 MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD);
    }

    // the last process can write the data
    if (CPU_Rank == CPU_Highest) {
	// write luminosities
	FILE *fd = 0;
	char *fd_filename;
	static bool fd_created_luminosity = false;

	if (asprintf(&fd_filename, "%s%s", outdir.c_str(), "monitor/luminosity.dat") == -1) {
	    logging::print_master(LOG_ERROR
				  "Not enough memory for string buffer.\n");
	    PersonalExit(1);
	}

	// check if file exists and we restarted
	if ((start_mode::mode == start_mode::mode_restart) &&
	    !(fd_created_luminosity)) {
	    fd = fopen(fd_filename, "r");
	    if (fd) {
		fd_created_luminosity = true;
	    }
	    fclose(fd);
	}
	// open logfile
	if (!fd_created_luminosity) {
	    fd = fopen(fd_filename, "w");
	} else {
	    fd = fopen(fd_filename, "a");
	}
	if (fd == NULL) {
	    logging::print_master(
		LOG_ERROR "Can't write 'luminosity.dat' file. Aborting.\n");
	    PersonalExit(1);
	}

	free(fd_filename);

	if (!fd_created_luminosity) {
	    // print header
	    fprintf(fd, "# PhysicalTime\tluminosities\n");
	    fd_created_luminosity = true;
	}

	fprintf(fd, "%.20e\t", sim::PhysicalTime);

	for (unsigned int i = 1; i < parameters::lightcurves_radii.size();
	     ++i) {
	    fprintf(fd, "%.20e\t", luminosity_values[i]);
	}
	fprintf(fd, "\n");

	// close file
	fclose(fd);

	// write dissipation
	static bool fd_created_dissipation = false;

	if (asprintf(&fd_filename, "%s%s", outdir.c_str(), "monitor/dissipation.dat") ==
	    -1) {
	    logging::print_master(LOG_ERROR
				  "Not enough memory for string buffer.\n");
	    PersonalExit(1);
	}

	// check if file exists and we restarted
	if ((start_mode::mode == start_mode::mode_restart) &&
	    !(fd_created_dissipation)) {
	    fd = fopen(fd_filename, "r");
	    if (fd) {
		fd_created_dissipation = true;
		fclose(fd);
	    }
	}
	// open logfile
	if (!fd_created_dissipation) {
	    fd = fopen(fd_filename, "w");
	} else {
	    fd = fopen(fd_filename, "a");
	}
	if (fd == NULL) {
	    logging::print_master(
		LOG_ERROR "Can't write 'dissipation.dat' file. Aborting.\n");
	    PersonalExit(1);
	}

	free(fd_filename);

	if (!fd_created_dissipation) {
	    // print header
	    fprintf(fd, "# PhysicalTime\tdissipation\n");
	    fd_created_dissipation = true;
	}

	fprintf(fd, "%.20e\t", sim::PhysicalTime);

	for (unsigned int i = 1; i < parameters::lightcurves_radii.size();
	     ++i) {
	    fprintf(fd, "%.20e\t", dissipation_values[i]);
	}
	fprintf(fd, "\n");

	// close file
	fclose(fd);
    }

    delete[] luminosity_values;
    delete[] dissipation_values;
}

/**
Write the corresponding fine grained output number
and the simulation time for each snapshot.
*/
void write_snapshot_time()
{
    FILE *fd = 0;
    static bool fd_created = false;

    if (CPU_Master) {

	const std::string filename = outdir + "snapshots/timeSnapshot.dat";
	auto fd_filename = filename.c_str();

	// check if file exists and we restarted
	if ((start_mode::mode == start_mode::mode_restart) && !(fd_created)) {
	    fd = fopen(fd_filename, "r");
	    if (fd) {
		fd_created = true;
		fclose(fd);
	    }
	}

	// open logfile
	if (!fd_created) {
	    fd = fopen(fd_filename, "w");
	} else {
	    fd = fopen(fd_filename, "a");
	}
	if (fd == NULL) {
	    logging::print_master(
		LOG_ERROR "Can't write 'timeSnapshot.dat' file. Aborting.\n");
	    PersonalExit(1);
	}

	if (!fd_created) {
	    // print header
	    fprintf(fd, "# Time log for course output.\n"
			"#version: 0.1\n"
			"#variable: 0 | time step | 1\n"
			"#variable: 1 | analysis time step | 1\n"
			"#variable: 2 | physical time | ");
	    fprintf(fd, "%s", units::time.get_cgs_factor_symbol().c_str());
	    fprintf(
		fd,
		"\n# One DT is %.18g (code) and %.18g (cgs).\n"
		"# Syntax: coarse output step <tab> fine output step <tab> physical time (code)\n",
		parameters::DT, parameters::DT * units::time.get_cgs_factor());
	    fd_created = true;
	}
    }

    if (CPU_Master) {
	fprintf(fd, "%u\t%u\t%#.16e\n", sim::N_snapshot, sim::N_monitor,
		sim::PhysicalTime);
	fclose(fd);
    }
}


/**
Write the corresponding fine grained output number
and the simulation time for each snapshot.
*/
void write_monitor_time()
{
    FILE *fd = 0;
    static bool fd_created = false;

    if (CPU_Master) {

	const std::string filename = outdir + "monitor/timeMonitor.dat";
	auto fd_filename = filename.c_str();

	// check if file exists and we restarted
	if ((start_mode::mode == start_mode::mode_restart) && !(fd_created)) {
	    fd = fopen(fd_filename, "r");
	    if (fd) {
		fd_created = true;
		fclose(fd);
	    }
	}

	// open logfile
	if (!fd_created) {
	    fd = fopen(fd_filename, "w");
	} else {
	    fd = fopen(fd_filename, "a");
	}
	if (fd == NULL) {
	    logging::print_master(
		LOG_ERROR "Can't write 'timeSnapshot.dat' file. Aborting.\n");
	    PersonalExit(1);
	}

	if (!fd_created) {
	    // print header
	    fprintf(fd, "# Time log for course output.\n"
			"#version: 0.1\n"
			"#variable: 0 | time step | 1\n"
			"#variable: 1 | analysis time step | 1\n"
			"#variable: 2 | physical time | ");
	    fprintf(fd, "%s", units::time.get_cgs_factor_symbol().c_str());
	    fprintf(
		fd,
		"\n# One DT is %.18g (code) and %.18g (cgs).\n"
		"# Syntax: coarse output step <tab> fine output step <tab> physical time (code)\n",
		parameters::DT, parameters::DT * units::time.get_cgs_factor());
	    fd_created = true;
	}
    }

    if (CPU_Master) {
	fprintf(fd, "%u\t%u\t%#.16e\n", sim::N_snapshot, sim::N_monitor,
		sim::PhysicalTime);
	fclose(fd);
    }
}

static std::istream &ignoreline(std::ifstream &in, std::ifstream::pos_type &pos)
{
    pos = in.tellg();
    return in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

static std::string getLastLine(std::ifstream &in)
{
    std::ifstream::pos_type pos = in.tellg();

    std::ifstream::pos_type lastPos;
    while (in >> std::ws && ignoreline(in, lastPos))
	pos = lastPos;

    in.clear();
    in.seekg(pos);

    std::string line;
    std::getline(in, line);
    return line;
}

std::string get_last_snapshot_id()
{
    const std::string filename = outdir + "snapshots/list.txt";
    std::ifstream file(filename);
    std::string last_id = getLastLine(file);
    return last_id;
}

std::int32_t get_latest_output_num(const std::string &snapshot_id)
{
	const std::string path = outdir + "snapshots/" + snapshot_id + "/misc.bin";

    logging::print_master(LOG_INFO "Getting output number of snapshot %s\n",
			  snapshot_id.c_str());

    std::ifstream misc_file(path, std::ios::in | std::ios::binary);

    if (!misc_file.is_open()) {
	logging::print_master(
	    LOG_INFO
	    "Can't read '%s' file in \"get_latest_output_num.\nAttempting to start fresh simulation.\n",
	    path.c_str());
	return -1;
    }

    output::misc_entry entry{0, 0, 0.0, 0.0, 0.0, 0.0, 0};

    misc_file.read((char *)&entry, sizeof(output::misc_entry));

    misc_file.close();

    return entry.timestep;
}

/**
	Checks the conservation of angular momentum over time.
*/
void CheckAngularMomentumConservation(t_data &data)
{
    static double totalStartAngularMomentum, gasStartAngularMomentum,
	planetsStartAngularMomentum;
    static int firstStart = 1;

    unsigned int nPlanet;

    double totalAngularMomentum, gasAngularMomentum,
	planetsAngularMomentum = 0.0, planetAngularMomentum = 0.0;

    FILE *fd;

    double xplanet, yplanet, vxplanet, vyplanet;
    double rpl, thetapl, vazimpl, masspl;

    gasAngularMomentum = quantities::gas_angular_momentum(data, RMAX);

    // computate angular momentum for each planet and sum up
    for (nPlanet = 0;
	 nPlanet < data.get_planetary_system().get_number_of_planets();
	 ++nPlanet) {
	xplanet = data.get_planetary_system().get_planet(nPlanet).get_x();
	yplanet = data.get_planetary_system().get_planet(nPlanet).get_y();
	rpl = sqrt(xplanet * xplanet + yplanet * yplanet);
	thetapl = atan2(yplanet, xplanet);
	vxplanet = data.get_planetary_system().get_planet(nPlanet).get_vx();
	vyplanet = data.get_planetary_system().get_planet(nPlanet).get_vy();
	vazimpl = -vxplanet * std::sin(thetapl) + vyplanet * std::cos(thetapl);
	masspl = data.get_planetary_system().get_planet(nPlanet).get_mass();
	planetAngularMomentum = masspl * rpl * vazimpl;
	planetsAngularMomentum += planetAngularMomentum;
    }

    totalAngularMomentum = gasAngularMomentum + planetsAngularMomentum;
    if (firstStart) {
	firstStart = 0;

	// sim::PhysicalTime < 1e-10 was the "old" condition for saving start values
	if (sim::PhysicalTime > 1e-10) {
	    logging::print_master(
		LOG_INFO
		"CheckAngularMomentumConservation is called for the first time very late: t=%f\n",
		sim::PhysicalTime);
	}

	planetsStartAngularMomentum = planetsAngularMomentum;
	gasStartAngularMomentum = gasAngularMomentum;
	totalStartAngularMomentum = totalAngularMomentum;
	logging::print_master(
	    LOG_INFO "time = %lg, Hp0 = %lg, Hg0 = %lg et Ht0 = %lg\n",
	    sim::PhysicalTime, planetsStartAngularMomentum, gasStartAngularMomentum,
	    totalStartAngularMomentum);
    }

    if (!CPU_Master) {
		return;
	}

	const std::string filename = outdir + "monitor/Momentum.dat";

    // open logfile
    fd = fopen(filename.c_str(), "a");
    if (fd == NULL) {
	logging::print_master(LOG_ERROR
			      "Can't write 'Momentum.dat' file. Aborting.\n");
	PersonalExit(1);
    }


    // computate absolute deviation from start values
    planetsAngularMomentum =
	fabs(planetsAngularMomentum - planetsStartAngularMomentum);
    gasAngularMomentum = fabs(gasAngularMomentum - gasStartAngularMomentum);
    totalAngularMomentum =
	fabs(totalAngularMomentum - totalStartAngularMomentum);

    // print to logfile
    fprintf(fd, "%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\n", sim::PhysicalTime,
	    planetsAngularMomentum, gasAngularMomentum, totalAngularMomentum,
	    totalAngularMomentum / totalStartAngularMomentum);

    // close file
    fclose(fd);
}


void write_ecc_peri_changes(const unsigned int snapshot_number, const unsigned monitor_number)
{
	FILE *fd = 0;
	char *fd_filename;
	static bool fd_created = false;

	if (CPU_Master) {

	if (asprintf(&fd_filename, "%s%s", output::outdir.c_str(), "monitor/eccentricity_change.dat") == -1) {
		logging::print_master(LOG_ERROR
				  "Not enough memory for string buffer.\n");
		PersonalExit(1);
	}
	// check if file exists and we restarted
	if ((start_mode::mode == start_mode::mode_restart) && !(fd_created)) {
		fd = fopen(fd_filename, "r");
		if (fd) {
		fd_created = true;
		fclose(fd);
		}
	}

	// open logfile
	if (!fd_created) {
		fd = fopen(fd_filename, "w");
	} else {
		fd = fopen(fd_filename, "a");
	}
	if (fd == NULL) {
		logging::print_master(
		LOG_ERROR "Can't write 'eccentricity_change.dat' file. Aborting.\n");
		PersonalExit(1);
	}

	free(fd_filename);

	if (!fd_created) {
		// print header
		fprintf(
		fd,
		"# Different torques and eccentricity changes by update steps.\n"
		"# Time unit is %.18g (cgs).\n"
		"# Torque unit is %.18g (cgs).\n"
		"# Eccentricity unit is 1.\n"
		"# Syntax:\n"
		"# 0 : coarse output step\n"
		"# 1 : fine output step\n"
		"# 2 : PhysicalTime\n"

		"# 3 : ecc change from source terms\n"
		"# 4 : ecc change from artificial viscosity\n"
		"# 5 : ecc change from viscosity\n"
		"# 6 : ecc change from transport\n"
		"# 7 : ecc change from damping\n"

		"# 8 : Periastron change from source terms\n"
		"# 9 : Periastron change from artificial viscosity\n"
		"# 10: Periastron change from from viscosity\n"
		"# 11: Periastron change from transport\n"
		"# 12: Periastron change from damping\n",
		units::time.get_cgs_factor(),
		units::torque.get_cgs_factor());
		fd_created = true;
	}
	}

	if (CPU_Master) {
	fprintf(fd, "%u\t%u\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\n",
		snapshot_number,
		monitor_number,
		sim::PhysicalTime,

		delta_ecc_source,
		delta_ecc_art_visc,
		delta_ecc_visc,
		delta_ecc_transport,
		delta_ecc_damp,

		delta_peri_source,
		delta_peri_art_visc,
		delta_peri_visc,
		delta_peri_transport,
		delta_peri_damp);
	fclose(fd);
	}

	delta_ecc_source = 0.0;
	delta_ecc_art_visc = 0.0;
	delta_ecc_visc = 0.0;
	delta_ecc_transport = 0.0;
	delta_ecc_damp = 0.0;

	delta_peri_source = 0.0;
	delta_peri_art_visc = 0.0;
	delta_peri_visc = 0.0;
	delta_peri_transport = 0.0;
	delta_peri_damp = 0.0;
}


} // namespace output
