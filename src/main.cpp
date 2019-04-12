#include <mpi.h>
#include <string.h>
#include <fstream>

#include "units.h"
#include "constants.h"
#include "time.h"
#include "global.h"
#include "init.h"
#include "SourceEuler.h"
#include "SideEuler.h"
#include "commbound.h"
#include "LowTasks.h"
#include "boundary_conditions.h"
#include "Pframeforce.h"
#include "Stockholm.h"
#include "Force.h"
#include "output.h"
#include "split.h"
#include "selfgravity.h"
#include "fpe.h"
#include "Interpret.h"
#include "parameters.h"
#include "logging.h"
#include "options.h"
#include "data.h"
#include "quantities.h"
#include "viscosity.h"
#include "Theo.h"
#include "particles.h"
#include "util.h"

int TimeToWrite;
int dimfxy = 11, Restart = 0;
static int InnerOutputCounter=0, StillWriteOneOutput;
extern int Corotating;
extern int SelfGravity, SGZeroMode;
extern bool Adiabatic;

unsigned int nTimeStep;

// TODO: NRad darf nicht größer sein als MAX1D
// copy operator for t_polargrid
// Handle SIGTERM
// write all polargrids on error

int main(int argc, char* argv[])
{
	t_data data;

	N_iter= 0;

	resize_radialarrays(MAX1D);

	unsigned int timeStepStart = 0;

	int CPU_NameLength;
	Force *force;
	char CPU_Name[MPI_MAX_PROCESSOR_NAME+1];

	// initialize MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &CPU_Rank);
	MPI_Comm_size(MPI_COMM_WORLD, &CPU_Number);
	selfgravity::mpi_init();

	// are we master CPU?
	CPU_Master = (CPU_Rank == 0 ? 1 : 0);

	// print some information about program
	logging::print_master(LOG_INFO "fargo: This file was compiled on %s, %s.\n",__DATE__,__TIME__);
	#ifdef GIT_COMMIT
	logging::print_master(LOG_INFO "fargo: Last git commit: %s\n", GIT_COMMIT);
	#endif
	#ifdef GIT_CHANGED
	logging::print_master(LOG_INFO "fargo: Files changed since git commit: %s\n", GIT_CHANGED);
	#endif
	#ifdef _GNU_SOURCE
	logging::print_master(LOG_INFO "fargo: This version of FARGO used _GNU_SOURCE\n");
	#endif
	#ifdef NDEBUG
	logging::print_master(LOG_INFO "fargo: This version of FARGO used NDEBUG. So no assertion checks!\n");
	#else
	logging::print_master(LOG_INFO "fargo: This version of FARGO does assertion checks! Compile with NDEBUG to speed up!\n");
	#endif

	// print information on which processor we're running
	MPI_Get_processor_name(CPU_Name, &CPU_NameLength);
	CPU_Name[CPU_NameLength] = '\0';
	logging::print(LOG_INFO "fargo: running on %s\n",CPU_Name);

	// control behavoir for floating point exceptions trapping (default is not to do anything)
	setfpe();

	// handle command line parameters
    options::parse(argc,argv);

    ReadVariables(options::parameter_file, data, argc, argv);
    // check if there is enough free space for all outputs (check before any output are files created)
    output::check_free_space(data);

	parameters::summarize_parameters();
	parameters::write_grid_data_to_file();

	units::print_code_units();
	units::write_code_unit_file();
	constants::print_constants();

	SplitDomain();

	TellEverything();

	if (options::memory_usage) {
        data.print_memory_usage(NRadial, NAzimuthal);
		PersonalExit(0);
	}

	if (options::disable)
		PersonalExit(0);

	data.set_size(GlobalNRadial, NAzimuthal, NRadial, NAzimuthal);

	init_radialarrays();
	force = AllocateForce(dimfxy);

	// Here planets are initialized feeling star potential but they do not feel disk potential
	data.get_planetary_system().read_from_file(PLANETCONFIG);

	logging::print_master(LOG_INFO "planets initialised.\n");

	if ((data.get_planetary_system().get_number_of_planets() == 0) && (Corotating == YES)) {
		logging::print_master(LOG_ERROR "Error: Corotating frame is not possible with 0 planets.\n");
		PersonalExit(1);
	}

	init_physics(data);

	if (parameters::integrate_particles)
		particles::init();

	// save starting values (needed for damping)
	copy_polargrid(data[t_data::V_RADIAL0], data[t_data::V_RADIAL]);
	copy_polargrid(data[t_data::V_AZIMUTHAL0], data[t_data::V_AZIMUTHAL]);
	copy_polargrid(data[t_data::DENSITY0], data[t_data::DENSITY]);
	copy_polargrid(data[t_data::ENERGY0], data[t_data::ENERGY]);

	// Initial Density is used to compute the circumplanetary mass with initial density field
	mdcp0 = CircumPlanetaryMass(data);

	if (options::restart) {
		timeStepStart = options::restart_from * NINTERM;
		logging::print_master(LOG_INFO "Restarting planetary system...\n");
		// restart planetary system
		data.get_planetary_system().restart(options::restart_from);

		logging::print_master(LOG_INFO "Reading misc data...\n");
		if (CPU_Master) {
			PhysicalTime = output::get_misc(options::restart_from, "physical time");
			OmegaFrame = output::get_misc(options::restart_from, "omega frame");
			FrameAngle = output::get_misc(options::restart_from, "frame angle");
		}

		// load grids at t = 0
		logging::print_master(LOG_INFO "Loading polargrinds at t = 0...\n");
		data[t_data::DENSITY].read2D((unsigned int)0);
		data[t_data::V_RADIAL].read2D((unsigned int)0);
		data[t_data::V_AZIMUTHAL].read2D((unsigned int)0);
		if (Adiabatic)
			data[t_data::ENERGY].read2D((unsigned int)0);

		// save starting values (needed for damping)
		copy_polargrid(data[t_data::V_RADIAL0], data[t_data::V_RADIAL]);
		copy_polargrid(data[t_data::V_AZIMUTHAL0], data[t_data::V_AZIMUTHAL]);
		copy_polargrid(data[t_data::DENSITY0], data[t_data::DENSITY]);
		copy_polargrid(data[t_data::ENERGY0], data[t_data::ENERGY]);

		// recalculate SigmaMed/EnergyMed
		RefillSigma(&data[t_data::DENSITY]);
		if (Adiabatic)
			RefillEnergy(&data[t_data::ENERGY]);

		// load grids at t = restart_from
		logging::print_master(LOG_INFO "Loading polargrinds at t = %u...\n",options::restart_from);
		data[t_data::DENSITY].read2D(options::restart_from);
		data[t_data::V_RADIAL].read2D(options::restart_from);
		data[t_data::V_AZIMUTHAL].read2D(options::restart_from);
		if (Adiabatic)
			data[t_data::ENERGY].read2D(options::restart_from);

		CommunicateBoundaries(&data[t_data::DENSITY], &data[t_data::V_RADIAL], &data[t_data::V_AZIMUTHAL], &data[t_data::ENERGY]);

		// recalculate everyting
		if (Adiabatic) {
			compute_sound_speed(data, true);
			compute_aspect_ratio(data, true);
		}
		viscosity::update_viscosity(data);
		compute_pressure(data, true);
		compute_temperature(data, true);
	} else {
		// create planet files
		data.get_planetary_system().create_planet_files();
		// create mass flow file
		if (parameters::write_massflow) {
			if (CPU_Master) {
				const std::string filename_info = std::string(OUTPUTDIR) + "/gasMassFlow1D.info";
				std::ofstream info_ofs(filename_info);
                info_ofs << "# Mass flow 1d radial, first line radii, from second line on, values at time in Quantities.dat, Nr = " << GlobalNRadial+1 << " , unit = " << data[t_data::MASSFLOW_1D].get_unit()->get_cgs_symbol() << " , bigendian = " << is_big_endian() << std::endl;
				const std::string filename = std::string(OUTPUTDIR) + "/gasMassFlow1D.dat";
				std::ofstream ofs(filename,  std::ios::binary);
				ofs.write( (char*)Radii.array, sizeof(*Radii.array)*(GlobalNRadial+1) );
			}
		}

        if (data[t_data::DENSITY].get_write_1D()) {
            if (CPU_Master) {
                char *tmp;

                if (asprintf(&tmp, "%s/gas%s1D.info",OUTPUTDIR,data[t_data::DENSITY].get_name())<0) {
                    die("Not enough memory!");
                }
                const std::string filename_info = std::string(tmp);
                std::ofstream info_ofs(filename_info);
                info_ofs << "# Surface density 1d radial, in first line alternating: | radii | surface density | minimum surface density | maximum surface density | values at time in timestepCoarse.dat, Nr = " << GlobalNRadial << " , unit = " << data[t_data::DENSITY].get_unit()->get_cgs_symbol()  << " , bigendian = " << is_big_endian() << std::endl;
                info_ofs.close();
            }
        }
	}
	PhysicalTimeInitial = PhysicalTime;

	logging::start_timer();
  
	for (nTimeStep = timeStepStart; nTimeStep <= NTOT; ++nTimeStep) {
		InnerOutputCounter++;

		if (InnerOutputCounter == 1) {
			InnerOutputCounter = 0;
			data.get_planetary_system().write_planets(TimeStep, true);
			//WriteBigPlanetSystemFile(sys, TimeStep);
            UpdateLog(data, force, TimeStep, PhysicalTime);
			if (Stockholm)
                UpdateLogStockholm(data, PhysicalTime);
		}

		if (NINTERM * (TimeStep = (nTimeStep / NINTERM)) == nTimeStep) {
			// Outputs are done here
			TimeToWrite = YES;

			// write polar grids
			output::write_grids(data, TimeStep, N_iter, PhysicalTime);
			// write planet data
			data.get_planetary_system().write_planets(TimeStep, false);
			// write quantities like energy, mass, ...
            output::write_quantities(data, TimeStep, nTimeStep, false);
			// write misc stuff (important for resuming)
			output::write_misc(TimeStep);
			// write time info for coarse output
			output::write_coarse_time(TimeStep, nTimeStep);
			// write disk quantities like eccentricity, ...
			if (parameters::write_torques)
				output::write_torques(data, TimeStep, false);
			if (parameters::write_lightcurves)
				output::write_lightcurves(data, TimeStep, false);
			// write particles
			if (parameters::integrate_particles)
				particles::write(TimeStep);
			if ((OnlyInit) || ((GotoNextOutput) && (!StillWriteOneOutput))) {
				PersonalExit(0);
			}
			StillWriteOneOutput--;
		} else {
			// write disk quantities like eccentricity, ...
			if (parameters::write_at_every_timestep)
                output::write_quantities(data, TimeStep, nTimeStep, true);
			if (parameters::write_torques)
				output::write_torques(data, TimeStep, true);
			if (parameters::write_lightcurves && parameters::write_at_every_timestep)
				output::write_lightcurves(data, TimeStep, true);
			TimeToWrite = NO;
		}
		AlgoGas(nTimeStep, force, data);
		SolveOrbits(data);
		if (parameters::write_massflow) {
			const std::string filename = std::string(OUTPUTDIR) + "/gasMassFlow1D.dat";
			data[t_data::MASSFLOW_1D].write(filename, TimeStep, data, true, true);
		}
	}

	logging::print_runtime_final();

	// free up everything
	delete [] options::parameter_file;
	FreeEuler();
	FreeForce(force);

	selfgravity::mpi_finalize();
	MPI_Finalize();

	return 0;
}
