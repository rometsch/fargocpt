import sys
sys.path.append('../../bin')
from fargocpt import run_fargo



N_procs = 2
N_OMP_threads = 7
fargo_args = ["auto"]

def run_simulation(params = {}, planet_params = None, configfile = "setup.yml"):
    params = dict(params)

    import ruamel.yaml
    yaml = ruamel.yaml.YAML()
    with open(configfile, "r") as infile:
        config = yaml.load(infile)
        
    for key, value in params.items():
        config[key] = value

    if planet_params is not None:
        for planet, plp in planet_params.items():
            for key, value in plp.items():
                config["nbody"][planet][key] = value

    with open(configfile, "w") as outfile:
        yaml.dump(config, outfile)

    run_fargo(N_procs, N_OMP_threads, fargo_args + [configfile])


for e in [0.1, 0.2]:
    for idm in ["0", "1"]:
        outputdir = f"output/a0.6_idm{idm}_e{e}"
        params = {
            "IndirectTermMode" : idm,
            "OutputDir": outputdir
        }
        planet_params = {
            1 : {
                "eccentricity" : e,
                "semi-major axis" : 0.6
            }
        }
        run_simulation(params, planet_params=planet_params
        )


