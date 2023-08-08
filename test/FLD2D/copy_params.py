#!/usr/bin/env python3

try:
    import ruamel.yaml
except ImportError:
    raise ImportError(
        "Please install ruamel.yaml with `python3 -m pip install ruamel.yaml`")
yaml = ruamel.yaml.YAML()


def write_setup_parameters(configfile, params):


    with open(configfile, "r") as infile:
        config = yaml.load(infile)

    for key, value in params.items():
        config[key] = value

    with open(configfile, "w") as outfile:
        yaml.dump(config, outfile)


if __name__ == "__main__":

    # get initial time
    with open("test_settings.yml", "r") as infile:
        params = yaml.load(infile)

    t0 = params["t0"]
    tfinal = params["tfinal"]
    steps = params["Nsteps"]
    dt = (tfinal - t0)/steps

    new_params = {
        "DT": dt,
        "FirstDT": dt,
        "RadiativeDiffusionTest2DSteps": params["Nsteps"],
        "RadiativeDiffusionMaxIterations": params["Niter"],
        "RadiativeDiffusionTest2DK": params["K"],
        "RadialSpacing": params["gridspacing"],
        "RadiativeDiffusionTolerance": params["tolerance"]
    }

    write_setup_parameters("setup.yml", new_params)
