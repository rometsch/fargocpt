# FargoCPT Simulation movies

This notebook shows how to make movies from output data using `matplotlib` and `ffmpeg`

First we create a new directory and change to it.


```python
example_name = "500_Movie"
example_dir = f"example_dirs/{example_name}"
import os
repo_root = os.path.abspath(os.path.join(os.getcwd(), "../"))
if not os.path.basename(os.getcwd()) == example_name:
    !mkdir -p $example_dir
    os.chdir(example_dir)
print(f"Current working directory: {os.getcwd()}")
print(f"Repository root directory: {repo_root}")
```

    Current working directory: /home/rometsch/repo/fargocpt/examples/example_dirs/500_Movie
    Repository root directory: /home/rometsch/repo/fargocpt


## Preparing a setup file

We'll take the example setup file from the examples directory and modify it in python.
If you want to create setup files for a parameter study, just copy the code and make your own setup creator script.


```python
configfile = "setup.yml"
!cp $repo_root/examples/config.yml $configfile
```

We'll use the `ruamel.yaml` package to read and write the setup file. This can be set up to preserve comments which is very useful if you want to trace your decisions later on.


```python
import ruamel.yaml
yaml = ruamel.yaml.YAML()
with open(configfile, "r") as infile:
    config = yaml.load(infile)
```


```python
config["nbody"][1]["accretion efficiency"] = "2"
config["MonitorTimestep"] = 0.0314 # monitor scalar files around every half orbit
config["Nmonitor"] = 5 # write a snapshot after every 20 monitor timesteps = every orbit
config["Nsnapshots"] = 100
config["cps"] = 4 # use medium resolution, this is only for show, but it should still look nice
config["Frame"] = "F"
config["OmegaFrame"] = 0
config["nbody"][1]["ramp-up time"] = 10
config["nbody"][1]["accretion efficiency"] = 10
config["DiskFeedback"] = "yes"


with open(configfile, "w") as outfile:
    yaml.dump(config, outfile)
```

## Running the code


```python
from fargocpt import run
np = 2 # Number of mpi processes. Should be equal to the number of numa nodes on your machine, check your cluster docu or run `lscpu` or `./run_fargo --print-numa` if you're on linux.
nt = 8 # Number of threads per mpi process, set it to the number of cores you want to use / number of MPI processes
run(["start", configfile], np=np, nt=nt, exe=repo_root+"/bin/fargocpt_exe", detach=False)
```

    Running command: mpirun -np 2 --report-pid /tmp/tmpu58irt78 --map-by ppr:1:numa --bind-to numa -x OMP_WAIT_POLICY=active -x OMP_PROC_BIND=close -x OMP_PLACES=cores -x OMP_NUM_THREADS=8 /home/rometsch/repo/fargocpt/bin/fargocpt_exe start setup.yml
    fargo process pid 179663
    
    [0] MPI rank #  0 runs as process 179666
    [1] MPI rank #  1 runs as process 179667
    [0] MPI rank #  0 OpenMP thread #  0 of  8 on cpt-kamino
    [1] MPI rank #  1 OpenMP thread #  0 of  8 on cpt-kamino
    [0] MPI rank #  0 OpenMP thread #  2 of  8 on cpt-kamino
    [1] MPI rank #  1 OpenMP thread #  3 of  8 on cpt-kamino
    [1] MPI rank #  1 OpenMP thread #  2 of  8 on cpt-kamino
    [1] MPI rank #  1 OpenMP thread #  1 of  8 on cpt-kamino
    [1] MPI rank #  1 OpenMP thread #  7 of  8 on cpt-kamino
    [1] MPI rank #  1 OpenMP thread #  5 of  8 on cpt-kamino
    [1] MPI rank #  1 OpenMP thread #  4 of  8 on cpt-kamino
    [0] MPI rank #  0 OpenMP thread #  5 of  8 on cpt-kamino
    [0] MPI rank #  0 OpenMP thread #  7 of  8 on cpt-kamino
    [0] MPI rank #  0 OpenMP thread #  1 of  8 on cpt-kamino
    [1] MPI rank #  1 OpenMP thread #  6 of  8 on cpt-kamino
    [0] MPI rank #  0 OpenMP thread #  6 of  8 on cpt-kamino
    [0] MPI rank #  0 OpenMP thread #  4 of  8 on cpt-kamino
    [0] MPI rank #  0 OpenMP thread #  3 of  8 on cpt-kamino
    [0] fargo: This file was compiled on Nov 14 2023, 20:52:20.
    [0] fargo: This version of FARGO used _GNU_SOURCE
    [0] fargo: This version of FARGO used NDEBUG. So no assertion checks!
    [0] Using parameter file setup.yml
    [0] Computing disk quantities within 5.00000e+00 L0 from coordinate center
    [0] BC: Inner composite = reflecting
    [0] BC: Outer composite = reflecting
    [0] BC: Sigma inner = zerogradient
    [0] BC: Sigma outer = zerogradient
    [0] BC: Energy inner = zerogradient
    [0] BC: Energy outer = zerogradient
    [0] BC: Vrad inner = reflecting
    [0] BC: Vrad outer = reflecting
    [0] BC: Vaz inner = keplerian
    [0] BC: Vaz outer = keplerian
    [0] DampingTimeFactor: 1.00000e-01 Outer damping time is computed at radius of 2.50000e+00
    [0] Damping VRadial to reference value at inner boundary.
    [0] Damping VRadial to reference value at outer boundary.
    [0] Damping VAzimuthal to reference value at inner boundary.
    [0] Damping VAzimuthal to reference value at outer boundary.
    [0] Damping SurfaceDensity to reference value at inner boundary.
    [0] Damping SurfaceDensity to reference value at outer boundary.
    [0] Damping Energy to reference value at inner boundary.
    [0] Damping Energy to reference value at outer boundary.
    [0] Radiative diffusion is disabled. Using fixed omega = 1.500000 with a maximum 50000 interations.
    [0] Indirect Term computed as effective Hydro center acceleratrion with shifting the Nbody system to the center.
    [0] Body force on gas computed via potential.
    [0] Using FARGO algorithm for azimuthal advection.
    [0] Using standard forward euler scheme for source terms.
    [0] Cps is set, overwriting Nrad and Naz!
    [0] Grid resolution set using cps = 4.000000
    [0] The grid has (Nrad, Naz) = (148, 504) cells with (4.013071, 4.010705) cps.
    [0] Computing scale height with respect to primary object.
    [0] Using isothermal equation of state. AdiabaticIndex = 1.400.
    [0] Viscosity is of alpha type with alpha = 1.000e-03
    [0] Defaulting to VanLeer flux limiter
    [0] Output information:
    [0]    Output directory: output/out/
    [0]     Number of files: 800
    [0]   Total output size: 0.00 GB
    [0]     Space Available: 38.76 GB
    [0] Initializing 8 RNGs per MPI process.
    [0] Warning : no `radii.dat' file found. Using default.
    [0] The first 1 planets are used to calculate the hydro frame center.
    [0] The mass of the planets used as hydro frame center is 1.000000e+00.
    [0] 2 planet(s) initialized.
    [0] Planet overview:
    [0] 
    [0]  #   | name                    | mass [m0]  | x [l0]     | y [l0]     | vx         | vy         |
    [0] -----+-------------------------+------------+------------+------------+------------+------------+
    [0]    0 | Star                    |          1 |          0 |         -0 |          0 |          0 |
    [0]    1 | Jupiter                 |  0.0009546033 |          1 |          0 |         -0 |   1.000477 |
    [0] 
    [0]  #   | e          | a          | T [t0]     | T [a]      | accreting  | Accretion Type |
    [0] -----+------------+------------+------------+------------+------------+----------------+
    [0]    0 |  6.368246e-17 |          1 |   6.280188 |   0.999548 |          0 |   No Accretion |
    [0]    1 |  6.368246e-17 |          1 |   6.280188 |   0.999548 |         10 |   Kley Accret. |
    [0] 
    [0]  #   | Temp [K]   | R [l0]     | irradiates | rampuptime |
    [0] -----+------------+------------+------------+------------+
    [0]    0 |       5778 |  0.0046505 |        yes |          0 |
    [0]    1 |          0 |  4.6505e-05 |         no |         10 |
    [0] 
    [0] Using Tscharnuter-Winkler (1979) artificial viscosity with C = 1.410000.
    [0] Artificial viscosity is used for dissipation.
    [0] Surface density factor: 2.50663
    [0] Tau factor: 0.5
    [0] Tau min: 0.01
    [0] Kappa factor: 1
    [0] Minimum temperature: 2.81162e-05 K = 3.00000e+00
    [0] Maximum temperature: 9.37206e+94 K = 1.00000e+100
    [0] Heating from viscous dissipation is enabled. Using a total factor of 1.
    [0] Cooling (beta) is disabled and reference temperature is floor. Using beta = 10.
    [0] Cooling (radiative) is enabled. Using a total factor of 1.
    [0] S-curve cooling is disabled. 
    [0] CFL parameter: 0.5
    [0] Opacity uses tables from Lin & Papaloizou, 1985
    [0] Particles are disabled.
    [0] Initializing Sigma(r) = 2.25093e-05 = 200 g cm^-2 * [r/(1 AU)]^(-0.5)
    [0] Total disk is mass is 0.000348846 = 6.9367e+29 g.
    [0] Writing output output/out/snapshots/0, Snapshot Number 0, Time 0.000000.
    [0] Writing output output/out/snapshots/reference, Snapshot Number 0, Time 0.000000.
    [0] Writing output output/out/snapshots/1, Snapshot Number 1, Time 0.157000.
    [0] Writing output output/out/snapshots/2, Snapshot Number 2, Time 0.314000.
    [0] Writing output output/out/snapshots/3, Snapshot Number 3, Time 0.471000.
    [0] Writing output output/out/snapshots/4, Snapshot Number 4, Time 0.628000.
    [0] Writing output output/out/snapshots/5, Snapshot Number 5, Time 0.785000.
    [0] Writing output output/out/snapshots/6, Snapshot Number 6, Time 0.942000.
    [0] Writing output output/out/snapshots/7, Snapshot Number 7, Time 1.099000.
    [0] Writing output output/out/snapshots/8, Snapshot Number 8, Time 1.256000.
    [0] Writing output output/out/snapshots/9, Snapshot Number 9, Time 1.413000.
    [0] Writing output output/out/snapshots/10, Snapshot Number 10, Time 1.570000.
    [0] Writing output output/out/snapshots/11, Snapshot Number 11, Time 1.727000.
    [0] Writing output output/out/snapshots/12, Snapshot Number 12, Time 1.884000.
    [0] Writing output output/out/snapshots/13, Snapshot Number 13, Time 2.041000.
    [0] Writing output output/out/snapshots/14, Snapshot Number 14, Time 2.198000.
    [0] Writing output output/out/snapshots/15, Snapshot Number 15, Time 2.355000.
    [0] Writing output output/out/snapshots/16, Snapshot Number 16, Time 2.512000.
    [0] Writing output output/out/snapshots/17, Snapshot Number 17, Time 2.669000.
    [0] Writing output output/out/snapshots/18, Snapshot Number 18, Time 2.826000.
    [0] Writing output output/out/snapshots/19, Snapshot Number 19, Time 2.983000.
    [0] Writing output output/out/snapshots/20, Snapshot Number 20, Time 3.140000.
    [0] Writing output output/out/snapshots/21, Snapshot Number 21, Time 3.297000.
    [0] Writing output output/out/snapshots/22, Snapshot Number 22, Time 3.454000.
    [0] Writing output output/out/snapshots/23, Snapshot Number 23, Time 3.611000.
    [0] Writing output output/out/snapshots/24, Snapshot Number 24, Time 3.768000.
    [0] Writing output output/out/snapshots/25, Snapshot Number 25, Time 3.925000.
    [0] Writing output output/out/snapshots/26, Snapshot Number 26, Time 4.082000.
    [0] Writing output output/out/snapshots/27, Snapshot Number 27, Time 4.239000.
    [0] Writing output output/out/snapshots/28, Snapshot Number 28, Time 4.396000.
    [0] Writing output output/out/snapshots/29, Snapshot Number 29, Time 4.553000.
    [0] Writing output output/out/snapshots/30, Snapshot Number 30, Time 4.710000.
    [0] Writing output output/out/snapshots/31, Snapshot Number 31, Time 4.867000.
    [0] Writing output output/out/snapshots/32, Snapshot Number 32, Time 5.024000.
    [0] Writing output output/out/snapshots/33, Snapshot Number 33, Time 5.181000.
    [0] Writing output output/out/snapshots/34, Snapshot Number 34, Time 5.338000.
    [0] Writing output output/out/snapshots/35, Snapshot Number 35, Time 5.495000.
    [0] Writing output output/out/snapshots/36, Snapshot Number 36, Time 5.652000.
    [0] Writing output output/out/snapshots/37, Snapshot Number 37, Time 5.809000.
    [0] Writing output output/out/snapshots/38, Snapshot Number 38, Time 5.966000.
    [0] Writing output output/out/snapshots/39, Snapshot Number 39, Time 6.123000.
    [0] Writing output output/out/snapshots/40, Snapshot Number 40, Time 6.280000.
    [0] Writing output output/out/snapshots/41, Snapshot Number 41, Time 6.437000.
    [0] Writing output output/out/snapshots/42, Snapshot Number 42, Time 6.594000.
    [0] Writing output output/out/snapshots/43, Snapshot Number 43, Time 6.751000.
    [0] Writing output output/out/snapshots/44, Snapshot Number 44, Time 6.908000.
    [0] Writing output output/out/snapshots/45, Snapshot Number 45, Time 7.065000.
    [0] Writing output output/out/snapshots/46, Snapshot Number 46, Time 7.222000.
    [0] Writing output output/out/snapshots/47, Snapshot Number 47, Time 7.379000.
    [0] Writing output output/out/snapshots/48, Snapshot Number 48, Time 7.536000.
    [0] Writing output output/out/snapshots/49, Snapshot Number 49, Time 7.693000.
    [0] Writing output output/out/snapshots/50, Snapshot Number 50, Time 7.850000.
    [0] Writing output output/out/snapshots/51, Snapshot Number 51, Time 8.007000.
    [0] Writing output output/out/snapshots/52, Snapshot Number 52, Time 8.164000.
    [0] Writing output output/out/snapshots/53, Snapshot Number 53, Time 8.321000.
    [0] Writing output output/out/snapshots/54, Snapshot Number 54, Time 8.478000.
    [0] Writing output output/out/snapshots/55, Snapshot Number 55, Time 8.635000.
    [0] Writing output output/out/snapshots/56, Snapshot Number 56, Time 8.792000.
    [0] Writing output output/out/snapshots/57, Snapshot Number 57, Time 8.949000.
    [0] Writing output output/out/snapshots/58, Snapshot Number 58, Time 9.106000.
    [0] Writing output output/out/snapshots/59, Snapshot Number 59, Time 9.263000.
    [0] Writing output output/out/snapshots/60, Snapshot Number 60, Time 9.420000.
    [0] Writing output output/out/snapshots/61, Snapshot Number 61, Time 9.577000.
    [0] Writing output output/out/snapshots/62, Snapshot Number 62, Time 9.734000.
    [0] Writing output output/out/snapshots/63, Snapshot Number 63, Time 9.891000.
    [0] Writing output output/out/snapshots/64, Snapshot Number 64, Time 10.048000.
    [0] Writing output output/out/snapshots/65, Snapshot Number 65, Time 10.205000.
    [0] Writing output output/out/snapshots/66, Snapshot Number 66, Time 10.362000.
    [0] Writing output output/out/snapshots/67, Snapshot Number 67, Time 10.519000.
    [0] Writing output output/out/snapshots/68, Snapshot Number 68, Time 10.676000.
    [0] Writing output output/out/snapshots/69, Snapshot Number 69, Time 10.833000.
    [0] Writing output output/out/snapshots/70, Snapshot Number 70, Time 10.990000.
    [0] Writing output output/out/snapshots/71, Snapshot Number 71, Time 11.147000.
    [0] Writing output output/out/snapshots/72, Snapshot Number 72, Time 11.304000.
    [0] Writing output output/out/snapshots/73, Snapshot Number 73, Time 11.461000.
    [0] Writing output output/out/snapshots/74, Snapshot Number 74, Time 11.618000.
    [0] Writing output output/out/snapshots/75, Snapshot Number 75, Time 11.775000.
    [0] Writing output output/out/snapshots/76, Snapshot Number 76, Time 11.932000.
    [0] Writing output output/out/snapshots/77, Snapshot Number 77, Time 12.089000.
    [0] Writing output output/out/snapshots/78, Snapshot Number 78, Time 12.246000.
    [0] Writing output output/out/snapshots/79, Snapshot Number 79, Time 12.403000.
    [0] Writing output output/out/snapshots/80, Snapshot Number 80, Time 12.560000.
    [0] Writing output output/out/snapshots/81, Snapshot Number 81, Time 12.717000.
    [0] Writing output output/out/snapshots/82, Snapshot Number 82, Time 12.874000.
    [0] Writing output output/out/snapshots/83, Snapshot Number 83, Time 13.031000.
    [0] Writing output output/out/snapshots/84, Snapshot Number 84, Time 13.188000.
    [0] Writing output output/out/snapshots/85, Snapshot Number 85, Time 13.345000.
    [0] Writing output output/out/snapshots/86, Snapshot Number 86, Time 13.502000.
    [0] Writing output output/out/snapshots/87, Snapshot Number 87, Time 13.659000.
    [0] Writing output output/out/snapshots/88, Snapshot Number 88, Time 13.816000.
    [0] Writing output output/out/snapshots/89, Snapshot Number 89, Time 13.973000.
    [0] Writing output output/out/snapshots/90, Snapshot Number 90, Time 14.130000.
    [0] Writing output output/out/snapshots/91, Snapshot Number 91, Time 14.287000.
    [0] Writing output output/out/snapshots/92, Snapshot Number 92, Time 14.444000.
    [0] Writing output output/out/snapshots/93, Snapshot Number 93, Time 14.601000.
    [0] Writing output output/out/snapshots/94, Snapshot Number 94, Time 14.758000.
    [0] Writing output output/out/snapshots/95, Snapshot Number 95, Time 14.915000.
    [0] Writing output output/out/snapshots/96, Snapshot Number 96, Time 15.072000.
    [0] Writing output output/out/snapshots/97, Snapshot Number 97, Time 15.229000.
    [0] Writing output output/out/snapshots/98, Snapshot Number 98, Time 15.386000.
    [0] Writing output output/out/snapshots/99, Snapshot Number 99, Time 15.543000.
    [0] Writing output output/out/snapshots/100, Snapshot Number 100, Time 15.700000.
    [0] -- Final: Total Hydrosteps 500, Time 15.70, Walltime 2.34 seconds, Time per Step: 4.69 milliseconds





    0



## Plotting

First, we define a function to plot the top down view on the disk.

Then, we plot all images seperately and store them in a directory.

Finally, we combine the single images into a movie by leveraging ffmpeg.


```python
import numpy as np
import matplotlib.colors as mplcolors
import matplotlib.pyplot as plt

def plot_field(loader, name, N, ax=None, dataunit=None, vmin=None, vmax=None, cmap="viridis", title=None):
    R, PHI, vals = loader.gas.vars2D.get(name, N, grid_for_plot=True)
    if dataunit is None:
        dataunit = vals.unit
    C = vals.to_value(dataunit)

    X = R*np.cos(PHI)
    Y = R*np.sin(PHI)

    if ax is None:
        fig, ax = plt.subplots(dpi=150)
    else:
        fig = ax.get_figure()

    norm = mplcolors.Normalize(vmin=vmin, vmax=vmax)

    # Hacky way to support arrays that are defined on the radial interfaces
    if C.shape[0] == X.shape[0]:
        C = C[:-1,:]

    pcm = ax.pcolormesh(X,Y,C, norm=norm, cmap=cmap)
    ax.set_aspect("equal")

    t = loader.snapshot_time[N].to_value("kyr")
    if title is None:
        title = ""
    else:
        title += "\n"
    title += f" t={t:.2e}kyr, N={N}"
    ax.set_title(title)

    cbar = fig.colorbar(pcm, ax=ax)
    cbar.set_label(f"{name} [{dataunit}]")
    
    return fig
```


```python
from fargocpt import Loader
l = Loader("output/out/")
```


```python
folder = "imgs"
!mkdir -p $folder
```


```python
# have a progress bar
!pip install tqdm
```

    Requirement already satisfied: tqdm in /home/rometsch/repo/fargocpt/python-venv/lib/python3.10/site-packages (4.66.1)


We plot all snapshots, this might take a while.


```python
from tqdm import tqdm
for n in tqdm(l.snapshots):
    fig = plot_field(l, "Sigma", n, dataunit="g/cm2", vmin=0, vmax=600, cmap="magma", title="Sigma");
    fig.savefig(f"imgs/Sigma{n}.jpg", dpi=300)
    plt.close(fig)
```

      0%|          | 0/101 [00:00<?, ?it/s]100%|██████████| 101/101 [00:15<00:00,  6.61it/s]


Finally, we make the movie.

To change the file size, you can adjust framerate and the quality parameter.
Try a few settings with a small number of frames to find your sweet stop.


```python
framerate = 60
quality = 20
!ffmpeg -y -framerate $framerate -i $folder/Sigma%d.jpg -c:v libx264 -crf $quality output.mp4
```

    ffmpeg version 4.4.2-0ubuntu0.22.04.1 Copyright (c) 2000-2021 the FFmpeg developers
      built with gcc 11 (Ubuntu 11.2.0-19ubuntu1)
      configuration: --prefix=/usr --extra-version=0ubuntu0.22.04.1 --toolchain=hardened --libdir=/usr/lib/x86_64-linux-gnu --incdir=/usr/include/x86_64-linux-gnu --arch=amd64 --enable-gpl --disable-stripping --enable-gnutls --enable-ladspa --enable-libaom --enable-libass --enable-libbluray --enable-libbs2b --enable-libcaca --enable-libcdio --enable-libcodec2 --enable-libdav1d --enable-libflite --enable-libfontconfig --enable-libfreetype --enable-libfribidi --enable-libgme --enable-libgsm --enable-libjack --enable-libmp3lame --enable-libmysofa --enable-libopenjpeg --enable-libopenmpt --enable-libopus --enable-libpulse --enable-librabbitmq --enable-librubberband --enable-libshine --enable-libsnappy --enable-libsoxr --enable-libspeex --enable-libsrt --enable-libssh --enable-libtheora --enable-libtwolame --enable-libvidstab --enable-libvorbis --enable-libvpx --enable-libwebp --enable-libx265 --enable-libxml2 --enable-libxvid --enable-libzimg --enable-libzmq --enable-libzvbi --enable-lv2 --enable-omx --enable-openal --enable-opencl --enable-opengl --enable-sdl2 --enable-pocketsphinx --enable-librsvg --enable-libmfx --enable-libdc1394 --enable-libdrm --enable-libiec61883 --enable-chromaprint --enable-frei0r --enable-libx264 --enable-shared
      libavutil      56. 70.100 / 56. 70.100
      libavcodec     58.134.100 / 58.134.100
      libavformat    58. 76.100 / 58. 76.100
      libavdevice    58. 13.100 / 58. 13.100
      libavfilter     7.110.100 /  7.110.100
      libswscale      5.  9.100 /  5.  9.100
      libswresample   3.  9.100 /  3.  9.100
      libpostproc    55.  9.100 / 55.  9.100
    Input #0, image2, from 'imgs/Sigma%d.jpg':
      Duration: 00:00:01.68, start: 0.000000, bitrate: N/A
      Stream #0:0: Video: mjpeg (Baseline), yuvj420p(pc, bt470bg/unknown/unknown), 1920x1440 [SAR 300:300 DAR 4:3], 60 fps, 60 tbr, 60 tbn, 60 tbc
    Stream mapping:
      Stream #0:0 -> #0:0 (mjpeg (native) -> h264 (libx264))
    Press [q] to stop, [?] for help
    [1;36m[libx264 @ 0x557fecd5efc0] [0musing SAR=1/1
    [1;36m[libx264 @ 0x557fecd5efc0] [0musing cpu capabilities: MMX2 SSE2Fast SSSE3 SSE4.2 AVX FMA3 BMI2 AVX2
    [1;36m[libx264 @ 0x557fecd5efc0] [0mprofile High, level 5.1, 4:2:0, 8-bit
    [1;36m[libx264 @ 0x557fecd5efc0] [0m264 - core 163 r3060 5db6aa6 - H.264/MPEG-4 AVC codec - Copyleft 2003-2021 - http://www.videolan.org/x264.html - options: cabac=1 ref=3 deblock=1:0:0 analyse=0x3:0x113 me=hex subme=7 psy=1 psy_rd=1.00:0.00 mixed_ref=1 me_range=16 chroma_me=1 trellis=1 8x8dct=1 cqm=0 deadzone=21,11 fast_pskip=1 chroma_qp_offset=-2 threads=24 lookahead_threads=4 sliced_threads=0 nr=0 decimate=1 interlaced=0 bluray_compat=0 constrained_intra=0 bframes=3 b_pyramid=2 b_adapt=1 b_bias=0 direct=1 weightb=1 open_gop=0 weightp=2 keyint=250 keyint_min=25 scenecut=40 intra_refresh=0 rc_lookahead=40 rc=crf mbtree=1 crf=20.0 qcomp=0.60 qpmin=0 qpmax=69 qpstep=4 ip_ratio=1.40 aq=1:1.00
    Output #0, mp4, to 'output.mp4':
      Metadata:
        encoder         : Lavf58.76.100
      Stream #0:0: Video: h264 (avc1 / 0x31637661), yuvj420p(pc, bt470bg/unknown/unknown, progressive), 1920x1440 [SAR 300:300 DAR 4:3], q=2-31, 60 fps, 15360 tbn
        Metadata:
          encoder         : Lavc58.134.100 libx264
        Side data:
          cpb: bitrate max/min/avg: 0/0/0 buffer size: 0 vbv_delay: N/A
    frame=  101 fps=0.0 q=-1.0 Lsize=     436kB time=00:00:01.63 bitrate=2187.8kbits/s speed=1.75x    
    video:434kB audio:0kB subtitle:0kB other streams:0kB global headers:0kB muxing overhead: 0.470529%
    [1;36m[libx264 @ 0x557fecd5efc0] [0mframe I:1     Avg QP:17.13  size: 47798
    [1;36m[libx264 @ 0x557fecd5efc0] [0mframe P:25    Avg QP:18.13  size:  7060
    [1;36m[libx264 @ 0x557fecd5efc0] [0mframe B:75    Avg QP:18.98  size:  2928
    [1;36m[libx264 @ 0x557fecd5efc0] [0mconsecutive B-frames:  1.0%  0.0%  0.0% 99.0%
    [1;36m[libx264 @ 0x557fecd5efc0] [0mmb I  I16..4: 20.8% 72.9%  6.4%
    [1;36m[libx264 @ 0x557fecd5efc0] [0mmb P  I16..4:  0.7% 12.8%  0.1%  P16..4:  4.6%  1.0%  0.6%  0.0%  0.0%    skip:80.3%
    [1;36m[libx264 @ 0x557fecd5efc0] [0mmb B  I16..4:  0.4%  1.3%  0.0%  B16..8:  8.4%  1.2%  0.1%  direct: 1.1%  skip:87.6%  L0:50.6% L1:47.7% BI: 1.7%
    [1;36m[libx264 @ 0x557fecd5efc0] [0m8x8 transform intra:86.4% inter:93.8%
    [1;36m[libx264 @ 0x557fecd5efc0] [0mcoded y,uvDC,uvAC intra: 20.9% 29.4% 2.5% inter: 0.7% 2.8% 0.0%
    [1;36m[libx264 @ 0x557fecd5efc0] [0mi16 v,h,dc,p: 48% 30% 12% 10%
    [1;36m[libx264 @ 0x557fecd5efc0] [0mi8 v,h,dc,ddl,ddr,vr,hd,vl,hu: 41% 23% 33%  1%  0%  0%  0%  0%  0%
    [1;36m[libx264 @ 0x557fecd5efc0] [0mi4 v,h,dc,ddl,ddr,vr,hd,vl,hu: 40% 23% 15%  4%  3%  4%  3%  4%  3%
    [1;36m[libx264 @ 0x557fecd5efc0] [0mi8c dc,h,v,p: 46% 25% 26%  3%
    [1;36m[libx264 @ 0x557fecd5efc0] [0mWeighted P-Frames: Y:0.0% UV:0.0%
    [1;36m[libx264 @ 0x557fecd5efc0] [0mref P L0: 60.4%  3.6% 20.1% 15.9%
    [1;36m[libx264 @ 0x557fecd5efc0] [0mref B L0: 82.4% 14.6%  3.0%
    [1;36m[libx264 @ 0x557fecd5efc0] [0mref B L1: 96.1%  3.9%
    [1;36m[libx264 @ 0x557fecd5efc0] [0mkb/s:2109.70


## The result

If you have the IPython package installed, the following should show you the video.


```python
from IPython.display import Video

# Display the video
Video(os.getcwd() + "/output.mp4")
```




<video src="500_Movie_files/output.mp4" controls  >
      Your browser does not support the <code>video</code> element.
    </video>




```python

```
