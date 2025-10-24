import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import diffusion

traj_340 = diffusion.traj_organiser(
    directory="/storage/chem/phuqdw/breaking_bad/diffusion/diff_340k_config1",
    is_dry=False,
)
msd_340, lagtimes_340 = diffusion.msd_calculator(
    topology_file="wet_cooling_ramp.tpr",
    sorted_trajectory_array=traj_340,
    api_residue_name="NAP",
    msd_cache_file_name="wet_340.npz",
)

diffusion.msd_plotter(msd_array=msd_340, lagtime_array=lagtimes_340, is_dry=False)

start_times = [
    200,
    250,
    500,
    600,
    600,
    600,
    600,
    600,
    1000,
    1000,
    1000,
    1500,
    2000,
    2000,
    2000,
    3000,
    3500,
    3500,
    3500,
    4000,
    4000,
    4000,
    4000,
    4000,
    5000,
    5000,
    5000,
    5000,
    5000,
    5000,
    5000,
    5000,
    5000,
    5000,
    5000,
    6000,
    6000,
    6000,
    7000,
    7000,
    7000,
]

end_times = [
    500,
    1000,
    1500,
    2000,
    2000,
    3000,
    3500,
    4000,
    4500,
    5000,
    5500,
    6000,
    6500,
    7000,
    7500,
    8000,
    8500,
    9000,
    9500,
    10000,
    10500,
    11000,
    11500,
    12000,
    12500,
    13000,
    13500,
    14000,
    14500,
    15000,
    15500,
    16000,
    16500,
    17000,
    17500,
    18000,
    19000,
    19500,
    20000,
]

d_340, lagtimes_end_340 = diffusion.diffusion_coefficients(
    msd_array=msd_340,
    lagtime_array=lagtimes_340,
    start_array_ps=start_times,
    end_array_ps=end_times,
)

diffusion.plot_diff_time(
    diff_array=d_340, lagtime_end_array=lagtimes_end_340, is_dry=False
)
