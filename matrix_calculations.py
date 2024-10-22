# matrix_calculations.py ---
#
# Filename: matrix_calculations.py
# Description: This code make the calculation of the connectivity matrix for the GBR
# Author: Javier Porobic
# Maintainer: Javier Porobic concat "javier.porobicgarate" "@" "csiro" ".au")
# Created: Thu Feb 16 13:05:47 2023 (+1100)
# Version: 0.0.1
# Package-Requires: ()
# Last-Updated:
#           By:Javier Porobic
#     Update #: 25
# URL:
# Doc URL:
# Keywords:
# Compatibility:
# This piece of code has been created to work on standard versions of python with
# publicly available libraries. The only configuration that has been hardwired is the
# location of the files and the commands related to the work on the high-performance
# computer of CSIRO.
#

# Commentary:
# This code utilizes particle dispersal tracks generated by OceanParcels to evaluate
# the level of connectivity among the Great Barrier Reef (GBR) reefs. The
# calculations for total connectivity take into account decay and competence, and are
# based on the functions established by Moneghetti et al. in 2019. These functions
# allow for a comprehensive assessment of the dispersal dynamics of particles,
# providing insights into the connectivity patterns among the GBR reefs. The code
# serves as a valuable tool for understanding the interconnectivity of reef systems
# and can contribute to informed decision-making in reef management and conservation
# efforts.
#

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
import xarray as xr
from scipy.integrate import quad
import geopandas as gpd
import numpy as np
import sys
from shapely.geometry import Point
import shapely
import math
from tqdm import tqdm
import pandas as pd
from numba import jit, njit
import numba #Numba translates Python functions to optimized machine code at runtime.
import os
from joblib import Parallel, delayed, parallel_backend
import time

## ~~~~~~~~~~~~~~~~~~ ##
## ~     Functions  ~ ##
## ~~~~~~~~~~~~~~~~~~ ##
def bathtub_curve(lmbda, v, sigma) :
    """
    Transforms the bathtub equation into a function that can be integrated
    using the trapezoidal function.

    Parameters
    ----------
    lmbda : float
        The scale parameter of the Weibull distribution used in the bathtub equation.
        Must be a positive integer or float.
    v : float
        The shape parameter of the Weibull distribution used in the bathtub equation.
        Must be a positive integer or float.
    sigma : float
        A constant that determines the shape of the bathtub curve. Must be a positive
        integer or float. It can be zero, becoming an exponential function.

    Returns
    -------
    function
        A lambda function that represents the bathtub curve equation. The function takes
        a single argument, t, which represents the age of the coral. The function returns
        the value of the bathtub curve equation at that age.
    """

    u = lambda t:(lmbda * v * pow((lmbda * t), v - 1)) / (1 - sigma * pow((lmbda * t), v))
    return(u)

def piecewise_decay(ages, Tcp, lmbda1, lmbda2, v1, v2, sigma1, sigma2):
    """
    Calculates the probability of survival of corals larvaes at different ages, using the piecewise
    Weibull-Weibull survival model described by Moneghetti et al. (2019).

    Parameters
    ----------
    ages : list
        A list of ages of the corals. Each age must be a positive integer or float.
    Tcp : float
        The age (inflection point) at which the corals transition from a Weibull survival curve to another Weibull
        survival curve. Must be greater than 0 and less than the maximum age in `ages`.
    lmbda1 : float
        The scale parameter of the Weibull survival curve in the first phase. Must be a positive
        integer or float.
    lmbda2 : float
        The scale parameter of the Weibull survival curve in the second phase. Must be a positive
        integer or float.
    v1 : float
        The shape parameter of the Weibull survival curve in the first phase. Must be a positive
        integer or float.
    v2 : float
        The shape parameter of the Weibull survival curve in the second phase. Must be a positive
        integer or float.
    sigma1 : float
        The standard deviation of the Gaussian noise added to the survival curve in the first phase.
        Must be a positive integer or float. It can be zero, becoming an exponential function.
    sigma2 : float
        The standard deviation of the Gaussian noise added to the survival curve in the second phase.
        Must be a positive integer or float. It can be zero, becoming an exponential function.

    Returns
    -------
    list
        A list of survival probabilities for the corals, calculated using the piecewise
        Weibull-Weibull survival model. Each survival probability corresponds to the age
        in the input list `ages`.
    """
    fx1   = bathtub_curve(lmbda1, v1, sigma1)
    fx2   = bathtub_curve(lmbda2, v2, sigma2)
    decay =[]
    for age in range(0, len(ages)):
        if(ages[age] < Tcp):
            area = quad(fx1, 0, ages[age])[0]
        else:
            area = quad(fx1, 0, Tcp)[0] + quad(fx2, Tcp, ages[age])[0]
        decay.append(math.exp(-area))
    return decay



def single_bathtub_decay(t, lambda_single, v_single, sigma_single):
    """
    Calculate decay values for given time points using a single bathtub curve.

    Parameters
    ----------
    t : array-like
        Time points to calculate decay for.
        lambda_single : float
        Lambda parameter for the bathtub curve.
    v_single : float
        Shape parameter for the bathtub curve.
    sigma_single : float
        Sigma parameter for the bathtub curve.

    Returns
    -------
    list
        A list of survival probabilities for the corals, calculated using a single
        bathtub curve. Each survival probability corresponds to the age
        in the input list `t`.
    """
    fx = bathtub_curve(lambda_single, v_single, sigma_single)
    decay = []      
    # Calculate decay for each time point
    for age in range(len(t)):
        area = quad(fx, 0, t[age])[0]
        decay.append(math.exp(-area))
    
    return decay

def get_decay_parameters(dsst):
    """
    Get decay parameters based on the given temperature.

    Parameters:
    -----------
    dsst : float
        The temperature value to look up parameters for.

    Returns:
    --------
    tuple
        A tuple containing (lmbda1, v1, sigma1) for the given temperature.
        Returns exact values if temperature matches, otherwise None.
    """
    # Create DataFrame inside the function
    decay_para = pd.DataFrame({
        'dSST':[-3,-2,-1,0,1,2,3,4],
        'temp': [25, 26, 27, 28, 29, 30, 31, 32],
        'lmbda1': [2.954000e-02, 1.677537e-03, 1.380000e-04, 3.450817e-05, 1.250000e-05, 3.762058e-06, 1.132246e-06, 3.407661e-07],
        'v1': [0.46120000, 0.30126014, 0.20690000, 0.16515002, 0.13860000, 0.11343958, 0.09284659, 0.07599191],
        'sigma1': [0.0000001275, 1.2623016324, 2.1545000000, 2.3369825925, 2.3833000000, 2.4480132916, 2.5127265832, 2.5774398748]
    })
    
    if dsst in decay_para['dSST'].values:
        row = decay_para[decay_para['dSST'] == dsst].iloc[0]
        return (row['lmbda1'], row['v1'], row['sigma1'])
    else:
        return None
    
def get_competence_parameters(dsst):
    """
    Get competence parameters based on the given temperature.

    Parameters:
    -----------
    dsst : float
        The temperature value to look up parameters for.

    Returns:
    --------
    tuple
        A tuple containing (tc, alpha, beta1, v) for the given temperature.
        Returns exact values if temperature matches, otherwise None.
    """
    compet_para = pd.DataFrame({
        'dsst':[-3,-2,-1,0,1,2,3,4],
        'temp':[25,26,27,28,29,30,31,32], 
        'tc':[5.380000, 5.168125, 4.890000, 4.413125, 3.870000, 3.360000, 2.850000, 2.340000],
        'alpha': [0.4497] * 8,
        'beta1': [0.01623000, 0.02190937, 0.02669000, 0.02877437, 0.02996000, 0.03159500, 0.03323000, 0.03486500],
        'v': [0.3981] * 8
    })
    
    if dsst in compet_para['dsst'].values:
        row = compet_para[compet_para['dsst'] == dsst].iloc[0]
        return (row['tc'], row['alpha'], row['beta1'], row['v'])
    else:
        return None

def piecewise_competence(ages, tc, Tcp, alpha, beta1, beta2, v):
    """
    Calculates the larval competence values at different ages (days), using the piecewise
    Weibull-exponential competence model. This function is a replica of the R code used by
    Moneghetti et al. (2019) to calculate competence.

    Parameters
    ----------
    ages : list
        A list of larvaes ages. Each age must be a positive integer or float.
    tc : float
        The age at which the larvaes reaches their maximum competence level. Must be a
        positive integer or float.
    Tcp : float
        The age at which the larvaes starts to experience a decline in competence. Must be
        greater than tc and a positive integer or float.
    alpha : float
        The scale parameter of the Weibull distribution. Must be a positive integer or float.
    beta1 : float
        The shape parameter of the Weibull distribution in the early decline phase. Must be
        a positive integer or float.
    beta2 : float
        The shape parameter of the Weibull distribution in the late decline phase. Must be
        a positive integer or float.
    v : float
        The exponential decay parameter in the early decline phase. Must be a positive
        integer or float.

    Returns
    -------
    list
        A list of competence values for larvaes, calculated using the piecewise
        Weibull-exponential competence model. Each competence value corresponds to the age
        in the input list `ages`.
    """
    competence = []
    for age in range(0, len(ages)):
        if(ages[age] < tc):
            area = 0
        if(ages[age] >= tc and ages[age] <= Tcp):
            fxtau_early = lambda tau: alpha * math.exp(-alpha * (tau - tc))* math.exp(- pow(( beta1 *(ages[age]-tau)), v))
            area  = quad(fxtau_early, tc, ages[age])[0]
        if(ages[age] > Tcp):
            fxtau_late_first  = lambda tau: alpha * math.exp(-alpha * (tau - tc))* math.exp(- pow(( beta1 *(Tcp-tau)), v)) * math.exp(-beta2 * (ages[age] -Tcp))
            fxtau_late_second = lambda tau: alpha * math.exp(-alpha * (tau - tc))* math.exp(- beta2 *(ages[age]-tau))
            area  = quad(fxtau_late_first, tc, Tcp)[0] + quad(fxtau_late_second, Tcp, ages[age])[0]
        competence.append(area)
    return(competence)


def single_competence(ages, tc, alpha, beta, v):
    """
    Calculates the larval competence values at different ages (days), using the piecewise
    Weibull-exponential competence model. This function is a replica of the R code used by
    Moneghetti et al. (2019) to calculate competence.

    Parameters
    ----------
    ages : list
        A list of larval ages. Each age must be a positive integer or float.
    tc : float
        The age at which the larvae reach their maximum competence level. Must be a
        positive integer or float.
    alpha : float
        The scale parameter of the Weibull distribution. Must be a positive integer or float.
    beta : float
        The shape parameter of the Weibull distribution in the decline phase. Must be
        a positive integer or float.
    v : float
        The exponential decay parameter in the decline phase. Must be a positive
        integer or float.

    Returns
    -------
    list
        A list of competence values for larvae, calculated using the piecewise
        Weibull-exponential competence model. Each competence value corresponds to the age
        in the input list `ages`.
    """
    competence = []
    for age in range(len(ages)):
        if ages[age] < tc:
            area = 0
        else:
            fxtau = lambda tau: alpha * math.exp(-alpha * (tau - tc)) * math.exp(-pow((beta * (ages[age] - tau)), v))
            area = quad(fxtau, tc, ages[age])[0]
        competence.append(area)
    return competence

def get_sst_values(reef_id, scenario):
    """
    Retrieve SST values for a specific reef ID across all years from the appropriate file.

    Parameters:
    -----------
    reef_id : int
        The ID of the reef.
    scenario : str
        The climate scenario ('2p6', '4p5', '7p0', or '8p5').

    Returns:
    --------
    pandas.Series
        A series containing SST values for the specified reef ID across all years.
    """
    # Define file paths for each scenario
    file_paths = {
        '2p6': '/datasets/work/oa-coconet/work/gbr_connectivity/Mandy_Cheung/input_data/4_04_CS126_round_deltaSST.csv',
        '4p5': '/datasets/work/oa-coconet/work/gbr_connectivity/Mandy_Cheung/input_data/4_04_CS245_round_deltaSST.csv',
        '7p0': '/datasets/work/oa-coconet/work/gbr_connectivity/Mandy_Cheung/input_data/4_04_CS370_round_deltaSST.csv',
        '8p5': '/datasets/work/oa-coconet/work/gbr_connectivity/Mandy_Cheung/input_data/4_04_CS585_round_deltaSST.csv'
    }
    
    # Check if the scenario is valid
    if scenario not in file_paths:
        raise ValueError(f"Invalid scenario: {scenario}. Must be one of: {', '.join(file_paths.keys())}")
    
    # Read only the specific row from the CSV file, starting from the 4th column, skipping the header
    sst_values = pd.read_csv(file_paths[scenario], index_col=0, header=0, usecols=range(3, None)).iloc[reef_id]
    
    return sst_values

def calculate_decay_and_competence(dsst_values, settled_age, settled_traj, ntraj):
    """
    Calculate decay, competence, and connectivity vectors for a given scenario and dsst values.
    Parameters:
    -----------
    dsst_values : pandas.Series
        A series containing dsst values for each year of simulation.
    settled_age : numpy.array
        Array of settled ages for particles.
    settled_traj : numpy.array
        Array of trajectory IDs for settled particles.
    ntraj : int
        Total number of trajectories.

    Returns:
    --------
    tuple
    Sum of max connectivity values for each year of simulation.
    """
    max_connectivity_vector = []
    for sst_value in dsst_values:
        tc, alpha, beta, v = get_competence_parameters(sst_value)
        lmbda, v_decay, sigma = get_decay_parameters(sst_value)
        
        decay = single_bathtub_decay(settled_age, lmbda, v_decay, sigma)
        competence = single_competence(settled_age, tc, alpha, beta, v)
        
        # Calculate connectivity and create DataFrame
        connect = pd.DataFrame({
            'connect': np.array(decay) * np.array(competence),
            'traj': settled_traj,
        })
        
        # Group by trajectory and calculate max
        connect_grouped = connect.groupby('traj')
        connect_max = connect_grouped['connect'].max() / ntraj
        total_connect_max = connect_max.sum()

        max_connectivity_vector.append(total_connect_max.to_numpy())

    return max_connectivity_vector

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~ In polygon algorithm and optimizers
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
@jit(nopython=True)
def point_in_polygon(x, y, polygon):
    """
    Determines whether a point is inside a polygon using the ray-casting algorithm.

    Parameters
    ----------
    x : float
        The x-coordinate of the point to test.
    y : float
        The y-coordinate of the point to test.
    polygon : list of tuples
        A list of tuples representing the vertices of the polygon, in order.

    Returns
    -------
    bool
        True if the point is inside the polygon, False otherwise.
    """
    num_vertices = len(polygon)
    is_inside = False
    previous_x, previous_y = polygon[0]
    ## start outside the polygon
    intersection_x = 0.0
    intersection_y = 0.0
    for i in numba.prange(num_vertices + 1):
        current_x, current_y = polygon[i % num_vertices]
        if y > min(previous_y, current_y):
            if y <= max(previous_y, current_y):
                if x <= max(previous_x, current_x):
                    if previous_y != current_y:
                        intersection_x = (y - previous_y) * (current_x - previous_x) / (current_y - previous_y) + previous_x
                    if previous_x == current_x or x <= intersection_x:
                        is_inside = not is_inside
        previous_x, previous_y = current_x, current_y

    return is_inside


@njit(parallel=False)
def points_in_polygon(xs, ys, miny, maxy, polygon):
    """
    Test whether a point is inside a given polygon using the point-in-polygon algorithm.

    This function tests each point in the `points` array to determine if it is inside the polygon
    defined by the vertices in the `polygon` array. The function uses the point-in-polygon Ray casting
    algorithm to perform this test. The algorithm performs the even-odd-rule algorithm to find out
    whether a point is in a given polygon. This runs in O(n) time where n is the number of edges of the polygon.

    Parameters:
    -----------
    points: numpy.ndarray of shape (n, 2)
        Array of n points to test. Each point is defined by its x and y coordinates in columns 0 and 1 respectively.
    polygon: numpy.ndarray of shape (m, 2)
        Array of m vertices defining the polygon. Each vertex is defined by its x and y coordinates in columns 0 and 1 respectively.

    Returns:
    --------
    D: numpy.ndarray of shape (n,)
        Boolean array indicating whether each point in `points` is inside the polygon (`True`) or not (`False`).
    """
    D = np.empty(len(ys), dtype=numba.boolean)
    for i in range(len(D)):
        if ys[i] >= miny and ys[i] <= maxy:
            D[i] = point_in_polygon(xs[i], ys[i], polygon)
        else:
            D[i] = False
    return D


def calc(source_reef):
    try:
        print(f"Processing job {source_reef}") 
        file_name = path + "/GBR1_H2p0_Coral_Release_" + release_start_day + "_Polygon_" +  str(source_reef) + '_Wind_3_percent_displacement_field.nc'
        if not os.path.exists(file_name):
            print('file missing - ' + str(source_reef))
        else:
            output_nc = xr.open_dataset(file_name)
            ntraj     = output_nc.dims['traj']
            particles = pd.DataFrame({
                'latitudes' : output_nc['lat'].values.ravel(),
                'longitudes' : output_nc['lon'].values.ravel(),
                'trajectories' : output_nc['trajectory'].values.ravel(),
                'age' : output_nc['age'].values.ravel() / 86400 ## Seconds to days
            })
            output_nc.close()
            
            print(f"Source reef {source_reef}: Initial particle count: {len(particles)}")
            
            # Cleaning the nans
            particles = particles.dropna()
            print(f"Source reef {source_reef}: Particle count after dropping NaNs: {len(particles)}")
            
            ## remove particles bellow minimum age
            particles = particles[particles['age'] > tc]
            print(f"Source reef {source_reef}: Particle count after age filter: {len(particles)}")
            
            ## set particles boundaries in model domain
            particle_max_lat = np.nanmax(particles['latitudes'].values)
            particle_min_lat = np.nanmin(particles['latitudes'].values)
            print(f"Source reef {source_reef}: Particle latitude range: {particle_min_lat} to {particle_max_lat}")
            
            # making boolean series
            upper_bound = data_shape['min_lat'] <= particle_max_lat
            mmax = upper_bound.to_numpy()
            inf_bound   = data_shape['max_lat'] >= particle_min_lat
            minf = inf_bound.to_numpy()
            boundary_reefs = np.where(np.multiply(minf, mmax))[0]
            print(f"Source reef {source_reef}: Number of potential sink reefs: {len(boundary_reefs)}")
            
            ## Get the dSST for each source reef
            dsst_2p6 = get_sst_values(source_reef, "2p6")
            print(f"Source reef {source_reef}: dSST for 2p6 scenario: {dsst_2p6}")
            dsst_4p5 = get_sst_values(source_reef, "4p5")
            print(f"Source reef {source_reef}: dSST for 4p5 scenario: {dsst_4p5}")
            dsst_7p0 = get_sst_values(source_reef, "7p0")
            print(f"Source reef {source_reef}: dSST for 7p0 scenario: {dsst_7p0}")
            dsst_8p5 = get_sst_values(source_reef, "8p5")
            print(f"Source reef {source_reef}: dSST for 8p5 scenario: {dsst_8p5}")
            
            year_simulations = len(dsst_2p6)
            num_scenarios = 4
            ## Creating empty arrays    
            connectivity_matrix_max   = np.zeros((year_simulations, num_scenarios, num_reefs))
            for sink_reef in boundary_reefs:
                reef_index = data_shape['FID'][sink_reef]
                polygon = np.array(list(data_shape['geometry'][sink_reef].exterior.coords))
                if(particles.size == 0):
                    break ## breaking the loop if not more particles
                # Are these points inside the polygon?
                p = points_in_polygon(
                    particles['longitudes'].values,
                    particles['latitudes'].values,
                    data_shape.min_lat[sink_reef],
                    data_shape.max_lat[sink_reef],
                    polygon)
                m = np.where(p)[0]
                if(len(m > 0)):
                    settled_age  = particles['age'].iloc[m].values
                    settled_traj = particles['trajectories'].iloc[m].values
                    ## get the model decay and competence for the entired projection
                    # Calculate decay and competence for the 2p6 scenario
                    connectivity_2p6 = calculate_decay_and_competence(dsst_2p6, settled_age, settled_traj, ntraj)
                    # Calculate decay and competence for the 4p5 scenario
                    connectivity_4p5 = calculate_decay_and_competence(dsst_4p5, settled_age, settled_traj, ntraj)
                    # Calculate decay and competence for the 7p0 scenario
                    connectivity_7p0 = calculate_decay_and_competence(dsst_7p0, settled_age, settled_traj, ntraj)
                    # Calculate decay and competence for the 8p5 scenario
                    connectivity_8p5 = calculate_decay_and_competence(dsst_8p5, settled_age, settled_traj, ntraj)
                    ## Save the connectivity for all years of simulation in the scenarios number and the sink_reef (reef index)
                    connectivity_matrix_max[:, 0, reef_index] = connectivity_2p6
                    connectivity_matrix_max[:, 1, reef_index] = connectivity_4p5
                    connectivity_matrix_max[:, 2, reef_index] = connectivity_7p0
                    connectivity_matrix_max[:, 3, reef_index] = connectivity_8p5
                
                    particles.drop(index = particles.iloc[m].index, inplace=True)
        return connectivity_matrix_max
    except Exception as e:
        print(f"Error processing job {source_reef}: {e}")
        return None

# ## ~~~~~~~~~~~~~~~~~~~~ ##
# ## ~      Parameters  ~ ##
# ## ~~~~~~~~~~~~~~~~~~~~ ##
# ## Parameters used in the decay and competence function
# ## Decay
# Tcp_decay = 2.583
# lmbda1    = 0.4
# lmbda2    = 0.019
# v1        = 2.892
# v2        = 1.716
# sigma1    = 0
# sigma2    = 0

# ## Competence
# tc       = 3.333
# Tcp_comp = 69.91245
# alpha    = 1.295
# beta1    = 0.001878001
# beta2    = 0.3968972
# v        = 0.364

## standar minimum tc = 2.340000. I used this value because is the minimum value, so I have all option cover
tc=2.340000

## ~~~~~~~~~~~~~~~ ##
## ~   Main Code ~ ##
## ~~~~~~~~~~~~~~~ ##
shapefile = '/datasets/work/oa-coconet/work/oceanparcels_gbr_Coral/Shape_files/gbr1_coral_1m_merged_buffer0p001.shp'
data_shape = gpd.read_file(shapefile)

## getting the boundaries of each reefs for
## entired GBR
num_reefs = data_shape.shape[0]
min_lat = []; max_lat = []
for i_polygon in range(0, num_reefs):
    min_lat.append(np.nanmin(np.array(data_shape['geometry'][i_polygon].bounds)[[1,3]]))
    max_lat.append(np.nanmax(np.array(data_shape['geometry'][i_polygon].bounds)[[1,3]]))
data_shape['min_lat']= min_lat
data_shape['max_lat']= max_lat

release_start_day = sys.argv[1]
print(f"release_start_day: {release_start_day}")
scenarios = ['2p6', '4p5', '7p0', '8p5']
len_scenarios = len(scenarios)
years_projection = list(range(2025, 2101))
len_years_projection = len(years_projection)
path='/datasets/work/oa-coconet/work/OceanParcels_outputs/Coral/' + release_start_day
print(f"path: {path}")
jobs = range(num_reefs)
n_jobs = int(os.getenv('SLURM_CPUS_ON_NODE', 10))
## print the parameters
print(f"Number of reefs: {num_reefs}")
print(f"Release start day: {release_start_day}")
print(f"Scenarios: {scenarios}")
print(f"Years projection: {years_projection}")
print(f"Number of jobs: {n_jobs}")
with parallel_backend(backend='loky', n_jobs=n_jobs):
    results_list = Parallel()(delayed(calc)(k) for k in jobs)

print('calculations done', time.strftime("%H:%M:%S"), flush = True)

## Creating empty arrays
connectivity_matrix_max = np.zeros((len_years_projection, len_scenarios, num_reefs, num_reefs))

for k in jobs:
    connectivity_matrix_max[:, :, k, :] = results_list[k]


## generate a netcdf file with the connectivity data
# Define chunk sizes

ds = xr.Dataset(
    {
        "connectivity": (["year", "scenario", "source_reef", "target_reef"], 
                         connectivity_matrix_max.astype(np.float32))
    },
    coords={
        "year": years_projection,
        "scenario": scenarios,
        "source_reef": np.arange(num_reefs),
        "target_reef": np.arange(num_reefs),
    },
)

# Add metadata to the dataset
ds.attrs["description"] = "Connectivity matrix for coral reefs under different climate scenarios (SIMIP 6 scenarios: 2p6, 4p5, 7p0, 8p5) for the years 2025 to 2100."
ds.attrs["units"] = "connectivity values"
ds.attrs["release_start_day"] = release_start_day
ds.attrs["created_on"] = time.strftime("%Y-%m-%d %H:%M:%S")
ds.attrs["created_by"] = "Javier Porobic, email: javier.porobicgarate@csiro.au"

# Save the dataset to a NetCDF file
output_file = f"connectivity_matrix_{release_start_day}.nc"
encoding = {
    'connectivity': {
        'zlib': True,
        'complevel': 9,
        'dtype': 'float32'
    }
}
ds.to_netcdf(output_file, encoding=encoding)
print(f"NetCDF file created: {output_file}")

