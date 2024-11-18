# README

## Overview
This repository contains code to simulate and analyze spatio-temporal data using the Spatio-Temporal Jump Model (ST-JM). The model efficiently handles spatio-temporal clustering with missing data and incorporates spatial and temporal dependencies between variables. 
For detailed information on the model specification, please refer to the associated paper available [here](https://arxiv.org/pdf/2411.09726).
The repository includes two main files:

1. **`Utils_.R`**: Contains utility functions for simulating spatio-temporal data, managing missing data, and running the ST-JM.
2. **`Example.R`**: Provides an example of how to use the ST-JM with simulated data and visualize the results.

---

## Files

### `Utils_.R`

This file includes the following key functions:

1. **`simulate_observations()`**:  
   Simulates observations from a multivariate normal distribution based on a latent state sequence. The function allows for continuous and categorical variables, along with missing data imputation capabilities.

2. **`generate_spatio_temporal_data()`**:  
   Generates spatio-temporal data with spatial and temporal persistence. The function returns the simulated states, complete and incomplete data matrices, and spatial points information.

3. **`generate_spatial_points()`**:  
   Generates random spatial coordinates for a given number of points.

4. **`order_states_condMean()`**:  
   Orders states based on their conditional means to ensure consistent state labeling across time and space.

5. **`STjumpDist()`**:  
   Implements the Spatio-Temporal Jump Model estimation algorithm for clustering spatio-temporal data. It accounts for both temporal transitions and spatial relationships between data points using penalty terms. Input data Y must be a data.frame in long format with columns m (identifier for spatial point) and t (time).

6. **`get_cat()`**:  
   Simulates categorical outcomes based on the latent states, incorporating state-conditional probabilities.

7. **`punct()`**:  
   Adds missing values to a dataset to simulate real-world data conditions.

8. **`initialize_states()`**:  
   Initializes the state sequence for the spatio-temporal clustering algorithm using a distance-based approach.

9. **`BAC()`**:  
   Computes the balanced accuracy between two clustering solutions, allowing for the evaluation of clustering performance.

---

### `Example.R`

This file demonstrates how to use the functions in `Utils_.R` to:

- Simulate spatio-temporal data.
- Run the Spatio-Temporal Jump Model (ST-JM).
- Visualize the clustering results across time and space.
  
It includes examples of generating spatial data points, running the clustering algorithm, and producing visualizations that illustrate the identified clusters.

---

## How to Use

1. Clone the repository and install the required packages (`MASS`, `cluster`, `StatMatch`, `caret`).
   
2. Load the utility functions from `Utils_.R`.
   
3. Run the example in `Example.R` to simulate spatio-temporal data and perform clustering using the ST-JM.

---

## Dependencies
This project requires the following R packages:
- `MASS`
- `cluster`
- `StatMatch`
- `caret`
