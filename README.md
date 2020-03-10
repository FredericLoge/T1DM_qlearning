# T1DM_qlearning
# R File descriptions

## Simulation

@_simulate_data.R
In this R file, we simulate data to be used for Reinforcement Learning, relying on the simglucose project on GitHub. The controller is randomized.

@_simulate_data_for_grid_search.R
Same as @_simulate_data.R for several baseline bolus advisors.

@_simulate_data_spec_policy.R
Same as @_simulate_data.R for a deterministic policy.

## Q-Learning w/ function approximation

@_parametrized_qlearning_V3.R
We apply our function approximation Q-Learning algorithm to simulated data.

@_pql_V3_foo.R
Functions used in @_parametrized_qlearning_V3.R

## Visualization

@_viz_results.R
Visualization of the policy found + daily blood glucose profiles
