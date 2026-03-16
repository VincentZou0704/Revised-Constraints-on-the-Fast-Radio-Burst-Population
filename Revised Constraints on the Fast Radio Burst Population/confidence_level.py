import numpy as np
from FRBpopulation import FuncZou as fz
from FRBpopulation.setup import *
from FRBpopulation.prepare import dm_model
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.interpolate import interp2d
from scipy import optimize, stats
import FRBpopulation.prepare.CHIME_Lin as CL
import FRBpopulation.MCMC.pzE as pzE
import FRBpopulation.prepare.Fluence as fce


# Simulated data
x = np.array([1, 2, 3, 4, 5])
y = np.array([1.2, 2.3, 3.8, 4.7, 5.9])

# Model function
def model_func(x, a, b):
    return a * x + b

# Residuals function
def residuals(params):
    return y - model_func(x, *params)

# Fit the model and estimate parameter values
params_fit, cov_matrix = optimize.curve_fit(model_func, x, y)

# Calculate the residuals
residuals_fit = residuals(params_fit)

# Degrees of freedom (number of data points minus number of model parameters)
df = len(x) - len(params_fit)

# Set the confidence level (e.g., 95%)
confidence_level = 0.95

# Calculate the critical value from the chi-square distribution
critical_value = stats.chi2.ppf(confidence_level, df)

# Calculate the standard errors of the parameter estimates
param_errors = np.sqrt(np.diag(cov_matrix))

# Calculate the confidence interval
lower_bound = params_fit - critical_value * param_errors
upper_bound = params_fit + critical_value * param_errors

confidence_interval = list(zip(lower_bound, upper_bound))

print("Parameter estimates:", params_fit)
print("Confidence interval:", confidence_interval)

