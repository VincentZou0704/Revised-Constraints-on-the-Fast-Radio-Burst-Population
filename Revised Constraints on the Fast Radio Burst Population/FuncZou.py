import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as IUS
from FRBpopulation.prepare import dm_model
import matplotlib.pyplot as plt
from FRBpopulation.setup import *


precision = 300    # the precision of integration/ decide the unique point
precision_vec = 300     # the precision of integration_vec/ decide the smooth
# In M1 corner plot, p = 300, pv = 50
func_precision = 100    # the precision of building function of CDF at (0,1) to sample


def integrate(Fun, lower_limit, upper_limit, *args, **kwargs):
    x = np.linspace(lower_limit, upper_limit, precision)
    y = Fun(x, *args,**kwargs)
    dx = (upper_limit - lower_limit) / precision
    result = (y * dx).sum(axis = 0)
    return result


def integrate_vec(Fun, lower_limit, upper_limit):
    x = np.linspace(lower_limit, upper_limit, precision_vec)
    arr = Fun(x[0])
    X, Y= np.meshgrid(x, arr)
    y = Fun(X.T)
    dx = (upper_limit - lower_limit) / precision_vec
    result = (y * dx).sum(axis = 0)
    return result


def stay_increase(arr):
    index, new_arr = [0], [arr[0]]
    for i in range(len(arr)-1):
        if arr[i] < arr[i+1]:
            index.append(i+1)
            new_arr.append(arr[i+1])
        else:
            arr[i+1] = arr[i]

    return np.array(new_arr), index


def save_by_index(arr, index):
    new_arr = []
    for index_i in index:
        new_arr.append(arr[index_i])

    return new_arr


def func_build(x, arr_a, arr_b, mode = 'sline'):
    arr_a_increase = stay_increase(arr_a)[0]
    index = stay_increase(arr_a)[1]
    arr_b_increase = save_by_index(arr_b, index)
    if mode == 'curve':
        spl = IUS(arr_a_increase, arr_b_increase)
        return spl(x)
    if mode == 'sline':
        return np.interp(x, arr_a_increase, arr_b_increase)


def get_sample(pdf, lower_limit, upper_limit, nun, *args, **kwargs):
    normalization = 1 / integrate(pdf, lower_limit, upper_limit, *args, **kwargs)
    x_segment = np.linspace(lower_limit, upper_limit, func_precision)
    cdf = normalization * integrate(pdf, lower_limit, x_segment, *args, **kwargs)
    data = np.random.rand(nun)
    sample = func_build(data, cdf, x_segment)
    return sample


if __name__ == '__main__':

    args = [sigma_host0, emu0, f_IGM_p, alpha0, F0]
    dm0 = 1000
    print(get_sample(dm_model.p_z_at_dm, 0.001, 3, 1, dm0, sigma_host0, emu0, f_IGM_p, alpha0, F0))

