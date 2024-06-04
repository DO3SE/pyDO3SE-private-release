"""Calculate the R Sq values of pod data."""
from pyDO3SE.tools.get_pod_values import get_pod_data, get_relative_yield, get_target_pody
import numpy as np
import sys
import json
from matplotlib import pyplot as plt
import csv
import os
# %%

# Polynomial Regression
def polyfit(x, y, degree):
    results = {}

    coeffs = np.polyfit(x, y, degree)

    # Polynomial Coefficients
    results['polynomial'] = coeffs.tolist()

    # r-squared
    p = np.poly1d(coeffs)
    # fit values, and mean
    yhat = p(x)                         # or [p(z) for z in x]
    ybar = np.sum(y) / len(y)          # or sum(y)/len(y)
    ssreg = np.sum((yhat - ybar)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
    sstot = np.sum((y - ybar)**2)    # or sum([ (yi - ybar)**2 for yi in y])
    results['determination'] = ssreg / sstot

    return results


# %%
r_vals = []
x = [1,2,3, 4, 5, 6]
y = [8, 6, 5, 4, 4, 4]
fit = polyfit(x, y, 1)
m, b = fit['polynomial']
r = fit['determination']
x = np.array(x)
plt.scatter(x, y)

plt.plot(x, m * x + b)
r
# ==============================================================================================================
# %%
dirs = os.listdir('../icp_runs/14_An_sweden_vary_vcmax/outputs')
pod_data = [get_pod_data(f'../icp_runs/14_An_sweden_vary_vcmax/outputs/{d}') for d in dirs]
pod_data
# %%
multi_out = {
    "SE87_CF_thermal": 0.00,
    "SE87_NF_thermal": 0.06,
    "SE87_NF+_thermal": 2.08,
    "SE88_CF_thermal": 0.02,
    "SE88_NF_thermal": 0.54,
    "SE88_NF+_thermal": 3.81,
    "SE88_NF++_thermal": 5.32,
    "SE94_NF23_thermal": 2.50,
    "SE94_NF+27_thermal": 4.06,
    "SE94_NF++26_thermal": 5.14,
    "SE94_NF+++15_thermal": 6.29,
    "SE95_NF25_thermal": 1.65,
    "SE95_NF+27_thermal": 2.95,
    "SE95_NF+post_thermal": 5.65,
    "SE97_CF23_thermal": 0.06,
    "SE97_NF30_thermal": 0.53,
    "SE97_NF+27_thermal": 4.01,
    "SE97_NF++25_thermal": 6.20,
    "SE99_CF11_thermal": 0.01,
    "SE99_NF+14_thermal": 6.74,
}

a = []
a_names = [d for d in dirs]
m = [float(multi_out.get(k, None)) for k in pod_data[0].keys()]
xtik = [k for k in pod_data[0].keys()]
for d in pod_data:
    a_in = []
    for k, v in d.items():
        a_in.append(float(v))
    a.append(a_in)

# %%
plt.plot(m, label="m")
for p, name in zip(a, a_names):
    plt.plot(p, label=name)
plt.xticks(range(len(xtik)), xtik, rotation='vertical')
plt.legend()
# plt.scatter(m, t)
# %%


# %%

# ============ GET R squared
yield_data = dict(get_relative_yield('../ga_model/relative_yield.csv'))
y = [float(yield_data.get(k, None)) for k in pod_data[0].keys()]
y
# %%
for x in a:
    print(len(x), len(y))
    plt.scatter(x, y)

# %%



# %%
vcmax_vals = [float(n.split('_')[-1]) for n in a_names]
vcmax_vals

# %%
r_vals = []
for x, v in zip(a, vcmax_vals):
    plt.scatter(x, y)
    fit = polyfit(x, y, 1)
    m, b = fit['polynomial']
    r = fit['determination']
    x = np.array(x)
    plt.plot(x, m * x + b, label=v)
    r_vals.append(r)
plt.xlabel("pody")
plt.ylabel("yield")
plt.legend()

# %%
r_vals_dict = dict(zip(a_names, r_vals))
r_vals_dict


# %%
vals = [x for x in sorted(zip(vcmax_vals, r_vals), key=lambda pair: pair[0])]
vals
# %%
x = [x for x, y in vals]
y = [y for x, y in vals]
plt.scatter(x, y)

plt.xlabel("vcmax")
plt.ylabel("Rsq")
# %%

m = [float(multi_out.get(k, None)) for k in pod_data[0].keys()]
y = [float(yield_data.get(k, None)) for k in pod_data[0].keys()]
plt.scatter(m, y)
fit = polyfit(m, y, 1)
fit
