"""Get POD values from folder of outputs."""
# %%
import numpy as np
import sys
import json
from matplotlib import pyplot as plt
import csv
import os


# # %%

# files = os.listdir('../icp_runs/10_An_match_params/inputs')
# files

# for f in files:
#     get_podY(f)

# %%


def get_headers(file):
    with open(file) as f:
        return f.readline().replace('\n', '').replace('\r\n', '').split(',')


def get_last_line_in_file(filename, headers):
    try:
        with open(filename, 'rb') as f:
            f.seek(-2, os.SEEK_END)
            while f.read(1) != b'\n':
                f.seek(-2, os.SEEK_CUR)
            last_line = f.readline().decode().replace(
                '\n', '').replace('\r\n', '').replace('\r', '').split(',')
            last_line_data = {h: d for h, d in zip(headers, last_line)}
            return last_line_data
    except OSError:
        Warning(f'Os Error while getting last line in file {filename}')
        return {h: None for h in headers}


def get_last_line_of_files(dir, files, headers):
    last_lines = [get_last_line_in_file(dir + '/' + f, headers) for f in files]
    return last_lines


def get_target_pody(target_pod_file):
    with open(target_pod_file, encoding='utf-8-sig', mode='r') as t_pod_data:
        spamreader = csv.reader(t_pod_data, delimiter=',')
        data = dict(row for row in spamreader)
        return data


def get_relative_yield(yield_file):
    with open(yield_file, encoding='utf-8-sig', mode='r') as yield_data:
        spamreader = csv.reader(yield_data, delimiter=',')
        data = [row for row in spamreader]
        return data


def extract_pod_from_last_line(line_data):
    return line_data['pody']


def get_pod_data(input_file_directory):
    files = [f for f in os.listdir(input_file_directory) if f.split('.')[-1] == 'csv']
    headers = get_headers(input_file_directory + "/" + files[0])
    last_lines = get_last_line_of_files(input_file_directory, files, headers)
    pod_data = [extract_pod_from_last_line(line_data) for line_data in last_lines]
    file_names = [f.split('.')[0] for f in files]
    return dict(zip(file_names, pod_data))


def get_and_save_pod_data(input_file_directory, output_directory):
    dirs = [d for d in os.listdir(input_file_directory)
            if os.path.isdir(f"{input_file_directory}/{d}")]
    print(dirs)
    pod_data = [get_pod_data(f'{input_file_directory}/{d}') for d in dirs]
    for dir, pods in zip(dirs, pod_data):
        name = os.path.basename(dir)
        print(f"Getting pod for{name}")
        os.makedirs(output_directory, exist_ok=True)
        with open(f"{output_directory}/{name}.json", 'w') as out_file:
            json.dump(pods, out_file)
        with open(f"{output_directory}/{name}.csv", 'w', newline='') as out_file:
            spamwriter = csv.writer(out_file, delimiter=',',
                                    quotechar='|', quoting=csv.QUOTE_MINIMAL)
            for k, v in pods.items():
                spamwriter.writerow([k, v])


if __name__ == "__main__":
    args = sys.argv[1:]
    input_file_directory, output_directory = args
    get_and_save_pod_data(input_file_directory, output_directory)

# # %%
# # ==============================================================================================================
# # %%
# dirs = os.listdir('../icp_runs/14_An_sweden_vary_vcmax/outputs')
# pod_data = [get_pod_data(f'../icp_runs/14_An_sweden_vary_vcmax/outputs/{d}') for d in dirs]
# pod_data
# # %%

# target_pod_data = get_target_pody('../icp_runs/14_An_sweden_vary_vcmax/ideal_pod_values.csv')
# target_pod_data

# # %%
# target_pod_data['SE87_CF_thermal']
# # %%
# multi_out = {
#     "SE87_CF_thermal": 0.00,
#     "SE87_NF_thermal": 0.06,
#     "SE87_NF+_thermal": 2.08,
#     "SE88_CF_thermal": 0.02,
#     "SE88_NF_thermal": 0.54,
#     "SE88_NF+_thermal": 3.81,
#     "SE88_NF++_thermal": 5.32,
#     "SE94_NF23_thermal": 2.50,
#     "SE94_NF+27_thermal": 4.06,
#     "SE94_NF++26_thermal": 5.14,
#     "SE94_NF+++15_thermal": 6.29,
#     "SE95_NF25_thermal": 1.65,
#     "SE95_NF+27_thermal": 2.95,
#     "SE95_NF+post_thermal": 5.65,
#     "SE97_CF23_thermal": 0.06,
#     "SE97_NF30_thermal": 0.53,
#     "SE97_NF+27_thermal": 4.01,
#     "SE97_NF++25_thermal": 6.20,
#     "SE99_CF11_thermal": 0.01,
#     "SE99_NF+14_thermal": 6.74,
# }

# a = []
# a_names = [d for d in dirs]
# m = [float(multi_out.get(k, None)) for k in pod_data[0].keys()]
# t = [float(target_pod_data.get(k, None)) for k in pod_data[0].keys()]
# xtik = [k for k in pod_data[0].keys()]
# for d in pod_data:
#     a_in = []
#     for k, v in d.items():
#         a_in.append(float(v))
#     a.append(a_in)

# # %%
# plt.plot(m, label="m")
# plt.plot(t, label="t")
# for p, name in zip(a, a_names):
#     plt.plot(p, label=name)
# plt.xticks(range(len(xtik)), xtik, rotation='vertical')
# plt.legend()
# # plt.scatter(m, t)
# # %%


# # %%

# # ============ GET R squared

# y = [float(yield_data.get(k, None)) for k in pod_data[0].keys()]
# y
# # %%
# for x in a:
#     print(len(x), len(y))
#     plt.scatter(x, y)

# # %%


# # Polynomial Regression
# def polyfit(x, y, degree):
#     results = {}

#     coeffs = np.polyfit(x, y, degree)

#     # Polynomial Coefficients
#     results['polynomial'] = coeffs.tolist()

#     # r-squared
#     p = np.poly1d(coeffs)
#     # fit values, and mean
#     yhat = p(x)                         # or [p(z) for z in x]
#     ybar = np.sum(y) / len(y)          # or sum(y)/len(y)
#     ssreg = np.sum((yhat - ybar)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
#     sstot = np.sum((y - ybar)**2)    # or sum([ (yi - ybar)**2 for yi in y])
#     results['determination'] = ssreg / sstot

#     return results


# # %%
# vcmax_vals = [float(n.split('_')[-1]) for n in a_names]
# vcmax_vals

# # %%
# r_vals = []
# for x, v in zip(a, vcmax_vals):
#     plt.scatter(x, y)
#     fit = polyfit(x, y, 1)
#     m, b = fit['polynomial']
#     r = fit['determination']
#     x = np.array(x)
#     plt.plot(x, m * x + b, label=v)
#     r_vals.append(r)
# plt.xlabel("pody")
# plt.ylabel("yield")
# plt.legend()

# # %%
# r_vals_dict = dict(zip(a_names, r_vals))
# r_vals_dict


# # %%
# vals = [x for x in sorted(zip(vcmax_vals, r_vals), key=lambda pair: pair[0])]
# vals
# # %%
# x = [x for x, y in vals]
# y = [y for x, y in vals]
# plt.scatter(x, y)

# plt.xlabel("vcmax")
# plt.ylabel("Rsq")
# # %%
