from multiprocessing import Queue, Process
import numpy as np
from copy import deepcopy
from pyDO3SE.tools.get_pod_values import get_relative_yield, get_target_pody
from pyDO3SE.run_model import run_model
from pyDO3SE.setup_model import Main_Overrides, setup_config, setup_external_state, setup_initial_state, setup_model, setup_model_processes
from pyDO3SE import main
from multiprocessing import Pool
import warnings
import os
from collections import namedtuple
from paramga.run import run_parallel as run_ga

# %%

# output_directory = 'temp'
# config_dir = 'examples/multirun/configs'
# input_dir = 'examples/multirun/inputs'
# config_file_type = 'json'


RunOutput = namedtuple('RunOutput', 'result error')


# def run_from_args(args):
#     log_level = 0
#     try:
#         main.main(
#             **args,
#             log_level=log_level,
#             output_results_only=True,
#         )
#         return RunOutput(0, None)
#     except Exception as e:
#         config_file = args['config_location']
#         input_file = args['data_location']
#         warnings.warn(f"Failed to run config {config_file} on data {input_file}")
#         return RunOutput(1, e)


# %%
config_location = '../icp_runs/16_An_sweden_vary_vcmax/configs_an_gsto/Sweden_An_gsto_90.json'
overrides = Main_Overrides()
base_model_config = setup_config(config_location, overrides)
input_dir = '../icp_runs/16_An_sweden_vary_vcmax/inputs'
input_files = os.listdir(input_dir)[0:7]
input_file_names = [f.split('.')[0] for f in input_files]
external_states = [setup_external_state(
    base_model_config, f'{input_dir}/{data_location}', overrides) for data_location in input_files]
# %%
start_days = [external_state.dd[0] - 1 for external_state in external_states]
end_days = [external_state.dd[-1] for external_state in external_states]
overrides_in = [overrides._replace(start_day=start_day, end_day=end_day)
                for start_day, end_day in zip(start_days, end_days)]
# %%
model_processes_all = [setup_model_processes(
    base_model_config, overrides) for overrides in overrides_in]
# %%

# %%
initial_states: List[ModelState] = [setup_initial_state(base_model_config, external_state, overrides)
                                    for external_state, overrides in zip(external_states, overrides_in)]


# target_pod_vals = get_target_pody('../icp_runs/10_An_match_params/ideal_pod_values.csv')
yield_data = dict(get_relative_yield('../ga_model/relative_yield.csv'))
yield_data = [float(yield_data.get(k, None)) for k in input_file_names]
yield_data

# %%
base_model_config
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


def run():
    def model(params, data):
        global initial_states, model_processes_all, external_states
        print("v_cmax_25 ", params['v_cmax_25'])
        args_to_run = []
        print("Getting args")
        diff = 0

        Q = Queue()

        def q_wrap(q, args):
            final_state, output_logs = run_model(*args)
            q.put(final_state.canopy_layer_component[0][0].POD_Y)

        procs = []
        outputs = []

        for name, external_state, initial_state, model_processes in zip(input_file_names, external_states, initial_states, model_processes_all):
            print(f"Getting args for {name}")
            config = deepcopy(base_model_config)
            config.Land_Cover.parameters[0].pn_gsto.V_cmax_25 = params['v_cmax_25']
            config.Land_Cover.parameters[0].pn_gsto.J_max_25 = params['v_cmax_25'] * 2
            # args_to_run.append((initial_state, config, external_state, model_processes))

            p = Process(target=q_wrap, args=(
                [Q, [initial_state, config, external_state, model_processes]]))
            procs.append(p)
            p.start()
            # Local run
            # final_state, output_logs = run_model(
            #     initial_state, config, external_state, model_processes)
            # final_pod = final_state.canopy_layer_component[0][0].POD_Y
            # target_pod = target_pod_vals[name]
            # print(final_pod, target_pod)
            # diff += abs(float(final_pod) - float(target_pod))

        for p in procs:
            res = Q.get()
            outputs.append(res)
        for p in procs:
            p.join()
        # TODO: Get the Rsq value here

        fit = polyfit(outputs, yield_data, 1)
        rsq = fit['determination']
        # diff = sum(abs(float(final_pod) - float(target_pod))
        #            for final_pod, target_pod in zip(outputs, target_pod_vals_list))

        print(f"---diff {rsq}")
        return 1 - rsq

    return model


# %%
model = run()
model({"v_cmax_25": 100}, [])
# %%


def loss_function(output, params):
    return output


# %%
base_config = {
    "v_cmax_25": 100,
}
ga_conf = {
    "v_cmax_25": {
        "type": "number",
        "min": 40,
        "max": 200,
        "step": 8,
    }
}
# %%
out = run_ga(
    base_config,
    ga_conf,
    run(),
    loss_function,
    [],
    max_iterations=10,
    verbose=True,
)
# %%
out
# %%
