"""Run comparisons of outputs."""
import os
from typing import List, Tuple
import pandas as pd
import warnings

from pyDO3SE.version import version
from pyDO3SE.Output.Output_Shape import output_fields_map
from pyDO3SE.Analysis.charts import multi_series_annual_graph


def pad_data(data, start, end, tstart, tend, pad_val=None):
    MAX_LENGTH = 600  # Maximum span for data from 0
    pad_value_start = pad_val if pad_val is not None else data[0]
    pad_value_end = pad_val if pad_val is not None else data[-1]
    data_out = list(pad_value_start for _ in range(start)) + data + \
        list(pad_value_end for _ in range(MAX_LENGTH))
    data_out = data_out[tstart:tend]
    return data_out


def is_output_file(fl):
    return fl[-4:] == '.csv' and fl != "external_data.csv"
    # return fl[-4:] == '.csv' and fl != "external_data.csv" and "latest" not in fl


def log(log_level):
    def _inner(*args, **kwargs):
        if log_level:
            print(*args, **kwargs)
    return _inner


def get_output_files(output_data_dir, directory_filter="") -> List[List[str]]:
    """Get each output file in directory.

    This filters out empty directories and directories with no csv files
    It also selects the first csv file in the directory that is not the external data file

    """
    out_dir = os.path.normpath(output_data_dir)
    return [[
        os.path.basename(dr),  # Dir name
        dr,  # Directory path
        list(sorted([fl for fl in f
                     if is_output_file(fl)])),  # Target file in dir
    ]
        for dr, d, f in sorted(os.walk(out_dir))
        if f and len([fv for fv in f if fv[-4:] == '.csv']) > 0 \
        and os.path.isdir(dr) \
        and directory_filter == "" or directory_filter in os.path.basename(dr)]
    # and os.path.basename(dr) != "latest"]


def get_series_titles(data_files):
    return [f"{name}-{f.split('.')[0]}"
            for name, _, files in data_files for f in files]


def create_comparison_graphs(
    data_files: List[List[str]],
    output_data_dir: str,
    fields_to_graph: List[str],
    output_dir: str,
    start_day: int = 0,
    end_day: int = 365,
    use_versioned_outdir: bool = True,
    use_versioned_outfile: bool = True,
    use_daily_average: bool = True,
    run_id: str = "",
    log_level: int = 1,
    *args,
    **kwargs,
):
    """Compare outputs in examples/outputs directory against DO3SE-UI.

    Parameters
    ----------
    output_data_dir : str
        The directory containing model output data within sub folders
    fields_to_graph : List[str]
        Fields to create graphs for
    output_dir : str, optional
        The location to save graphs
    use_versioned_outdir : bool, optional
        If true then save outputs in a sub directory with the version, by default True

    *args, **kwargs passed to plt.plot
    """
    logger = log(log_level)
    output_data_dir = os.path.normpath(output_data_dir)
    logger(
        f"\n\n====== Collecting files in {output_data_dir} ======\n===========================\n")

    # make sure outputdir exists
    output_dir_versioned = f'{output_dir}/{version}' if use_versioned_outdir else output_dir
    os.makedirs(output_dir_versioned, exist_ok=True)

    series_titles = get_series_titles(data_files)

    data_sets = {f: [] for f in fields_to_graph}
    logger(f"Found the following outputs \n {series_titles}")

    logger("\n\n====== Preparing data ======\n===========================\n")

    # TODO: Use this to auto set x lims
    full_data_start_day = None
    full_data_end_day = None


    for name, dir, files in data_files:
        for file in files:
            logger(name, " | File: ", file)
            data = pd.read_csv(os.path.join(dir, file))

            if "dd" in data and data['dd'] is not None:
                data_start_day = data['dd'].astype(
                    int).iloc[0] - 1 if "dd" in data and data['dd'] is not None else start_day
                data_end_day = data['dd'].astype(
                    int).iloc[-1] if data['dd'] is not None else end_day
            else:
                data_start_day = start_day
                data_end_day = end_day

            assert data_start_day is not None
            assert data_end_day is not None

            full_data_start_day = min(full_data_start_day, data_start_day) if full_data_start_day is not None else data_start_day
            full_data_end_day = max(full_data_end_day, data_end_day) if full_data_end_day is not None else data_end_day

            for f in fields_to_graph:
                # Note if data does not exist we set to 0
                if f not in data.columns:
                    warnings.warn(f'{f} not in {name} data')
                    field_data = [0 for _ in range(start_day, end_day)]
                    data_sets[f].append(field_data)
                elif type(data[f][0]) == str:
                    warnings.warn(f'{f} is str in {name} data')
                    print(f'{f} is str in {name} data')
                    # raise ValueError("str in data!")
                    field_data = [0 for _ in range(start_day, end_day)]
                    data_sets[f].append(field_data)
                else:
                    field_data = pad_data(data[f].values.tolist(), data_start_day * 24,
                                          data_end_day * 24, start_day * 24, end_day * 24)
                    data_sets[f].append(field_data)
    logger(f'\n\nStart day: {full_data_start_day}')
    logger(f'end day: {full_data_end_day}')
    logger("\n\n===== Printing Graphs ===== \n=========================== \n")
    failed_graphs = []
    for f in fields_to_graph:
        try:
            logger(f"== Generating graph for {f}\n")
            if len(data_sets[f]) == 0:
                warnings.warn(f"\nData set {f} is empty!\nIt must not exist in any output.")
                continue
            # data_sets_filtered = [d for d in data_sets[f] if not isinstance(d[0], str)]
            min_y = 0
            try:
                min_y = min([min(d) for d in data_sets[f]])
            except TypeError as e:
                print(f"Failed to get min value of {f}")
                continue
                # raise e
            except ValueError as e:
                print(f'Failed to get min value of {f}')
            bottom_y_lim = min_y if min_y < 0 else 0
            try:
                multi_series_annual_graph(
                    f,
                    data_sets[f],
                    series_titles,
                    output_fields_map[f],
                    output_dir=output_dir_versioned,
                    label_x_days=30 if use_daily_average else 1,
                    average_step=24 if use_daily_average else 1,
                    chart_id=version if use_versioned_outfile else run_id,
                    linestyles=[{
                        'style': 'dashed' if s == 'DO3SE_ui' else 'solid',
                        'width': 0.3,
                        'zorder': 100 if s == 'DO3SE_ui' else i,
                    } for i, s in enumerate(series_titles)],
                    bottom_y_lim=bottom_y_lim,
                    start_day=start_day,
                    end_day=end_day,
                    figsize=(8, 3),
                    dpi=150,
                    log_level=log_level,
                    *args,
                    **kwargs,
                )
            except Exception as e:
                raise e
        except Exception as e:
            print(f'Failed to create graph for {f}')
            failed_graphs.append((f, e))
    if len(failed_graphs) > 0:
        logger(f'The following charts failed:\n')
        logger('\n'.join([
            f'- {f}' for f, e in failed_graphs
        ]))
    logger(f'Outputs saved to {output_dir}')


pn_output_fields = [
    'gsto_l',
    'gsto_bulk',
    'gsto_canopy',
    'gsto',
    'rsto',
    'rsto_c',
    'rsto_l',
    "td_dd",
    # 'ts_c',
    'dvi',
    'PARsun',
    'PARshade',
    'fO3_d',
    'fO3_l',
    'fst',
    'fst_acc',
    'A_n',
    'A_n_canopy',
    'lai',
    'sai',
    'o3_ppb',
    'o3_ppb_i',
    'o3_nmol_m3',
    'canopy_height',
    'swp',
    'pody',
    'pod0',
    'ustar',
    'ra',
    'rb',
    'rsur',
    'rinc',
    'f_LS',
    'f_LA',
    'f_VPD',
    'npp',
    'npp_acc',
    'c_root',
    'c_stem',
    'c_leaf',
    'c_harv',
    'c_resv',
    'f_phen',
    'leaf_f_phen',
]

multip_output_fields = [
    'gsto',
    'gsto_l',
    'gsto_bulk',
    'rsto_c',
    'rsto_l',
    # "td_dd",
    # 'ts_c',
    # 'dvi',
    'PARsun',
    'PARshade',
    'fst',
    # 'A_n',
    'lai',
    'sai',
    'o3_ppb',
    'o3_ppb_i',
    'o3_nmol_m3',
    'canopy_height',
    'pody',
    'pod0',
    'swp',
    'ei',
    'et',
    'es',
    'f_phen',
    'leaf_f_phen',
    'f_light',
    'leaf_f_light',
    'f_temp',
    'f_VPD',
    'f_SW',
    'f_O3',
    'ustar',
    'u50',
    'ra',
    'rb',
    'rsur',
    'rinc',
]
