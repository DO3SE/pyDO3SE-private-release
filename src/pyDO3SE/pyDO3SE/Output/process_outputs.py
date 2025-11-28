import warnings
import numpy as np
import csv
from datetime import datetime, date
import datetime as dt
import json
import os
import pickle
from pathlib import Path
from typing import Any, Callable, List, Tuple, Optional
import pandas as pd
from dataclasses import asdict
from data_helpers.cls_parsing import unpack
from data_helpers.encoders import AdvancedJsonEncoder
from data_helpers.list_helpers import flatten_list
from proflow.Objects.Process import Process

from do3se_phenology.plots import plot_phenology_from_config
from do3se_phenology.units import TimeTypes
from pyDO3SE.Model_State.model_state_loader import dump_state_to_file

from pyDO3SE.optional_dependencies import xarray as xr
from pyDO3SE.External_State.External_State_Config import FileTypes
from pyDO3SE.Analysis.util import (
    output_log_to_field_data,
    day_of_year_to_month,
    output_log_to_data_with_mm,
)
from pyDO3SE.Analysis.charts import annual_graph, monthly_diurnal_graph
from pyDO3SE.Output.Output_Shape import output_fields_map, output_fields
from pyDO3SE.External_State import External_State_Shape
from pyDO3SE.Model_State import Model_State_Shape
from pyDO3SE.Config.Config_Shape import Config_Shape
from pyDO3SE.util.error_handling import OutputError
from pyDO3SE.util.logger import Logger, get_git_revision_hash
from pyDO3SE.util.Objects import Field
from pyDO3SE.version import version as model_version
from .OutputConfig import OutputOptions

type OutputsFile = None


def dump_config_to_string(
    config: Config_Shape,
) -> str:
    return json.dumps(unpack(config), indent=4, cls=AdvancedJsonEncoder, sort_keys=True)


def dump_config_to_file_json(
    config: Config_Shape,
    target_path: Path,
) -> OutputsFile:
    with open(target_path, "w") as cf:
        json.dump(unpack(config), cf, indent=4, cls=AdvancedJsonEncoder, sort_keys=True)


def dump_config_to_file_binary(
    config: Config_Shape,
    target_path: Path,
) -> OutputsFile:
    """Pickle the model config.

    Parameters
    ----------
    state : Config_Shape
        Config Shape Object
    target_path : Path
        Location to save config

    """
    with open(target_path, "wb") as configfile:
        pickle.dump(config, configfile)


def dump_output_to_file_csv(
    output_data: List[dict[str, Any]],
    target_path: Path,
) -> OutputsFile:
    # TODO: Could filter fields here
    full_logs = pd.DataFrame(output_data)
    full_logs.to_csv(target_path)


def dump_output_to_file_netcdf(
    output_data: List[dict],
    target_path: Path,
) -> OutputsFile:
    # TODO: Could filter fields here
    full_logs = pd.DataFrame(output_data)
    ds = full_logs.to_xarray()
    # ds = xr.DataSet.from_dataframe(full_logs)
    xr.save_mfdataset([ds], paths=[target_path])


def merge_netcdf_grid_data(
    output_data: np.ndarray,  # data with shape(x, y ,t) where each element is a dict
    lat_data: List[List[float]],
    lon_data: List[List[float]],
    time_data: List[np.datetime64],
) -> xr.Dataset:
    variables = output_data[0][0][0].keys()
    X = len(output_data)
    Y = len(output_data[0])
    T = len(output_data[0][0])
    data_vars = {
        k: (
            ["x", "y", "time"],
            [[[output_data[y][x][t][k] for t in range(T)] for x in range(X)] for y in range(Y)],
        )
        for k in variables
    }

    # TODO: Could filter fields here
    return xr.Dataset(
        data_vars=data_vars,
        coords=dict(
            lon=(["x", "y"], lon_data),
            lat=(["x", "y"], lat_data),
            time=time_data,
            # reference_time=reference_time,
        ),
        attrs=dict(description="DO3SE outputs"),
    )


def dump_output_to_netcdf_grid(
    full_output_data: dict,
    output_shape: Tuple[int, int, int],
    lat_data: np.ndarray,
    lon_data: np.ndarray,
    time_data: np.ndarray,
    output_fields: List[str],
) -> xr.Dataset:
    """Dump output to netcdf file.


    Output data should be a dictionary where each value has shape (x, y, t)
    The lat and lon data should have shape (x, y)
    Time length should be t

    Parameters
    ----------
    full_output_data : dict
        dictionary of key: str and values: data with shape(x, y ,t) where each element is a dict
    output_shape : Tuple[int, int, int]
        Shape of output data (x, y, t)
    lat_data : List[List[float]]
        matrix of lon data for each grid point, shape(x, y)
    lon_data : List[List[float]]
        matrix of lon data for each grid point, shape(x, y)
    time_data : List[np.datetime64]
        matrix of time data with shape(t)
    output_fields : List[str]
        Fields to save

    Returns
    -------
    xr.Dataset
        xr.Dataset

    """
    variables = output_fields or list(full_output_data.keys())

    try:
        data_vars = {k: (["x", "y", "time"], full_output_data[k]) for k in variables}
        return xr.Dataset(
            data_vars=data_vars,
            coords=dict(
                lon=(["x", "y"], lon_data),
                lat=(["x", "y"], lat_data),
                time=time_data,
                # reference_time=reference_time,
            ),
            attrs=dict(description="DO3SE outputs"),
        )
    except KeyError as e:
        print(e)
        raise OutputError(
            f"Error saving to netcdf. Check output shape: {output_shape}, lon_shape: {lon_data.shape}, lat_shape: {lat_data.shape}"
        ) from e
    except ValueError as e:
        raise OutputError(
            f"Error saving to netcdf. Check output shape: {output_shape}, lon_shape: {lon_data.shape}, lat_shape: {lat_data.shape}"
        ) from e


def dump_output_to_file_netcdf_grid(
    full_output_data: dict,
    output_shape: Tuple[int, int, int],
    lat_data: np.ndarray,
    lon_data: np.ndarray,
    time_data: np.ndarray,
    output_fields: List[str],
    target_path: Path,
) -> OutputsFile:
    """Dump the output data to a NETCDF file.

    Output data should be a dictionary where each value has shape (x, y, t)
    The lat and lon data should have shape (x, y)
    Time length should be t


    Parameters
    ----------
    full_output_data : dict
        dictionary of key: str and values: data with shape(x, y ,t) where each element is a dict
    output_shape : Tuple[int, int, int]
        Shape of output data (x, y, t)
    lat_data : List[List[float]]
        matrix of lon data for each grid point, shape(x, y)
    lon_data : List[List[float]]
        matrix of lon data for each grid point, shape(x, y)
    time_data : List[np.datetime64]
        matrix of time data with shape(t)
    target_path : Path
        Location to save nc file

    Returns
    -------
    IO
        [description]
    """
    ds = dump_output_to_netcdf_grid(
        full_output_data=full_output_data,
        output_shape=output_shape,
        lat_data=lat_data,
        lon_data=lon_data,
        time_data=time_data,
        output_fields=output_fields,
    )
    xr.save_mfdataset([ds], paths=[target_path])


def dump_model_processes_info_to_string(
    model_processes: List[Process],
    detailed: bool = True,
    allow_errors: bool = False,
    silent: bool = False,
) -> list[str]:
    flattened_processes: List[Process] = [p for p in flatten_list(model_processes)]
    with warnings.catch_warnings():
        if silent:
            warnings.simplefilter("ignore")
        # p.human creates a lot of warnings about code that it cannot parse.
        flattened_processes_human = [p.human(allow_errors, silent) for p in flattened_processes]
    if detailed:
        flattened_process_comments = [
            "\n".join(
                filter(
                    lambda a: a,
                    [
                        # "________________________________________",
                        " ",
                        (
                            "# "
                            if "=====" in p.get("comment", "")
                            else "## "
                            if "==" in p.get("comment", "")
                            else "### "
                        )
                        + f"{p.get('comment', None) or p.get('func')}",
                        f"\tGroup: '{p.get('group', 'No group')}'",
                        "\n```json",
                        json.dumps(
                            {
                                "config_inputs": p.get("config_inputs", None),
                                "external_state_inputs": p.get("external_state_inputs", None),
                                "state_inputs": p.get("state_inputs", None),
                                "state_outputs": p.get("state_outputs", None),
                            },
                            indent=2,
                        ),
                        "\n```",
                    ],
                )
            )
            for p in flattened_processes_human
        ]
    else:
        flattened_process_comments: list[str] = [
            p.get("comment", "") or p.get("func", "") for p in flattened_processes_human
        ]
    return flattened_process_comments


def dump_model_processes_info_to_file(
    model_processes: List[Process],
    output_file: Path,
    detailed: bool = True,
    allow_errors: bool = False,
    silent: bool = False,
):
    output = dump_model_processes_info_to_string(
        model_processes,
        detailed,
        allow_errors,
        silent,
    )
    with open(output_file, "w") as f:
        f.write("\n".join(output))


def create_diurnal_plots(
    start_day: int,
    end_day: int,
    output_data,
    field: Field,
    runid: str | int,
    output_dir: Path,
    observed_diurnal_data=None,
):
    """Create a diurnal plot from model output data."""
    print(f"Creating diurnal plots for {field}")
    try:
        start_month = day_of_year_to_month(start_day)
        end_month = day_of_year_to_month(end_day)

        data_with_mm = output_log_to_data_with_mm(output_data, field)
        # match max/min to the annual plot
        ylim = (min([d[field.id] for d in data_with_mm]), max([d[field.id] for d in data_with_mm]))
        if observed_diurnal_data is not None and field.id in observed_diurnal_data.columns:
            ylim = (
                min([ylim[0], observed_diurnal_data[field.id].min()]),
                max([ylim[1], observed_diurnal_data[field.id].max()]),
            )

        for m in range(start_month, end_month + 1):
            df = pd.DataFrame([d for d in data_with_mm if d["mm"] == m])
            _odf = None
            if observed_diurnal_data is not None and field.id in observed_diurnal_data.columns:
                _odf = observed_diurnal_data.query(f"mm == {m}")[["mm", "hr", field.id]]
                # rename column to observed
                _odf.rename(columns={field.id: "observed"}, inplace=True)

            monthly_diurnal_graph(
                df,
                observed_diurnal_data=_odf,
                field=field,
                month=m,
                ylim=ylim,
                output_dir=output_dir,
                chart_id=f"{field.id}-{runid}",
            )
    except Exception as e:
        raise Exception(f"""Failed to create diurnal plot for field: \"{field}\" for run {runid}\
            start_day: {start_day}
            end_day: {end_day}
            """) from e

def plot_diurnal_charts(
    runid: str | int,
    fields_to_graph: List[str],
    observed_diurnal_path: Optional[Path],
    output_dir: Path,
    output_data,
    start_day: int,
    end_day: int,
):
    """Creates diurnal plots of the given fields.

    Parameters
    ----------
    runid : str
        [description]
    fields_to_graph : List[str]
        [description]
    observed_diurnal_path : Path
        [description]
    output_dir : Path
        [description]
    output_data : [type]
        [description]
    start_day : int
        [description]
    end_day : int
        [description]

    """
    # if observed diurnal data provided
    observed_diurnal_data = None
    if observed_diurnal_path:
        # read data
        observed_diurnal_data = pd.read_csv(observed_diurnal_path)

    for f in fields_to_graph:
        field = output_fields_map[f]
        # monthly diurnal plot
        try:
            create_diurnal_plots(
                start_day=start_day,
                end_day=end_day,
                output_data=output_data,
                field=field,
                runid=runid,
                observed_diurnal_data=observed_diurnal_data,
                output_dir=output_dir,
            )
        except Exception as e:
            print(e)
            warnings.warn(f"Failed to create diurnal plot for {field}")


def plot_annual_charts(
    runid: str | int,
    fields_to_graph: list[str],
    output_dir: Path,
    output_data: list[dict],
    start_day: int,
    end_day: int,
):
    """Create annual plots from all fields.

    Parameters
    ----------
    runid : str
        Run id unique identifier
    fields_to_graph : list[str]
        List of field ids to graph
    output_dir : Path
        Directory to save plots
    output_data : list[dict]
        List of output data dictionaries where each item is a time step
    start_day : int
        start day to plot
    end_day : int
        end day to plot
    """
    for f in fields_to_graph:
        field = output_fields_map[f]
        data = output_log_to_field_data(output_data, field)

        try:
            annual_graph(
                data,
                field,
                output_dir=output_dir,
                chart_id=f"{field.id}-{runid}",
                start_day=start_day,
                end_day=end_day,
                label_x_days=20,
                figsize=(20, 10),
            )
        except Exception as e:
            print(e)
            warnings.warn(f"Failed to create annual plot for {field}")


def generate_run_notes(
    runnotes: Optional[List[str]] = None,
    time_taken: Optional[date | str] = None,
    time_taken_setup: Optional[date | str] = None,
    config_version: Optional[str | int] = None,
    model_version: Optional[str] = None,
    errors: List[str] = [],
) -> List[str]:
    git_commit = get_git_revision_hash()
    from do3se_met.version import version as do3se_met_version
    from do3se_phenology.version import version as do3se_phenology_version
    from proflow.version import version as proflow_version

    pip_list_pretty = []
    try:
        # NOTE: Uses legacy api
        from pip._internal.utils.misc import get_installed_distributions  # type: ignore

        pip_list = get_installed_distributions()
        pip_list_pretty = [f"{p.project_name}\t{p._version}" for p in pip_list]
    except Exception as e:
        print(e)
        pass
    dependency_version_notes = []
    dependency_version_notes.append("===== dependency versions ======\n")

    dependency_version_notes.append(f"do3se_met: {do3se_met_version}\n")
    dependency_version_notes.append(f"do3se_phenology: {do3se_phenology_version}\n")
    dependency_version_notes.append(f"proflow: {proflow_version}\n")

    dependency_version_notes.append("=====================\n")

    return [
        str(runnotes),
        f"Date:\t {datetime.now().strftime('%d/%m/%Y %H:%M:%S')}",
        f"Model took:\t {time_taken}",
        f"Setup took:\t {time_taken_setup}",
        f"Using model version:\t {model_version}",
        f"Using git commit:\t {git_commit}",
        f"Using config version:\t {config_version}",
        str(dependency_version_notes),
        "Full dependencies: \n",
        str("\n".join(pip_list_pretty)),
        *["Errors: ", *[f"{e}" for e in errors]],
    ]


def export_failed_run_output(
    external_state: External_State_Shape,
    output_directory: Path,
    final_config: Config_Shape,
    initial_state: Model_State_Shape,
    model_processes: List[Process],
    additional_log_notes: Optional[List[str]] = None,
    logger: Callable[[str], None] = Logger(),
    verbose: bool = False,
    **kwargs,
):
    """Output results ofa failed pyDO3SE run.

    Parameters
    ----------
    external_state : External_State_Shape
        Processed external state that the model has run with
    output_directory : Path
        Location to save all outputs
    final_config : Config_Shape
        Processed config
    initial_state : Model_State_Shape
        Initital model state
    model_processes : List[Process]
        Model processes that are ran
    additional_log_notes : List[str], optional
        Additional notes to output, by default None
    logger: Callable[[str], None], optional
        Log function, by default print
    **kwargs: OutputOptions
        Options on what to output

    """

    options = OutputOptions(**kwargs)

    os.makedirs(output_directory, exist_ok=True)
    try:
        if options.save_external_processed_data:
            logger("Saving external processed data")
            row_count = len(external_state.dd)  # type: ignore
            external_state_out = pd.DataFrame(
                {
                    k: v[0:row_count] if v is not None else None
                    for k, v in asdict(external_state).items()
                }
            )
            external_state_out.to_csv(f"{output_directory}/external_data.csv")

        if options.save_processed_config:
            # Save out processed config
            logger("Saving processed config")
            dump_config_to_file_json(
                final_config, Path(f"{output_directory}/processed_config.json")
            )

        if options.save_initial_state:
            # Save out processed initial state
            logger("Saving initial state")
            dump_state_to_file(initial_state, Path(f"{output_directory}/initial_state.json"))

        if options.save_model_processes:
            # Save model processes info to file
            logger("save model processes")
            model_processes_out = (
                model_processes
                if type(model_processes) is type([])
                else model_processes.values() # type: ignore
                if type(model_processes) is type({})
                else None
            )
            if model_processes_out is None:
                raise ValueError(f"Invalid model processes type: {type(model_processes)}")

            dump_model_processes_info_to_file(
                model_processes_out,
                Path(f"{output_directory}/model_processes.md"),
                detailed=False,
                allow_errors=True,
                silent=not verbose,
            )
        if options.save_model_processes_detailed:
            logger("save model processes detailed")
            model_processes_out = (
                model_processes
                if type(model_processes) is type([])
                else model_processes.values() # type: ignore
                if type(model_processes) is type({})
                else None
            )
            if model_processes_out is None:
                raise ValueError(f"Invalid model processes type: {type(model_processes)}")

            dump_model_processes_info_to_file(
                model_processes_out,
                Path(f"{output_directory}/model_processes_detailed.md"),
                allow_errors=True,
                silent=not verbose,
            )

    except Exception as e:
        warnings.warn(
            f"There was an error processing the outputs! \nError\n========\n{e}\n========\n"
        )
        raise e from e
    finally:
        # Save notes
        if options.save_logs:
            log_notes = generate_run_notes(
                additional_log_notes,
                None,
                None,
                final_config.VERSION,
                model_version,
            )

            if log_notes:
                with open(f"{output_directory}/notes.log", "w") as notesfile:
                    notesfile.write("\n".join(log_notes))


def export_output(
    output_data: List[dict[str, Any]],
    final_state: Model_State_Shape,
    external_state: External_State_Shape,
    output_directory: Path,
    final_config: Config_Shape,
    initial_state: Model_State_Shape,
    model_processes: List[Process],
    output_filename: str = "pyDO3SE_output.csv",
    fields_to_graph: List[str] = [],
    observed_diurnal_path: Optional[Path] = None,
    time_taken: Optional[dt.date] = None,
    time_taken_setup: Optional[dt.date] = None,
    runid: Optional[str | int] = None,
    additional_log_notes: Optional[List[str]] = None,
    dump_output_to_file: Callable[
        [List[dict[str, Any]], Path], OutputsFile
    ] = dump_output_to_file_csv,
    logger: Callable[[str], None] = Logger(),
    verbose: bool = False,
    **kwargs,
):
    """Output results of pyDO3SE.

    Parameters
    ----------
    output_data : List[dict]
        Output data from pyDO3SE run_model
    final_state : Model_State_Shape
        The final model state
    external_state : External_State_Shape
        Processed external state that the model has run with
    output_directory : Path
        Location to save all outputs
    final_config : Config_Shape
        Processed config
    initial_state : Model_State_Shape
        Initital model state
    model_processes : List[Process]
        Model processes that are ran
    output_filename : str, optional
        File name for saved model hourly outputs, by default "pyDO3SE_output.csv"
    fields_to_graph : List[str], optional
        Fields to plot(Should match fields in Output_Shape), by default []
    observed_diurnal_path : Path, optional
        Path to observed diurnal data, by default None
    time_taken: Date
        Time taken to run model
    runid : str, optional
        run id, by default None
    additional_log_notes : List[str], optional
        Additional notes to output, by default None
    dump_output_to_file : Callable[[List[dict], Path], IO], optional
        function to dump output to file, by default dump_output_to_file_csv
    logger: Callable[[str], None], optional
        Log function, by default print
    **kwargs: OutputOptions
        Options on what to output

    """
    options = OutputOptions(**kwargs)

    os.makedirs(output_directory, exist_ok=True)
    try:
        if options.save_hourly_output_data:
            output_filename = "pyDO3SE_output.csv" if output_filename is None else output_filename
            output_location = Path(f"{output_directory}/{output_filename}")
            logger(f"Saving hourly output to {output_location}")
            dump_output_to_file(output_data, output_location)

        if options.save_external_processed_data:
            logger("Saving external processed data")
            row_count = len(external_state.dd or [])
            external_state_out = pd.DataFrame(
                {
                    k: v[0:row_count] if v is not None else None
                    for k, v in asdict(external_state).items()
                }
            )
            external_state_out.to_csv(f"{output_directory}/external_data.csv", index=False)

        if options.save_processed_config:
            # Save out processed config
            logger("Saving processed config")
            dump_config_to_file_json(
                final_config, Path(f"{output_directory}/processed_config.json")
            )

        if options.save_initial_state:
            # Save out processed initial state
            logger("Saving initial state")
            dump_state_to_file(initial_state, Path(f"{output_directory}/initial_state.json"))

        if options.save_model_processes:
            # Save model processes info to file
            logger("save model processes")
            model_processes_out = (
                model_processes
                if type(model_processes) is type([])
                else model_processes.values() # type: ignore
                if type(model_processes) is type({})
                else None
            )
            if model_processes_out is None:
                raise ValueError(f"Invalid model processes type: {type(model_processes)}")
            dump_model_processes_info_to_file(
                model_processes_out,
                Path(f"{output_directory}/model_processes.md"),
                detailed=False,
                allow_errors=True,
                silent=not verbose,
            )
        if options.save_model_processes_detailed:
            logger("save model processes detailed")
            model_processes_out = (
                model_processes
                if type(model_processes) is type([])
                else model_processes.values() # type: ignore
                if type(model_processes) is type({})
                else None
            )
            if model_processes_out is None:
                raise ValueError(f"Invalid model processes type: {type(model_processes)}")

            dump_model_processes_info_to_file(
                model_processes_out,
                Path(f"{output_directory}/model_processes_detailed.md"),
                allow_errors=True,
                silent=not verbose,
            )

        if options.save_output_heading_info:
            df = pd.DataFrame(
                [field._asdict() for field in output_fields if field.id in output_data[0].keys()]
            )
            df.to_csv(f"{output_directory}/output_fields_info.csv", index=False)

        if (options.plot_annual_charts or options.plot_diurnal_charts) and fields_to_graph:
            # Validate fields to graph
            for f in fields_to_graph:
                try:
                    output_fields_map[f]
                except KeyError:
                    warnings.warn(
                        f"{f} not in output fields map. Run cli command `pyDO3SE_cli available-outputs"
                    )
                    continue

        if options.plot_annual_charts and fields_to_graph and len(fields_to_graph):
            logger(f"Plotting annual charts. Fields: {fields_to_graph}")
            logger(f"Found the following fields to plot: {fields_to_graph}")
            plot_dir_annual = Path(f"{output_directory}/annual_plots")
            os.makedirs(plot_dir_annual, exist_ok=True)
            start_day = final_config.Location.start_day
            assert external_state.dd is not None, "external_state.dd is None"
            start_day = start_day if start_day is not None else int(external_state.dd[0])
            end_day = final_config.Location.end_day
            end_day = end_day if end_day is not None else int(external_state.dd[-1])
            plot_annual_charts(
                runid=runid or 0,
                fields_to_graph=fields_to_graph,
                output_dir=plot_dir_annual,
                output_data=output_data,
                start_day=start_day,
                end_day=end_day,
            )

        if options.plot_diurnal_charts and fields_to_graph and len(fields_to_graph):
            logger(f"Plotting diurnal charts. Fields: {fields_to_graph}")
            logger(f"Found the following fields to plot: {fields_to_graph}")
            diurnal_plot_dir = Path(f"{output_directory}/diurnal_plots")
            os.makedirs(diurnal_plot_dir, exist_ok=True)

            start_day = final_config.Location.start_day
            assert external_state.dd is not None, "external_state.dd is None"
            start_day = start_day if start_day is not None else int(external_state.dd[0])
            end_day = final_config.Location.end_day
            end_day = end_day if end_day is not None else int(external_state.dd[-1])
            plot_diurnal_charts(
                runid=runid or 0,
                fields_to_graph=fields_to_graph,
                observed_diurnal_path=observed_diurnal_path,
                output_dir=diurnal_plot_dir,
                output_data=output_data,
                start_day=start_day,
                end_day=end_day,
            )

        if options.plot_phenology_charts:
            logger("Plotting phenology charts")
            phenology_plot_dir = Path(f"{output_directory}/phenology_plots")
            os.makedirs(phenology_plot_dir, exist_ok=True)
            day_count = int(
                (final_config.Location.end_day or 365) - (final_config.Location.start_day or 0)
                if final_config.Location.start_day is not None
                and final_config.Location.end_day is not None
                else 365
            )
            plot_phenology_from_config(
                final_config.Land_Cover.parameters[0].phenology,
                final_config.Land_Cover.phenology_options,
                nP=final_config.Land_Cover.nP,
                output_location=Path(f"{phenology_plot_dir}/phenology.png"),
                day_count=day_count,
                plot_dd=final_config.Land_Cover.phenology_options.time_type == TimeTypes.JULIAN_DAY,
                plot_td=final_config.Land_Cover.phenology_options.time_type
                == TimeTypes.THERMAL_TIME,
                plot_f_phen=True,
                plot_lengths=True,
                plot_carbon=final_config.Land_Cover.phenology_options.time_type
                == TimeTypes.THERMAL_TIME,
                plot_growing=final_config.Land_Cover.phenology_options.time_type
                == TimeTypes.THERMAL_TIME,
            )
        if options.save_final_state:
            dump_state_to_file(final_state, Path(f"{output_directory}/final_state.json"))
    except Exception as e:
        warnings.warn(
            f"There was an error processing the outputs! \nError\n========\n{e}\n========\n"
        )
        raise e from e
    finally:
        # Save notes
        if options.save_logs:
            log_notes = generate_run_notes(
                additional_log_notes,
                time_taken,
                time_taken_setup,
                final_config.VERSION,
                model_version,
            )

            if log_notes:
                with open(f"{output_directory}/notes.log", "w") as notesfile:
                    notesfile.write("\n".join(log_notes))


def export_single_hour_output(
    output_data_dir: Path,
    output_logs: List[dict],
    final_state: Model_State_Shape,
    output_type: FileTypes = FileTypes.CSV,
):
    os.makedirs(output_data_dir, exist_ok=True)
    if output_type == FileTypes.NETCDF:
        output_data_file = Path(f"{output_data_dir}/output_data.nc")
        dump_output_to_file_netcdf(output_logs, output_data_file)
    else:
        output_data_file = Path(f"{output_data_dir}/output_data.csv")
        dump_output_to_file_csv(output_logs, output_data_file)
    output_state_file = Path(f"{output_data_dir}/output_state.nc")
    dump_state_to_file(final_state, output_state_file)


def get_headers(file):
    with open(file) as f:
        return f.readline().replace("\n", "").replace("\r\n", "").split(",")


def extract_state_val(line_data, field):
    return line_data[field]


def get_final_state_values(input_file_directory, field):
    files = [f for f in os.listdir(input_file_directory) if f.split(".")[-1] == "csv"]
    headers = get_headers(input_file_directory + "/" + files[0])
    # print(headers)
    last_lines = get_last_line_of_files(input_file_directory, files, headers)
    data = [line_data.get(field, None) for line_data in last_lines]
    file_names = [f.split(".")[0] for f in files]
    return dict(zip(file_names, data))


def get_last_line_in_file(filename, headers):
    try:
        with open(filename, "rb") as f:
            f.seek(-2, os.SEEK_END)
            while f.read(1) != b"\n":
                f.seek(-2, os.SEEK_CUR)
            last_line = list(csv.reader([f.readline().decode()]))[0]
            last_line_data = {h: d for h, d in zip(headers, last_line)}
            return last_line_data
    except OSError:
        Warning(f"Os Error while getting last line in file {filename}")
        return {h: None for h in headers}


def get_last_line_of_files(dir, files, headers):
    last_lines = [get_last_line_in_file(dir + "/" + f, headers) for f in files]
    return last_lines


def get_and_save_final_state_data(input_file_directory, output_directory, field_name):
    data = get_final_state_values(input_file_directory, field_name)
    name = os.path.basename(input_file_directory)
    print(f"Getting pod for{name}")
    os.makedirs(output_directory, exist_ok=True)
    with open(f"{output_directory}/{name}.json", "w") as out_file:
        json.dump(data, out_file)
    with open(f"{output_directory}/{name}.csv", "w", newline="") as out_file:
        spamwriter = csv.writer(out_file, delimiter=",", quotechar="|", quoting=csv.QUOTE_MINIMAL)
        for k, v in data.items():
            spamwriter.writerow([k, v])


def get_and_save_final_state_data_multi_dir(input_file_directory, output_directory, field_name):
    dirs = [
        d for d in os.listdir(input_file_directory) if os.path.isdir(f"{input_file_directory}/{d}")
    ]
    data = [get_final_state_values(f"{input_file_directory}/{d}", field_name) for d in dirs]
    print(dirs)
    for dir, data_group in zip(dirs, data):
        name = os.path.basename(dir)
        print(f"Getting pod for{name}")
        os.makedirs(output_directory, exist_ok=True)
        with open(f"{output_directory}/{name}.json", "w") as out_file:
            json.dump(data_group, out_file)
        with open(f"{output_directory}/{name}.csv", "w", newline="") as out_file:
            spamwriter = csv.writer(
                out_file, delimiter=",", quotechar="|", quoting=csv.QUOTE_MINIMAL
            )
            for k, v in data_group.items():
                spamwriter.writerow([k, v])


def extract_final_results(
    input_directory: Path,
    output_directory: Path,
    field_name: str,
):
    get_and_save_final_state_data(input_directory, output_directory, field_name)
