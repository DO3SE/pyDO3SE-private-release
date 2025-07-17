from dataclasses import dataclass
from typing import List


@dataclass
class OutputConfig:
    """The output config defines what the pyDO3Se model outputs.

    The Fields should be a list of string using the field ids from here: :data:`pyDO3SE.Output.Output_Shape.output_fields`.

    """

    fields: List[str] = None
    log_multilayer: bool = False
    include_spacers: bool = True


@dataclass
class OutputOptions:
    """Model output options.

    """

    save_hourly_output_data: bool = True
    save_external_processed_data: bool = True
    plot_annual_charts: bool = True
    plot_diurnal_charts: bool = True
    save_processed_config: bool = True
    save_initial_state: bool = True
    save_final_state: bool = True
    save_model_processes: bool = True
    save_model_processes_detailed: bool = True
    save_output_heading_info: bool = True
    save_logs: bool = True


def output_results_only_options():
    return OutputOptions(
        save_hourly_output_data=True,
        save_external_processed_data=False,
        plot_annual_charts=False,
        plot_diurnal_charts=False,
        save_processed_config=False,
        save_initial_state=False,
        save_final_state=False,
        save_model_processes=False,
        save_model_processes_detailed=False,
        save_output_heading_info=False,
        save_logs=False,
    )


def output_options_none():
    return OutputOptions(
        save_hourly_output_data=False,
        save_external_processed_data=False,
        plot_annual_charts=False,
        plot_diurnal_charts=False,
        save_processed_config=False,
        save_initial_state=False,
        save_final_state=False,
        save_model_processes=False,
        save_model_processes_detailed=False,
        save_output_heading_info=False,
        save_logs=False,
    )
