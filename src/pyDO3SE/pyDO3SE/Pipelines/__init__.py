"""The Defaults module contains the pipelines for running the model.

Pipelines run a list of process on the state using Proflow the functional process running library.

The key pipelines are:

- :mod:`Default Processes <pyDO3SE.Pipelines.default_processes>` - Main model processes
- :mod:`Config Init Processes <pyDO3SE.Pipelines.config_init_processes>` - Processes for setting default config settings
- :mod:`External State Init Processes <pyDO3SE.Pipelines.es_init_processes>` - Processes for calculating external state values
- :mod:`State Init Processes <pyDO3SE.Pipelines.state_init_processes>` - Processes for calculating initial internal state values
- :mod:`Validation Processes <pyDO3SE.Pipelines.validation_processes>` - Processes for validating initial state

"""
