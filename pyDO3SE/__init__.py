"""pyDO3SE Package


This is the main pyDO3SE package.


How the model works
-------------------

This is a brief overview on how the model is linked together.



Model Entry Points
^^^^^^^^^^^^^^^^^^

There are a number of entrypoints to the model. Each of these leads to the :mod:`main <pyDO3SE.main>`
module

- :mod:`CLI <pyDO3SE.tools.cli>` - Entrypoint for running from the command line
- :mod:`main <pyDO3SE.main>` - Entrypoint for running the model

Model Inputs
^^^^^^^^^^^^^^^^^^^

**Configuration**:

The model configuration is a json file. For more info on how this is setup see the link below:

- :mod:`Config_Shape <pyDO3SE.Config.Config_Shape>`


**External Data**:

External data is provided in a .csv file where each row represents an hour.
For more information see the link below:

- :mod:`External State Shape <pyDO3SE.External_State.External_State_Shape>`

Model State
^^^^^^^^^^^

**Internal State**:

This is the state inside the model. Outputs are extracted from this state.
For more information see the link below:

- :mod:`Internal State Shape <pyDO3SE.Model_State.Model_State>`


Outputs
^^^^^^^

These are the output variables that can be logged and saved as a csv file and charts.
For more information see the link below:

- :mod:`Output Shape <pyDO3SE.Outputs.Output_Shape>`

Pipelines
^^^^^^^^^

The model runs a number of pipelines. These process the inputs then run the model processes.
For more information see the link below:

- :mod:`Pipelines <pyDO3SE.Pipelines>`


Constants
^^^^^^^^^

These are constants used by the model that are not altered by inputs configurations.
For more information see the link below:

- :mod:`Constants <pyDO3SE.Constants>`


"""
