##########
Model Data
##########

The internal model data is split into Config, External Data and State

Config
--------------
This is the initial config of the model. It is set at the start and should be considered imutable within the model. It is normally a .yml or .json file


External_State
--------------
External Data(State) is inputed at the start of the model and is also considered imutable within the model. Weather data is an example of external data. It is often a csv file.

Model_State
--------------
State is modified as the model runs.

