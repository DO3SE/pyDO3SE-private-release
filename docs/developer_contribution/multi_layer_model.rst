=================
Multi Layer model
=================

The DO3SE model can run with multiple canopy layers, components, and leaf populations.
It can also outputs values at custom heights.


---------
Key words
---------
The following keyword/variable name suffixes are used to identify the layer that the
variable or heading refers to.

 - "nLC" - Number of Land Cover types
 - "nL" - Number of Layers
 - "nP" - Number of leaf populations
 - "iLC" - Land Cover type index
 - "iL" - Layer index
 - "iCH" - Custom height index
 - "iP" - Leaf population index


------------------
Layer count limits
------------------

The architecture of the model state means that we have to set an internal limit
to the number of layers created in the model state. You can override this value
by defining the following environment variables:

 - MAX_NUM_OF_CANOPY_LAYERS: int
 - MAX_NUM_OF_CUSTOM_LAYERS: int
 - MAX_NUM_OF_CANOPY_COMPONENTS: int
 - MAX_NUM_OF_LEAF_POPULATIONS: int

 To set an environment variable in linux use the following command:
 `export MAX_NUM_OF_CANOPY_LAYERS=10`