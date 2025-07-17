"""do3se phenology

The DO3SE phenology module estimates key plant life stages as thermal time intervals.
It aims to align many different phenology terms into a single thermal time function.
This will allow multiple phenology parameters used by different parts of the model
to all be defined simply from the thermal time data and a small number of reference points.

The three key reference points are **Sowing date**, **Mid anthesis** and **Season end Date**


Glossary
----------
.. csv-table:: Phenology Variables
   :file: phenology_variables.csv
   :widths: 10,10,10,10,50,10
   :header-rows: 1

"""


import os
import sys
sys.path.append(os.path.abspath(__file__))
