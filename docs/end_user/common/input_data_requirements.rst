=============================
Model Input Data Requirements
=============================

Input Data
==========

 - Input data should currently be provided with each row representing an hour and row 1 representing hour 0 day 1 (i.e Jan 1st 12am)

Old Documentation
=================


.. rubric:: 2.1.                  DOSE input data requirements
   :name: dose-input-data-requirements

 

 

+--------------------------+--------------------------+--------------------------+
| **Input data**           | **Height (m)**           | **Units**                |
+--------------------------+--------------------------+--------------------------+
| Ozone concentration (O)  | Reference height and     | ppb                      |
|                          | surface canopy height    |                          |
+--------------------------+--------------------------+--------------------------+
| Horizontal windspeed (u) | Reference height and     | m s\ :sup:`-1`           |
|                          | surface canopy height    |                          |
+--------------------------+--------------------------+--------------------------+
| Turbulent stress (t)     | surface                  | kg m\ :sup:`-1`          |
|                          |                          | s\ :sup:`-2`             |
+--------------------------+--------------------------+--------------------------+
| Heat flux density (Hd)   | surface                  | W m\ :sup:`-2`           |
+--------------------------+--------------------------+--------------------------+
| Air pressure (p)         | surface                  | N m\ :sup:`-2`           |
+--------------------------+--------------------------+--------------------------+
| Global radiation         | Top of canopy            | W m\ :sup:`-2`           |
|                          |                          |                          |
| **OR**                   |                          |                          |
|                          |                          |                          |
| Photon Photosynthetic    | Top of canopy            | mmol m\ :sup:`-2`        |
| Flux Density (PPFD)      |                          | s\ :sup:`-1`             |
+--------------------------+--------------------------+--------------------------+
| Air temperature (T) /    | surface                  | K                        |
| Leaf temperature (T)     |                          |                          |
+--------------------------+--------------------------+--------------------------+
| Vapour Pressure Deficit  | Surface                  | kPa                      |
| (VPD)                    |                          |                          |
+--------------------------+--------------------------+--------------------------+
| Precipitation (Pr)       | Ground                   | mm                       |
+--------------------------+--------------------------+--------------------------+
| **Site/Grid specific     | **Character**            | **Ideal units**          |
| variables**              |                          |                          |
+--------------------------+--------------------------+--------------------------+
| Latitude and Longitude   | -                        | :sup:`o`\  , ‘           |
| (lat & long)             |                          |                          |
+--------------------------+--------------------------+--------------------------+
| Elevation (e)            | -                        | m a.s.l.                 |
+--------------------------+--------------------------+--------------------------+
| Target canopy height     | -                        | m                        |
| (tgt)                    |                          |                          |
+--------------------------+--------------------------+--------------------------+
| Soil texture             | coarse / medium / fine   | -                        |
+--------------------------+--------------------------+--------------------------+

 

N.B. Input data greyed out are only needed to estimate atmospheric
stability (and are not required by the DOSE model interface since here
we assume a stable atmosphere); the DOSE model will default to a medium
soil texture if this information is not provided

 

Update table; make it very clear which parameters are required for which
version of model!

 

The DOSE model’s total dry deposition scheme is shown in Fig 1

.. _msoanchor-1:

[b1]

 . The model assumes the key resistances to O deposition are the
aerodynamic resistance (Ra), the quasi-laminar sub-layer resistance (Rb)
above the canopy, and the surface resistance (Rsur). Rsur comprises two
resistance paths in series; the stomatal and non-stomatal resistance.
The latter represents a) within canopy aerodynamic resistance (Rinc) and
subsequent b) soil resistance to decomposition at the soil surface (Rgs)
which encompasses features such as leaf litter and ground vegetation
under forest canopies, as well as c) resistance to adsorption to the
external plant parts (Rext) including cuticle, bark etc… a comment on
planned revision of estimation of Rext here?
|image1|

 

As such, the loss rate of ozone to the ground surface, within a volume
unit area and height Δ\ *z* is given by the product of the deposition
velocity (*Vg*) at height *z* and the concentration (C) at that height
as described in eq. 1.

 

|image2|

.. _ref393790837:

1

 

Where *Vg* is estimated as given in eq. 2

|image3|

.. _ref393790933:

2

 

.. _toc36708831:

-
