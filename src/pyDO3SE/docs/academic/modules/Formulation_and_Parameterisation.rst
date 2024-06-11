Model Inputs and Parameterisation
=================================

Download the original word file  :download:`Formulation_and_Parameterisation.docx <Formulation_and_Parameterisation.docx>`

.. container:: WordSection1

   .. _Toc36708829:

   Model Inputs and Parameterisation

   .. _Toc36708830:

   Describe input requirements; describe landcover and species parameterisations (splitting out between multiplicative and photosynthetic as appropriate

   Contents

   `DO3SE input data requirements.1 <#toc49265269>`__

    

   .. rubric::  
      :name: section

   .. _Toc49265269:

    

   .. rubric:: Model Data
      :name: model-data

   Model data is split into the following categories:

   .. rubric:: Configuration
      :name: configuration

   This is derived from the configuration file provided by the user.

    

   .. rubric:: External State
      :name: external-state

   This is derived from the data provided by the user. The data is
   cleaned and processed ready for use by the DO3SE model (See Input
   requirements below).

    

   .. rubric:: Model State
      :name: model-state

   | Model state is the internal state of the model. It defines the
     current state of growth, photosynthesis etc. The output is derived
     from the model state.

   .. rubric:: DO3SE input data requirements(CHECK THESE)
      :name: do3se-input-data-requirementscheck-these

    

    

   +-----------------------+-----------------------+-----------------------+
   | **Input data**        | **Height (m)**        | **Units**             |
   +-----------------------+-----------------------+-----------------------+
   | Ozone concentration   | Reference height and  | ppb                   |
   | (O3)                  | surface canopy height |                       |
   +-----------------------+-----------------------+-----------------------+
   | Horizontal windspeed  | Reference height and  | m s\ -1               |
   | (u)                   | surface canopy height |                       |
   +-----------------------+-----------------------+-----------------------+
   | Turbulent stress (t)  | surface               | kg m\ -1 s\ -2        |
   +-----------------------+-----------------------+-----------------------+
   | Heat flux density     | surface               | W m\ -2               |
   | (Hd)                  |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Air pressure (p)      | surface               | N m\ -2               |
   +-----------------------+-----------------------+-----------------------+
   | Global radiation      | Top of canopy         | W m\ -2               |
   |                       |                       |                       |
   | **OR**                |                       |                       |
   |                       |                       |                       |
   | Photon Photosynthetic | Top of canopy         | mmol m\ -2 s\ -1      |
   | Flux Density (PPFD)   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Air temperature       | surface               | K                     |
   | (Tair) / Leaf         |                       |                       |
   | temperature (Tleaf)   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Vapour Pressure       | Surface               | kPa                   |
   | Deficit (VPD)         |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Precipitation (Pr)    | Ground                | mm                    |
   +-----------------------+-----------------------+-----------------------+
   | **Site/Grid specific  | **Character**         | **Ideal units**       |
   | variables**           |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Latitude and          | -                     | :sup:`o`\ , ‘         |
   | Longitude (lat &      |                       |                       |
   | long)                 |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Elevation (e)         | -                     | m a.s.l.              |
   +-----------------------+-----------------------+-----------------------+
   | Target canopy height  | -                     | m                     |
   | (tgt)                 |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Soil texture          | coarse / medium /     | -                     |
   |                       | fine                  |                       |
   +-----------------------+-----------------------+-----------------------+

    

   N.B. Input data greyed out are only needed to estimate atmospheric
   stability (and are not required by the DO\ 3\ SE model interface
   since here we assume a stable atmosphere); the DO\ 3\ SE model will
   default to a medium soil texture if this information is not provided

    

   Update table; make it very clear which parameters are required for
   which version of model!

    
