Multilayer Model
================

Download the original word file  :download:`MultiLayerModel.docx <MultiLayerModel.docx>`

.. container:: WordSection1

   Multilayer Model

   Upscaling from leaf level to canopy and overview of deposition
   methods. This documentation gives an overview of where the multilayer
   approach has been implemented into the DO3SE model. More detailed
   explanations can be found in specific documentation for each module.

    

    

   Contents

   `Methods.2 <#toc67397480>`__

   `Model Flow..2 <#toc67397481>`__

   `Glossary.3 <#toc67397482>`__

   `References.4 <#toc67397483>`__

    

   .. rubric:: General Overview
      :name: general-overview

   .. rubric:: Layer number
      :name: layer-number

   We count the layer numbers from the top. I.e the top layer is 0.

   .. _Toc67397481:

   \_Toc67397481

    

   .. rubric:: Upscaling
      :name: upscaling

   Upscaling in the model is performed after we calculate gsto at the
   leaf level. More details on this can be found in the surface
   resistance documentation.

    

   .. rubric:: Per layer environment
      :name: per-layer-environment

   Local environment is calculated for each layer including solar
   energy, gas concentration and moisture. The values at each layer take
   into account per layer resistances calculated by the model.

   .. rubric:: Solar Energy
      :name: solar-energy

   A multilayer model is setup to calculate the PAR at each layer in
   respect to shadowing and refraction from the canopy above. More
   details on this can be found in the Meteorological documentation.

   .. rubric:: Ozone Concentration
      :name: ozone-concentration

   Ozone deposition is calculated per layer in the Ozone Deposition
   module.

    

   .. rubric:: CO2 Concentration
      :name: co2-concentration

   WIP: This still needs implementing in the model.

    

   .. rubric:: Soil moisture
      :name: soil-moisture

   WIP this still needs implementing.

    

   .. rubric:: Phenology
      :name: phenology

   WIP: The theory for this is being developed.

   .. _Toc49267197:

   \_Toc49267197

   .. _Toc49348697:

    

   .. _Toc67397482:

   \_

   .. rubric:: Glossary
      :name: glossary

   +-----------------+-----------------+-----------------+-----------------+
   | **Term**        | **Description** | **Short name**  | ** **           |
   +-----------------+-----------------+-----------------+-----------------+
   |                 |                 |                 |                 |
   +-----------------+-----------------+-----------------+-----------------+
   |                 |                 |                 |                 |
   +-----------------+-----------------+-----------------+-----------------+
   |                 |                 |                 |                 |
   +-----------------+-----------------+-----------------+-----------------+
   |                 |                 |                 |                 |
   +-----------------+-----------------+-----------------+-----------------+

   .. rubric::  
      :name: section

   .. _Toc67397483:

   \_

   .. rubric:: References
      :name: references

   ·        

    

    

.. container::

   --------------
