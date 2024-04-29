=============
Soil Moisture
=============

Soil moisture in the pyDO3SE model is setup as below:


Config Requirements
===================

Here we explain the configuration Requirements.


.. literalinclude:: ../../../examples/soil_moisture/configs/penman_monteith.json
   :language: json
   :linenos:
   :caption: Penman Monteith Config Overrides


External Data Requirements
==========================

Here we explain the external data Requirements.


+--------------------------+--------------------------+--------------------------+
| **Id**                   | **Description**          | **Units**                |
+--------------------------+--------------------------+--------------------------+
| precip                   | Precipitation            | mm/hour                  |
+--------------------------+--------------------------+--------------------------+


Model Code
==========

Please find below links to the relevant model code:

- Soil Moisture module code: :mod:`pyDO3SE.plugins.soil_moisture`