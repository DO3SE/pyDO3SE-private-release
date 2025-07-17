===========
Legacy Info
===========

Please find below information on the legacy of the DO3SE model code


Existing and Past Versions
==========================

**Table 1.** Evolution of the DOSE model (update)

 

+--------------------+--------------------+--------------------+--------------------+
| **Version**        | **Key modules /    | **Parameterisation | **Reference**      |
|                    | additions**        | **                 |                    |
+--------------------+--------------------+--------------------+--------------------+
| DOSE Ver1          | Ra, Rb, gsto, LAI  | X forest species   | Emberson et al     |
|                    | scaling            |                    | (2000)             |
|                    |                    | X crops species    |                    |
|                    |                    |                    |                    |
|                    |                    | X grassland        |                    |
|                    |                    | species            |                    |
+--------------------+--------------------+--------------------+--------------------+
| DO3SE Ver2         |                    |                    | EMEP Unified Model |
|                    |                    |                    | (2003)             |
|                    |                    |                    |                    |
|                    |                    |                    |                    |
+--------------------+--------------------+--------------------+--------------------+
| DO3SE Ver3         |                    |                    | EMEP status report |
|                    |                    |                    | (2006)             |
+--------------------+--------------------+--------------------+--------------------+
| DOSE Ver2          | Latitude derived   |                    | Emberson et al     |
|                    | GS, SMD (VPD       |                    | (2007)             |
|                    | driven)            |                    |                    |
+--------------------+--------------------+--------------------+--------------------+
| DOSE Ver3.1 &      | SMD (Penman        |                    | Büker et al.       |
| interfaced version | Monteith)          |                    | (2012)             |
+--------------------+--------------------+--------------------+--------------------+
|                    |                    |                    |                    |
+--------------------+--------------------+--------------------+--------------------+
| DOSE Ver4          | Multi-layer &      | Grassland          | Emberson et al.    |
|                    | multi-PFT model    |                    | (in prep)          |
+--------------------+--------------------+--------------------+--------------------+
|                    |                    |                    |                    |
+--------------------+--------------------+--------------------+--------------------+
|                    |                    |                    |                    |
+--------------------+--------------------+--------------------+--------------------+

 


DO3SE UI (Fortran & Python GUI)
-------------------------------
https://github.com/SEI-DO3SE/DO3SE-UI
This is the original DO3SE model written in Fortran that has a python GUI using wxWidgets.
This model can be compiled and run on any OS.

DO3SE model (Fortran & Python CLI)
----------------------------------
https://github.com/SEI-DO3SE/DO3SE-model
The DO3SE model is an updated version of the DO3SE UI model also written in Fortran.
It has a python wrapper that exposes the model to python environments such as CLI and Jupter Notebooks
It includes an initial attempt at upscalling and multi layer modelling.

pyDO3SE (Python)
----------------
https://github.com/SEI-DO3SE/pyDO3SE
This is the most recent update to the DO3SE model and is a complete rewrite in python.
