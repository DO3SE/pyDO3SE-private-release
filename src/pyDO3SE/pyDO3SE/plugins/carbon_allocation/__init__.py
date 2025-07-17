"""Carbon allocation module.

The carbon allocation manages carbon pool allocation in the plant.
It convert net photosynthesis into accumulated daily carbon then distributes
these between each part of the plant.
From the carbon pools we can also then calculate the LAI and plant height.

Setup
-----

To setup the carbon allocation module a number of config changes need to be made.

 - set `config.carbon_allocation.use_carbon_allocation` to `true`
 - set the carbon allocation parameters in `config.carbon_allocation`
 - set `config.Land_Cover.height_method` to `carbon`
 - set `config.Land_Cover.LAI_method` to `carbon`

"""
