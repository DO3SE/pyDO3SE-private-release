================
Ozone Deposition
================

This documentation relates to running the model purely to tranfer ozone from
measured height to an arbitary canopy height.

WIP

Requirements for ozone deposition runs

Model Processes
===============

 - (a) Reset and Calculate the resistance model for O3 over the target canopy
 - (b) Calculate Top Layer Canopy ozone
 - (c) Calculate estimated windspeed at canopy
 - (d) Calculate ustar_ref and monin obukhov length
        or
        Set ustar_ref to external input
        Calculate ustar_ref using simple method
 - (e) calculate measured wind canopy displacement parameters
 - (f) Calculate gsto - multiplicative
 - (g) Calculate f_light


Config Requirements
===================

a)
I(config.Location.z_O3, as_='z_O3'),
I(config.Location.O3_d, as_='O3_d'),
I(config.Location.O3_z0, as_='O3_z0'),
I(config.Location.izr, as_='izr'),
I(config.resistance.ra_calc_method, as_='ra_method'),

b)
I(config.Land_Cover.nL, as_='nL'),
I(config.Land_Cover.nLC, as_='nLC'),
I(config.Location.Rsoil, as_='Rsoil'),
I(ra_calc_method, as_='ra_calc_method'),
I(rsur_calc_method, as_='rsur_calc_method'),
I(rext_calc_method, as_='rext_calc_method'),

c)
I(config.Location.OTC, as_='o_top_chamber'),
I(config.Location.z_u, as_='z_u'),
I(config.Location.u_z0, as_='u_z0'),
I(config.Location.u_d, as_='u_d'),
I(config.Location.izr, as_='izr'),

d)
I(config.Location.z_u, as_='z_u'),
I(config.Location.u_d, as_='u_d'),
I(config.Location.u_z0, as_='u_z0'),

e)
I(config.Land_Cover.land_cover_type == LandCoverType.FOREST, as_="is_forest"),

f)
I(config.Land_Cover.parameters[iLC].gsto.VPD_crit, as_='VPD_crit'),
I(config.Land_Cover.parameters[iLC].multip_gsto.gmax, as_='gmax'),
I(config.Land_Cover.parameters[iLC].multip_gsto.gmorph, as_='gmorph'),
I(config.Land_Cover.parameters[iLC].gsto.fmin, as_='fmin'),

g)
I(config.Land_Cover.parameters[iLC].multip_gsto.f_lightfac, as_='f_lightfac'),

External_State
--------------

a)
I(lget(e_state.O3, row_index), as_='O3'),

b)
I(e_state.Ts_C[row_index], as_='Ts_C'),
    I(lget(e_state.Hd, row_index), as_='Hd'),
    I(e_state.snow_depth[row_index], as_='snow_depth'),
    I(lget(e_state.P, row_index), as_='P'),
] if ra_calc_method == 'heat_flux' else [],

c)
I(lget(e_state.u, row_index), as_="u"),

d)
I(lget(e_state.Hd, row_index), as_='Hd'),
I(lget(e_state.P, row_index), as_='P'),
I(e_state.Ts_C[row_index] + T0, as_='Tk'),
I(e_state.ustar_ref[row_index], as_='ustar_ref_in'),

e)
NONE

f)
I(e_state.VPD_dd[row_index], as_='VPD_dd')

g)
I(lget(e_state.sinB, row_index), as_='sinB'),
I(lget(e_state.PAR, row_index), as_='PAR'),

combined:
^^^^^^^^^
O3
Ts_C
u
VPD_dd

optional:
^^^^^^^^^
ustar_ref
Hd
snow_depth
P
PAR

State
-----

a)
[LINKED - b]    I(state.canopy.rmodel_O3.Rsur[top_layer_index], as_='Rsur_top_layer'),
[LINKED - b]    I(state.canopy.rmodel_O3.Rb, as_='Rb_top_layer'),
[LINKED - b]    I(state.canopy.rmodel_O3.Ra, as_='Ra_top_layer'),
[LINKED - c]    I(state.met.ustar, as_='ustar'),
[LINKED - d]    I(state.met.ustar_ref, as_='ustar_ref'),  # TODO: Should be ozone ustar_ref
[LINKED - d]    I(state.met.L, as_='L'),  # TODO: Should be ozone L
[LINKED - e]    I(state.canopy.d, as_='d'),
[LINKED - e]    I(state.canopy.z0, as_='z0'),

b)
[LINKED - c]    I(state.met.ustar, as_='ustar'),
[INPUT]         I(state.canopy.canopy_height, as_='canopy_height'), # from config
[INPUT]         I([[state.canopy_layer_component[iL][iLC].SAI for iLC in range(nLC)]
    for iL in range(nL)], as_='SAI_values'),
[INPUT]         I([[state.canopy_layer_component[iL][iLC].LAI for iLC in range(nLC)]
    for iL in range(nL)], as_='LAI_values'),
I([[state.canopy_layer_component[iL][iLC].mean_gsto for iLC in range(nLC)]
    for iL in range(nL)], as_='mean_gsto_values'),

c)
[INPUT] I(state.canopy.canopy_height, as_='h'),
[LINKED - e] I(state.canopy.d, as_='d'),
[LINKED - e] I(state.canopy.z0, as_='z0'),
[LINKED - d] I(state.met.L, as_='L'),
[LINKED - d] I(state.met.ustar_ref, as_='ustar_ref'),

d)
I(state.met.u_i, as_='u'),

e)
[INPUT] I(state.canopy.canopy_height, as_="h"),

f)
I(state.canopy_component_population[iLC]
              [iP].mean_gsto_per_layer[iL], as_='initial_leaf_gsto'),
I(state.canopy_layer_component[iL][iLC].mean_gsto, as_='initial_mean_gsto'),

[INPUT] I(state.canopy_layer_component[iL][iLC].gsto_params.f_phen, as_='f_phen'),
[INPUT] I(state.canopy_layer_component[iL][iLC].gsto_params.leaf_f_phen, as_='leaf_f_phen'),
I(state.canopy_layer_component[iL][iLC].gsto_params.f_light, as_='f_light'),
I(state.canopy_layer_component[iL][iLC].gsto_params.leaf_f_light, as_='leaf_f_light'),
I(state.canopy_layer_component[iL][iLC].gsto_params.f_temp, as_='f_temp'),
I(state.canopy_layer_component[iL][iLC].gsto_params.f_VPD, as_='f_VPD'),
I(state.canopy_layer_component[iL][iLC].gsto_params.f_SW, as_='f_SW'),
I(state.canopy_layer_component[iL][iLC].gsto_params.f_O3, as_='f_O3'),

g)
I(state.canopy.LAI_total, as_='LAI'),
I(state.canopy_layers[iL].micro_met.PARsun, as_='PARsun'),
I(state.canopy_layers[iL].micro_met.PARshade, as_='PARshade'),
I(state.canopy_layer_component[iL][iLC].LAIsunfrac, as_='LAIsunfrac'),
