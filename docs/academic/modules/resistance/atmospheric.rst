Atmospheric Resistance (Ra)
===========================

Download the original word file  :download:`atmospheric.docx <atmospheric.docx>`

.. container:: WordSection1

   .. _Toc36708831:

   Atmospheric Resistance (Ra)

   The atmospheric resistance (Ra) term describes the resistance to O\ 3
   transfer due to mechanical and thermal atmospheric mixing processes
   between a specific reference height in the atmosphere at which the
   O\ 3 concentration is modelled or measured (zR) and the top of the
   plant canopy (z1) for which O\ 3 deposition or uptake is being
   calculated.

   Contents

   `Model Flow..2 <#toc50042651>`__

   `Atmospheric resistance (Ra)3 <#toc50042652>`__

   `2.2.1        Rasimp3 <#toc50042653>`__

   `2.2.2        Racomp3 <#toc50042654>`__

    

    

   .. _Toc50042651:

   \_

   .. rubric:: Model Flow
      :name: model-flow

   TODO:

    

   .. _Toc50042652:

   \_

   .. rubric:: Atmospheric resistance (Ra)
      :name: atmospheric-resistance-ra-1

    

   We have to explain when Ra is needed and when not (e.g. for
   leaf-level applications with available ozone data measured at top of
   plant canopy)

   The atmospheric resistance (Ra) term describes the resistance to O\ 3
   transfer due to mechanical and thermal atmospheric mixing processes
   between a specific reference height in the atmosphere at which the
   O\ 3 concentration is modelled or measured (zR) and the top of the
   plant canopy (z1) for which O\ 3 deposition or uptake is being
   calculated.

   Ra is estimated using two different methods, the selection of which
   depends primarily upon the complement of available input data. The
   more data intensive method incorporates both mechanical and thermal
   atmospheric mixing processes in the estimation of Ra (here termed
   Ra\ comp) and requires heat flux data (both sensible and latent)
   which are often unavailable at site-specific locations. The less data
   intensive method only incorporates mechanical atmospheric mixing
   processes (i.e. assumes neutral stability of the atmospheric profile)
   in the estimation of Ra (here termed Ra\ simp).

   These Ra formulations both incorporate canopy roughness
   characteristics which will in part determine atmospheric turbulence.
   The vegetation roughness length (*zo)* and displacement height *(d)*
   are approximately 10 % and 70 % of the vegetation height (*h*)
   respectively (references). Default values for *h* used in EMEP
   DO\ :sub:`3`\ SE are given in Table x\ `[p1] <#msocom-1>`__\ \ \  .
   Here, z1 is more accurately defined as the top of the quasi-laminar
   canopy boundary layer (zo+d).

    

   The interfaced version of the DO\ 3\ SE model uses the Rasimp model
   formulations.

    

   .. _Toc36708832:

   \_

   .. _Toc50042653:

   \_

   .. rubric:: 2.2.1     Rasimp
      :name: rasimp

   Assuming neutral stability, Ra\ simp and friction velocity (u\*) can
   be estimated as follows using standard methods consistent with those
   described in the UNECE (2004).

             Ra\ simp =

    

   The friction velocity, u\ \* (m s\ -1) is derived from :-

              u\*\ `[B2] <#msocom-2>`__\ \ \    =

   Where k is the von Karman constant (0.41), u\* is the friction
   velocity, u\ (z) is the windspeed at height z, d is the displacement
   height and zo is the roughness length.

   .. _Toc36708833:

   \_

   .. _Toc50042654:

   \_

   .. rubric:: 2.2.2     Racomp
      :name: racomp

    

   Where data are available to determine the actual stability profile of
   the atmosphere, Ra\ comp can be approximated following the procedure
   given by Garland (1978), using the Monin-Obukhov similarity relations
   for the dimensionless wind shear and potential temperature gradient
   suggested by Businger et al. (1971):

   Ra\ comp =\ `[p3] <#msocom-3>`__\ \ \  

   Where k is the von Karman constant (0.41), u\* is the friction
   velocity, d is the displacement height, zo is the vegetation
   roughness length and L is the Monin-Obukhov length. ψh is the
   integrated stability function for heat.

   According to Arya (1988), integration of simplified forms of the
   similarity functions with respect to height yields the formulations
   used for the stability function for heat (ψh) and momentum (ψm)

    

   ψh(ξ) = ψm (ξ) = -5
   ξ                                                                          if
   ξ ≥ 0 (stable)

   ψ\ m(ζ ) = ln                if ζ < 0 (unstable)

   ψ\ h(ζ ) = 2
   ln                                                                        if
   ζ < 0 (unstable)

   where x\ m = (1-16 ζ ) 1/4 and x\ h = (1-16 ζ
   )\ `[SB4] <#msocom-4>`__\ \ \ \  1/2 and ζ = (z-d)/L or zo/L.

    

   The friction velocity, u\ \* (m s\ -1) is derived from :-

   u\* = ,

    

   where τ is turbulent surface stress (kg / m / s\ -2) and ρ is surface
   density of dry air (kg/m3) derived from :-

   ρ=

    

   where P is surface air pressure (Pa), Rmass is the mass gas constant
   for dry air (J / Kg / K) and T\ 2 is the surface air temperature
   (Kelvin).

   The Monin-Obukhov length,L, is derived from

              L = -

   `[p5] <#msocom-5>`__\ \ \  

   where c\ p is the specific heat capacity of dry air; k is von
   Karman’s constant (0.41), g is the gravitational acceleration (9.81 m
   s\ -2) and H is the surface flux of sensible heat (W m\ -2).

.. container::

   --------------

   .. container::

      .. container:: msocomtxt

         .. _msocom-1:

         \_msocom-1

          \ \ \ `[p1] <#msoanchor-1>`__\ This table has still to be
         included

   .. container::

      .. container:: msocomtxt

         .. _msocom-2:

         \_msocom-2

          \ \ \ `[B2] <#msoanchor-2>`__\ Not sure we have always
         modelled this correctly (though in most code it has not
         mattered as we have not used the resulting Ra…need to check!)

   .. container::

      .. container:: msocomtxt

         .. _msocom-3:

         \_msocom-3

          \ \ \ `[p3] <#msoanchor-3>`__\ We might want to stress the
         difference between zR and z1 here

   .. container::

      .. container:: msocomtxt

         .. _msocom-4:

         \_msocom-4

          \ \ \ `[SB4] <#msoanchor-4>`__\ Emep model update:

         We have x_h and x_m

          

   .. container::

      .. container:: msocomtxt

         .. _msocom-5:

         \_msocom-5

          \ \ \ `[p5] <#msoanchor-5>`__\ Double check this?
