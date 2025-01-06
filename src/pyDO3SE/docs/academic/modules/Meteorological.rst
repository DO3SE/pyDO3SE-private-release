Meteorological
==============

Download the original word file  :download:`Meteorological.docx <Meteorological.docx>`

.. container:: WordSection1

   .. _Toc36708868:

   Meteorological Derivations

   The meteorological module provides methods for converting and
   cleaning meteorological data into a format required by DO3SE.

   Contents

   `Transfer functions for wind speed and ozone
   concentration.2 <#toc50039730>`__

   `10.1 Wind speed.2 <#toc50039731>`__

   `10.2 Ozone concentration.4 <#toc50039732>`__

    

    

   .. _Toc50039730:

   \_

   .. _msoanchor-1:

   \_

   .. rubric:: Transfer functions for wind speed and ozone
      concentration[p1] 
      :name: transfer-functions-for-wind-speed-and-ozone-concentrationp1

   It is common for wind speed and O3 concentrations to be measured
   above a surface vegetation type that is different from that for which
   O3 deposition and stomatal O3 flux estimates are being performed
   (e.g. measurements are often made over a grassland surface where the
   “target” deposition or flux canopy is a nearby forest). When such a
   situation occurs, transfer methods need to be applied to allow the
   derivation of the necessary variables (namely O3 concentration and
   wind speed) at the “target” canopy height from variables that are
   measured above the alternative surface vegetation type, here referred
   to as the “reference” canopy.

   To be able to estimate wind speed and O3 concentrations for those
   situations where either one or both of the variables is measured
   above a “reference” surface vegetation cover requires derivation of
   these variables at a height sufficiently removed from the surface
   vegetation that the atmosphere will be decoupled from the underlying
   vegetation, and hence the O3 concentration or wind speed can be
   assumed to be independent of the underlying surface vegetation. In
   the following calculations, a height of 50 m is assumed appropriate
   to meet these criteria. The O3 concentration and wind speed variables
   can then be “bought back down” to the “target” canopy and hence
   re-coupled with the appropriate land surface type.

   The coupling of the vegetation to the atmosphere is described by the
   friction velocity (u*) which is, in part, determined by the surface
   vegetation, as such it is necessary to estimate u\* for canopies
   above which the variables wind speed and O3 concentration are
   measured, as well as for the “target” canopy; these will be referred
   to as u*Ref,w u*Ref,o and u*,tgt, respectively in the following text.

   .. _Toc36708869:

   \_

   .. _Toc50039731:

   \_

   .. rubric:: 10.1 Wind speed
      :name: wind-speed

   The wind speed at the “target” canopy height is necessary for the
   estimation of the boundary layer resistance of representative upper
   canopy leaves for leaf O\ 3 flux estimates.

   If the original wind speed was NOT measured over the “target” canopy
   but over a canopy of a different land-cover type (here referred to as
   the “reference” canopy), then the friction velocity for the reference
   canopy (u\*ref,w) has to be estimated as follows :-

    

   u\ \*ref,w =

    

   where u\ ref ,w is the wind speed measured at height h\ ref,w\ above
   the reference canopy over which the wind speed measurements are made,
   k is the von Karman constant (0.41), d\ ref,w is the canopy
   displacement height (where d\ ref,w equals h\ ref,w \* 0.7) and
   zo\ ref,w\ is the canopy roughness length (where zo\ ref,w\ equals
   h\ ref,w \* 0.1), with h\ ref,w being the height of the reference
   canopy above which the wind speed is measured. N.B. Note the wind
   speed has to be measured ABOVE the “reference” canopy with height
   h\ ref\ C\ w i.e. h\ ref,w > h\ ref\ C\ w.

   The u\ \*ref,w\ can then be used to estimate wind speed at an upper
   level height u\ ref,w\ (hup), here h\ up is assumed to be 50 m; where
   the atmosphere would be expected to be largely decoupled from any
   underlying surface vegetation layer.

    

   u\ ref,w\ (hup) =

    

   This 50 m wind speed u\ ref,w\ (hup) can then be used to estimate
   u\ \*tgt for the target canopy for which O\ 3 flux is to be estimated
   as follows :- 

    

   u\ \*tgt =

    

   where d\ tgt and zo\ tgt are the “target” canopy displacement and
   roughness length values calculated as 0.7 \* h\ tgt and 0.1 \* h\ tgt
   where h\ tgt is the “target” canopy height.

   The wind speed at the target canopy height u\ tgt\ (htgt) can then be
   estimated according to the following equation assuming a stable
   atmosphere :-

    

   u\ tgt\ (htgt) =

    

   If the original wind speed WAS measured over the target canopy
   u\ w\ (hw) at height h\ w then the wind speed at the top of the
   canopy, u\ tgt\ (htgt), can then be calculated as shown below
   assuming a stable atmosphere :-

    

   u\ tgt\ (htgt) =

    

   Figure x describes the variables necessary to estimate wind speed at
   the target canopy height u\ tgt\ (htgt) of a forest canopy from
   measurements of wind speed over a reference canopy (which in this
   example is represented by a grassland canopy).

    

   **Figure x**\ Variables necessary to estimate wind speed at a target
   canopy height u\ tgt\ (htgt) represented as a forest canopy from
   measurements of wind speed made over a reference canopy represented
   by a grassland canopy u\ ref,w\ (href,w).

   .. _Toc36708870:

   \_

   .. _Toc50039732:

   \_

   .. rubric:: 10.2 Ozone concentration
      :name: ozone-concentration

   To estimate the O\ 3 concentration at the canopy height it is
   necessary to estimate the atmospheric resistance to O\ 3 transfer
   from its measurement height to the top of the “target” canopy. As for
   wind speed, Ra is a function of the coupling between the canopy and
   the atmosphere (i.e. u*). If the O\ 3 concentration has not been
   measured above the “target” canopy it is necessary to estimate u\*
   for the “reference” canopy above which the O\ 3 concentration was
   measured (u\*ref,o). This can be achieved using the previous
   estimates of wind speed at the height (hup) (i.e. a height where the
   atmosphere is decoupled from the underlying surface vegetation)
   assuming that at h\ up wind speed is independent of local surface
   features so that u\ ref,o\ (hup) = u\ ref,w\ (hup) where u\ ref,o is
   the wind speed for the “ozone reference” canopy. The friction
   velocity for the “ozone reference” canopy can then be calculated as
   :-

    

   u\ \*ref,o=

    

   where d\ ref,o is the canopy displacement height (where d\ ref,o
   equals h\ ref,o \* 0.7) and zo\ ref,o is the canopy roughness length
   (where zo\ ref,o equals h\ ref,o \* 0.1) where h\ ref,o is the canopy
   height above which the O\ 3 concentration is measured, here referred
   to as the “reference” canopy.  

   The O\ 3 concentration at h\ up can then be estimated from the
   atmospheric resistance (Raref,o\ (ho,h\ up)) to O\ 3 between the
   measurement height over the reference canopy at which the O\ 3
   concentration is measured C\ ref,o\ (href,o) and h\ up:-

              Ra\ :sub:`ref,o`\ (h:sub:`o`,h\ :sub:`up`)=\ 

   Where h\ o is the O\ 3 concentration measurement height.

    

   The deposition velocity Vg\ ref,o at h\ up for the “ozone reference”
   canopy is estimated by :-

              Vg\ :sub:`ref,o`\ =

    

   Where

   Ra\ ref,o\ (zoref,o\ +d\ ref,o,h\ up) =

    

   Rb\ ref,o = quasi laminar resistance for „ozone reference” canopy

   Rb\ ref,o =

   Rsur\ :sub:`ref,o` = bulk surface resistance for “ozone reference”
   canopy.

   Rsur\ :sub:`ref,o` is estimated by the DO\ :sub:`3`\ SE model as
   described previously in this documentation.

   When using the DO\ :sub:`3`\ SE interface Ver.1 it is assumed that
   the R\ surof both the “target” and “reference” canopy are the same to
   avoid additional computational time. However, the user could perform
   two model runs to avoid uncertainties in the modelling that might
   result from differences in the Rsurof the reference and target
   canopies for both wind speed and ozone conversions. These model runs
   would be used to estimate the wind speed and/or ozone concentration
   at h\ :sub:`up` over the “reference” canopy (i.e. using
   DO\ :sub:`3`\ SE model parameterisation for this reference canopy)
   which can then be used in the final model run as inputs over the
   “target” canopy.

    

   The O\ :sub:`3` concentration at h\ :sub:`up` for the O\ :sub:`3`
   “reference” canopy is estimated by :-

   C\ :sub:`ref,o`\ (h:sub:`up`) =

   Here we assume that at h\ :sub:`up` the O\ :sub:`3` concentration is
   independent of local surface features so that
   C\ :sub:`tgt`\ (h:sub:`up`) = C\ :sub:`ref,o`\ (h:sub:`up`) where
   C\ :sub:`tgt` is the O\ :sub:`3` concentration at the top of the
   target canopy.

   Aerodynamic resistance between h\ :sub:`up` and h\ :sub:`tgt` for the
   “target” canopy is calculated by :-

   Ra\ :sub:`tgt` (h:sub:`up,` h\ :sub:`tgt`) =

   Deposition velocity at h\ up for the “target” canopy is estimated as
   :-

               Vg\ tgt\ (hup) =

   Where

   Ra\ tgt\ (zotgt\ +d\ tgt,h\ up) =

    

   Rb\ :sub:`tgt` = quasi laminar resistance for “target” canopy

   Rb\ :sub:`tgt` =

   Rsur\ :sub:`tgt` = bulk surface resistance for “target” canopy.

   Rsur\ :sub:`tgt` is estimated by the DO\ :sub:`3`\ SE model as
   described previously in this documentation.

   Finally, O\ :sub:`3` concentration at the “target” canopy is
   estimated as:-

   C\ :sub:`tgt` (h:sub:`tgt`) =
   C\ :sub:`ref,o`\ (h:sub:`up`).[1-Ra\ :sub:`tgt`
   (h:sub:`tgt`,h\ :sub:`up`).Vg\ :sub:`tgt` (h:sub:`up`)]

    

   If the O\ :sub:`3` concentration is measured above the “target”
   canopy, it is only necessary to estimate the aerodynamic resistance
   (Ra:sub:`tgt`\ (h:sub:`o`,h\ :sub:`tgt`))  between the O\ :sub:`3`
   measurement height (h:sub:`o`) and the target canopy height
   (h:sub:`tgt`), this can be achieved by substituting the h\ :sub:`up`
   for h\ :sub:`o` as follows :-

   Ra\ :sub:`tgt` (h:sub:`o,` h\ :sub:`tgt`) =

   Deposition velocity at h\ :sub:`up` for the “target” canopy is
   estimated as :-

                  Vg\ :sub:`tgt`\ (h:sub:`o`) =

   Where

   Ra\ :sub:`tgt`\ (zo:sub:`tgt`\ +d\ :sub:`tgt`,h\ :sub:`o`) =

    

   Rb\ :sub:`tgt` = quasi laminar resistance for “target” canopy

   Rb\ :sub:`tgt` =

   Rsur\ :sub:`tgt` = bulk surface resistance for “target” canopy.

   Rsur\ :sub:`tgt` is estimated by the DO\ :sub:`3`\ SE model as
   described previously in this documentation.

   Finally, O\ :sub:`3` concentration at the “target” canopy is
   estimated as:-

    

   C\ tgt (htgt) = C\ ref,o\ (ho).[1-Ra\ tgt (htgt,h\ o).Vg\ tgt (ho)]

    

   The O\ 3 concentration at the canopy height is then used in estimates
   of stomatal O\ 3 flux for assessments of O\ 3 impacts. Figure x
   describes the variables necessary to estimate O\ 3 concentration at
   the “target” canopy height C\ tgt\ (htgt) of a forest canopy from
   measurements of O\ 3 concentration made over a “reference” canopy
   (which in this example is represented by a grassland canopy).

    

   ** **

   ** **

   ** **

   ** **

   **Figure x**\ Variables necessary to estimate O\ 3 concentration at a
   “target” canopy height C\ tgt\ (htgt) represented as a forest canopy
   from measurements of O\ 3 concentration made over a “reference”
   canopy represented by a grassland canopy C\ ref,o\ (href,o).

    

    

.. container:: WordSection2

   .. rubric:: Above canopy PAR
      :name: above-canopy-par

   We calculate above canopy direct and indirect PAR components.

   Uses *calc_Idrctt_Idfuse* in code

    

   .. rubric:: Calculate potential PAR:
      :name: calculate-potential-par

            m = 1.0 / sinB

          # Potential direct PAR

           pPARdir = 600 \* exp(-0.185 \* (P / seaP) \* m) \* sinB

           # Potential diffuse PAR

           pPARdif = 0.4 \* (600 - pPARdir) \* sinB

           # Potential total PAR

           pPARtotal = pPARdir + pPARdif

   .

   .. rubric:: Sky Transmissivity
      :name: sky-transmissivity

   Can use either PAR input or cloudFrac to calculate the sky
   transmissivity.

                  ST = max(0.21, min(0.9, PAR / pPARtotal))

                  or

                  ST = 1.0 - 0.75 \* (cloudFrac \*\* 3.4)

   .. rubric:: Direct and Diffuse PAR component
      :name: direct-and-diffuse-par-component

    

           # Direct and diffuse fractions

           # A = 0.9

           # B = 0.7

           if cloudFrac is None or cloudFrac < 0.9:

               fPARdir = (pPARdir / pPARtotal) \* (1 - ((0.9 - ST) /
   0.7)**(2.0 / 3.0))

           else:

               fPARdir = (pPARdir / pPARtotal)

           fPARdif = 1 - fPARdir

    

           # Apply calculated direct and diffuse fractions to PARtotal

           Idrctt = fPARdir \* PAR

           Idfuse = fPARdif \* PAR

    

   .. rubric:: Sunlit and Shaded Leaf Calculations
      :name: sunlit-and-shaded-leaf-calculations

   To calculate sunlit and shaded leaves we use the equations from
   Farquhar 1997. We calculate the sun and shaded values per m\ :sup:`2`
   then multiply this by the leaf area index(LAI) to find the value for
   the layer/canopy.

   `Sunlit shaded equations diagram <media/sunlitshaded-diagram.pdf>`__

   **Constants used**

   Sigma = 0.15  # Leaf scattering coefficient of PAR (p_i + T_i)

   **Equations**

   All equations reference Farquhar 1997

       k_b = 0.5 / sinB  # Beam radiation extinction coefficient

       k_b_alt = 0.46 / sinB  # beam and scattered beam PAR extinction
   coefficient

   **Eq A20.**

       P_h = (1 - (1 - sigma) \*\* 0.5) / (1 + (1 - sigma) \*\* 0.5)

   **Eq A19.**

       P_cb = 1 - exp

   .. _msoanchor-2:

   [SB2] (-2 * P_h * k_b / (1 + k_b))Eq A21Eq A5.    Ir_diffuse = (1 - P_cd) * k_d_alt * Ir_dfuse_0 * exp(-k_d_alt * LAI_c)Eq A11.    Ir_beam_sun = (1 - sigma) * Ir_beam_0 * cosA / sinBEq A8.    Eq A12       PARsun = PARshade + Ir_beam_sunE1 A7    PARshade = Ir_diffuse + Ir_scattered_bMultilayer ApproachTo calculate the PAR sun and shade at each layer of the canopy we use the cumulative LAI(LAI_c) value at each layer of the canopy where LAI_c = 0 at the top of the canopy. We then calculate the PAR sun at each layer. See surface resistance documentation for how we then combine and upscale this back to the canopy level.resistance\surface_resistance\Up-scaling.html  ReferencesSunlit Shaded-        Farqhuar 1997 - DE PURY, D.G.G. and FARQUHAR, G.D. (1997), Simple scaling of photosynthesis from leaves to canopies without the errors of big‐leaf models. Plant, Cell & Environment, 20: 537-557. https://doi.org/10.1111/j.1365-3040.1997.00094.x

   .. container::

      --------------

      .. container::

         .. container:: msocomtxt

          \ \ `[p1] <#msoanchor-1>`__\ especially important for
         interface application

   .. container::

      .. container:: msocomtxt

         .. _msocom-2:

         \_msocom-2

          \ \ `[SB2] <#msoanchor-2>`__\ Note ‘-‘ was missing in Farquhar
         paper reference
