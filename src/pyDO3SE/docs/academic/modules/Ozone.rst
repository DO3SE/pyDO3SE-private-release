Estimating Ozone Deposition (Total & Stomatal)
==============================================

Download the original word file  :download:`Ozone.docx <Ozone.docx>`

.. container:: WordSection1

   Estimating Ozone Deposition (Total & Stomatal)

   Calculate the multi-layer ozone resistance model, deposition velocity
   and multi-layer ozone concentration. Calculated at leaf level
   separately.

   Contents

   `Ozone Deposition.1 <#toc50029803>`__

   `Total and Stomatal Ozone flux (Ftot and Fst)1 <#toc50029804>`__

   `Total O3 flux (Ftot)1 <#toc50029805>`__

   `Stomatal O3 flux (Fst[p5] )1 <#toc50029806>`__

   `Accumulated stomatal flux (AFstY)1 <#toc50029807>`__

    

   | 

    

   .. rubric:: Model Flow
      :name: model-flow

   TODO: Model Flow

   .. _Toc50029803:

    

   .. rubric:: Ozone Deposition
      :name: ozone-deposition

    

   The DO\ 3\ SE model’s total dry deposition scheme is shown in Fig
   1\ `[B1] <#msocom-1>`__\ \  . The model assumes the key resistances
   to O\ 3 deposition are the aerodynamic resistance (Ra), the
   quasi-laminar sub-layer resistance (Rb) above the canopy, and the
   surface resistance (Rsur). Rsurcomprises two resistance paths in
   series; the stomatal and non-stomatal resistance. The latter
   represents a) within canopy aerodynamic resistance (Rinc) and
   subsequent b) soil resistance to decomposition at the soil surface
   (Rgs) which encompasses features such as leaf litter and ground
   vegetation under forest canopies, as well as c) resistance to
   adsorption to the external plant parts (Rext) including cuticle, bark
   etc… a comment on planned revision of estimation of Rext here?

    

   As such, the loss rate of ozone to the ground surface, within a
   volume unit area and height Δ\ z is given by the product of the
   deposition velocity (Vg) at height z\ ref and the concentration (Ci)
   at that height as described in eq.1.

    

   .. _Ref393790837:

   1

    

   Where *Vg* is estimated as given in eq. 2

   .. _Ref393790933:

   2

   .. _Toc50029804:

   \_

   .. _Toc50027007:

   \_

   .. _Toc49267186:

   \_

   .. rubric:: Total and Stomatal Ozone flux (Ftot and Fst)
      :name: total-and-stomatal-ozone-flux-ftot-and-fst

    

   .. _Toc50029805:

   \_

   .. _Toc50027008:

   \_

   .. _Toc36708859:

   \_

   .. _Toc49267187:

   \_

   .. rubric:: Total O3 flux (Ftot)
      :name: total-o3-flux-ftot

   The loss rate of O\ 3 to the land surface (Ftot) is dependent upon
   the volume, defined by the unit area and height at which the O\ 3
   concentration is given and the deposition velocity between that
   height and the surface.

    

   In the interface O\ 3 concentrations can be provided either above a
   “reference” canopy or the “target” canopy (see section 3.1). As such,
   it is theoretically possible to calculate different Ftot values using
   the following formulations:

    

   The O\ 3 loss rate (Ftot50) from 50m above the canopy (h50) when O\ 3
   concentrations are provided above a “reference” canopy :-

    

              Vg\ :sub:`(50)` = |image0|

    

               Ftot\ :sub:`50` = -Vg(:sub:`50`) \* O\ :sub:`3`\ (*50*)

    

   The O\ 3 loss rate (FtotO3hz) from the height above the “target”
   canopy (O3\ h\ z) that the O\ 3 concentration (O3\ (z)) was provided
   :- add schematic here

    

               Vg(O3hz) = |image1|

    

               FtotO3hz = -Vg(O3hz) \* O3(z)

    

    

   Finally, it is also possible to estimate the O\ 3 loss rate (Ftoth)
   from the canopy top (h) when O\ 3 concentrations are calculated for
   the top of the “target” canopy as is done in the interface code :-

    

               Vg(h) = |image2|

    

               Ftoth = -Vg(h) \* O3(h)

    

    

   In the current version of the interface Ftot\ h is the O\ 3 loss rate
   provided as output.

   .. _Toc50029806:

   \_

   .. _Toc50027009:

   \_

   .. _Toc36708860:

   \_

   .. _Toc49267188:

   \_

   .. _msoanchor-5:

   \_

   .. rubric:: Stomatal O3 flux (Fst[p5] )
      :name: stomatal-o3-flux-fstp5

    

   The estimation of stomatal flux of ozone (Fst) is based on the
   assumption that the concentration of O\ 3 at the top of the canopy
   (O3\ (h)) represents a reasonable estimate of the O\ 3 concentration
   at the upper surface of the laminar layer of the sunlit upper canopy
   leaves (and the flag leaf in the case of wheat). If O\ 3\ (h) is the
   concentration of O\ 3 at canopy top (height h, unit: m), in nmol
   m\ -3, then Fst (nmol m\ -2 PLA s\ -1), is given by:

    

   Fst = O3(h) \* gsto \* |image3|

    

   Where rc is equal to 1/(g\ sto\ +g\ ext), here g\ sto and g\ ext are
   given in units of m/s. At normal temperatures and air pressure, the
   conversion is made by dividing the conductance values expressed in
   mmol m\ -2 s\ -1 by 41000 to given conductance in m/s. A value for
   g\ ext of 1/2500 is used to maintain consistency with the value of
   2500 s/m used in the canopy scale estimate of Rsur.   

    

   .. _Toc50029807:

   \_

   .. _Toc50027010:

   \_

   .. _Toc36708861:

   \_

   .. _Toc49267189:

   \_

   .. rubric:: Accumulated stomatal flux (AFstY)
      :name: accumulated-stomatal-flux-afsty

    

   The accumulated flux above an O\ 3 stomatal flux rate threshold of Y
   nmol m\ -2 s\ -1 (AFst\ Y) is calculated as described below with the
   accumulation estimated over the respective accumulation period.

    

   AFstY  =  |image4| for Fsti ³ Y nmol m-2 PLA
   s-1                                             

    

   where Fst\ i is the hourly O\ 3 mean flux in nmol m\ -2 PLA s\ -1,
   and n is the number of hours within the accumulation period.

    

.. container::

   --------------

   .. container::

      .. container:: msocomtxt

         .. _msocom-1:

         \_msocom-1

          \ \ `[B1] <#msoanchor-1>`__\ Can we identify the different
         versions in this figure? Maybe we need separate figures for
         separate model versions?

.. |image0| image:: ../../../../Documentation/SplitSections%20-%20Copy/Photosynthesis_files/image002.png
   :width: 114px
   :height: 41px
.. |image1| image:: ../../../../Documentation/SplitSections%20-%20Copy/Photosynthesis_files/image003.png
   :width: 130px
   :height: 41px
.. |image2| image:: ../../../../Documentation/SplitSections%20-%20Copy/Photosynthesis_files/image004.png
   :width: 67px
   :height: 39px
.. |image3| image:: ../../../../Documentation/SplitSections%20-%20Copy/Photosynthesis_files/image005.png
   :width: 48px
   :height: 39px
.. |image4| image:: ../../../../Documentation/SplitSections%20-%20Copy/Photosynthesis_files/image006.png
   :width: 86px
   :height: 46px
