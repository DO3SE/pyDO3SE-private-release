Multiplicative stomatal resistance (Rsto)
=========================================

Download the original word file  :download:`Multiplicative.docx <Multiplicative.docx>`

.. container:: WordSection1

   Multiplicative stomatal resistance (Rsto)

    

   Contents

   `3.1.1.       Maximum stomatal conductance ( <#toc49276282>`__\ )2

    

   Mean canopy level stomatal conductance to O\ 3 (G\ sto), the inverse
   of stomatal resistance (R\ sto), is calculated with the
   multiplicative model of Emberson et al. (2000a):

    

               Gsto = gmax \* Fphen \* Flight \* max{fmin, (ftemp \*
   fVPD \* fSWP)}

    

   where g\ max is the leaf/needle-level maximum stomatal conductance
   (mmol O\ 3 m\ -2 PLA s\ -1), F\ phen is the bulk canopy stomatal
   conductance relationship with seasonal canopy age, F\ light is the
   whole canopy stomatal conductance relationship with irradiance
   penetrating within the canopy, f\ min is the minimum daylight
   stomatal conductance, and f\ temp, f\ VPD and f\ SWP are the stomatal
   conductance relationships with temperature, vapour pressure deficit
   and soil water potential respectively. This calculation gives a mean
   g\ sto value for all leaves/needles in the canopy. As such, G\ sto
   can be scaled according to LAI to estimate whole canopy stomatal
   conductance (Gsto\ canopy\ \ `[p1] <#msocom-1>`__\ \  ).

   Leaf level stomatal conductance to O\ 3 (g\ sto) is calculated again
   with the multiplicative model of Emberson et al. (2000a) but
   following the formulations of UNECE (2004) :-

    

   g\ sto = g\ :sub:`max` \* {min (f:sub:`phen`, fO\ 3) \*
   f\ :sub:`light` \* {max (f:sub:`min`, f\ :sub:`temp` \* f\ :sub:`VPD`
   \* f\ :sub:`SWP`})

    

   where g\ max is the leaf/needle-level maximum stomatal conductance
   (mmol O\ 3 m\ -2 PLA s\ -1),\ `[p2] <#msocom-2>`__\  f\ phen is the
   leaf/needle level stomatal conductance relationship with leaf/needle
   age, fO3 is..... (see below), f\ light is the stomatal conductance
   relationship with irradiance at the top of the canopy, f\ min is the
   minimum daylight stomatal conductance, and f\ temp, f\ VPD and f\ SWP
   are the stomatal conductance relationships with temperature, vapour
   pressure deficit and soil water potential
   respectively\ \ `[p3] <#msocom-3>`__\ \  . This calculation gives a
   leaf/needle level g\ sto value for leaves/needles representative of
   those at the top of the canopy (or in the case of wheat, the flag
   leaf). As such, g\ sto can then be used to estimate upper canopy
   leaf/needle stomatal O\ 3 flux for risk assessment. 

   The term fO\ 3, which allows for O\ 3 exposure to cause early
   senescence, is currently only established for wheat and potato and
   therefore is set equal to 1 for all other cover-types and species.

   The flux-effect models developed by Pleijel et al. (2002) and
   Danielsson et al. (2003) include a function to allow for the
   influence of ozone concentrations on stomatal conductance (fO\ 3) on
   wheat and potato via the onset of early senescence. As such this
   function is used in association with the fphen function to estimate
   gsto. The fO\ 3 function typically operates over a one-month period
   and only comes into operation if it has a stronger
   senescence-promoting effect than normal senescence. The functions are
   given in Equations 3.21 and 3.22. The ozone function for spring wheat
   (based on Danielsson et al. (2003) but recalculated for PLA):

   fO\ 3= ((1+(AFst0/11.5)10)-1) [xxx]

   where AFst0 is accumulated from Astart

   The ozone function for potato (based on Pleijel et al. (2002)):

   fO\ 3 = ((1+(AOT0/40)5)-1) [xxx]

   where AOT0 is accumulated from Astart

    

   .. _Toc49276282:

   \_

   .. _Toc36708843:

   \_

   .. rubric:: Maximum stomatal conductance ()
      :name: maximum-stomatal-conductance

   Table 1lists absolute maximum stomatal conductance values (in mmol
   O\ 3 m\ -2 PLA s\ -1, denoted, where the exponent m stands for the
   median of g\ s values taken from the literature) and f\ min values
   (provided as a fraction of).  These values are given for the
   vegetation cover types (as necessary for regional deposition
   modelling) and by key species (from which the cover typeand f\ min
   values have been derived).

   The values for are provided often provided for gases other than O\ 3
   (i.e. H\ 2\ O vapour or CO\ 2), for total or projected leaf/needle
   areas and in conductance units of e.g. mm s\ -1 rather than molar
   units.

    

   For pressure (P) and temperature (T), g\ max in m s-1 units is given
   by:

    

   g\ max= RT/P

    

   R is here the gas-constant (8.314 J/mol/K). At
   normal\ `[p4] <#msocom-4>`__\ \  temperature and pressure, g\ max ≈
   /14000.

   `[p5] <#msocom-5>`__\  

    

    

    

    

.. container:: WordSection2

   .. _Ref395016154:

   \_Ref395016154

   .. _Ref395016168:

   Table

   1 Default deposition land-cover and species class values for  and
   f\ :sub:`min` parameters.
    

   +-------------+-------------+-------------+-------------+-------------+
   | **Land-cove | **Climate   | **g\ max**\ | **f\ min**  | **Reference |
   | r           | region**    |  mmol       |             | **          |
   | type &      |             | O\ 3 m\ -2  |             |             |
   | Species**   |             | PLA s\ -1   |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   | ** **       |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   | **Coniferou |             | 160         | 0.1         | Simpson et  |
   | s           |             |             |             | al. (2003)  |
   | Forests     |             |             |             |             |
   | (CF)**      |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   | Norway      | Northern    | 112         | 0.1         | UNECE       |
   | spruce      | Europe      | (111-118(11 |             | (2004);     |
   |             |             | 9))         |             | Zimmerman   |
   | (Picea      |             |             |             | et al.      |
   | abies)      |             |             |             | (1988)      |
   |             |             |             |             | [112];      |
   |             |             |             |             | Sellin.     |
   |             |             |             |             | (2001)      |
   |             |             |             |             | [119];      |
   |             |             |             |             | Hansson et  |
   |             |             |             |             | al.(in      |
   |             |             |             |             | prep) [111] |
   +-------------+-------------+-------------+-------------+-------------+
   | Scots Pine  | Atlantic    | 180         | 0.1         | UNECE       |
   |             | Central     | (171-188)   |             | (2004);     |
   | (Pinus      | Europe      |             |             | Whitehead   |
   | sylvestris) |             |             |             | et al.      |
   |             |             |             |             | (1984)      |
   |             |             |             |             | [188];      |
   |             |             |             |             | Beadle et   |
   |             |             |             |             | al. (1985)  |
   |             |             |             |             | [175];      |
   |             |             |             |             | Sturm et    |
   |             |             |             |             | al. (1998)  |
   |             |             |             |             | [171]       |
   +-------------+-------------+-------------+-------------+-------------+
   | Norway      | Continental | 125         | 0.16        | UNECE       |
   | Spruce      | Central     | (87-140)    |             | (2004) cf.  |
   |             | Europe      |             |             | Körner et   |
   | (Picea      |             |             |             | al. (1979)  |
   | abies)      |             |             |             | [87]; Dixon |
   |             |             |             |             | et al.      |
   |             |             |             |             | (1995)      |
   |             |             |             |             | [121];      |
   |             |             |             |             | Emberson et |
   |             |             |             |             | al. (2000)  |
   |             |             |             |             | [130];      |
   |             |             |             |             | Zweifel et  |
   |             |             |             |             | al. (2000,  |
   |             |             |             |             | 2001, 2002) |
   |             |             |             |             | [140];      |
   |             |             |             |             | Korner.     |
   |             |             |             |             | (1979)      |
   |             |             |             |             |             |
   |             |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   | **Deciduous |             | 134         | 0.13        | Simpson et  |
   | Forests**   |             |             |             | al. (2003)  |
   |             |             |             |             |             |
   | **(DF)**    |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   | **Generic   | **All       | **150       | **0.1**     | **UNECE     |
   | Deciduous** | Europe**    | (100-180)** |             | (2004)**    |
   +-------------+-------------+-------------+-------------+-------------+
   | Silver      | Northern    | 196         | 0.1         | UNECE       |
   | birch       | Europe      | (180-211)   |             | (2004);     |
   |             |             |             |             | Uddling et  |
   | (Betula     |             |             |             | al. (2005a) |
   | pendula)    |             |             |             | [180];      |
   |             |             |             |             | Sellin et   |
   |             |             |             |             | al.(2005)   |
   |             |             |             |             | [211]       |
   +-------------+-------------+-------------+-------------+-------------+
   | Beech       | Atlantic    | 150         | 0.1         | UNECE       |
   |             | Central     | (100-180)   |             | (2004);     |
   | (Fagus      | Europe      |             |             |             |
   | sylvatica)  |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   | Oak         | Atlantic    | 230         | 0.06        | UNECE       |
   |             | Central     | (177-325)   |             | (2004);     |
   | (Quercus    | Europe      |             |             | Breda et    |
   | petraea &   |             |             |             | al. (1995)  |
   | robur)      |             |             |             | [228];      |
   |             |             |             |             | Epron &     |
   |             |             |             |             | Dreyer.     |
   |             |             |             |             | (1993)      |
   |             |             |             |             | [177];      |
   |             |             |             |             | Breda et    |
   |             |             |             |             | al. (1993a) |
   |             |             |             |             | [233];      |
   |             |             |             |             | Breda et    |
   |             |             |             |             | al. (1993b) |
   |             |             |             |             | [275];      |
   |             |             |             |             | Q.robur     |
   |             |             |             |             | from Epron  |
   |             |             |             |             | & Dreyer.   |
   |             |             |             |             | (1993)      |
   |             |             |             |             | [198];      |
   |             |             |             |             | Dolman &    |
   |             |             |             |             | Van den     |
   |             |             |             |             | Burg.(1988) |
   |             |             |             |             | [264]       |
   +-------------+-------------+-------------+-------------+-------------+
   | Beech       | Continental | 150         | 0.13        | UNECE       |
   |             | Central     | (132-300)   |             | (2004);     |
   | (Fagus      | European    |             |             | Nunn et al. |
   | sylvatica)  |             |             |             | (2005)      |
   |             |             |             |             | [147];      |
   |             |             |             |             | Matyssek et |
   |             |             |             |             | al. (2004)  |
   |             |             |             |             | [132]; Keel |
   |             |             |             |             | et al.      |
   |             |             |             |             | (2007)      |
   |             |             |             |             | [180];      |
   |             |             |             |             | Kutsch et   |
   |             |             |             |             | al. (2001)  |
   |             |             |             |             | [300];      |
   |             |             |             |             | Freeman     |
   |             |             |             |             | (1998) cf.  |
   |             |             |             |             | Medlyn et   |
   |             |             |             |             | al.(2001)   |
   |             |             |             |             | [180]; cf.  |
   |             |             |             |             | Körner et   |
   |             |             |             |             | al.(1979)   |
   |             |             |             |             | [150];      |
   |             |             |             |             | Schaub      |
   |             |             |             |             | (pers.      |
   |             |             |             |             | comm.)      |
   |             |             |             |             | [137];      |
   |             |             |             |             | Korner.     |
   |             |             |             |             | (1994)      |
   +-------------+-------------+-------------+-------------+-------------+
   | Beech       | Mediterrane | 145         | 0.02        | UNECE       |
   |             | an          | (100-183)   |             | (2004);     |
   | (Fagus      | Europe      |             |             | Raftoyannis |
   | sylvatica)  |             |             |             | & Radoglou  |
   |             |             |             |             | (2002)      |
   |             |             |             |             | [156]; cf.  |
   |             |             |             |             | Körner et   |
   |             |             |             |             | al.(1979)   |
   |             |             |             |             | [100; 140]; |
   |             |             |             |             | Nunn et al. |
   |             |             |             |             | (2005)      |
   |             |             |             |             | [147];      |
   |             |             |             |             | Matyssek et |
   |             |             |             |             | al.(2004)   |
   |             |             |             |             | [132];      |
   |             |             |             |             | Aranda et   |
   |             |             |             |             | al. (2000)  |
   |             |             |             |             | [183]       |
   +-------------+-------------+-------------+-------------+-------------+
   | **Needlelea |             | 180         | 0.13        | Simpson et  |
   | f           |             |             |             | al. (2003)  |
   | Forests**   |             |             |             |             |
   |             |             |             |             |             |
   | **(NF)**    |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   | Aleppo Pine | Mediterrane | 215         | 0.15        | UNECE       |
   |             | an          |             |             | (2004);     |
   | *(Pinus     | Europe      |             |             | Elvira et   |
   | halepensis) |             |             |             | al (2007)   |
   | *           |             |             |             | [215]       |
   +-------------+-------------+-------------+-------------+-------------+
   | **Broadleaf |             | 200         | 0.03        | Simpson et  |
   | Forests**   |             |             |             | al. (2003)  |
   |             |             |             |             |             |
   | **(BF)**    |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   | **Generic   | **All       | **175       | **0.02**    | **UNECE     |
   | Evergreen   | Europe**    | (70-365)**  |             | (2004)**    |
   | Mediterrane |             |             |             |             |
   | an**        |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   | Holm Oak    | Mediterrane | 180         | 0.02        | Rhizopoulos |
   |             | an          | (134-365\ \ |             | & Mitrakos  |
   | (Quercus    | Europe      |  `[B6] <#ms |             | (1990)      |
   | ilex)       |             | ocom-6>`__\ |             | [250];      |
   |             |             |  \  )       |             | Manes et    |
   |             |             |             |             | al.(1997)   |
   |             |             |             |             | [366];      |
   |             |             |             |             | Filho et    |
   |             |             |             |             | al. (1998)  |
   |             |             |             |             | [225];      |
   |             |             |             |             | Tognetti et |
   |             |             |             |             | al.(1998)   |
   |             |             |             |             | [195]; Sala |
   |             |             |             |             | & Tenhunen  |
   |             |             |             |             | (1994)      |
   |             |             |             |             | [165];      |
   |             |             |             |             | Alonso et   |
   |             |             |             |             | al. (2007)  |
   |             |             |             |             | [183, 191]; |
   |             |             |             |             | Infante et  |
   |             |             |             |             | al. (1999)  |
   |             |             |             |             | [323];      |
   |             |             |             |             | Castell et  |
   |             |             |             |             | al.(1994)   |
   |             |             |             |             | [177];      |
   |             |             |             |             | Damesin et  |
   |             |             |             |             | al. (1998)  |
   |             |             |             |             | [171];      |
   |             |             |             |             | Mediavilla  |
   |             |             |             |             | & Escudero  |
   |             |             |             |             | (2003)      |
   |             |             |             |             | [122];      |
   |             |             |             |             | Corcuera et |
   |             |             |             |             | al. (2005)  |
   |             |             |             |             | [134];      |
   |             |             |             |             | Gratani et  |
   |             |             |             |             | al.(2000)   |
   |             |             |             |             | [159];      |
   |             |             |             |             | Bussotti &  |
   |             |             |             |             | Ferretti    |
   |             |             |             |             | (2007)      |
   |             |             |             |             | [166, 188,  |
   |             |             |             |             | 156]        |
   +-------------+-------------+-------------+-------------+-------------+
   | **Temperate |             | 300         | 0.01        | Simpson et  |
   | crops**     |             |             |             | al. (2003)  |
   |             |             |             |             |             |
   | **(TC)**    |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   | **Generic   | **All       | **450**     | **0.01**    | **UNECE     |
   | crop**      | Europe**    |             |             | (2004)**    |
   +-------------+-------------+-------------+-------------+-------------+
   | Wheat       | All Europe  | 450         | 0.01        | UNECE       |
   |             |             |             |             | (2004);Arau |
   | (Triticum   |             |             |             | s           |
   | aestivum)   |             |             |             | et          |
   |             |             |             |             | al.(1989) ; |
   |             |             |             |             | Ali et al.  |
   |             |             |             |             | (1999) ;    |
   |             |             |             |             | Gruters et  |
   |             |             |             |             | al. (1995); |
   |             |             |             |             |             |
   |             |             |             |             | Korner et   |
   |             |             |             |             | al. (1979); |
   |             |             |             |             | Danielsson  |
   |             |             |             |             | et al.      |
   |             |             |             |             | (2003);De   |
   |             |             |             |             | la Torre    |
   |             |             |             |             | (2004):     |
   |             |             |             |             |             |
   |             |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   | **Mediterra |             | 156         | 0.019       | Simpson et  |
   | nean        |             |             |             | al. (2003)  |
   | crops**     |             |             |             |             |
   |             |             |             |             |             |
   | **(MC)**    |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   | Maize       | All Europe  | 305 (27     | 0.05        | ICP         |
   |             |             | s.d.)       |             | Vegetation  |
   | *(Zea       |             |             |             | contract    |
   | mays)*      |             |             |             | report      |
   |             |             |             |             | (2006);     |
   |             |             |             |             | Körner et   |
   |             |             |             |             | al. (1979)  |
   |             |             |             |             | [315];      |
   |             |             |             |             | Sinclair et |
   |             |             |             |             | al.         |
   |             |             |             |             | (1975) [355 |
   |             |             |             |             | ];          |
   |             |             |             |             | Stigter &   |
   |             |             |             |             | Lammers     |
   |             |             |             |             | (1974)      |
   |             |             |             |             | [295]; Tard |
   |             |             |             |             | ieu         |
   |             |             |             |             | et al.      |
   |             |             |             |             | (1991);     |
   |             |             |             |             | Ozier-Lafon |
   |             |             |             |             | tein        |
   |             |             |             |             | et          |
   |             |             |             |             | al.(1998)   |
   |             |             |             |             | [300]       |
   +-------------+-------------+-------------+-------------+-------------+
   | Sunflower   | All Europe  | 370 (230    | 0.05        | ICP         |
   |             |             | s.d.)       |             | Vegetation  |
   | *(Helianthu |             |             |             | contract    |
   | s           |             |             |             | report      |
   | annuus)*    |             |             |             | (2006);     |
   |             |             |             |             | Ward &      |
   |             |             |             |             | Bunce       |
   |             |             |             |             | (1986)      |
   |             |             |             |             | [390];      |
   |             |             |             |             | Connor &    |
   |             |             |             |             | Jones       |
   |             |             |             |             | (1985)      |
   |             |             |             |             | [325];      |
   |             |             |             |             | Wookey et   |
   |             |             |             |             | al.(1991)   |
   |             |             |             |             | [153];      |
   |             |             |             |             | Schurr et   |
   |             |             |             |             | al. (1992)  |
   |             |             |             |             | [397];      |
   |             |             |             |             | Rivelli et  |
   |             |             |             |             | al. (2002)  |
   |             |             |             |             | [586];      |
   |             |             |             |             | Hirasawa et |
   |             |             |             |             | al. (1995)  |
   |             |             |             |             | [473];      |
   |             |             |             |             | Wample &    |
   |             |             |             |             | Thornton    |
   |             |             |             |             | (1984)      |
   |             |             |             |             | [233]; Fay  |
   |             |             |             |             | & Knapp     |
   |             |             |             |             | (1996)      |
   |             |             |             |             | [1104];     |
   |             |             |             |             | Steduto et  |
   |             |             |             |             | al. (2000)  |
   |             |             |             |             | [732];      |
   |             |             |             |             | Turner et   |
   |             |             |             |             | al. (1984)  |
   |             |             |             |             | [323];      |
   |             |             |             |             | Turner et   |
   |             |             |             |             | al. (1985)  |
   |             |             |             |             | [372];      |
   |             |             |             |             | Quick et    |
   |             |             |             |             | al. (1992)  |
   |             |             |             |             | [350];      |
   |             |             |             |             | Körner et   |
   |             |             |             |             | al.(1979)   |
   |             |             |             |             | [385; 355;  |
   |             |             |             |             | 355]        |
   +-------------+-------------+-------------+-------------+-------------+
   | Tomato      | All Europe  | 285 (74     | 0.01        | ICP         |
   |             |             | s.d.)       |             | Vegetation  |
   | *(Solanum   |             |             |             | contract    |
   | lycopersicu |             |             |             | report      |
   | m)*         |             |             |             | (2006);     |
   |             |             |             |             | Duniway     |
   |             |             |             |             | (1971)      |
   |             |             |             |             | [385],      |
   |             |             |             |             | Moreshet &  |
   |             |             |             |             | Yocum       |
   |             |             |             |             | (1972)      |
   |             |             |             |             | [200];Kater |
   |             |             |             |             | ji          |
   |             |             |             |             | et          |
   |             |             |             |             | al.(1998)   |
   |             |             |             |             | [216];      |
   |             |             |             |             | Bakker      |
   |             |             |             |             | (1991)      |
   |             |             |             |             | [283];      |
   |             |             |             |             | Boulard et  |
   |             |             |             |             | al. (1991)  |
   |             |             |             |             | [208];      |
   |             |             |             |             | Pirker et   |
   |             |             |             |             | al. (2003)  |
   |             |             |             |             | [307, 384,  |
   |             |             |             |             | 288]        |
   +-------------+-------------+-------------+-------------+-------------+
   | Grape vine  | All Europe  | 215 (51     | 0.01        | ICP         |
   |             |             | s.d.)       |             | Vegetation  |
   | *(Vitis     |             |             |             | contract    |
   | vinifera)*  |             |             |             | report      |
   |             |             |             |             | (2006);     |
   |             |             |             |             | Schultz     |
   |             |             |             |             | (2003)      |
   |             |             |             |             | [225], Naor |
   |             |             |             |             | & Wample    |
   |             |             |             |             | (1995)      |
   |             |             |             |             | [210];      |
   |             |             |             |             | Schultz     |
   |             |             |             |             | (2003a)     |
   |             |             |             |             | [134,       |
   |             |             |             |             | 150],       |
   |             |             |             |             | Medrano et  |
   |             |             |             |             | al. (2003)  |
   |             |             |             |             | [279, 204]; |
   |             |             |             |             | Patakas et  |
   |             |             |             |             | al. (2003)  |
   |             |             |             |             | [211, 216,  |
   |             |             |             |             | 201];       |
   |             |             |             |             | Winkel &    |
   |             |             |             |             | Rambal      |
   |             |             |             |             | (1993)      |
   |             |             |             |             | [267, 193]; |
   |             |             |             |             | Winkel &    |
   |             |             |             |             | Rambal      |
   |             |             |             |             | (1990)      |
   |             |             |             |             | [264, 336,  |
   |             |             |             |             | 216];       |
   |             |             |             |             | Correia et  |
   |             |             |             |             | al.(1995)   |
   |             |             |             |             | [276];      |
   |             |             |             |             | Jacobs et   |
   |             |             |             |             | al.(1996)   |
   |             |             |             |             | [188];      |
   |             |             |             |             | Massman et  |
   |             |             |             |             | al.(1994)   |
   |             |             |             |             | [315]       |
   +-------------+-------------+-------------+-------------+-------------+
   | **Root      |             | 360         | 0.02        | Simpson et  |
   | crops**     |             |             |             | al. (2003)  |
   |             |             |             |             |             |
   | **(RC)**    |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   | Potato      | All Europe  | 750         | 0.01        | UNECE       |
   |             |             |             |             | (2004);Jeff |
   | (Solanuum   |             |             |             | ries        |
   | tuberosum)  |             |             |             | (1994); Vos |
   |             |             |             |             | & Groenwald |
   |             |             |             |             | (1989);     |
   |             |             |             |             | Marshall &  |
   |             |             |             |             | Vos (1991); |
   |             |             |             |             | Pleijel et  |
   |             |             |             |             | al. (2002); |
   |             |             |             |             | Danielsson  |
   |             |             |             |             | (2003)      |
   |             |             |             |             |             |
   |             |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   | **Semi-Natu |             | 60          | 0.01        | Simpson et  |
   | ral         |             |             |             | al. (2003)  |
   | /           |             |             |             |             |
   | Moorland**  |             |             |             |             |
   |             |             |             |             |             |
   | **(SNL)**   |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   | **Grassland |             | 270         | 0.01        | Simpson et  |
   | **          |             |             |             | al. (2003)  |
   |             |             |             |             |             |
   | **(GR)**    |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   | Perennial   | All Europe  | 295         | 0.02        | ICP         |
   | rye grass   |             |             |             | Vegetation  |
   |             |             |             |             | contract    |
   | (Lolium     |             |             |             | report      |
   | perenne)    |             |             |             | (2009);     |
   |             |             |             |             | Sheehy et   |
   |             |             |             |             | al. (1975); |
   |             |             |             |             | Gay (1986); |
   |             |             |             |             | Nijs et al. |
   |             |             |             |             | (1989);     |
   |             |             |             |             | Ferris et   |
   |             |             |             |             | al. (1996); |
   |             |             |             |             | Jones et    |
   |             |             |             |             | al.         |
   |             |             |             |             | (1996) ;Coy |
   |             |             |             |             | le          |
   |             |             |             |             | (personal   |
   |             |             |             |             | communicati |
   |             |             |             |             | on);        |
   |             |             |             |             | Mills &     |
   |             |             |             |             | Hayes       |
   |             |             |             |             | (pers.      |
   |             |             |             |             | comm.)      |
   |             |             |             |             |             |
   |             |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   | Clover      | All Europe  | 360         | 0.02        | ICP         |
   |             |             |             |             | Vegetation  |
   | (Trifolium  |             |             |             | contract    |
   | repens)     |             |             |             | report      |
   |             |             |             |             | (2009);     |
   |             |             |             |             | Degl'Innoce |
   |             |             |             |             | nti         |
   |             |             |             |             | et          |
   |             |             |             |             | al.(2003);  |
   |             |             |             |             | Nussbaum,   |
   |             |             |             |             | and Fuhrer  |
   |             |             |             |             | (2000)      |
   +-------------+-------------+-------------+-------------+-------------+
   | **Mediterra |             | 213         | 0.014       | Simpson et  |
   | nean        |             |             |             | al. (2003)  |
   | scrub**     |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+

    

    is given in mmol O\ 3 m\ -2 PLA s\ -1. PLA is projected leaf area
   (m2). Species g\ max is the median of all values rounded to the
   nearest 5 mmol O\ 3 m\ -2 PLA s\ -1. A/H describes whether the
   stomata are found over the entire leaf or (A=amphistomatous) or one
   side (H=hypostomatous). The conversion of total to projected leaf
   area for needles uses a factor of 2.6. The median of species values
   are rounded to the nearest multiple of 5; values given in square
   brackets are the study specific maximum g\ s values in mmol O\ 3
   m\ -2 PLA s\ -1. The mean cover type  values are rounded to the
   nearest multiple of 10. \*gs data given for the flag leaf. † Q. ilex
   is chosen to represent the climax vegetation type of this broad cover
   type. Cover type f\ min values are assigned according to the lowest
   species f\ min value within that group, key reference used to
   establish this value is provided though other references may exist
   that support the value selected.

.. container:: WordSection3

   .. rubric:: Irradiance (flight)
      :name: irradiance-flight

   The ‘big leaf’ DO\ 3\ SE model uses a one-layer “two big-leaf”
   approach to calculate the stomatal conductance. This assumes that the
   canopy can be represented as a single ‘big-leaf’ leaf split into
   shaded and sun lit fractions. This model follows the accepted
   principle that radiation attenuation through canopies follows Beer’s
   law and that such radiation penetration must separately consider
   direct and diffuse radiation, due to their different attenuation
   through canopies, and visible and near infra-red wavelengths due to
   differential absorptance by leaves(Goudriaan, 1977). The DO\ 3\ SE
   model uses the radiation attenuation method ofWeiss & Norman (1985)to
   estimate the quantity and quality of radiation distribution through
   the canopy. Application of this modelrequires that the
   Photosynthetically Active Radiation (PAR), which has a wavelength of
   400 to 700 nm, reaching the top of the canopy has to be
   differentiated into direct and diffuse irradiance. This is dependent
   upon transmission through the atmosphere which is a function of
   optical air mass (om), atmospheric pressure, (P) and hence altitude.
   If atmospheric pressure is not recorded, it can be estimated
   according to eq.6.

   .. _Ref393799528:

   6

   where P\ o is the atmospheric pressure at sea level (101.325 kPa), g
   is the acceleration due to gravity (9.81 m s\ -2), h\ site is the
   site altitude in m and H\ atm is the atmospheric scale height which
   is equal to 7400 m.

   The value of om can be calculated according toWeiss & Norman,
   (1985)as in eq.7

   .. _Ref393799709:

   7

    

   where sinβ is the solar elevation. This parameter varies over the
   course of the day as a function of latitude and day length as
   described in eq.8, this eq. and the other solar geometry equations
   required for its calculation are taken fromCampbell & Norman, (1998).

   .. _Ref393808678:

   8

    

   where β is the solar elevation above the horizontal, λ is the
   latitude, δ is the angle between the sun’s rays and the equatorial
   plane of the earth (solar declination), hr is the hour angle of the
   sun and is given by  wheret is time and t\ o is the time at solar
   noon.

   The solar declination (δ) is calculated according to eq.9.

    

   .. _Ref393811007:

   9

                                                                                                  

   where t\ d is the year day.

   The time, t is in hours (standard local time), ranging from 0 to 23.
   Solar noon (to) varies during the year by an amount that is given by
   the equation of time (e, in min) and calculated by:-

   10

    

   where LC is the longitude correction. LC is + 4 or –4 minutes for
   each degree you are either east or west of the standard meridian. e
   is a 15 to 20 minute correction which depends on year day according
   to eq.11.

    

   .. _Ref393813337:

   11

    

   where f = 279.575 + 0.9856 t\ d in degrees.

   It is also necessary to calculate the day length so that the hour
   angle of the sun can be calculated throughout the day. Day length is
   defined as the number of hours that the sun is above the horizon and
   requires the hour angle of the sun, hr, at sunrise or sunset to be
   calculated with eq.12.

    

   .. _Ref393813427:

   12

    

    

   so that day length in hours equals 2\ hr/15.

   The formulations ofWeiss & Norman, (1985)can then be used to estimate
   the potential direct (pPARdir) (eq.13) and diffuse (pPARdiff) (eq.14)
   irradiances in W/m\ 2 which are necessary to estimate irradiance
   penetration and quality (whether direct or diffuse) into the canopy.

    

   .. _Ref393813597:

   13

    

   where the 600 (W m\ -2) represents the average amount of PAR
   radiation available at the top of the atmosphere, estimated according
   to the solar constant (1320 W m\ -2), of which 0.45 is the PAR
   fraction(Jones, 1992)and 0.185 represents the extinction coefficient.

    

   .. _Ref393813868:

   14

    

   where the term in brackets represents the total available PAR diffuse
   radiation. 0.4 is the fraction of intercepted PAR beam radiation that
   is converted to downward diffuse radiation at the surface. The
   potential total PAR beam radiation (pPARtotal) is then the sum of
   pPAR\ dir and pPAR\ diff. The actual total PAR (PAR\ total) is
   measured at each site (or provided as modelled values). To allow for
   variability in the calibration of the measurement apparatus or
   modelled data we estimate the sky transmissivity (ST) during the
   daylight period as in eq.15where the pPAR\ total is not allowed to
   exceed the PAR\ total.

   .. _Ref393814129:

   15

    

   *ST*\ is confined within 0.9 and 0.21 to deal with situations when
   the zenith angle may be greater than 80\ o (i.e. at sunrise and
   sunset). Estimation of this value allows the fraction of the direct
   and diffuse components of the radiation beam to be defined during the
   daylight periods using eqs.16and17.

    

   .. _Ref393814315:

   16

    

   .. _Ref393814370:

   17

    

   The actual PAR\ dir and PAR\ diff can then simply be calculated by
   multiplying the respective fPAR with the actual total PAR (PARtotal).

   Estimations of the diffuse and direct irradiance fractions are
   necessary to calculate the PAR incident on the sunlit (LAIsun) (see
   eq.18) and shaded (LAIshade) (see eq.19) portions of the canopy.

    

   .. _Ref393814592:

   18

    

   .. _Ref393814676:

   19

    

              

   PAR\ sun, which is dependent on the mean angle between leaves and the
   sun, is calculated using a modified “big leaf” version of the canopy
   radiation transfer model of(Zhang et al., 2001)where the flux density
   of PAR on sunlit leaves is calculated as described in eq.20

    

   .. _Ref393814942:

   20

    

    

   Where PAR\ dir is the actual direct PAR above the canopy (as
   calculated previously) and\ *q*\ is the angle between a leaf and the
   sun. For these calculations it is assumed that the canopy has a
   spherical leaf inclination distribution (\ *q*\ ) constant at 60
   degrees.

    

   PAR\ shadeis calculated semi-empirically using eq.21

    

   .. _Ref393815151:

   21

    

    

   where PAR\ diff is the actual diffuse PAR above the canopy.

    

   These values are then used to estimate the mean canopy relative
   stomatal conductance (Flight) as a function of irradiance using the
   f\ light function in eq.22 according to the proportions of sunlit and
   shaded leaf area.

   .. _Ref393815299:

   22

    

   Where PPFD represents the photosynthetic photon flux density in units
   of μmol m\ -2 s\ -1 (conversion from PAR in units of W m\ -2 to PPFD
   in units of μmol m\ -2 s\ -1 can be achieved using a conversion
   factor of 4.57 afterJones (1992), this provides estimates of PPFDsun
   and PPFDshade from values provided in eqs.20and21respectively. The
   above calculations give estimates off\ lightfor individual leaves or
   needles at the top of the canopy.

   To estimate canopyG\ sto, it is necessary to allow for the variable
   penetration of irradiance into the canopy which can be achieved by
   scalingf\ lightaccording to the sunlit and shaded\ *LAI*\ fractions
   of the canopy (see eqs23and24). Sunlit (flight\ sun) and shaded
   (flight\ shade) leaves are then weighted according to the fraction of
   sunlit and shaded LAI and summed to give the total
   irradiance-dependant canopy conductance (Flight) as described in
   eq.25. This method has been simplified so that rather than
   integrating F\ light over the canopy, a “big leaf” approach has been
   employed. This means that F\ light will be underestimated since the
   model does not allow for the variability in direct and diffuse
   irradiance within the canopy. However, since these parameters vary
   only slightly the model under-estimations are relatively small.

    

   .. _Ref393815894:

   23

    

    

   .. _Ref393815899:

   24

    

   F\ lightis calculated as : -

    

   .. _Ref393816051:

   25

    

    

   This mean canopy F\ light value is used to estimate mean canopy
   leaf/needle G\ sto. The canopy G\ sto is then calculated in equation
   using the cover-type specific LAI to upscale from the leaf/needle to
   the canopy level.

   ** **

   .. rubric:: Temperature (ftemp)
      :name: temperature-ftemp

    

   Add text

    

   f\ T = max{f\ min, (T-Tmin) / (Topt – Tmin)*[(Tmax-T) /
   (Tmax-Topt)]bt}

   where bt = (Tmax-Topt)/(Topt-Tmin)

    

.. container:: WordSection4

   Table (x) Default deposition land-cover and species class parameters
   for flight (\ **a) and ftemp (Tmin, T\ opt and T\ max)
   calculations.**\ `[p7] <#msocom-7>`__\  

    

   +---------+---------+---------+---------+---------+---------+---------+
   | **Land- | **Clima | **Fligh | **T\ mi | **T\ op | **T\ ma | **Refer |
   | cover   | te      | t       | n**     | t**     | x**     | ence**  |
   | type &  | region* | factor  |         |         |         |         |
   | Species | *       | (a)**   |         |         |         |         |
   | **      |         |         |         |         |         |         |
   +---------+---------+---------+---------+---------+---------+---------+
   | ** **   |         |         |         |         |         |         |
   +---------+---------+---------+---------+---------+---------+---------+
   | **Conif |         | 0.0083  | 1       | 18      | 36      | Simpson |
   | erous   |         |         |         |         |         | et al.  |
   | Forests |         |         |         |         |         | (2003)  |
   | (CF)**  |         |         |         |         |         |         |
   +---------+---------+---------+---------+---------+---------+---------+
   | Norway  | Norther | 0.006   | 0       | 20      | 200*\*  | UNECE   |
   | spruce  | n       |         |         |         |         | (2004); |
   |         | Europe  |         |         |         |         | Karlsso |
   | (Picea  |         |         |         |         |         | n       |
   | abies)  |         |         |         |         |         | et al.  |
   |         |         |         |         |         |         | (2000); |
   |         |         |         |         |         |         | Hansson |
   |         |         |         |         |         |         | et      |
   |         |         |         |         |         |         | al.(in  |
   |         |         |         |         |         |         | prep)   |
   +---------+---------+---------+---------+---------+---------+---------+
   | Scots   | Atlanti | 0.006   | 0       | 20      | 36      | UNECE   |
   | Pine    | c       |         |         |         |         | (2004); |
   |         | Central |         |         |         |         | Beadle  |
   | (Pinus  | Europe  |         |         |         |         | et al.  |
   | sylvest |         |         |         |         |         | (1985); |
   | ris)    |         |         |         |         |         | Sturm   |
   |         |         |         |         |         |         | et      |
   |         |         |         |         |         |         | al.(199 |
   |         |         |         |         |         |         | 8);     |
   |         |         |         |         |         |         | Ng.     |
   |         |         |         |         |         |         | (1979)  |
   +---------+---------+---------+---------+---------+---------+---------+
   | Norway  | Contine | 0.01    | 0       | 14      | 35      | UNECE   |
   | Spruce  | ntal    |         |         |         |         | (2004); |
   |         | Central |         |         |         |         | Thoene  |
   | (Picea  | Europe  |         |         |         |         | et al.  |
   | abies)  |         |         |         |         |         | (1991), |
   |         |         |         |         |         |         | Korner  |
   |         |         |         |         |         |         | et al.  |
   |         |         |         |         |         |         | (1995), |
   |         |         |         |         |         |         | Zweifel |
   |         |         |         |         |         |         | et      |
   |         |         |         |         |         |         | al.(200 |
   |         |         |         |         |         |         | 0,      |
   |         |         |         |         |         |         | 2001,20 |
   |         |         |         |         |         |         | 02)     |
   |         |         |         |         |         |         |         |
   |         |         |         |         |         |         |         |
   +---------+---------+---------+---------+---------+---------+---------+
   | **Decid |         | 0.006   | 6       | 20      | 34      | Simpson |
   | uous    |         |         |         |         |         | et al.  |
   | Forests |         |         |         |         |         | (2003); |
   | **      |         |         |         |         |         |         |
   |         |         |         |         |         |         |         |
   | **(DF)* |         |         |         |         |         |         |
   | *       |         |         |         |         |         |         |
   +---------+---------+---------+---------+---------+---------+---------+
   | **Gener | **All   | **0.006 | **0**   | **21**  | **35**  | **UNECE |
   | ic      | Europe* | **      |         |         |         | (2004)* |
   | Deciduo | *       |         |         |         |         | *       |
   | us**    |         |         |         |         |         |         |
   +---------+---------+---------+---------+---------+---------+---------+
   | Silver  | Norther | 0.0042  | 5       | 20      | 200*\*  | UNECE   |
   | birch   | n       |         |         |         |         | (2004); |
   |         | Europe  | (0.006) |         |         |         | Uddling |
   | (Betula |         |         |         |         |         | et al.  |
   | pendula |         |         |         |         |         | (2005a) |
   | )       |         |         |         |         |         |         |
   |         |         |         |         |         |         | Osonubi |
   |         |         |         |         |         |         | &       |
   |         |         |         |         |         |         | Davies. |
   |         |         |         |         |         |         | (1980); |
   |         |         |         |         |         |         | Oksanen |
   |         |         |         |         |         |         | .       |
   |         |         |         |         |         |         | (pers.  |
   |         |         |         |         |         |         | comm.)( |
   |         |         |         |         |         |         | 2003)   |
   +---------+---------+---------+---------+---------+---------+---------+
   | Beech   | Atlanti | 0.006   | 0       | 21      | 35      | UNECE   |
   |         | c       |         |         |         |         | (2004)  |
   | (Fagus  | Central |         |         |         |         |         |
   | sylvati | Europe  |         |         |         |         |         |
   | ca)     |         |         |         |         |         |         |
   +---------+---------+---------+---------+---------+---------+---------+
   | Oak     | Atlanti | 0.003   | 0       | 20      | 35      | UNECE   |
   |         | c       |         |         |         |         | (2004); |
   | (Quercu | Central |         |         |         |         | Breda   |
   | s       | Europe  |         |         |         |         | et al.  |
   | petraea |         |         |         |         |         | (1993b) |
   | &       |         |         |         |         |         | ;       |
   | robur)  |         |         |         |         |         | Dolman. |
   |         |         |         |         |         |         | (1988)  |
   +---------+---------+---------+---------+---------+---------+---------+
   | Beech   | Contine | 0.006   | 5       | 16      | 33      | UNECE   |
   |         | ntal    |         |         |         |         | (2004); |
   | (Fagus  | Central |         |         |         |         | Braun   |
   | sylvati | Europea |         |         |         |         | et al.  |
   | ca)     | n       |         |         |         |         | (in     |
   |         |         |         |         |         |         | prep)   |
   +---------+---------+---------+---------+---------+---------+---------+
   | Beech   | Mediter | 0.006   | 4       | 21      | 37      | UNECE   |
   |         | ranean  |         |         |         |         | (2004); |
   | (Fagus  | Europe  |         |         |         |         | Damesin |
   | sylvati |         |         |         |         |         | et al.  |
   | ca)     |         |         |         |         |         | (1998); |
   |         |         |         |         |         |         | Rico et |
   |         |         |         |         |         |         | al.     |
   |         |         |         |         |         |         | (1996)  |
   +---------+---------+---------+---------+---------+---------+---------+
   | **Needl |         | 0.013   | 4       | 20      | 37      | Simpson |
   | eleaf   |         |         |         |         |         | et al.  |
   | Forests |         |         |         |         |         | (2003)  |
   | **      |         |         |         |         |         |         |
   |         |         |         |         |         |         |         |
   | **(NF)* |         |         |         |         |         |         |
   | *       |         |         |         |         |         |         |
   +---------+---------+---------+---------+---------+---------+---------+
   | Aleppo  | Mediter | 0.013   | 10      | 27      | 38      | UNECE   |
   | Pine    | ranean  |         |         |         |         | (2004); |
   |         | Europe  |         |         |         |         | Elvira  |
   | *(Pinus |         |         |         |         |         | et al.  |
   | halepen |         |         |         |         |         | (2007)  |
   | sis)*   |         |         |         |         |         |         |
   +---------+---------+---------+---------+---------+---------+---------+
   | **Broad |         | 0.009   | 4       | 20      | 37      | Simpson |
   | leaf    |         |         |         |         |         | et al.  |
   | Forests |         |         |         |         |         | (2003)  |
   | **      |         |         |         |         |         |         |
   |         |         |         |         |         |         |         |
   | **(BF)* |         |         |         |         |         |         |
   | *       |         |         |         |         |         |         |
   +---------+---------+---------+---------+---------+---------+---------+
   | **Gener | **All   | **0.009 | **2**   | **23**  | **38**  | **UNECE |
   | ic      | Europe* | **      |         |         |         | (2004)* |
   | Evergre | *       |         |         |         |         | *       |
   | en      |         |         |         |         |         |         |
   | Mediter |         |         |         |         |         |         |
   | ranean* |         |         |         |         |         |         |
   | *       |         |         |         |         |         |         |
   +---------+---------+---------+---------+---------+---------+---------+
   | Holm    | Mediter | 0.012   | 1       | 23      | 39      | These   |
   | Oak     | ranean  |         |         |         |         | refs    |
   |         | Europe  |         |         |         |         | are in  |
   | (Quercu |         |         |         |         |         | the     |
   | s       |         |         |         |         |         | hard    |
   | ilex)   |         |         |         |         |         | copy    |
   |         |         |         |         |         |         | version |
   +---------+---------+---------+---------+---------+---------+---------+
   | **Tempe |         | 0.009   | 12      | 26      | 40      | Simpson |
   | rate    |         |         |         |         |         | et al.  |
   | crops** |         |         |         |         |         | (2003)  |
   |         |         |         |         |         |         |         |
   | **(TC)* |         |         |         |         |         |         |
   | *       |         |         |         |         |         |         |
   +---------+---------+---------+---------+---------+---------+---------+
   | **Gener | **All   | **0.010 | **12**  | **26**  | **40**  | **UNECE |
   | ic      | Europe* | 5**     |         |         |         | (2004)* |
   | crop**  | *       |         |         |         |         | *       |
   +---------+---------+---------+---------+---------+---------+---------+
   | Wheat   | All     | 0.0105  | 12      | 26      | 40      | UNECE   |
   |         | Europe  |         |         |         |         | (2004); |
   | (Tritic |         |         |         |         |         | Gruters |
   | um      |         |         |         |         |         | et al.  |
   | aestivu |         |         |         |         |         | (1995); |
   | m)      |         |         |         |         |         | Bunce   |
   |         |         |         |         |         |         | (2000)  |
   |         |         |         |         |         |         |         |
   |         |         |         |         |         |         |         |
   +---------+---------+---------+---------+---------+---------+---------+
   | **Medit |         | 0.0048  | 0       | 25      | 51      | Simpson |
   | erranea |         |         |         |         |         | et al.  |
   | n       |         |         |         |         |         | (2003)  |
   | crops** |         |         |         |         |         |         |
   |         |         |         |         |         |         |         |
   | **(MC)* |         |         |         |         |         |         |
   | *       |         |         |         |         |         |         |
   +---------+---------+---------+---------+---------+---------+---------+
   | Maize   | All     | 0.0048  | 2       | 25      | 48      | ICP     |
   |         | Europe  |         |         |         |         | Vegetat |
   | *(Zea   |         |         |         |         |         | ion     |
   | mays)*  |         |         |         |         |         | contrac |
   |         |         |         |         |         |         | t       |
   |         |         |         |         |         |         | report  |
   |         |         |         |         |         |         | (2006); |
   |         |         |         |         |         |         | Betheno |
   |         |         |         |         |         |         | d       |
   |         |         |         |         |         |         | &       |
   |         |         |         |         |         |         | Tardieu |
   |         |         |         |         |         |         | (1990); |
   |         |         |         |         |         |         | Turner  |
   |         |         |         |         |         |         | & Begg  |
   |         |         |         |         |         |         | (1973); |
   |         |         |         |         |         |         | Rochett |
   |         |         |         |         |         |         | e       |
   |         |         |         |         |         |         | et al.  |
   |         |         |         |         |         |         | (1991); |
   |         |         |         |         |         |         | Machado |
   |         |         |         |         |         |         | & Lagoa |
   |         |         |         |         |         |         | (1994); |
   |         |         |         |         |         |         | Guilion |
   |         |         |         |         |         |         | i       |
   |         |         |         |         |         |         | et al.  |
   |         |         |         |         |         |         | (2000); |
   |         |         |         |         |         |         | Olioso  |
   |         |         |         |         |         |         | et al.  |
   |         |         |         |         |         |         | (1995); |
   |         |         |         |         |         |         | Ozier-L |
   |         |         |         |         |         |         | afontei |
   |         |         |         |         |         |         | n       |
   |         |         |         |         |         |         | et al.  |
   |         |         |         |         |         |         | (1998); |
   |         |         |         |         |         |         | Rodrigu |
   |         |         |         |         |         |         | ez      |
   |         |         |         |         |         |         | &       |
   |         |         |         |         |         |         | Davies  |
   |         |         |         |         |         |         | (1982)  |
   |         |         |         |         |         |         |         |
   |         |         |         |         |         |         |         |
   +---------+---------+---------+---------+---------+---------+---------+
   | Sunflow | All     | 0.002   | 2       | 25      | 48      | ICP     |
   | er      | Europe  |         |         |         |         | Vegetat |
   |         |         |         |         |         |         | ion     |
   | *(Helia |         |         |         |         |         | contrac |
   | nthus   |         |         |         |         |         | t       |
   | annuus) |         |         |         |         |         | report  |
   | *       |         |         |         |         |         | (2006); |
   |         |         |         |         |         |         | Turner  |
   |         |         |         |         |         |         | (1970); |
   |         |         |         |         |         |         | Fay &   |
   |         |         |         |         |         |         | Knapp   |
   |         |         |         |         |         |         | (1996). |
   |         |         |         |         |         |         |         |
   |         |         |         |         |         |         | For     |
   |         |         |         |         |         |         | t_opt,  |
   |         |         |         |         |         |         | t_min   |
   |         |         |         |         |         |         | and     |
   |         |         |         |         |         |         | t_maxma |
   |         |         |         |         |         |         | ize     |
   |         |         |         |         |         |         | paramet |
   |         |         |         |         |         |         | erisati |
   |         |         |         |         |         |         | on      |
   |         |         |         |         |         |         | used as |
   |         |         |         |         |         |         | default |
   |         |         |         |         |         |         |         |
   |         |         |         |         |         |         |         |
   +---------+---------+---------+---------+---------+---------+---------+
   | Tomato  | All     | 0.0175  | 0       | 21      | 35      | ICP     |
   |         | Europe  |         |         |         |         | Vegetat |
   | *(Solan |         |         |         |         |         | ion     |
   | um      |         |         |         |         |         | contrac |
   | lycoper |         |         |         |         |         | t       |
   | sicum)* |         |         |         |         |         | report  |
   |         |         |         |         |         |         | (2006); |
   |         |         |         |         |         |         | Boulard |
   |         |         |         |         |         |         | et al.  |
   |         |         |         |         |         |         | (1991); |
   |         |         |         |         |         |         | Bakker  |
   |         |         |         |         |         |         | (1991); |
   |         |         |         |         |         |         | Starck  |
   |         |         |         |         |         |         | et al.  |
   |         |         |         |         |         |         | (2000). |
   |         |         |         |         |         |         |         |
   |         |         |         |         |         |         |         |
   +---------+---------+---------+---------+---------+---------+---------+
   | Grape   | All     | 0.0076  | 9       | 30      | 43      | ICP     |
   | vine    | Europe  |         |         |         |         | Vegetat |
   |         |         |         |         |         |         | ion     |
   | *(Vitis |         |         |         |         |         | contrac |
   | vinifer |         |         |         |         |         | t       |
   | a)*     |         |         |         |         |         | report  |
   |         |         |         |         |         |         | (2006); |
   |         |         |         |         |         |         | Schultz |
   |         |         |         |         |         |         | (2003a) |
   |         |         |         |         |         |         | ;       |
   |         |         |         |         |         |         | Winkel  |
   |         |         |         |         |         |         | &       |
   |         |         |         |         |         |         | Rambal  |
   |         |         |         |         |         |         | (1993); |
   |         |         |         |         |         |         | Winkel  |
   |         |         |         |         |         |         | &       |
   |         |         |         |         |         |         | Rambal  |
   |         |         |         |         |         |         | (1990); |
   |         |         |         |         |         |         | Massman |
   |         |         |         |         |         |         | et al.  |
   |         |         |         |         |         |         | (1994); |
   |         |         |         |         |         |         | Jacobs  |
   |         |         |         |         |         |         | et al.  |
   |         |         |         |         |         |         | (1996); |
   |         |         |         |         |         |         | Correia |
   |         |         |         |         |         |         | et al   |
   |         |         |         |         |         |         | (1995); |
   |         |         |         |         |         |         | Flexas  |
   |         |         |         |         |         |         | et      |
   |         |         |         |         |         |         | al.(199 |
   |         |         |         |         |         |         | 9);     |
   |         |         |         |         |         |         | Schultz |
   |         |         |         |         |         |         | (2003); |
   |         |         |         |         |         |         | Shultz  |
   |         |         |         |         |         |         | (2003a) |
   |         |         |         |         |         |         | ;       |
   |         |         |         |         |         |         | Massman |
   |         |         |         |         |         |         | et al   |
   |         |         |         |         |         |         | (2003); |
   |         |         |         |         |         |         | Schultz |
   |         |         |         |         |         |         | (2003)  |
   |         |         |         |         |         |         |         |
   |         |         |         |         |         |         |         |
   |         |         |         |         |         |         |         |
   |         |         |         |         |         |         |         |
   +---------+---------+---------+---------+---------+---------+---------+
   | **Root  |         | 0.0023  | 8       | 24      | 50      | Simpson |
   | crops** |         |         |         |         |         | et al.  |
   |         |         |         |         |         |         | (2003)  |
   | **(RC)* |         |         |         |         |         |         |
   | *       |         |         |         |         |         |         |
   +---------+---------+---------+---------+---------+---------+---------+
   | Potato  | All     | 0.005   | 13      | 28      | 39      | UNECE   |
   |         | Europe  |         |         |         |         | (2004); |
   | (Solanu |         |         |         |         |         | Ku et   |
   | um      |         |         |         |         |         | al.     |
   | tuberos |         |         |         |         |         | (1977)  |
   | um)     |         |         |         |         |         | ;       |
   |         |         |         |         |         |         | Dwelle  |
   |         |         |         |         |         |         | et al   |
   |         |         |         |         |         |         | (1981)  |
   +---------+---------+---------+---------+---------+---------+---------+
   | **Semi- |         | 0.009   | 1       | 18      | 36      | Simpson |
   | Natural |         |         |         |         |         | et al.  |
   | /       |         |         |         |         |         | (2003)  |
   | Moorlan |         |         |         |         |         |         |
   | d**     |         |         |         |         |         |         |
   |         |         |         |         |         |         |         |
   | **(SNL) |         |         |         |         |         |         |
   | **      |         |         |         |         |         |         |
   +---------+---------+---------+---------+---------+---------+---------+
   | **Grass |         | 0.009   | 12      | 26      | 40      | Simpson |
   | land**  |         |         |         |         |         | et al.  |
   |         |         |         |         |         |         | (2003)  |
   | **(GR)* |         |         |         |         |         |         |
   | *       |         |         |         |         |         |         |
   +---------+---------+---------+---------+---------+---------+---------+
   | Perenni | All     | 0.007   | 10      | 25      | 40      | ICP     |
   | al      | Europe  |         |         |         |         | Vegetat |
   | rye     |         |         |         |         |         | ion     |
   | grass   |         |         |         |         |         | contrac |
   |         |         |         |         |         |         | t       |
   | (Lolium |         |         |         |         |         | report  |
   | perenne |         |         |         |         |         | (2009)  |
   | )       |         |         |         |         |         |         |
   +---------+---------+---------+---------+---------+---------+---------+
   | Clover  | All     | 0.008   | 10      | 27      | 43      | ICP     |
   |         | Europe  |         |         |         |         | Vegetat |
   | (Trifol |         |         |         |         |         | ion     |
   | ium     |         |         |         |         |         | contrac |
   | repens) |         |         |         |         |         | t       |
   |         |         |         |         |         |         | report  |
   |         |         |         |         |         |         | (2009)  |
   +---------+---------+---------+---------+---------+---------+---------+
   | **Medit |         | 0.012   | 4       | 20      | 37      | Simpson |
   | erranea |         |         |         |         |         | et al.  |
   | n       |         |         |         |         |         | (2003)  |
   | scrub** |         |         |         |         |         |         |
   +---------+---------+---------+---------+---------+---------+---------+

.. container:: WordSection5

   .. rubric:: Vapour pressure deficit (fVPD)
      :name: vapour-pressure-deficit-fvpd

    

    

   Add text

   f\ VPD=max{f\ min, min{1, (1-fmin)*(VPDmin-VPD)/(VPDmin-VPDmax) +
   f\ min}}

    

    

   For wheat, potato and the generic crop
   there\ \ `[p8] <#msocom-8>`__\ \   is another effect on stomata by
   water relations which can be modelled using VPD. During the
   afternoon, the air temperature typically decreases, which is
   normally, but not always, followed (if the absolute humidity of the
   air remains constant or increases) by declining VPD. According to the
   fVPD function this would allow the stomata to re-open if there had
   been a limitation by fVPD earlier during the day. Most commonly this
   does not happen. This is related to the fact that during the day the
   plant loses water through transpiration at a faster rate than it is
   replaced by root uptake. This results in a reduction of the plant
   water potential during the course of the day and prevents stomata
   re-opening in the afternoon. The plant water potential then recovers
   during the following night when the rate of transpiration is low. A
   simple way to model the extent of water loss by the plant is to use
   the sum of hourly VPD values during the daylight hours (as suggested
   by Uddling et al., 2004). If there is a large sum it is likely to be
   related to a larger amount of transpiration, and if the accumulated
   amount of transpiration during the course of the day (as represented
   by a VPD sum) exceeds a certain value, then stomatal re-opening in
   the afternoon does not occur. This is represented by the VPDsum
   function (ΣVPD) which is calculated in the following  manner:

    

   If ΣVPD ≥ ΣVPD_crit, then gsto_hour_n+1 ≤ gsto_hour_n

    

   Where g\ sto_hour_n and g\ sto_hour_n+1 are the g\ sto values for
   hour n and hour n+1 respectively calculated according to the g\ sto
   equation.

    

    

    

   .. rubric:: Soil Water Potential (fSWP)
      :name: soil-water-potential-fswp

    

    

.. container::

   --------------

   .. container::

      .. container:: msocomtxt

         .. _msocom-1:

         \_msocom-1

          \ \ `[p1] <#msoanchor-1>`__\ I wonder whether a separate text
         box or schematic would be good here to show the difference
         between gsto, Gsto and Gstocanopy

   .. container::

      .. container:: msocomtxt

         .. _msocom-2:

         \_msocom-2

          \ \ `[p2] <#msoanchor-2>`__\ repetition

   .. container::

      .. container:: msocomtxt

         .. _msocom-3:

         \_msocom-3

          \ \ `[p3] <#msoanchor-3>`__\ repetition

   .. container::

      .. container:: msocomtxt

         .. _msocom-4:

         \_msocom-4

          \ \ `[p4] <#msoanchor-4>`__\ what is “normal”?

   .. container::

      .. container:: msocomtxt

         .. _msocom-5:

         \_msocom-5

          \ \ `[p5] <#msoanchor-5>`__\ update – where is this from? Do
         we need this?

   .. container::

      .. container:: msocomtxt

         .. _msocom-6:

         \_msocom-6

          \ \ `[B6] <#msoanchor-6>`__\ Give stdev. rather than range
         here?

   .. container::

      .. container:: msocomtxt

         .. _msocom-7:

         \_msocom-7

          \ \ `[p7] <#msoanchor-7>`__\ Check and update (add units) –
         task for Richard?

   .. container::

      .. container:: msocomtxt

         .. _msocom-8:

         \_msocom-8

          \ \ `[p8] <#msoanchor-8>`__\ How about trees? This links to
         fSWP part of documentation
