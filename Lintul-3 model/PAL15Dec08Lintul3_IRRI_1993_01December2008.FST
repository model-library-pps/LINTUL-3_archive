*-------------------------------------------------------------------------*
* Copyright 2013. Wageningen University, Plant Production Systems group,  *
* P.O. Box 430, 6700 AK Wageningen, The Netherlands.                      *
* You may not use this work except in compliance with the Licence.        *
* You may obtain a copy of the Licence at:                                *
*                                                                         *
* http://models.pps.wur.nl/content/licence-agreement                      *
*                                                                         *
* Unless required by applicable law or agreed to in writing, software     *
* distributed under the Licence is distributed on an "AS IS" basis,       *
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
*-------------------------------------------------------------------------*

************************************************************************
*   LINTUL3 is an extended version of LINTUL1 (the version of LINTUL   *
*          for optimal growth conditions) and LINTUL2 includes a simple*
*          water balance for studying effects of drought. LINTUL2-N    *
*          includes N-limitation on crop growth. The latter program is *
*          called LINTUL3.                                             *
***********************Defining Subroutines*****************************

DEFINE_CALL ASTRO(INPUT,INPUT,                  OUTPUT)

DEFINE_CALL SUBDVS(INPUT,INPUT,INPUT,INPUT,INPUT,  OUTPUT)

DEFINE_CALL PENMAN(INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,   ...
                   OUTPUT,OUTPUT,OUTPUT,OUTPUT,OUTPUT,OUTPUT)

DEFINE_CALL EVAPTR(INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT, ...
                   INPUT,INPUT,INPUT,INPUT,...
                   OUTPUT,OUTPUT,OUTPUT,OUTPUT,OUTPUT)

DEFINE_CALL DRUNIR(INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT, ...
                   INPUT,INPUT,INPUT,INPUT,  OUTPUT,OUTPUT,OUTPUT)

DEFINE_CALL GLA(INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,...
                INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,      OUTPUT)

DEFINE_CALL  SUBPAR (INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,...
                     OUTPUT,OUTPUT,OUTPUT,OUTPUT,OUTPUT,OUTPUT)

DEFINE_CALL GROWTH(INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,...
                    INPUT,                         OUTPUT,OUTPUT)

DEFINE_CALL RELGR(INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,...
                  INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,...
                                        OUTPUT,OUTPUT,OUTPUT,OUTPUT)
                  
DEFINE_CALL RNUSUB(INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,...
                   INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,...
                   OUTPUT,OUTPUT,OUTPUT)

DEFINE_CALL NOPTM(INPUT,INPUT,INPUT,    OUTPUT,OUTPUT)

DEFINE_CALL NDEMND(INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT, ...
                   INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,...
                   INPUT,INPUT,INPUT,INPUT,...
                   OUTPUT,OUTPUT,OUTPUT,OUTPUT)

DEFINE_CALL NTRLOC(INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,...
                   INPUT,INPUT,         OUTPUT,OUTPUT,OUTPUT,OUTPUT)

DEFINE_CALL NTRANS(INPUT,INPUT,INPUT,INPUT,INPUT,  OUTPUT,OUTPUT,OUTPUT)

DEFINE_CALL RNLD(INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,...
                                         OUTPUT,OUTPUT,OUTPUT)
                                         
DEFINE_CALL NNINDX(INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,    OUTPUT)

DEFINE_CALL DEATHL(INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,...
                   INPUT,INPUT,INPUT,INPUT,          OUTPUT,OUTPUT,...
                      OUTPUT,OUTPUT,OUTPUT,OUTPUT,OUTPUT,OUTPUT,OUTPUT)

************************************************************************
TITLE LINTUL2_N with adaptation to N-Limitation
*----------------------------------------------------------------------*
***   Initial conditions and run control.
INITIAL

*Initial conditions.
      INCON ZERO  = 0.

* Initial Temperature sum (TSUM) and initial developmental stage (DVS).
      INCON TSUMI = 362.
      DVSI        = TSUMI/TSUMAN

* Temperature sum for anthesis, and maturity
     PARAM TSUMAN = 1420.; TSUMMT= 580.

* In case of transplanted rice, TSUMI accounts for the TSUM from
* emergence to transplanting and the initial DVS is calculated using
* TSUMI. But in case of wheat TSUMI and DVSI are zero because it is
* directly sown in the field.

* Initial amount of water present in the rooted depth at the start of
* the calculations, based on the initial water content.
      WAI  = 1000. * ROOTDI * WCI

* Initial weight (g/m2) of leaves, stem, roots, and storage organs at
* transplanting.
      INCON WLVGI = 0.86
      INCON WSTI  = 0.71
      INCON WRTLI = 1.58
      INCON WSOI  = 0.00

*     Initial LAI.
      LAII   = WLVGI * ISLA
      ISLA   = SLAC * SLACFI
      SLACFI = AFGEN (SLACF,DVSI)
      
* Initial amount of N (g/m2) in leaves, stem, roots, and storage organs.
      ANLVI = NFRLVI* WLVGI
      ANSTI = NFRSTI* WSTI
      ANRTI = NFRRTI* WRTLI
      ANSOI = 0.

* Initial water content in cm3 of water/(cm3 of soil).
      PARAM WCI    = 0.60

* Initial root depth (m).
      INCON ROOTDI = 0.05
      
* Initial fraction of N (g N g-1 DM) in leaves, stem, and roots.
      PARAM NFRLVI = 0.027;  NFRSTI = 0.018; NFRRTI= 0.016
*---------------------------------------------------------------------*
*Run control.
FINISH DVS > 2.01
FINISH TSUM> TTSUM
TIMER STTIME = 25.; FINTIM = 200.; DELT = 1.0; PRDEL = 1.
TRANSLATION_GENERAL DRIVER='EUDRIV'

*----------------------------------------------------------------------*
PRINT DVS,TSUM,TAGBM,WST,WLVG,WLVD,WSO,LAI,NTAC,WRT,WSTOB,...
TAGBMO,WLVGOB,WSOOB,WLVDOB,LAIOB,GTSUM,CBALAN,...
TRANRF,NNI,SLA,FRACT,FRTWET,FLVT,FSTT,FSOT,RWLVG,RWST,RWRT,RWSO,DLV,...
CUMPAR,LUECAL
*----------------------------------------------------------------------*

DYNAMIC
WEATHER WTRDIR='C:\SYS\Weather\'; CNTR='brazil';ISTN=3; IYEAR=1975

*----------Emergence,Temperature sum and Developmental stages----------*

      DTR    = RDD/1.E+6
      DAVTMP = 0.5 * (TMMN + TMMX)
      DTEFF  = MAX ( 0., DAVTMP-TBASE )
      EMERG  = REAAND(TIME-DOYEM+1.,WC-WCWP)* INSW(-LAI,1.,0.)

      TSUM   = INTGRL(TSUMI, RTSUMP)
      RTSUM  = DTEFF * EMERG
      TTSUM  = TSUMAN + TSUMMT
      
* Photoperiodic effect.
      PHOTPF  = INSW (TSUM-TSUMAN, AFGEN (PHOTTB,DAYL),1.)
      RTSUMP  = RTSUM * PHOTPF

***********************Calling Subroutines******************************
* Calling the subroutine for calculating the astrological daylength.
CALL ASTRO(DOY, LAT,                    DAYL)

* Calling the subroutine for converting TSUM to the developmental stage.
CALL SUBDVS(TIME,DOYEM,TSUM,TSUMAN,TSUMMT,         DVS)

* Calling the subroutine for Potential evaporation and transpiration.
CALL PENMAN(DAVTMP,VP,DTR,LAI,WN,RNINTC,...
            RLWN,NRADC,PENMRC,PENMD,PEVAP,PTRAN)

* Calling the subroutine for actual rates of evaporation and
* transpiration.
CALL EVAPTR(PEVAP,PTRAN,ROOTD,WA,WCAD,WCWP,WCFC,WCWET,WCST,TRANCO,...
            DELT,WMFAC,          WCCR,EVAP,TRAN,FR,AVAILF)

* Calling the subroutine for rates of drainage, runoff and irrigation.
CALL DRUNIR(RAIN,RNINTC,EVAP,TRAN,IRRIGF,DRATE,DELT,WA,ROOTD,WCFC,...
            WCST,WMFAC,         DRAIN,RUNOFF,IRRIG)

* Calling the subroutine for daily increase of leaf area index.
CALL GLA(TIME,DOYEM,DTEFF,LAII,RGRL,DELT,SLA,LAI,GLV,NLAI,...
         WC,WCWP,DVS,TRANRF,NNI,               GLAI)

* Calling the subroutine for relative modification for root and shoot
* allocation.
CALL SUBPAR (NPART,TRANRF,NNI,FRTWET,FLVT,FSTT,FSOT,...
                     FSHMOD,FLVMOD,FRT,FLV,FST,FSO)

* Calling the subroutine for total growth rate.
CALL GROWTH(TIME,DOYEM,DTR,K,NLUE,LAI,LUE,TRANRF,NNI,   PARINT,GTOTAL)

* Calling the subroutine for relative growth rate of roots, leaves, stem
* and storage organs.
CALL RELGR(TIME,DOYEM,EMERG,WLVGI,WRTLI,WSTI,WSOI,GTOTAL,FLV,FRT,...
           FST,FSO,DLV,DRRT,DELT,    RWLVG,RWRT,RWST,RWSO)
             
* Calling the subroutine to compute the partitioning of the total
* N uptake rate (NUPTR) over the leaves, stem and roots.
CALL RNUSUB(TIME,DOYEM,EMERG,NDEML,NDEMS,NDEMR,NUPTR,...
            NDEMTO,ANLVI,ANSTI,ANRTI,DELT,...
            RNULV,RNUST,RNURT)

* Calling the subroutine for maximum nitrogen concentration of leaves
* and stem.
CALL NOPTM(FRNX,NMAXLV,NMAXST, NOPTLV,NOPTST)

* Calling the subroutine for N demand of leaves, roots and stem storage
* organs.
CALL NDEMND(NMAXLV,NMAXST,NMAXRT,NMAXSO,WLVG,WST,WRT,WSO,...
            RWLVG,RWST,RWRT,RWSO, ANLV,ANST,...
            ANRT,ANSO,TCNT,DELT,            NDEML,NDEMS,NDEMR,NDEMSO)

* Calling the subroutine for translocatable N in leaves, stem, roots and
* storage organs.
CALL NTRLOC(ANLV,ANST,ANRT,WLVG,WST,WRT,RNFLV,RNFST,RNFRT,...
            FNTRT,...
            ATNLV,ATNST,ATNRT,ATN)

* Calling the subroutine for N translocated from leaves, stem, and roots.
CALL NTRANS(RNSO,ATNLV,ATNST,ATNRT,ATN,  RNTLV,RNTST,RNTRT)

* Calling the subroutine for N loss due to death of leaves and roots.
CALL RNLD(DVS,WRT,RDRRT,RNFLV,DLV,RNFRT,DVSDR,...
          DRRT,RNLDLV,RNLDRT)

* Calling the subroutine for Nitrogen Nutrition Index.
CALL NNINDX(TIME,DOYEM,EMERG,NFGMR,NRMR,NOPTMR,    NNI)

* Calling the subroutine for relative death rate of leaves.
CALL DEATHL(TIME,DOYEM,TSUM,TSUMAN,RDRTMP,RDRSHM,LAI,LAICR,WLVG,...
            RDRNS,NNI,SLA,      RDRDV,RDRSH,RDR,DLV,DLVS,DLVNS,...
            DLAIS,DLAINS,DLAI)

*** Leaf growth and LAI.
      GLV    = FLV * GTOTAL
      LAI    = INTGRL(LAII, RLAI)

*   Net rate of change of Leaf area.
      RLAI   = GLAI - DLAI

*  Relative death rate of leaves due to senescence/ageing.
      RDRTMP = AFGEN(RDRT, DVS)

*   Maximum N concentration in the leaves, from which the values of the
*   stem and roots are derived, as a function of development stage.

* From ORYZA_N
    NMAXLV = AFGEN (NMXLV, DVS)
    NMAXST = LSNR * NMAXLV
    NMAXRT = LRNR * NMAXLV

    PARAM NMAXSO = 0.0175
    PARAM LSNR = 0.50; LRNR = 0.37

**Rate of change of N in organs
    RNLV = RNULV-RNTLV-RNLDLV
    RNST = RNUST-RNTST
    RNRT = RNURT-RNTRT-RNLDRT

**Total vegetative biomass.
    TBGMR =WLVG+WST

** Total above ground biomass
    TAGBM = WLV + WST + WSO
     
*  Maximum N concentration as the fraction of maximum N concentration.
    PARAM FRNX = 0.5

* Maximum N content in the plant.
    NOPTS = NOPTST* WST
    NOPTL = NOPTLV* WLVG
    NOPTMR = (NOPTL+ NOPTS)/NOTNUL(TBGMR)

**Actual N content in various organs.
    ANLV =  INTGRL (ANLVI,RNLV)
    ANST =  INTGRL (ANSTI,RNST)
    ANRT =  INTGRL (ANRTI,RNRT)
    ANSO =  INTGRL (ANSOI,RNSO)

* Rate of N uptake in grains.
    RNSO =  AMIN1 (NDEMSO,NSUPSO)
    
* N supply to the storage organs.
     NSUPSO = INSW (DVS-DVSNT,0.,ATN/TCNT)
     PARAM DVSNT = 0.8
     
* N concentration (g N g-1 DM) of the leaves, stem, roots and storage
* organs
     NFLV = ANLV/ NOTNUL(WLVG)
     NFST = ANST/ NOTNUL(WST)
     NFRT = ANRT/ NOTNUL(WRT)
     NFSO = ANSO/ NOTNUL(WSO)

* Total N in green matter of the plant.
     NUPGMR = ANLV + ANST
     NFGMR  = NUPGMR/NOTNUL(TBGMR)

**Average residual N concentration.
     NRMR   = (WLVG*RNFLV+WST*RNFST)/NOTNUL(TBGMR)

*****************SOIL NITROGEN SUPPLY***********************************
*    Total Nitrogen demand.
     NDEMTO = MAX (0.0,(NDEML + NDEMS + NDEMR))

*    Nitrogen uptake limiting factor at low moisture conditions in the
*    rooted soil layer before anthesis. After anthesis there is no
*    uptake from the soil anymore.
     NLIMIT = INSW(DVS-DVSNLT , INSW(WC-WCWP,0.,1.) , 0.0)
     PARAM DVSNLT = 1.0

*----------------------------------------------------------------
*    Total uptake of N.
     NTAKE  = INTGRL (ZERO,NUPTR)

*    Soil N supply (g N m-2 d-1) through mineralization.
     RTMIN  = 0.05

*--------------------Fertilizer N---------------------------------------
*    Total nutrient uptake.
     NUPTR = (MAX (0., MIN (NDEMTO, TNSOIL))* NLIMIT ) / DELT * ...
             INSW(TIME-DOYEM, 0.,1.)
     
*    Amount of N available each day for crop uptake.
     TNSOIL = MAX(0.,FERTNP)+ RTMIN * DELT

*    Fertilizer N pool present in the soil.
     FERTNP = INTGRL (ZERO, RFERNP)
     
*    Rate of fertlizer N addition to the pool, as function of fertilizer
*    input and crop uptake.
     RFERNP = FERTNS/DELT -MAX(0.,(NUPTR-RTMIN))

*---------------Fertilizer application---------------------------------*
     FERTN  = AFGEN  (FERTAB, TIME)
     FERTNS = FERTN * NRF
     NRF    = AFGEN (NRFTAB, DVS)

**------------------Nitrogen Balance---------------------------------***
     NLOSSL = INTGRL (ZERO,RNLDLV)
     NLOSSR = INTGRL (ZERO,RNLDRT)

     NBALAN = NTAKE+(ANLVI+ANSTI+ANRTI+ANSOI)-(ANLV+ANST+ANRT+ANSO...
              +NLOSSL+NLOSSR)

*-------- Growth rates and dry matter production of plant organs-------*
* Biomass partitioning functions under (water and nitrogen)non-stressed
* situations
     FRTWET = AFGEN(FRTTB, DVS )
     FLVT   = AFGEN(FLVTB, DVS )
     FSTT   = AFGEN(FSTTB, DVS)
     FSOT   = AFGEN(FSOTB, DVS )
     FRACT  = FLV+ FST+ FSO + FRT

*  Weight of green leaves, dead leaves, stem, storage organs and roots
     WLVG   = INTGRL( WLVGI,RWLVG)
     WLVD   = INTGRL( ZERO, DLV  )
     WST    = INTGRL( WSTI, RWST )
     WSO    = INTGRL( WSOI, RWSO )
     WRT    = INTGRL( WRTLI,RWRT )

*    Relative death rate of leaves due to N stress.
     RDRNS  = 0.03

*    Development stage above which death of leaves and roots start.
     PARAM DVSDR = 0.8
      
*  Total leaf weight.
     WLV    = WLVG + WLVD

*  Weight of storage organs in terms of tons/ha.
     WSOTHA = WSO / 100.
                       
*    Root growth.
     ROOTD  = INTGRL( ROOTDI, RROOTD)
     RROOTD = RRDMAX * INSW( WC-WCWP, 0., 1. ) * ...
               REAAND(ROOTDM-ROOTD,TSUMAN-TSUM ) * EMERG

*    Water-Nitrogen stress factor
     RNW = MIN(TRANRF,NNI)

**------------------Biomass Carbon Balance---------------------------***

     GTSUM  = INTGRL(ZERO, GTOTAL)
     WDRT   = INTGRL(ZERO, DRRT)

     CBALAN = GTSUM+(WRTLI+WLVGI+WSTI+WSOI)-(WLV+WST+WSO+WRT+WDRT)

** For calculation of LUE (it can be removed once after checking, however,
*  I have put these new additions in these two (c and N) balances sections
     PAR   = DTR*0.50
     LUECAL = GTOTAL/PAR
     CUMPAR= INTGRL(ZERO, PAR)

***--------------------Soil moisture balance---------------------------*

*    Water management
*    We will represent two water management situations in the model
*    as irrigated up to the field capacity:                  WMFAC = 0.0
*    as irrigated up to saturation, thus mimicking flooding; WMFAC = 1.0
*    In both cases parameter IRRIGF must be taken 1.0. "Irrigation" thus
*    means that the moisture content is kept at optimal conditions. It
*    has not to do with irrigation that is done in reality.
*    Here, the crop being flooded is rice, and the soil water content
*    is always saturated.

     PARAM WMFAC = 1.0
      
     WA  = INTGRL(WAI,RWA)
     RWA = (RAIN+EXPLOR+IRRIG)-(RNINTC+RUNOFF+TRAN+EVAP+DRAIN)

*    Rainfall interception by the crop
     RNINTC = MIN( RAIN, 0.25*LAI )

*    Exploration of water in soil when roots grow downward.
     EXPLOR = 1000. * RROOTD * INSW(WMFAC-1.,WCFC,WCST)

*    Water content in the rootzone
     WC  = 0.001* WA/NOTNUL(ROOTD)

**   Growth reduction function for water stress(actual trans/potential)
     TRANRF = TRAN / NOTNUL(PTRAN)

**   Total Exploration
     TEXPLO = INTGRL(ZERO, EXPLOR)
      
***-------------Functions and parameters for rice----------------------*

*    Specific Leaf area(m2/g).
     SLA = SLAC * AFGEN (SLACF,DVS)*EXP(-NSLA * (1.-NNI))

*    Specific leaf area constant.
     PARAM  SLAC = 0.02

*    Base temperature for rice crop.
     PARAM TBASE = 8.
     
*    The day of the year (July 13 = DOY 194) on which crop emerges.
     PARAM DOYEM = 25.

*    Relative growth rate of LAI at the exponential growth phase
*    (oC d)-1, critical LAI above which mutual shading of leaves occurs,
*    and the maximum relative death rate of leaves due to shading.

     PARAM RGRL = 0.0085    ; LAICR = 4.0    ; RDRSHM = 0.03

*    Light use efficiency and extinction coefficient.
     PARAM LUE = 3.0;  K = 0.6

*    Maximum root depth and maximum rate of increase in rooting depth
*    (m d-1) for a rice crop.
     PARAM ROOTDM = 1.0; RRDMAX = 0.010

*    Soil hydraulic properties
     PARAM WCAD  = 0.01 ; WCWP = 0.21; WCFC = 0.47
     PARAM WCWET = 0.50; WCST = 0.55
     
*    Transpiration constant (mm/day) indicating the level of drought
*    tolerance of the rice crop.
     PARAM TRANCO = 3.

*    Maximum drainage rate of the soil (mm/day) and Irrigation factor
*    (if 1.0, the soil will be always at saturation, if WMFAC=1.0    )
*    (if 1.0, the soil will be always at field capacity, if WMFAC=0.0)
     PARAM DRATE = 30.; IRRIGF = 1.0
     
*    Time coefficient(days) for N translocation.
     PARAM TCNT = 10.

*    Relative death rate of roots.
     RDRRT = 0.03

*    Residual N concentration in leaves, stem and roots.
     PARAM RNFLV= 0.004 ; RNFST= 0.0015; RNFRT= 0.002

*    Nitrogen translocation from roots as a fraction of the total amount
*    of nitrogen translocated from leaves and stem.
     PARAM FNTRT = 0.15

*    Extinction coefficient for  Nitrogen distribution down the canopy
     PARAM NLUE = 0.20

*    Coefficient for the effect of N stress on LAI reduction(during juvenile
*    phase)
     PARAM NLAI = 1.0

*    Coefficient for the effect of N stress on SLA reduction
     PARAM NSLA = 1.0

*    Coefficient for the effect of N stress on leaf biomass reduction
     PARAM NPART= 1.0

************************** Functions ***********************************
*    Function to include the effect of photoperiodicity.
     FUNCTION PHOTTB = 0.0,0.0,       8.,0.0,    10.,1.0, ...
                       12.,1.0,      13.,0.8,    14.,0.6, ...
                       18.,0.0

*    Relative death rate of leaves as a function of Developmental stage
*    (adopted from ORYZA2000)
     FUNCTION RDRT  = 0.,0.00, 0.6,0.00, 1.0,.015, 1.6,0.025, 2.1,0.05

*    Leaf area correction function as a function of development stage,
*    DVS. Reference: Drenth, H.,ten Berge, H.F.M. and Riethoven, J.J.M.
*    1994, p.10. (Complete reference under Observed data.)
     FUNCTION SLACF     =   0.0 ,1.72,   0.21,1.72,    0.24,1.72, ...
                            0.33,1.32,   0.7 ,1.20,    1.01,1.00, ...
                            2.0 ,0.75,   2.1 ,0.75

*    Maximum N concentration in the leaves, from which the values of the
*    stem and roots are derived, as a function of development stage.
     FUNCTION NMXLV = 0.0,0.05,    0.4,0.05,   0.7,0.04, ...
                      1.0,0.03,    2.0,0.02,   2.1,0.02

********** Partitioning coefficients ***********************************

    FUNCTION FRTTB =  0.0 ,0.30,   0.48,0.30,  0.62,0.12, 0.69,0.11 ,...
                      0.84,0.11,   0.92,0.10,  1.00,0.08, ...
                      1.38,0.00,   2.10,0.0

    FUNCTION FLVTB =  0.0 ,0.315,  0.48,0.35,  0.62,0.44, 0.69,0.463,...
                      0.84,0.463,  0.92,0.45,  1.00,0.00, ...
                      1.38,0.00 ,  2.10, 0.0

    FUNCTION FSTTB =  0.0 ,0.385,  0.48,0.35,  0.62,0.44, 0.69,0.427,...
                      0.84,0.427,  0.92,0.27,  1.00,0.00, ...
                      1.38,0.00 ,  2.1, 0.0

    FUNCTION FSOTB =  0.0 ,0.00,   0.48,0.00,  0.62,0.00, 0.69,0.00,...
                      0.84,0.00,   0.92,0.18,  1.00,0.92, ...
                      1.38,1.00,   2.1, 1.00

*************************VALIDATION************************************
***OBSERVED DATA, LOS BANOS,IRRI,1993, IR72; UNDER 0,100,400 KG N/HA****
*
* Source:  Drenth, H.,ten Berge, H.F.M. and Riethoven, J.J.M.(Editors).
* ORYZA simulation modules for potential and nitrogen limited rice
* production SARP Research Proceedings-December, 1994.
* DLO-Research Institute for Agrobiology and Soil fertility, Wageningen,
* WAU-Department of Theoretical Production Ecology, Wageningen,
* IRRI-International Rice Research Institute, Los Banos,
* Pages: 197-210.

************************************************************************
* Fertilizer application @ 0 KgN/ha*
*----------------------------------------------------------------------*
* Leaf area index.
LAIOB  = AFGEN (LAIOBT, TIME)

* Weight of Leaves, stems, roots and storage organs.
WSTOB  = AFGEN (WSTOBT, TIME)
WLVGOB = AFGEN (WLVGOT, TIME)
WSOOB  = AFGEN (WSOOBT, TIME)
WRTOB  = AFGEN (WRTOBT, TIME)
WLVDOB = AFGEN (WLVDOT, TIME)

WLVOB  = WLVGOB + WLVDOB
TAGBMO = WLVOB + WSTOB + WSOOB

* Amount N in leaves, stems, storage organs and roots.
ANLVOB = AFGEN (ANLVOT, TIME)
ANSTOB = AFGEN (ANSTOT, TIME)
ANSOOB = AFGEN (ANSOOT, TIME)
ANRTOB = AFGEN (ANRTOT, TIME)

* Nitrogen concentration in leaves, stem and storage organs.
NLFCOB = ANLVOB/NOTNUL(WLVOB)
NSTCOB = ANSTOB/NOTNUL(WSTOB)
NSOCOB = ANSOOB/NOTNUL(WSOOB)

NTAG   = ANLV   + ANST   + ANSO
NTAGOB = ANLVOB + ANSTOB + ANSOOB

NTAC   = NTAG/NOTNUL(TAGBM)
NTACOB = NTAGOB/NOTNUL(TAGBMO)

RTSH   = WRT/TAGBM
RTSHOB = WRTOB/TAGBMO

* Fertilizer application as a function of TIME.
FUNCTION FERTAB =  25.,0.,   26.,0.,   27.,0.,   33.,0.,...
                   34.,0.,   35.,0.,   57., 0.,  58.,0.,...
                   59.,0.,   69.,0.,   70.,0.,   71.,0.,...
                   83.,0.,   84.,0.,   85.,0.,   97.,0.,...
                   98.,0.,   99.,0.,   300.,0.,  365.,0.


*Fertilizer nitrogen recovery fraction
FUNCTION NRFTAB =  0.0,0.30,     0.2, 0.35,    0.4,0.50,    ...
                   0.8,0.75,     1.0, 0.75,    2.5,0.75

FUNCTION LAIOBT = 25.,  0.118 ,   36.0,  0.118,   50.0,  0.414,...
                  64.0,  1.388,   76.0,  1.920,   90.0,  1.388,...
                  98.0,  1.506,   112.,  0.768

FUNCTION WSTOBT = 25., 0.71,     36., 21.8,     50., 25.2,...
                  64., 67.2,     76., 134.5,    90., 277.33,...
                  98., 302.5,    112.,277.33

FUNCTION WLVGOT = 25., 0.86,    36., 3.4,     50., 33.61,...
                  64., 58.82,   76., 100.84,  90., 84.03,...
                  98., 109.24,  112.,58.82

FUNCTION WSOOBT = 25., 0.0,     36., 0.0,     50., 0.0,...
                  64., 0.0,     76., 0.0,     90., 84.03,...
                  98., 193.28,  112.,403.36

FUNCTION TAGBOT = 25., 1.58,     36., 25.2,     50., 58.8,...
                  64., 126.1,   76., 235.3,    90., 453.8,...
                  98., 605.0,  112., 739.5


FUNCTION WRTOBT = 25.,0.0,    128.,0.

FUNCTION WLVDOT = 25.,0.0,    128.,0.

FUNCTION ANLVOT = 25.,0.,     128.,0.

FUNCTION ANSTOT = 25.,0.,     128.,0.

FUNCTION ANSOOT = 25.,0.,     128.,0.

FUNCTION ANRTOT = 25.,0.,     128.,0.

END

*---Fertilizer application @ 100 KgN /ha-------------------------------*
FUNCTION FERTAB =  25.,0.,   26.,1.43,    27.,0.,...
                   32.,0.,   33.,1.43,    34.,0.,...
                   39.,0.,   40.,1.43,    41.,0.,...
                   46.,0.,   47.,1.43,    48.,0.,...
                   53.,0.,   54.,1.43,    55.,0.,...
                   60.,0.,   61.,1.43,    62.,0.,...
                   67.,0.,   68.,1.43,    69.,0.,...
                   98.,0.,   99.,0.0,    300.,0.,...
                   365.,0.

FUNCTION LAIOBT =  25.,  0.118 ,  36.0,  0.118,   50.0,  0.5,...
                   64.0, 2.118,   76.0,  2.176,   90.0,  3.118,...
                   98.0, 2.324,   116.,  1.912

FUNCTION WSTOBT =  25., 0.71,     36., 21.8,     50., 25.2,...
                   64., 84.0,     76., 193.3,    90., 352.9,...
                   98., 378.2,    116.,521.0

FUNCTION WLVGOT =  25., 0.86,     36., 3.4,     50., 33.61,...
                   64., 75.63,    76., 134.45,  90., 142.86,...
                   98., 126.05,   116.,126.05

FUNCTION WSOOBT =  25., 0.0,      36., 25.4,     50., 0.0,...
                  64., 0.0,       76., 0.0,     90., 100.84,...
                  98., 235.29,    116.,764.706

FUNCTION TAGBOT =  25.,1.58,      36., 0.0,     50., 42.0,...
                  64., 159.7,     76.,327.7,     90.,596.6, ...
                  98., 739.5,     116.,1411.8

FUNCTION WRTOBT = 25.,0.0,   128.,0.
FUNCTION WLVDOT = 25.,0.0,    128.,0.
FUNCTION ANLVOT = 25.,0.,     128.,0.
FUNCTION ANSTOT = 25.,0.,     128.,0.
FUNCTION ANSOOT = 25.,0.,     128.,0.
FUNCTION ANRTOT = 25.,0.,     128.,0.

END

*---Fertilizer application @ 400 KgN /ha-------------------------------*
FUNCTION FERTAB =  25.,0.,   26.,5.71,    27.,0.,   ...
                   32., 0.,  33.,5.71,    34.,0.,    ...
                   39.,0.,   40.,5.71,    41.,0.,    ...
                   46., 0.,  47.,5.71,    48.,0.,    ...
                   53.,0.,   54.,5.71,    55.,0.,    ...
                   60.,0.,   61.,5.71,    62.,0.,    ...
                   67.,0.,   68.,5.71,    69.,0.,    ...
                   98.,0.,   99.,0.0,    300.,0.,   ...
                   365.,0.

FUNCTION LAIOBT =  25.,  0.118 ,  36.0,  0.119,   50.0,  0.475,...
                   64.0, 3.381,   76.0,  5.784,   90.0,  5.992,...
                   98.0, 4.479,   116.,  3.174


FUNCTION WSTOBT = 25., 0.71,      36.,21.8 ,     50., 34.3,...
                  64., 120.2,     76., 248.9,    90., 438.3,...
                  98., 473.5,     116.,604.9

FUNCTION WLVGOT = 25., 0.86,      36., 8.58,     50., 33.61,...
                  64., 103.0,     76., 248.93,  90., 300.43,...
                  98., 240.34,    116.,197.43

FUNCTION WSOOBT = 25., 0.0,       36., 0.0,     50., 0.0,...
                  64., 0.0,       76., 0.0,     90., 128.,...
                  98., 333.,      116.,949.

FUNCTION TAGBOT  = 25., 1.58,      36.,30.38,      50., 52.,...
                  64., 223.,      76., 498.,     90., 867.,...
                  98., 1047.,     116.,1751.

FUNCTION WRTOBT = 25., 0.0,    128.,0.

FUNCTION WLVDOT = 25.,0.0,      128.,0.

FUNCTION ANLVOT = 25.,0.,     128.,0.

FUNCTION ANSTOT = 25.,0.,     128.,0.

FUNCTION ANSOOT = 25.,0.,     128.,0.

FUNCTION ANRTOT = 25.,0.,      128.,0.

END
STOP

*************************Subroutines************************************
* ---------------------------------------------------------------------*
*  SUBROUTINE ASTRO                                                    *
*  Purpose: This subroutine calculates the astrological daylength,     *
*           based on Goudriaan and van Laar 1994, around p.30.         *
* ---------------------------------------------------------------------*

      SUBROUTINE ASTRO(DOY, LAT,
     $                 DAYL)
      IMPLICIT REAL (A-Z)
      SAVE
      
*-----PI VALUE
      PI     = 3.1415926

* SINE AND COSINE OF LATITUDE
      SINLAT = SIN (PI * LAT/180.)
      COSLAT = COS (PI * LAT/180.)

* MAXIMAL SINE OF DECLINATION
      SINDCM = SIN (PI * 23.45/180.)

* SINE AND COSINE OF DECLINATION  (EQUATIONS 3.4, 3.5)

      SINDEC = -SINDCM * COS (2.* PI * (DOY+10.)/365.)
      COSDEC = SQRT (1.-SINDEC * SINDEC)

* THE TERMS A AND B ACCORDING TO EQUATION 3.3

      A     = SINLAT * SINDEC
      B     = COSLAT * COSDEC

* DAYLENGTH ACCORDING TO EQUATION 3.6

      DAYL  = 12. * (1. + (2./PI)* ASIN (A/B))

      RETURN
      END
      
* ---------------------------------------------------------------------*
*  SUBROUTINE SUBDVS                                                   *
*  Purpose: To compute the developmental stage corresponding to the    *
*  temperature sum.                                                    *
* ---------------------------------------------------------------------*

      SUBROUTINE SUBDVS(TIME,DOYEM,TSUM,TSUMAN,TSUMMT,    DVS)
      IMPLICIT REAL (A-Z)
      SAVE

      DVS = 0.

      IF( TIME .LE. DOYEM ) THEN
      DVS1 = TSUM/TSUMAN
      DVS2 = 0.0

      ELSEIF( TIME .GE. DOYEM .AND. TSUM .LE. TSUMAN ) THEN
            DVS1   = TSUM/TSUMAN

      ELSE
            DVS2   = (TSUM-TSUMAN)/TSUMMT

      ENDIF

      DVS = DVS1+DVS2

      RETURN
      END

* ---------------------------------------------------------------------*
*  SUBROUTINE PENMAN                                                   *
*  Purpose: Computation of the PENMAN EQUATION.                        *
* ---------------------------------------------------------------------*

      SUBROUTINE PENMAN(DAVTMP,VP,DTR,LAI,WN,RNINTC,
     $                  RLWN,NRADC,PENMRC,PENMD,PEVAP,PTRAN)
      IMPLICIT REAL (A-Z)
      SAVE

      DTRJM2 = DTR * 1.E6
      BOLTZM = 5.668E-8
      LHVAP  = 2.4E6
      PSYCH  = 0.067

      BBRAD  = BOLTZM * (DAVTMP+273.)**4 * 86400.
      SVP    = 0.611 * EXP(17.4 * DAVTMP / (DAVTMP + 239.))
      SLOPE  = 4158.6 * SVP / (DAVTMP + 239.)**2
      RLWN   = BBRAD * MAX(0.,0.55*(1.-VP/SVP))
      NRADS  = DTRJM2 * (1.-0.15) - RLWN
      NRADC  = DTRJM2 * (1.-0.25) - RLWN
      PENMRS = NRADS * SLOPE/(SLOPE+PSYCH)
      PENMRC = NRADC * SLOPE/(SLOPE+PSYCH)

      WDF    = 2.63 * (1.0 + 0.54 * WN)
      PENMD  = LHVAP * WDF * (SVP-VP) * PSYCH/(SLOPE+PSYCH)

      PEVAP  = EXP(-0.5*LAI)  * (PENMRS + PENMD) / LHVAP
      PTRAN  = (1.-EXP(-0.5*LAI)) * (PENMRC + PENMD) / LHVAP
      PTRAN  = MAX( 0., PTRAN-0.5*RNINTC )

      RETURN
      END

* ---------------------------------------------------------------------*
*  SUBROUTINE EVAPTR                                                   *
*  Purpose: To compute actual rates of evaporation and transpiration.  *
* ---------------------------------------------------------------------*

      SUBROUTINE EVAPTR(PEVAP,PTRAN,ROOTD,WA,WCAD,WCWP,WCFC,WCWET,WCST,
     $                  TRANCO,DELT,WMFAC,
     $                  WCCR,EVAP,TRAN,FR,AVAILF)
      IMPLICIT REAL (A-Z)
      SAVE

      WC   = 0.001 * WA   / NOTNUL(ROOTD)
      WAAD = 1000. * WCAD * ROOTD
      WAFC = 1000. * WCFC * ROOTD
      EVAP = PEVAP * LIMIT( 0., 1., (WC-WCAD)/(WCFC-WCAD) )
      WCCR = WCWP + MAX( 0.01, PTRAN/(PTRAN+TRANCO) * (WCFC-WCWP) )
      
*     If the soil is flooded for rice, then the soil is assumed to be
*     permanently saturated and there is no effect of the high water
*     content on crop growth, because rice has aerenchyma.
*     Thus FR is formulated as below:

      IF (WMFAC .GE. 1.0 ) THEN
         IF (WC.GT.WCCR) THEN
             FR = 1.
         ELSE
             FR = LIMIT( 0., 1., (WC-WCWP)/(WCCR-WCWP)  )
         ENDIF

*     If soil is irrigated but not flooded, then soil water content is
*     assumed to be at field capacity and the critical water content
*     that affects crop growth (FR) is formulated as below:

      ELSE
         IF (WC.GT.WCCR) THEN
             FR = LIMIT( 0., 1., (WCST-WC)/(WCST-WCWET) )
         ELSE
             FR = LIMIT( 0., 1., (WC-WCWP)/(WCCR-WCWP)  )
         ENDIF
      ENDIF

      TRAN = PTRAN * FR

         AVAILF = MIN( 1., ((WA-WAAD)/DELT)/NOTNUL(EVAP+TRAN) )

      EVAP = EVAP * AVAILF
      TRAN = TRAN * AVAILF

      RETURN
      END

* ---------------------------------------------------------------------*
*  SUBROUTINE DRUNIR                                                   *
*  Purpose: To compute rates of drainage, runoff and irrigation.       *
* ---------------------------------------------------------------------*

      SUBROUTINE DRUNIR(RAIN,RNINTC,EVAP,TRAN,IRRIGF,
     $                  DRATE,DELT,WA,ROOTD,WCFC,WCST,WMFAC,
     $                  DRAIN,RUNOFF,IRRIG)
      IMPLICIT REAL (A-Z)
      SAVE

      WC   = 0.001 * WA   / NOTNUL(ROOTD)
      WAFC = 1000. * WCFC * ROOTD
      WAST = 1000. * WCST * ROOTD

      DRAIN  = LIMIT( 0., DRATE, (WA-WAFC)/DELT +
     $               (RAIN - RNINTC - EVAP - TRAN)                  )
      RUNOFF =          MAX( 0., (WA-WAST)/DELT +
     $               (RAIN - RNINTC - EVAP - TRAN - DRAIN)          )

      IF (WMFAC .GE. 1.0) THEN
*        If a soil is irrigated by flooding, then soil water content is
*        kept at saturation via "irrigation events".
         IRRIG  = IRRIGF * MAX( 0., (WAST-WA)/DELT -
     $                  (RAIN - RNINTC - EVAP - TRAN - DRAIN - RUNOFF) )
      ELSE

*        If soil is irrigated but not flooded, then soil water content
*        is kept at field capacity via "irrigation events".
         IRRIG  = IRRIGF * MAX( 0., (WAFC-WA)/DELT -
     $                  (RAIN - RNINTC - EVAP - TRAN - DRAIN - RUNOFF) )
      ENDIF
      
      RETURN
      END

* ---------------------------------------------------------------------*
*  SUBROUTINE GLA                                                      *
*  Purpose: This subroutine computes daily increase of leaf area index *
*           (ha leaf/ ha ground/ d).                                   *
* ---------------------------------------------------------------------*

      SUBROUTINE GLA(TIME,DOYEM,DTEFF,LAII,RGRL,DELT,SLA,LAI,GLV,NLAI,
     $               WC,WCWP,DVS,TRANRF,NNI,
     $               GLAI)
      IMPLICIT REAL (A-Z)
      SAVE

*---- Growth during maturation stage:
      GLAI = SLA * GLV

*---- Growth during juvenile stage:
      IF ((DVS .LT. 0.2).AND.(LAI .LT. 0.75))
     $   GLAI = (LAI * (EXP(RGRL * DTEFF * DELT) - 1.)/ DELT )* TRANRF*
     $           EXP(-NLAI* (1.0 - NNI))

*---- Growth at day of seedling emergence:
      IF ((TIME.GE.DOYEM).AND.(LAI.EQ.0.).AND.(WC.GT.WCWP))
     $   GLAI = LAII / DELT

*---- Growth before seedling emergence:
      IF (TIME.LT.DOYEM) GLAI = 0.

      RETURN
      END

*----------------------------------------------------------------------*
*  SUBROUTINE SUBPAR                                                  *
*  Purpose: Dry matter partitioning fractions: leaves, stem and storage organs.
*----------------------------------------------------------------------*
      SUBROUTINE SUBPAR (NPART,TRANRF,NNI,FRTWET,FLVT,FSTT,FSOT,
     $FSHMOD,FLVMOD,FRT,FLV,FST,FSO)
      IMPLICIT REAL (A-Z)
      SAVE

      IF(TRANRF .LT. NNI) THEN
*     Water stress is more severe as compared to nitrogen stress and
*     partitioning will follow the original assumptions of LINTUL2*

      FRTMOD = MAX( 1., 1./(TRANRF+0.5))
      FRT    = FRTWET * FRTMOD
      FSHMOD = (1.-FRT) / (1.-FRT/FRTMOD)
      FLV    = FLVT * FSHMOD
      FST    = FSTT * FSHMOD
      FSO    = FSOT * FSHMOD

      ELSE

*    Nitrogen stress is more severe as compared to water stress and the
*    less partitioning to leaves will go to the roots, THE STEMS AND THE
*    STORAGE ORGANS.

      FLVMOD = EXP(-NPART* (1.0-NNI))
      FLV    = FLVT * FLVMOD
      MODIF  = (1.-FLV)/(1.-(FLV/FLVMOD))
      FST    = FSTT *  MODIF
      FRT    = FRTWET* MODIF
      FSO    = FSOT *  MODIF

      ENDIF

      RETURN
      END

* ---------------------------------------------------------------------*
*  SUBROUTINE GROWTH                                                   *
*  Purpose: To compute the total growth rate.                          *
* ---------------------------------------------------------------------*

      SUBROUTINE GROWTH(TIME,DOYEM,DTR,K,NLUE,LAI,LUE,TRANRF,NNI,
     $                                                    PARINT,GTOTAL)
      IMPLICIT REAL (A-Z)
      SAVE

      PARINT = 0.5 * DTR * (1.- EXP(-K*LAI)) * INSW(TIME-DOYEM,0.,1.)

      IF(TRANRF .LE. NNI) THEN
*     Water stress is more severe as compared to nitrogen stress and
*     partitioning will follow the original assumptions of LINTUL2*

          GTOTAL = LUE * PARINT * TRANRF

      ELSE

*     Nitrogen stress is more severe as compared to water stress and the
*     less partitioning to leaves will go to the roots*

          GTOTAL = LUE * PARINT * EXP(-NLUE*(1.0-NNI))

      ENDIF

      RETURN
      END

* ---------------------------------------------------------------------*
*  SUBROUTINE RELGR                                                    *
*  Purpose: To compute the relative growth rate of roots, leaves, stem *
*           and storage organs.                                        *
* ---------------------------------------------------------------------*

      SUBROUTINE RELGR(TIME,DOYEM,EMERG,WLVGI,WRTLI,WSTI,WSOI,GTOTAL,
     $                 FLV,FRT,FST,FSO,DLV,DRRT,DELT,
     $                 RWLVG,RWRT,RWST,RWSO)
      IMPLICIT REAL (A-Z)
      SAVE

      IF (TIME .GE. DOYEM .AND. EMERG .EQ. 1.0) THEN
          RWLVG = GTOTAL * FLV - DLV
          RWRT  = GTOTAL * FRT - DRRT
          RWST  = GTOTAL * FST
          RWSO  = GTOTAL * FSO

      ELSE
          RWLVG = 0.
          RWRT  = 0.
          RWST  = 0.
          RWSO  = 0.

      ENDIF

      RETURN
      END
* ---------------------------------------------------------------------*
*  SUBROUTINE RNUSUB                                                   *
*  Purpose: To compute the partitioning of the total N uptake rate     *
*           (NUPTR) over the leaves, stem, and roots.                  *
* ---------------------------------------------------------------------*

      SUBROUTINE RNUSUB(TIME,DOYEM,EMERG,NDEML,NDEMS,NDEMR,NUPTR,
     $                  NDEMTO,ANLVI,ANSTI,ANRTI,DELT,
     $                  RNULV,RNUST,RNURT)
      IMPLICIT REAL (A-Z)
      SAVE
      IF (TIME .GE. DOYEM .AND. EMERG .EQ.1)THEN
          RNULV = (NDEML / NOTNUL(NDEMTO))* NUPTR
          RNUST = (NDEMS / NOTNUL(NDEMTO))* NUPTR
          RNURT = (NDEMR / NOTNUL(NDEMTO))* NUPTR

      ELSE
          RNULV = 0.
          RNUST = 0.
          RNURT = 0.

      ENDIF

      RETURN
      END

* ---------------------------------------------------------------------*
*  SUBROUTINE NOPTM                                                    *
*  Purpose: To compute the maximum nitrogen concentration of the organs*
* ---------------------------------------------------------------------*

      SUBROUTINE NOPTM(FRNX,NMAXLV,NMAXST, NOPTLV,NOPTST)
      IMPLICIT REAL (A-Z)
      SAVE

      NOPTLV= FRNX * NMAXLV
      NOPTST= FRNX * NMAXST

      RETURN
      END

* ---------------------------------------------------------------------*
*  SUBROUTINE NDEMND                                                   *
*  Purpose: To compute the nitrogen demand of the vegetative organs.   *
* ---------------------------------------------------------------------*

      SUBROUTINE NDEMND(NMAXLV,NMAXST,NMAXRT,NMAXSO,WLVG,WST,WRT,WSO,
     $                  RWLVG,RWST,RWRT,RWSO,
     $                  ANLV,ANST,ANRT,ANSO,TCNT,DELT,
     $                  NDEML,NDEMS,NDEMR,NDEMSO)
      IMPLICIT REAL (A-Z)
      SAVE

      NDEML  =  AMAX1 (NMAXLV*WLVG -  ANLV, 0.)
      NDEMS  =  AMAX1 (NMAXST*WST  - ANST, 0.)
      NDEMR  =  AMAX1 (NMAXRT*WRT  - ANRT, 0.)
      NDEMSO =  AMAX1 (NMAXSO*WSO  - ANSO, 0.)/TCNT

      RETURN
      END

* ---------------------------------------------------------------------*
*  SUBROUTINE NTRLOC                                                   *
*  Purpose: To compute the translocatable N in the organs.             *
* ---------------------------------------------------------------------*

      SUBROUTINE NTRLOC(ANLV,ANST,ANRT,WLVG,WST,WRT,RNFLV,RNFST,RNFRT,
     $                  FNTRT,
     $                  ATNLV,ATNST,ATNRT,ATN)
      IMPLICIT REAL (A-Z)
      SAVE

      ATNLV = AMAX1 (0. , ANLV-WLVG*RNFLV)
      ATNST = AMAX1 (0. , ANST-WST*RNFST)
      ATNRT = AMIN1((ATNLV + ATNST) * FNTRT, ANRT-WRT*RNFRT)
      ATN   = ATNLV +  ATNST + ATNRT

      RETURN
      END
* ---------------------------------------------------------------------*
*  SUBROUTINE NTRANS                                                   *
*  Purpose: To compute the N translocated from different organs.       *
* ---------------------------------------------------------------------*

      SUBROUTINE NTRANS(RNSO,ATNLV,ATNST,ATNRT,ATN,  RNTLV,RNTST,RNTRT)
      IMPLICIT REAL (A-Z)
      SAVE

      RNTLV= RNSO* ATNLV/ NOTNUL(ATN)
      RNTST= RNSO* ATNST/ NOTNUL(ATN)
      RNTRT= RNSO* ATNRT/ NOTNUL(ATN)

      RETURN
      END
* ---------------------------------------------------------------------*
*  SUBROUTINE RNLD                                                     *
*  Purpose: To compute the N loss due to death of leaves and roots.    *
* ---------------------------------------------------------------------*

      SUBROUTINE RNLD(DVS,WRT,RDRRT,RNFLV,DLV,RNFRT,DVSDR,
     $                                DRRT,RNLDLV,RNLDRT)
      IMPLICIT REAL (A-Z)
      SAVE

      IF (DVS .LT. DVSDR)THEN
         DRRT= 0.
      ELSE
         DRRT= WRT * RDRRT
      ENDIF

      RNLDLV= RNFLV* DLV
      RNLDRT= RNFRT* DRRT

      RETURN
      END
* ---------------------------------------------------------------------*
*  SUBROUTINE NNINDX                                                   *
*  Purpose: To compute Nitrogen Nutrition Index                        *
* ---------------------------------------------------------------------*

      SUBROUTINE NNINDX(TIME,DOYEM,EMERG,NFGMR,NRMR,NOPTMR,    NNI)
      IMPLICIT REAL (A-Z)
      SAVE

      TINY=0.001

      IF (TIME .GE. DOYEM .AND. EMERG .EQ. 1.) THEN

         NNI= LIMIT (TINY,1.0, ((NFGMR-NRMR)/NOTNUL(NOPTMR-NRMR)))
      ELSE
         NNI = 0.0

      ENDIF

      RETURN
      END

* ---------------------------------------------------------------------*
*  SUBROUTINE DEATHL                                                   *
*  Purpose: To compute the relative death rate of leaves due to age,*
*           shading amd due to nitrogen stress.                        *
* ---------------------------------------------------------------------*
      SUBROUTINE DEATHL(TIME,DOYEM,TSUM,TSUMAN,RDRTMP,RDRSHM,LAI,
     $                  LAICR,WLVG,RDRNS,NNI,SLA,
     $                 RDRDV,RDRSH,RDR,DLV,DLVS,DLVNS,DLAIS,DLAINS,DLAI)
      IMPLICIT REAL (A-Z)
      SAVE

      IF (TSUM .LT. TSUMAN)THEN
         RDRDV = 0.
      ELSE
         RDRDV = RDRTMP
      ENDIF

      RDRSH = MAX(0.,RDRSHM * (LAI-LAICR) / LAICR)
      RDR   = MAX(RDRDV, RDRSH)

      IF (TIME .GE. DOYEM .AND. NNI .LT. 1.) THEN
          DLVNS   = WLVG *RDRNS * (1.-NNI)
          DLAINS  = DLVNS * SLA
      ELSE
          DLVNS   = 0.
          DLAINS  = 0.
      ENDIF

      DLVS  = WLVG * RDR
      DLAIS = LAI *  RDR

      DLV   = DLVS + DLVNS
      DLAI  = DLAIS + DLAINS

      RETURN
      END
* ---------------------------------------------------------------------*

ENDJOB

