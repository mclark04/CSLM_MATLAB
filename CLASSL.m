function [HLAK, LLAK, BLAK, NLAK, TLAK, T0, HDPTH, LKICEH, SNICEH,ROFICEH,...
    SNO, RHOSNO, TSNOW, ALBSNO, WSNOW, CDH, CDM, QSENS, TFLUX, QEVAP, EVAP,...
    QFLUX, EVPPOT, EVAPB, GT, QSURF, DRAG, ST, SU, SV, SQ, SH, QLWAVG,...
    ALVS, ALIR,  FSGL, FLGL, HFSL, HEVL, HMFL, HTCL, FSGS, FLGS, HFSS,...
    HEVS, HMFN, HTCS, PCPL, PCPN, QFL, QFN, ROFN, FICE, FLS, GZEROSL,...
    EXPW, DTEMP, TKE, DELU, GRED, RHOMIX,...QSWINV, QSWINI,...    
    QLWIN, UWIND,...
    VWIND, TA, QA, RHOAIR, PADRY, PRES, CSZ, ZREFM, ZREFH, ZDIAGM, ZDIAGH,...
    TR, TS, RHOSNI, RADJ, ASVDAT, ASIDAT, FSDB, FSFB, FSSB, REFSN,...R, S,  
    BCSN, JL, NLAKMAX, ISLFD, IZREF, ITG, IALS, NBS,... ILG, IL1, IL2,
    ISNOALB, IGL,  TSED] = ...IRSTRT, NSTEP, IYEAR, IDAY, IHOUR, IMIN,
 CLASSL(HLAK, LLAK, BLAK, NLAK, TLAK,...
    T0, HDPTH, LKICEH, SNICEH, ROFICEH, ...%CLASSL
    SNO, RHOSNO, TSNOW, ALBSNO, WSNOW,...%CDH, CDM, QSENS,TFLUX, QEVAP, EVAP, QFLUX, EVPPOT, EVAPB, GT, QSURF, DRAG,...
    ST, SU,SV, SQ, SH,...%QLWAVG, ALVS, ALIR, FSGL, FLGL, HFSL, HEVL, HMFL, HTCL, FSGS, FLGS, HFSS, HEVS, HMFN, HTCS, PCPL, PCPN,...
    QFL,...%QFN, ROFN, FICE,FLS, GZEROSL,...    
    EXPW, DTEMP, TKE, DELU, GRED, RHOMIX,...
    QSWINV, QSWINI,QLWIN, UWIND, VWIND, TA, QA,...
    RHOAIR, PADRY, PRES, CSZ, ZREFM, ZREFH,...
    ZDIAGM, ZDIAGH, R, TR, S, TS, RHOSNI, RADJ,...
    ASVDAT, ASIDAT, FSDB, FSFB, FSSB, REFSN, BCSN,...
    ILG, IL1, IL2, JL, NLAKMAX, ISLFD, IZREF,...
    ITG, IALS, NBS, ISNOALB, IGL, IRSTRT, NSTEP, IYEAR, IDAY, IHOUR, IMIN,...
    TSED,TFREZ,RHOICE,RHOW,HCPW,SPHICE,CLHMLT,VMIN,...)%Added because Globals (and below
    TKEMIN,DELZLK,GRAV,DELSKIN,VKC,CLHVAP,SPHAIR,EMSW,SBC,TCW,DELT,...
    TKECN,TKECF,TKECE,TKECS,TKECL,SPHW,DHMAX,HDPTHMIN,DUMAX,DELMAX,...
    DELMIN,TCICE,HCPICE,CGRAV,CPD,CKARM,...
    fID_71,fID_72,fID_73,fID_74,fID_75,fID_76,fID_77,fID_82)

%CLASSL
% !-------------------------------------- LICENCE BEGIN ------------------------------------
% !Environment Canada - Atmospheric Science and Technology License/Disclaimer,
% !                     version 3; Last Modified: May 7, 2008.
% !This is free but copyrighted software; you can use/redistribute/modify it under the terms
% !of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
% !version 3 or (at your option) any later version that should be found at:
% !http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
% !
% !This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
% !without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
% !See the above mentioned License/Disclaimer for more details.
% !You should have received a copy of the License/Disclaimer along with this software;
% !if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
% !CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
% !
% !-------------------------------------- LICENCE END --------------------------------------
% =====================================================================
% C     * FEB 15/18 - M.MACKAY    ESSENTIALLY IDENTICAL TO CLIMATE MODEL VERSION.
% C     *                         INPUT PARAMETERS IN CALL SLIGHTLY DIFFERENT
% C     * NOV 09/17 - M.LAZARE.   DEFINE DUMMY ARRAY "ZPONDL", INITIALIZE
% C     *                         IT TO ZERO AND PASS TO TSOLVE IN
% C     *                         APPROPRIATE PLACE.
% C     *                         THIS IS TO SUPPORT COMPATABILITY WITH
% C     *                         THE CHANGES ON THE LAND SIDE FOR INCLUDING
% C     *                         PONDING.                    
% C     * JAN 21/17 - D.VERSEGHY. INCLUDE EFFECTS OF SNOW ICE PRODUCTION
% C     *                         IN ROFN AND HTCS DIAGNOSTICS.
% C     * DEC 19/16 - M.LAZARE/   - DEFINE INTERNAL VARIABLES {THLIQ,THLMIN,DELZW}
% C     *                           TO SATISFY REVISED CALL TO TSOLVE SUPPORTING
% C     *                           "WLOST" BUGFIX.
% C     *                         - INITIALIZE SNOW-RELATED VARIABLES TO ZERO
% C     *                           AT THE END OF LOOP 120, TO AVOID GETTING NAN'S
% C     *                           LATER ON.
% C     *                         - DEFINE {ALVSG,ALIRG} BASED ON ICE ALBEDOS
% C     *                           REQUIRED AS INPUT TO SNOALBA.
% C     *                         - CALCULATE {ALBW,ALBI} HERE (ALBI NEEDED
% C     *                           TO DEFINE {ALVSG,ALIRG}) AND PASS IN TO
% C     *                           TSOLVL, INSTEAD OF HAVING THEM LOCAL IN
% C     *                           TSOLVL, FOR CONSISTENCY SAKE.
% C     *                         - BUGFIX!!! RMELT->RALB IN CALL SNOALBW
% C     *                                     AND RMELT REMOVED AS ARRAY.
% C     * MAY 05/16 - D.VERSEGHY. CLIMATE MODEL VERSION; RESTORE LINEAR VARIATION
% C     *                         OF FRACTIONAL ICE COVER
% C     * APR 11/16 - M.MACKAY.   THERMAL EXPANSION OF ICE BUG CORRECTED
% C     *                         ICE DRAFT,FREEBOARD NOW INCLUDED
% C     *                         SW ATTENUATION THROUGH LEADS INCLUDED
% C     * JAN  25/16 - M.MACKAY.  FRACTIONAL ICE COVER INCLUDED BUT SET
% C     *                         TO 0 OR 100% FOR NOW
% C     * SEP 09/15 - D.VERSEGHY. ADD CALLS TO CLASS SUBROUTINES AND
% C     *                         SUPPLEMENTARY CODE FOR MODELLING OF
% C     *                         SNOW ON LAKE ICE; ADDITIONAL DIAGNOSTIC 
% C     *                         CALCULATIONS.
% C     * APR 13/15 - D.VERSEGHY. RE-ORDERING OF SUBROUTINE CALL;
% C     *                         COSMETIC CHANGES TO CODE AND NAMING.
% C     * FEB 1/2013 - M.MACKAY.  SOLVES FOR LAKE TEMPERATURE PROFILES,
% C     *                         MIXED LAYER DEPTH, AND FLUXES.
% C     *                        	SFC ENERGY BALANCE COMPUTED OVER
% C     *                         SKIN OF FINITE WIDTH
%     C=======================================================================
%       IMPLICIT NONE
% C
% C ----* LAKE MODEL VARIABLES *----------------------------------------
% C
%       INTEGER NLAKMAX
%       INTEGER,DIMENSION(ILG) :: NLAK
%       REAL,DIMENSION(ILG) :: HLAK, LLAK, BLAK, T0, HDPTH,LKICEH,SNICEH,
%      1                      EXPW,DTEMP,TKE,DELU,GRED,RHOMIX,TSED,ROFICEH
%       REAL,DIMENSION(ILG) :: SNO, RHOSNO, TSNOW, ALBSNO, WSNOW
%       REAL,DIMENSION(ILG,NLAKMAX) :: TLAK, QFLX, FFLX
% C
% C ----* INTEGER CONSTANTS *-------------------------------------------
% C
%       INTEGER ISNOW,ISLFD,IZREF,ITG,ILG,IL1,IL2,JL,IALS,NBS,ISNOALB,IGL 
%       INTEGER IRSTRT,NSTEP,IYEAR,IDAY,IHOUR,IMIN
% C
% C ----* INPUT FIELDS *------------------------------------------------
% C
%       REAL,DIMENSION(ILG) :: QSWINV, QSWINI, QLWIN, UWIND, VWIND, TA, 
%      1                       QA, RHOAIR, PADRY, PRES, CSZ, ZREFM, ZREFH,
%      2                       ZDIAGM, ZDIAGH, R, TR, S, TS, RHOSNI, RADJ
% C
%       REAL   ASVDAT(ILG),  ASIDAT(ILG),  REFSN (ILG),  BCSN  (ILG)     
% C
%       REAL   FSDB(ILG,NBS), FSFB(ILG,NBS), FSSB(ILG,NBS)
% C
% C ----* DIAGNOSTIC OUTPUT FIELDS *-------------------------------------
% C
%       REAL,DIMENSION(ILG) :: FSGL, FLGL, HFSL, HEVL, HMFL, HTCL, DRAGL,
%      1                       FSGS, FLGS, HFSS, HEVS, HMFN, HTCS, DRAGS,
%      2                       PCPL, QFL, PCPN, QFN, ROFN,
%      3                       CDH, CDM, QSENS, TFLUX, QEVAP, EVAP, QFLUX,
%      4                       EVPPOT, EVAPB, GT, QSURF, DRAG,
%      5                       ST, SU, SV, SQ, SH, QLWAVG, ALVS, ALIR,
%      6                       ILMO, HBL, UE, FTEMP, FVAP, OSW 
% C
% C ----* INTERNAL ARRAYS *----------------------------------------------
% C
%       REAL,DIMENSION(ILG) :: Q0,DISS,BFLX,FQU,FSHEAR,FENTRA,TRAN,
%      1                       CQ1A, CQ1B, CQ2A, CQ2B, CQ3A, CQ3B,
%      2                       CQ1BI, CQ2BI, CQ3BI, FICE, CDHL, CDML,
%      3                       HRPCPL,HSPCPL,HCONV, ZPONDL
%       REAL,DIMENSION(ILG) :: TCMIX, VA, QSWIN, G0, F0, USTAR, LKICEH0
%       REAL,DIMENSION(NLAKMAX) :: TLAK0
% C
%       REAL   FLS   (ILG),  HCPSNO(ILG),  ZSNOW (ILG),
%      1       ALBW  (ILG),  ALBI  (ILG)
% C
%       REAL   ALVSSN(ILG),  ALIRSN(ILG),  ALVSSC(ILG),  ALIRSC(ILG),     
%      1       ALVSG (ILG),  ALIRG (ILG),  TRSNOWC(ILG),
%      2       ALVSL (ILG),  ALIRL (ILG)
% C                                                                       
%       REAL   GCONSTS(ILG), GCOEFFS(ILG), CPHCHS(ILG),  TCSNOW(ILG),
%      1       ZRSLDM(ILG),  ZRSLDH(ILG),  ZRSLFM(ILG),  ZRSLFH(ILG),
%      2       ZDSLM (ILG),  ZDSLH (ILG),  ZOSCLM(ILG),  ZOSCLH(ILG),
%      3       ZOMLNS(ILG),  ZOELNS(ILG),  ZOM   (ILG),  ZOH   (ILG),
%      4       TVIRTA(ILG),  TPOTA (ILG),  CRIBS (ILG),  CEVAP (ILG), 
%      5       TSTART(ILG),  FCOR  (ILG),  PCP   (ILG),  FTOT  (ILG),
%      6       TSURF (ILG),  QSURFL(ILG),  TSNBOT(ILG),  FSCALL(ILG),
%      7       RPCN  (ILG),  TRPCN (ILG),  SPCN  (ILG),  TSPCN (ILG)
% C
%       REAL   QMELT (ILG),  RHOMAX(ILG),
%      1       RALB  (ILG),  WLOST (ILG),  SN    (ILG),  SL    (ILG),
%      2       RN    (ILG),  RL    (ILG),
%      2       TRUNOF(ILG),  RUNOFF(ILG),  TOVRFL(ILG),  OVRFLW(ILG)
% C
%       REAL   HTC   (ILG,IGL), THLIQ (ILG,IGL), THLMIN (ILG,IGL),
%      1       DELZW (ILG,IGL)
% C
%       INTEGER             IWATER(ILG),   ISAND (ILG,IGL)
% C                                                                                 
% C     * TSOLVE OUTPUT ARRAYS.                                                  
% C                                                                       
%       REAL QSWNS (ILG),    QLWOS (ILG),    QTRANSL(ILG),   QSENSS(ILG), 
%      1     QEVAPS(ILG),    EVAPS (ILG),    TZEROS(ILG),    QZEROS(ILG), 
%      2     GSNOW (ILG),    GZEROSL(ILG),   QMELTS(ILG),    CDHS  (ILG), 
%      3     CDMS  (ILG),    RIBS  (ILG),    CFLUXS(ILG),    FTEMPS(ILG), 
%      4     FVAPS (ILG),    ILMOS (ILG),    UES   (ILG),    HS    (ILG)  
% C                                                                       
%       INTEGER          ITERCT(ILG,6,50),   IEVAP (ILG)                  
% C                                                                       
% C     * RADIATION BAND-DEPENDENT ARRAYS.                                          
% C                                                                       
%       REAL TRSNOWG(ILG,NBS), ALSNO(ILG,NBS), 
%      1     TRTOP  (ILG,NBS)                                             
% C                                                                       
% C     * INTERNAL WORK ARRAYS FOR TSOLVE.                                           
% C                                                                       
%       REAL TSTEP (ILG),    TVIRTS(ILG),    EVBETA(ILG),    Q0SAT (ILG), 
%      1     RESID (ILG),    DCFLXM(ILG),    CFLUXM(ILG),    WZERO (ILG),
%      2     A     (ILG),    B     (ILG),                                 
%      3     LZZ0  (ILG),    LZZ0T (ILG),    FM    (ILG),    FH    (ILG) 
% C                                                                       
%       INTEGER              ITER  (ILG),    NITER (ILG),    JEVAP (ILG), 
%      1                     KF    (ILG)                                  
% C                                                                       
% C ----* COMMON BLOCK PARAMETERS *--------------------------------------
% C
%       REAL DELT,TFREZ,RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN,TCW,TCICE,
%      1     TCSAND,TCCLAY,TCOM,TCDRYS,RHOSOL,RHOOM,HCPW,HCPICE,HCPSOL,
%      2     HCPOM,HCPSND,HCPCLY,SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,
%      3     TCGLAC,CLHMLT,CLHVAP,DELTA,CGRAV,CKARM,CPD,AS,ASX,CI,BS,
%      4     BETA,FACTN,HMIN,ANGMAX                          
%       REAL TKECN,TKECF,TKECE,TKECS,TKECL,HDPTHMIN,
%      1     TKEMIN,DELMAX,DELMIN,EMSW,DELZLK,DELSKIN,DHMAX,DUMAX
% [   DELT,TFREZ,RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN,TCW,TCICE,...
%     TCSAND,TCCLAY,TCOM,TCDRYS,RHOSOL,RHOOM,HCPW,HCPICE,HCPSOL,...
%     HCPOM,HCPSND,HCPCLY,SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,...
%     TCGLAC,CLHMLT,CLHVAP,DELTA,CGRAV,CKARM,CPD,AS,ASX,CI,BS,...
%     BETA,FACTN,HMIN,ANGMAX,TKECN,TKECF,TKECE,TKECS,TKECL,HDPTHMIN,...
%     TKEMIN,DELMAX,DELMIN,EMSW,DELZLK,DELSKIN,DHMAX,DUMAX]=...
%     makeBlanks(1,1); %cut this because it was overwriting passed vars
% C
% C ----* CLASS COMMON BLOCKS *------------------------------------------
% C
%       COMMON /CLASS1/ DELT,TFREZ                                       
%       COMMON /CLASS2/ RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN
%       COMMON /CLASS3/ TCW,TCICE,TCSAND,TCCLAY,TCOM,TCDRYS,
%      1                RHOSOL,RHOOM
%       COMMON /CLASS4/ HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,
%      1                SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,
%      2                TCGLAC,CLHMLT,CLHVAP
%       COMMON /PHYCON/ DELTA,CGRAV,CKARM,CPD
%       COMMON /CLASSD2/ AS,ASX,CI,BS,BETA,FACTN,HMIN,ANGMAX              
% 
%       COMMON /LAKECON/ TKECN,TKECF,TKECE,TKECS,HDPTHMIN,        
%      2                 TKEMIN,DELMAX,DELMIN,EMSW,DELZLK,DELSKIN,
%      3                 DHMAX,TKECL,DUMAX
% C
% C ----* LOCAL VARIABLES *---------------------------------------------
% C
%       INTEGER I,J,K,JBOT,JTOP,NLEV,JMIX,JMIXP1
%       REAL Z,DAY,DECL,HOUR,TKE_OLD,ZBOT,TBOT,WEDB,
%      1     Z2,TMIX,ZMIX,HTOP,HBOT,ZTOP,TTOP,
%      2     RHO1,RHO2,TC1,TC2,TBAR,ATTEN1,ATTEN2,ATTEN3,DELTHRM,HCAP,
%      3     EAVAIL,EHEAT,ECOOL,EFLUX,ZICE,NEWICE,TC,SNOLIM,ICELIM,
%      4     FACTM,FACTH,RATIOM,RATIOH
%       REAL TSEDT,TSEDB,TCSED,HCPSED,DSED,FFLXSED,ICEMASS
%       REAL ICECAP,PORE,MASSL,SNO0,ALPHA,ETA,ZSNOW0,ICE0,BBLFAC
%       REAL MASS,TC0,RHO0,RHOIW,ICETOP,ICEBOT,ICETOP0,ICEBOT0,ROFICE,
%      1     HFREZ,HWARM,HFLX,HFLXI,HFLXS,HCPRAT,QMLTS,QSED,XXX,SHELTER
%       INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(15,307)
% C
% C=======================================================================
%  SEDIMENT and WIND SHELTERING PROPERTIES  ------------------------------------------------
% C     DSED=10.0       !THICKNESS OF SEDIMENT LAYER (m)  (STD)
% C     TSEDB=TFREZ+6.0 !CONSTANT TEMP AT BOTTOM OF SEDIMENT (K)  (STD)
% C     HCPSED=2.13E6  !VOLUMETRIC HEAT CAPACITY OF SEDIMENT      (STD)
% C     TCSED=2.5      !THERMAL CONDUCTIVITY OF SEDIMENT (0 FOR ADIABATIC LOWER BC)
% C     QSED=0.06	     !FRACTION OF NET SW THAT REACHES SEDIMENTS (0 to turn off)
% C     SHELTER=1.0    !wind sheltering factor reduces surface drag (1.0 turns off)
% C----------------------------------------------------------------------

%C
%C ----* INTERNAL ARRAYS *-
%C
% [Q0,DISS,BFLX,FQU,FSHEAR,FENTRA,TRAN,...
%     CQ1A, CQ1B, CQ2A, CQ2B, CQ3A, CQ3B,...
%     CQ1BI, CQ2BI, CQ3BI, CDHL, CDML,...

[HRPCPL,HSPCPL,HCONV, ZPONDL,FICE] = makeZeros(ILG,1);%
[TCMIX, VA, QSWIN,USTAR, LKICEH0] = makeZeros(ILG,1);% G0, F0, 
[TLAK0] = makeZeros(NLAKMAX,1);
[HCPSNO, ZSNOW, ALBW, ALBI,FLS] = makeZeros(ILG,1);%FLS,
[ALVSSN, ALIRSN, ...ALVSSC, ALIRSC,
    ALVSG, ALIRG,... TRSNOWC,
    ALVSL, ALIRL] = makeZeros(ILG,1);
[GCONSTS, GCOEFFS, CPHCHS,  TCSNOW,...
    ZRSLDM,  ZRSLDH,  ZRSLFM,  ZRSLFH,...
    ZDSLM ,  ZDSLH ,  ZOSCLM,  ZOSCLH,...
    ZOMLNS,  ZOELNS,  ZOM   ,  ZOH   ,...
    TVIRTA,  TPOTA ,  CRIBS ,  ...CEVAP ,
    TSTART,  FCOR  ,  PCP   ,  FTOT  ,...
    TSURF ,  TSNBOT,  FSCALL,...QSURFL,  
    RPCN  ,  TRPCN ,  SPCN  ,  TSPCN ] =  makeZeros(ILG,1);
[RHOMAX,...QMELT ,
    RALB  ,  WLOST ,  SN    ,  SL    ,...
    RN    ,  RL    ,...
    TRUNOF,  RUNOFF,  TOVRFL,  OVRFLW] = makeZeros(ILG,1);
[HTC, THLIQ, THLMIN, DELZW] = makeZeros(ILG,IGL);
[IWATER]= makeZeros(ILG,1);
%[ISAND] = makeZeros(ILG,IGL);
%C
%C     * TSOLVE OUTPUT ARRAYS.
%C
[QLWOS ,    QSENSS,... QSWNS ,    QTRANSL,   
	QEVAPS,    EVAPS ,    TZEROS,    QZEROS,... 
	GSNOW ,    QMELTS,    CDHS  ,...
	CDMS  ,    RIBS  ,    CFLUXS,    FTEMPS,...
 	FVAPS ,    ILMOS ,    UES   ,    HS    ] = makeZeros(ILG,1); %   GZEROSL,     
ITERCT = single(nan(ILG,6,50));
%IEVAP = makeZeros(ILG,1);
%C
%C     * RADIATION BAND-DEPENDENT ARRAYS.
%C
[TRSNOWG, ALSNO, TRTOP] = makeZeros(ILG,NBS);  
%C
%C     * INTERNAL WORK ARRAYS FOR TSOLVE.   
%C
[   TSTEP ,    TVIRTS,    EVBETA,    Q0SAT ,... 
	RESID ,    DCFLXM,    CFLUXM,    WZERO ,...
	A     ,    B     ,...
	LZZ0  ,    LZZ0T ,    FM    ,    FH    ,... 
    ITER  ,    NITER ,    ... JEVAP ,
	KF    ] = makeZeros(ILG,1);   
%C
%C ----* LOCAL VARIABLES *---------------------------------------------
%C
% [I,J,K,JBOT,JTOP,NLEV,JMIX,JMIXP1,...
%     Z,DAY,DECL,HOUR,TKE_OLD,ZBOT,TBOT,WEDB,...    
%     Z2,TMIX,ZMIX,HTOP,HBOT,ZTOP,TTOP,...
%     RHO1,RHO2,TC1,TC2,TBAR,ATTEN1,ATTEN2,ATTEN3,DELTHRM,HCAP,...
%     EAVAIL,EHEAT,ECOOL,EFLUX,ZICE,NEWICE,TC,SNOLIM,ICELIM,...
%     FACTM,FACTH,RATIOM,RATIOHRHO1,...
%     TSEDT,TSEDB,TCSED,HCPSED,DSED,FFLXSED,ICEMASS...
%     ICECAP,PORE,MASSL,SNO0,ALPHA,ETA,ZSNOW0,ICE0,BBLFAC,...
%     MASS,TC0,RHO0,RHOIW,ICETOP,ICEBOT,ICETOP0,ICEBOT0,ROFICE,...
%     HFREZ,HWARM,HFLX,HFLXI,HFLXS,HCPRAT,QMLTS,QSED,XXX,SHELTER]= makeZeros(1,1);
% C
% C ----* DIAGNOSTIC OUTPUT FIELDS *-------------------------------------
% C
[ DRAGS, HTCL, PCPN, PCPL, HTCS, HMFN, QFN, ROFN,GZEROSL,...
    HMFL,OSW,CDM,CDH,DRAG,QSURF,ALVS,ALIR,FSGS,FLGS,HFSS,...
    HEVS,QLWAVG,QSENS,TFLUX,QEVAP,EVAP,QFLUX,EVPPOT,EVAPB,GT] = makeZeros(ILG,1);

% [FSGL, FLGL, HFSL, HEVL, HMFL, HTCL, DRAGL,...
%     FSGS, FLGS, HFSS, HEVS, HMFN, HTCS, DRAGS,...
% 	PCPL, QFL, PCPN, QFN, ROFN,...
% 	CDH, CDM, QSENS, TFLUX, QEVAP, EVAP, QFLUX,...
% 	EVPPOT, EVAPB, GT, QSURF, DRAG,...
% 	ST, SU, SV, SQ, SH, QLWAVG, ALVS, ALIR,...
% 	ILMO, HBL, UE, FTEMP, FVAP, OSW] = makeZeros(ILG,1); %This has not been checked to see if they overwrite anything
% C
%[DRAGL] = makeZeros(ILG,1);

%Other yet to be catagorized arrays to pre-allocate
[QFLX,FFLX] = makeZeros(ILG,NLAKMAX);


% Definitions starting on line 221 of the fortran code
DSED=single(10.0);%       !THICKNESS OF SEDIMENT LAYER (m)  (STD)
TSEDB=TFREZ+6.0;% !CONSTANT TEMP AT BOTTOM OF SEDIMENT (K)  (STD)
HCPSED=single(2.13E6);%   !VOLUMETRIC HEAT CAPACITY OF SEDIMENT     (STD)
TCSED=single(0.0);%       !THERMAL CONDUCTIVITY OF SEDIMENT (0 FOR ADIABATIC LOWER BC)
QSED=single(0.0);%        !FRACTION OF NET SW THAT REACHES SEDIMENTS (0 to turn off)
% C     QSED=0.06
SHELTER=single(1.0);%     !wind sheltering factor reduces surface drag (1.0 turns off)
% C ---------------------------------------------------------------------
ISNOW=single(1);
SNOLIM=single(0.10);
RHOIW=RHOICE/RHOW;
% C======================================================================
% C    * CLASS SMALL LAKE MODULE
% C-----------------------------------------------------------------------
% C    * INITIAL CONDITIONS (NSTEP=1)
% C    * EQUATION OF STATE FROM FARMER AND CARMACK (1981,JPO), BUT 
% C    * NEGLECTS SALINITY AND PRESSURE EFFECTS
% C
T1=T0;
if (NSTEP == 1 || IRSTRT == 1)     
  for I=IL1:IL2
      LKICEH(I)=0.0; %    !INITIAL ICE COVER SET TO ZERO
      SNICEH(I)=0.0; %    !INITIAL SNOW ICE  SET TO ZERO
      ROFICEH(I)=0.0; %    !INITIAL runoff ICE  SET TO ZERO
      T0(I)=TLAK(I,1);%   !INITIAL SKIN TEMP SET TO FIRST LAYER TEMP
      NLEV=NLAK(I);%
      TSED(I)=TLAK(I,NLEV);%   !INITIAL SEDIMENT TEMP SET TO LOWEST LAYER TEMP
      TKE(I)=TKEMIN;
      DELU(I)=0.0;
% C ---* 
% C ---* INITIAL MIXED LAYER DEPTH ESTIMATED BASED ON INITIAL T PROFILE
% C ---* 
      for J=1:NLAK(I)-1
          JMIX=J;
          TTOP=TLAK(I,J);
          TBOT=TLAK(I,J+1);
          if (TTOP-TBOT > 1.0) 
              %Exit loop
              break;
          end
      end
      HDPTH(I)=DELZLK*JMIX;
      DTEMP(I)=TLAK(I,JMIX)-TLAK(I,JMIX+1);% !see Spigel et al 1986
% C ---* 
% C ---* INITIAL MIXED LAYER TEMP (CELSIUS), DENSITY, EXPANSIVITY, 
% C ---* AND REDUCED GRAVITY
% C ---* 
     TCMIX(I)=( (TLAK(I,1)+TLAK(I,JMIX))/2.0) -TFREZ; 
     %CALL EQNST(EXPW(I),RHOMIX(I),TCMIX(I),0.0)
     [EXPW(I),RHOMIX(I)]=EQNST(TCMIX(I),0.0);
     GRED(I) = GRAV*abs(RHOW-RHOMIX(I))/RHOW;
     if (GRED(I) < 0) %this is impossible because of the abs, unless G is pointed upsidedown. 
         error('RUNLAKE -- In CLASSL GRED < 0');
     end
  end
end
% C-----------------------------------------------------------------------
% C     * ICE COVER AND SNOW PARAMETERS AT BEGINNING OF TIMESTEP
% C---------------------------------------
for I=IL1:IL2
      ZPONDL(I)=0.0;
      LKICEH0(I)=LKICEH(I);
      ICELIM=LLAK(I)*1.0E-5;        
      FICE(I)=MIN(ICELIM,LKICEH(I))/ICELIM;   
      FTOT(I)=1.0;
      if(SNO(I)>0.0)
          ZSNOW(I)=SNO(I)/RHOSNO(I);
          if(ZSNOW(I)>=(SNOLIM-0.00001))
              FLS(I)=1.0;
          else
              FLS(I)=ZSNOW(I)/SNOLIM;
              ZSNOW(I)=SNOLIM;
              WSNOW(I)=WSNOW(I)/FLS(I);
          end
      else
          ZSNOW(I)=0.0;
          FLS(I)=0.0;
      end
      if(S(I)>0.0)   
          if(LKICEH(I)>=DELSKIN)
              if(FLS(I)>0.0)    
                  SN(I)=S(I)*FICE(I)/FLS(I);
              else
                  SN(I)=S(I);
              end
              SL(I)=(1.0-FICE(I))*S(I);
          else
              SN(I)=0.0;
              SL(I)=S(I);
          end
      else
          SN(I)=0.0;
          SL(I)=0.0;
      end
      if(R(I)>0.0) 
          if(FLS(I)>0.0)
              RN(I)=R(I);
              RL(I)=(1.0-FLS(I))*R(I);
          else
              RN(I)=0.0;
              RL(I)=R(I);
          end
      else
          RN(I)=0.0;
          RL(I)=0.0;
      end
      if(FLS(I)>0)
          FSCALL(I)=FLS(I);
          TSTART(I)=TSNOW(I);
      else
          FSCALL(I)=FICE(I);
          TSTART(I)=TFREZ;
      end
      PCPN(I) =FLS(I)*RHOW*RN(I)+FSCALL(I)*RHOSNI(I)*SN(I);
      PCPL(I) =RHOW*RL(I)+RHOSNI(I)*SL(I);
      HTCS(I)=0.0;
      HMFN(I)=0.0;
      QFN(I)=0.0;
      ROFN(I)=0.0;
      GZEROSL(I)=0.0;
      QMELTS(I)=0.0;
      HRPCPL(I)=(TR(I)+TFREZ-T0(I))*RL(I)*HCPW;
      HSPCPL(I)=(TS(I)+TFREZ-T0(I))*SL(I)*SPHICE*RHOSNI(I);
      HCONV(I)=CLHMLT*SL(I)*RHOSNI(I);
      HTCL(I)=HRPCPL(I)+HSPCPL(I)-HCONV(I);
% c
% C         * INITIALIZE OTHER SNOW-RELATED VARIABLES.
% C
      SPCN  (I)=0;
      TSPCN (I)=0;
      RPCN  (I)=single(0);
      TRPCN (I)=0;
      CDMS  (I)=0;
      CDHS  (I)=0;
      DRAGS (I)=0;
      TZEROS(I)=0;
      QZEROS(I)=0;
      ALVSSN(I)=0;
      ALIRSN(I)=0;
      QLWOS (I)=0;
      QSENSS(I)=0;
      QEVAPS(I)=0;
      EVAPS (I)=0;
      HCPSNO(I)=0;
      TCSNOW(I)=0;
      GSNOW (I)=0;          
end
%     -----------------------------------------------------------------------
% C     * WIND SPEED AND PRECIP
% C
for I=IL1:IL2
  VA(I)=MAX(VMIN,sqrt(UWIND(I)*UWIND(I)+VWIND(I)*VWIND(I)));
  FCOR(I)=2.0*7.29E-5*sin(RADJ(I));
  PCP(I)=R(I)+S(I);
end
% C
% C     * DEFINE ICE AND WATER ALBEDOS, USED IN TSOLVL.
% C     * STARTING "GROUND ALBEDOS" FOR ROUTINE SNOALBA.
% C
for I=IL1:IL2
    ALBW(I)=0.045/MAX(CSZ(I),0.1);%          !std value for water
    %white space
    ALBI(I)=0.08+0.44*(LKICEH(I))^0.28;%    !thin ice albedo (Vavrus et al 1996)
    ALBI(I)=MIN(ALBI(I),0.44);
    %white space
    if(LKICEH(I)>0.0)
        ALVSG(I)=ALBI(I);
        ALIRG(I)=ALBI(I);
    else
        ALVSG(I)=0.0;
        ALIRG(I)=0.0;
    end
end
% C
% C-----------------------------------------------------------------------
% C     * SNOW PACK ENERGY AND WATER BALANCE
% C
% 1-Be careful transcribing, some varaibles change name in the subroutine
[ALVSSN,ALIRSN,ALVSSC,ALIRSC,TRSNOWC, ALSNO, TRSNOWG, ...ALBSNO, FSDB, FSFB,
        ...RHOSNO, REFSN,BCSN,SNO,CSZ,ZSNOW,FLS ,ASVDAT,ASIDAT,ALVSG, ALIRG, 
        ] = ...ILG,IGL,IL1,IL2,JL,IALS,NBS,ISNOALB
    SNOALBA(ALVSSN,ALIRSN,ALBSNO, ALSNO, TRSNOWG,...ALVSSC,ALIRSC,TRSNOWC, 
        ZSNOW,FLS,ASVDAT,ASIDAT, ...FSDB, FSFB, RHOSNO, REFSN,BCSN,SNO,CSZ,
        IL1,IL2,JL,IALS,NBS,ISNOALB);% ALVSG, ALIRG, ILG,IGL,

[GCOEFFS,GCONSTS,CPHCHS,TCSNOW,HCPSNO,IWATER,ZRSLDM,ZRSLDH,ZRSLFM,ZRSLFH,...
        ZDSLM,ZDSLH,ZOSCLM,ZOSCLH,ZOMLNS,ZOELNS,ZOM,ZOH,TVIRTA,TPOTA,CRIBS,...
        DRAGS,CEVAP,IEVAP,ISAND,...FLS,ZSNOW,TSNOW,RHOSNO,WSNOW,ZREFM,ZREFH,
        ]=...,JL,ILG,IL1,IL2,IGL,ZDIAGM,ZDIAGH,TA,QA,VA,IZREF
    TLSPREP(GCOEFFS,GCONSTS,CPHCHS,TCSNOW,HCPSNO,IWATER,ZRSLDM,ZRSLDH,ZRSLFM,...
        ZRSLFH,ZDSLM,ZDSLH,ZOSCLM,ZOSCLH,ZOMLNS,ZOELNS,ZOM,ZOH,TVIRTA,TPOTA,CRIBS,...
        DRAGS,FLS,ZSNOW,TSNOW,RHOSNO,WSNOW,ZREFM,ZREFH,ZDIAGM,...CEVAP,IEVAP,ISAND,
        ZDIAGH,TA,QA,VA,IZREF,IL1,IL2,IGL,...ILG,JL,
        CGRAV,CPD,CKARM,HCPICE,RHOICE,HCPW,RHOW,CLHVAP,CLHMLT);%added globals

[...ISNOW,FLS,
        QSWNS,QLWOS,QTRANSL,QSENSS,QEVAPS,EVAPS,...
        TZEROS,QZEROS,GSNOW,QMELTS,CDHS,CDMS,RIBS,CFLUXS,...
        FTEMPS,FVAPS,ILMOS,UES,HS,ZRSLFH,ZRSLFM,ZOH,ZOM,FCOR,ILG,...%passed through FLZSURFZ
        ...QLWIN,TPOTA,QA,VA,PADRY,RHOAIR,
        ...ALVSSN,ALIRSN,CRIBS,CPHCHS,CEVAP,TVIRTA,
        ...ZOSCLH,ZOSCLM,
        ...GCONSTS,GCOEFFS,TSTART,PCP,TRSNOWG,FSSB,ALSNO,
        ...THLIQ,THLMIN,DELZW,RHOSNO,ZSNOW,ZPONDL,
        ITERCT,...IWATER,IEVAP,ISAND,      
        ...ISLFD,ITG,IGL,IL1,IL2,JL,NBS,ISNOALB,
        TSTEP,TVIRTS,EVBETA,Q0SAT,RESID,...    
        DCFLXM,CFLUXM,WZERO,TRTOP,...A,B,
        LZZ0,LZZ0T,FM,FH,ITER,NITER,JEVAP,KF] = ...
    TSOLVE(ISNOW,FLS,...
        QLWOS,QSENSS,QEVAPS,EVAPS,...QSWNS,QTRANSL,
        TZEROS,QZEROS,GSNOW,QMELTS,CDHS,CDMS,RIBS,CFLUXS,...
        FTEMPS,FVAPS,ILMOS,UES,HS,ZRSLFH,ZRSLFM,ZOH,ZOM,FCOR,...%passed through FLZSURFZ
        QLWIN,TPOTA,QA,VA,PADRY,RHOAIR,...
        ALVSSN,ALIRSN,CRIBS,CPHCHS,CEVAP,TVIRTA,...
        ZOSCLH,ZOSCLM,...
        GCONSTS,GCOEFFS,TSTART,TRSNOWG,FSSB,ALSNO,...PCP,
        THLIQ,THLMIN,DELZW,RHOSNO,ZSNOW,ZPONDL,...
        IWATER,IEVAP,ITERCT,ISAND,...      
        ISLFD,ITG,ILG,IL1,IL2,NBS,ISNOALB,...IGL,JL,
        TSTEP,TVIRTS,EVBETA,Q0SAT,RESID,...    
        DCFLXM,CFLUXM,WZERO,TRTOP,...A,B,
        LZZ0,LZZ0T,FM,FH,ITER,NITER,KF,...JEVAP,
        DELT,TFREZ,RHOW,SBC,SPHAIR,GRAV,VKC);   
    
[GSNOW,TSNOW,WSNOW,RHOSNO,QMELTS,GZEROSL,...
        TSNBOT,HTCS,HMFN,QFN,EVAPS,RPCN,TRPCN,SPCN,TSPCN,...
        HCPSNO...GCONSTS,GCOEFFS,T0,ZSNOW,TCSNOW,QTRANSL,
        ...RN,TR,SN,TS,TZEROS,RHOSNI,
        ] = ...FLS,DELSKIN,ILG,IL1,IL2,JL
    TLSPOST(GSNOW,TSNOW,WSNOW,RHOSNO,QMELTS,GZEROSL,...
        TSNBOT,HTCS,HMFN,EVAPS,...QFN,RPCN,TRPCN,SPCN,TSPCN,
        T0,ZSNOW,TCSNOW,HCPSNO,...GCONSTS,GCOEFFS,QTRANSL,
        RN,TR,SN,TS,TZEROS,RHOSNI,...
        FLS,DELSKIN,IL1,IL2,...ILG,JL,
        RHOW,DELT,TFREZ,CLHMLT,HCPICE,RHOICE,HCPW);%added these in because globals
    
[RHOSNO,ZSNOW,HCPSNO,TSNOW,EVAPS,QFN,QFL,HTCS,...
        WLOST,TRUNOF,RUNOFF,TOVRFL,OVRFLW,...
        WSNOW] = ...FLS,RPCN,SPCN,RHOSNI,,ILG,IL1,IL2,JL
    SNOVAP(RHOSNO,ZSNOW,HCPSNO,TSNOW,EVAPS,QFN,QFL,HTCS,...
        WLOST,TRUNOF,RUNOFF,TOVRFL,OVRFLW,...
        FLS,RPCN,SPCN,RHOSNI,WSNOW,IL1,IL2,...ILG,JL,
        DELT,TFREZ,RHOW,HCPICE,RHOICE,HCPW);%globals
    
[ZSNOW,TSNOW,QMELTS,RPCN,TRPCN,GZEROSL,RALB,...
        HMFN,HTCS,HTC,HCPSNO,WSNOW...FLS,RHOSNO,
        ] =...ISAND,IGL,ILG,IL1,IL2,JL
    TMELT(ZSNOW,TSNOW,QMELTS,RPCN,TRPCN,GZEROSL,RALB,...
        HMFN,HTCS,HTC,FLS,HCPSNO,RHOSNO,WSNOW,...
        ISAND,IL1,IL2,...IGL,ILG,JL,
        DELT,TFREZ,RHOW,HCPICE,RHOICE,HCPW,CLHMLT);
[RPCN,TRPCN,ZSNOW,TSNOW,RHOSNO,HCPSNO,WSNOW,...
        HTCS,HMFN,PCPL,ROFN] = ...,FLS,ILG,IL1,IL2,JL
    SNINFL(RPCN,TRPCN,ZSNOW,TSNOW,RHOSNO,HCPSNO,WSNOW,...
        HTCS,HMFN,PCPL,ROFN,FLS,ILG,IL1,IL2,JL,...
        TFREZ,DELT,HCPICE,RHOICE,HCPW,RHOW,CLHMLT);%added globals
[ALBSNO,RHOSNO,ZSNOW,HCPSNO,...TSNOW,
        RHOMAX...FLS,SPCN,RALB,WSNOW,ISAND,
        ]= ...ILG,IGL,IL1,IL2,JL
    SNOALBW(ALBSNO,RHOSNO,ZSNOW,HCPSNO,TSNOW,...              
        FLS,SPCN,RALB,WSNOW,RHOMAX,...                 
        IL1,IL2,...ISAND,ILG,IGL,JL,
        DELT,RHOW,HCPICE,RHOICE,HCPW);%added globals            
[ALBSNO,TSNOW,RHOSNO,ZSNOW,HCPSNO,HTCS] = ...FSCALL,SPCN,TSPCN,RHOSNI,WSNOW,ILG,IL1,IL2,JL                                                                                  
    SNOADD(ALBSNO,TSNOW,RHOSNO,ZSNOW,HCPSNO,HTCS,...
        FSCALL,SPCN,TSPCN,RHOSNI,WSNOW,IL1,IL2,...ILG,JL
        TFREZ,DELT,HCPICE,RHOICE,HCPW,RHOW);%added globals

%    
% ----- set up file writing
%
%   This is very different than the way FORTRAN handles files.

format6010='%4d %3d %3d %3d %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %7.3f %7.3f %7.3f %10.3e %10.3e %10.3e %10.3e %10.3e'; 
% 6010  FORMAT(I4,1X,3(I3,1X),8F8.2,2F7.3,5E10.3)
format6020='%4d %3d %3d %3d';
for i=1:200
    format6020=[format6020,' %7.2f'];
end
%format6020=[format6020,'\n'];
% 6020  FORMAT(I4,1X,3(I3,1X),200F7.2)
%format6030='%4d %3d %3d %3d %10.2e %10.2e %10.2e %10.2e %5.2f %6.2f %3d %7.2f %7.2f %7.2f %10.2e %10.2e';
format6030='%4d %3d %3d %3d %10.2e %10.2e %10.2e %10.2e %5.2f %6.2f %3d %7.2f %7.2f %7.2f %10.2e %10.2e %10.2e';
% 6030  FORMAT(I4,1X,3(I3,1X),4E10.2,F5.2,F6.2,I3,3F7.2,2E10.2)
format6040='%4d %3d %3d %3d %10.3e %10.3e %10.3e %7.2f %7.2f %7.2f %8.3f %8.3f %10.3e %10.3e';
% 6040  FORMAT(I4,1X,3(I3,1X),3E10.3,3F7.2,F8.3,F8.3,E10.3)
format6050='%4d %3d %3d %3d %10.3e %10.3e %10.3e %10.3e %10.3e';
% 6050  FORMAT(I4,1X,3(I3,1X),5E10.3)
format6060='%4d %3d %3d %3d %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e';
% 6060  FORMAT(I4,1X,3(I3,1X),8E10.3)
format6070='%4d %3d %3d %3d %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f';
% 6070  FORMAT(I4,1X,3(I3,1X),6F7.1) 

% baseDir='D:\Dropbox\PostDoc\UofS_Course\Project\TestOutput';
% fID_71=fopen([baseDir,'\LAKE.of1'],'A');
% fID_72=fopen([baseDir,'\LAKE.of2'],'A');
% fID_73=fopen([baseDir,'\LAKE.of3'],'A');
% fID_74=fopen([baseDir,'\LAKE.of4'],'A');
% fID_75=fopen([baseDir,'\LAKE.of5'],'A');
% fID_76=fopen([baseDir,'\LAKE.of6'],'A');
% fID_77=fopen([baseDir,'\LAKE.of7'],'A');
%fID_78=fopen([baseDir,'\LAKE.of8'],'a');%Not written to here
%fID_79=fopen([baseDir,'\LAKE.of9'],'a');%Not written to here
   
    
%C
for I=IL1:IL2%DO I=IL1,IL2
    if(ZSNOW(I)>0.0)
        TSNOW(I)=TSNOW(I)+TFREZ;
        SNO(I)=ZSNOW(I)*RHOSNO(I)*FSCALL(I);
        WSNOW(I)=WSNOW(I)*FLS(I);
    else
        TSNOW(I)=0.0;
        RHOSNO(I)=0.0;
        SNO(I)=0.0;
        WSNOW(I)=0.0;
    end
    % C-----------------------------------------------------------
    % C Melt last bit of snow if below threshold or no ice exists
    % C
    if(SNO(I)>0.0 && (SNO(I)<1.0E-6 || LKICEH(I) <=0.01))
        ROFN(I)=ROFN(I)+(SNO(I)+WSNOW(I))/DELT;
        PCPL(I)=PCPL(I)+(SNO(I)+WSNOW(I))/DELT;
        SL(I)=SL(I)+SNO(I)/(RHOSNO(I)*DELT);
        RL(I)=RL(I)+WSNOW(I)/(RHOW*DELT);
        HTCS(I)=HTCS(I)-TSNOW(I)*(SPHICE*SNO(I)+SPHW*WSNOW(I))/DELT;
        TSNOW(I)=0.0;
        RHOSNO(I)=0.0;
        SNO(I)=0.0;
        WSNOW(I)=0.0;
    end
    % C---------------------------------------------------------------
    % C Heat capacity calculation for use if necessary by: (1) heat 
    % C added from snowmelt/rain runoff; and (2) snow ice production
    % C Any heat added goes into layer 1, not the skin layer.
    ICEBOT=RHOIW*LKICEH0(I);%   !bottom of floating ice
    ZTOP=DELSKIN;%             !top of first layer
    ZBOT=DELSKIN + DELZLK;%    !bottom of first layer
    if (ICEBOT >= ZBOT) 
        HCAP=HCPICE;
    elseif (ICEBOT <= ZTOP)
        HCAP=HCPW;
    else 
        Z=ICEBOT-ZTOP;
        HCAP=(Z*HCPICE + (DELZLK-Z)*HCPW)/DELZLK;
    end
    % C-----------------------------------------------------------
    % C Add heat of snow melt/rain reaching bottom of snowcover to ice
    % C Heat added to first layer below skin
    % C ROFN cools then runs off - ie does not add to LKICEH
    % C
    % C-----------------------------------------------------------
    % C RUNOFF ICE PRODUCTION/HEATING TURNED OFF: JUNE 2016, M.MACKAY
    % C
    %TLAKOld1=TLAK(1,:);
    if (ROFN(I)>0.0 && LKICEH0(I)>0.0) 
        ROFICE=0.0;%             !runoff ice production turned off
        % Cmdm         ROFICE=DELT*ROFN(I)/RHOICE    !ice production this timestep
        % Cmdm         ROFICEH(I)=ROFICEH(I)+ROFICE  !cumulative diagnostic
        % C Cool ROFN to TFREZ, then add heat to ice layer below skin
        HWARM=HCPW*(TRPCN(I)-0.0)*DELT*ROFN(I)/RHOW;
        HFREZ=RHOICE*CLHMLT*ROFICE;
        TLAK(I,1)=TLAK(I,1) + (HFREZ+HWARM)/(HCAP*DELZLK);
        % Cmdm         LKICEH(I)=LKICEH(I)+ROFICE
    end
    % C======================================================================
    % C SNOW ICE
    % c---------
    % C Produce snow-ice if weight of snow cover sufficient
    % C Pore volume in snow fills with water then freezes
    % C
    if(SNO(I)>0.0) 
        ROFN(I)=ROFN(I)+SNO(I)/DELT;
        HTCS(I)=HTCS(I)-HCPSNO(I)*TSNOW(I)*SNO(I)/(RHOSNO(I)*DELT);
    end
    BBLFAC=single(1.0);%                         !non-bubble fraction in snow ice
    ICECAP =LKICEH0(I)*(RHOW-RHOICE);%   !snow holding capacity of ice
    ZBOT = MIN(10.0, (NLAK(I)-1.)*DELZLK);%   !limit snow ice production to 10m or depth of lake (less 0.5m)
    %TLAKOld2=TLAK(1,:);
    if (SNO(I)>ICECAP && LKICEH0(I)>0.0 && LKICEH0(I)<=ZBOT )
        ICE0=LKICEH0(I);%            !initial ice thickness
        SNO0=SNO(I);%                     !initial snow mass
        ZSNOW0=SNO(I)/RHOSNO(I);%         !initial mean snow depth 
        PORE=(RHOICE-RHOSNO(I))/RHOICE;%  !pore volume in snow
        ALPHA=(RHOW*PORE*BBLFAC + RHOICE*(1.0-PORE))/RHOICE;
        ETA=(RHOW-RHOICE)/RHOSNO(I);
        %C Final snow cover mass equals new ice holding capacity
        SNO(I)=RHOSNO(I)*ETA*(ICE0+ALPHA*ZSNOW0)/(1.0+ALPHA*ETA);
        LKICEH(I)=SNO(I)/(RHOSNO(I)*ETA);
        %C Mass of liquid water that freezes in slush layer
        MASSL=RHOW*PORE*BBLFAC*(SNO0-SNO(I))/RHOSNO(I);
        %C Cumulative snow-ice diagnostic
        SNICEH(I)=SNICEH(I)+LKICEH(I)-ICE0;
%         C Add heat of fusion to snow pack
%         C First warm snow to TFREZ, then freeze lake water
%         C Lake water is assumed to be at TFREZ already (ie comes from layer that
%         C contains ICEBOT, which is at TFREZ)
        HWARM=HCPSNO(I)*(TFREZ-TSNOW(I))*(SNO0-SNO(I))/RHOSNO(I);
        HFREZ=CLHMLT*MASSL;
        HFLX=HFREZ-HWARM;
        %C Latent heat goes goes into snowpack unless snow too thin        
        if (SNO(I) >= 5.0E-3) 
        % C Partition latent heat between ice and snow based on heat capacities
        % Cmdm test   HCPRAT=HCPSNO(I)/HCPICE
        % Cmdm test   HFLXS=HFLX*(HCPRAT/(1.0+HCPRAT))
        % Cmdm test   HFLXS=HFLX/2.0
        % Cmdm test   HFLX=0.0   !mdm test
            HFLXS=HFLX;%           !all latent heat goes into snow cover
            HFLXI=HFLX-HFLXS;
            TLAK(I,1)=TLAK(I,1) + HFLXI/(HCAP*DELZLK);% !mdm test
            TSNOW(I)=TSNOW(I)+HFLXS/(HCPSNO(I)*SNO(I)/RHOSNO(I)); 
            if(TSNOW(I)>TFREZ) 
                QMLTS=(TSNOW(I)-TFREZ)*HCPSNO(I)*SNO(I)/RHOSNO(I);% !J/m2
                if (QMLTS > SNO(I)*CLHMLT)
                    %C All snow melts. Excess heat put in lake
                    %print*, "all snow melts during flooding"
                    fprintf('\nall snow melts during flooding');
                    HFLX=QMLTS-SNO(I)*CLHMLT;
                    TLAK(I,1)=TLAK(I,1) + HFLX/(HCAP*DELZLK);
                    SNO(I)=0.0;
                    TSNOW(I)=0.0;
                    RHOSNO(I)=0.0;
                    WSNOW (I)=0.0;
                    ZSNOW (I)=0.0;
                else
                    %C Some snow melts
                    SNO(I)=SNO(I)-QMLTS/CLHMLT;
                    TSNOW(I)=TFREZ;
                end
            end                                                  
        else
            %print*, "no snow: all heat into lake"
            fprintf('\nno snow: all heat into lake');
            TLAK(I,1)=TLAK(I,1) + (HFREZ-HWARM)/(HCAP*DELZLK);
        end
    end
    if(SNO(I)>0.0)
        ROFN(I)=ROFN(I)-SNO(I)/DELT;
        HTCS(I)=HTCS(I)+HCPSNO(I)*TSNOW(I)*SNO(I)/(RHOSNO(I)*DELT);
    end
end%ENDDO
% C
% C-----------------------------------------------------------------------
% C     * LAKE DRAG COEFFICIENTS 
% C
[CDML,CDHL,DRAGL] =...,VA,T0,TA,QA,PRES,ZREFM,ZREFH,FICE,FLS,ILG,IL1,IL2
    DRCOEFL(VA,T0,TA,QA,PRES,ZREFM,ZREFH,FICE,FLS,IL1,IL2,...CDML,CDHL,DRAGL,ILG,
        GRAV,VKC,TFREZ);%added these because gloabl
% dat=[NSTEP,CDML(1),CDHL(1),DRAGL(1),...
%             VA(1),T0(1),TA(1),QA(1),PRES(1),ZREFM(1),ZREFH(1),FICE(1),FLS(1),IL1(1),IL2(1),...
%             GZEROSL(1),QTRANSL(1),HTCL(1),FLS(1),FICE(1),ALBW(1),ALBI(1),...
%             GRAV(1),VKC(1),TFREZ(1)];
% dat=[NSTEP(1),CDML(1),CDHL(1),DRAGL(1),VA(1),T0(1),... 
%     TA(1),QA(1),PRES(1),ZREFM(1),ZREFH(1),FICE(1),FLS(1),...
%     ILG(1),IL1(1),IL2(1)];
% n=length(dat);
% forStr=['%5d'];
% for i=2:n
%      forStr=[forStr,' %12.8e'];
% end
% forStr=[forStr,'\n'];
% fprintf(fID_82,forStr,dat);    
    
% C-----------------------------------------------------------------------
% C     * TOTAL SOLAR RADIATION, AND WATER FRICTION VELOCITY
% C       Reduce momentum transfer due to sheltering
% C
for I=IL1:IL2%DO 150 I=IL1,IL2
    QSWIN(I) = QSWINV(I) + QSWINI(I);
    CDML(I)=CDML(I)*SHELTER;
    if (LKICEH(I) <= 0) 
        USTAR(I)=VA(I)*sqrt(CDML(I)*RHOAIR(I)/RHOW);
    else
        USTAR(I)=0.0;%            !no wind stress through ice
    end
end%150   CONTINUE
% C
% C------------------------
% C REMOVE STATIC INSTABILITY if PRESENT BY MIXING WITH ADJACENT LAYER
% C  - include skin layer if no ice is present
% C  -  Sept 2, 2011 M MacKay
T2=T0;
[T0,TLAK]=...LKICEH,,RHOIW,NLAK,NLAKMAX,ILG,IL1,IL2,IYEAR,IHOUR,IDAY,IMIN
    FREECONV(LKICEH,T0,TLAK,RHOIW,NLAK,IL1,IL2,...NLAKMAX,ILG,IYEAR,IHOUR,IDAY,IMIN,
        TFREZ,DELSKIN,DELZLK);%need this, no globals
    

% C-----------------------------------------------------------------------
% C --- * COMPUTE WATER EXTINCTION COEFFICIENTS
% C
[CQ1A,CQ1B,CQ2A,CQ2B,CQ3A,CQ3B,CQ1BI,CQ2BI,CQ3BI] =...BLAK,IL1,IL2,ILG,
    LKTRANS (BLAK,IL1,IL2);%CQ1A,CQ1B,CQ2A,CQ2B,CQ3A,CQ3B,ILG,CQ1BI,CQ2BI,CQ3BI
% C
% C-----------------------------------------------------------------------
% C --- * COMPUTE SURFACE ENERGY BALANCE FOR CURRENT TIMESTEP AND STEP
% C --- * FORWARD SKIN TEMP T0.  COMPUTE ICE COVER IN SKIN LAYER.
% [TLAK,T0,LKICEH,Q0,F0,QSURFL,... 
%             FSGL,FLGL,HFSL,HEVL,QFL,ALVSL,ALIRL,...
%             QSWIN,QLWIN,CSZ,TA,QA,VA,PRES,RHOAIR,CDHL,...
%             GZEROSL,QTRANSL,HTCL,FLS,FICE,ALBW,ALBI,...
%             CQ1A,CQ1B,CQ2A,CQ2B,CQ3A,CQ3B,...
%             CQ1BI,CQ2BI,CQ3BI,G0,...
%             NLAKMAX,ILG,IL1,IL2] = ...
%     TSOLVL(TLAK,T0,LKICEH,Q0,F0,QSURFL,... 
%             FSGL,FLGL,HFSL,HEVL,QFL,ALVSL,ALIRL,...
%             QSWIN,QLWIN,CSZ,TA,QA,VA,PRES,RHOAIR,CDHL,...
%             GZEROSL,QTRANSL,HTCL,FLS,FICE,ALBW,ALBI,...
%             CQ1A,CQ1B,CQ2A,CQ2B,CQ3A,CQ3B,...
%             CQ1BI,CQ2BI,CQ3BI,G0,...
%             NLAKMAX,ILG,IL1,IL2,...
%             DELZLK,DELSKIN,RHOICE,RHOW,TFREZ,...%needed for globals
%             CLHVAP,SPHAIR,EMSW,SBC,TCW,DELT,HCPW,...%needed for class I think
%             CLHMLT,TCICE,HCPICE);%needed for ice stuff
%Reduced to:
T0old=T0;
LKICEHold=LKICEH;
QFLold=QFL;
T3=T0;
[T0,LKICEH,Q0,F0,QSURFL,FSGL,FLGL,HFSL,...
           HEVL,QFL,ALVSL,ALIRL,G0] = ...
    TSOLVL(TLAK,T0,LKICEH,QFL,ALVSL,ALIRL,...
           QSWIN,QLWIN,TA,QA,VA,PRES,RHOAIR,CDHL,...
           GZEROSL,QTRANSL,HTCL,FLS,FICE,ALBW,ALBI,...
           CQ1A,CQ1B,CQ2A,CQ2B,CQ3A,CQ3B,...
           CQ1BI,CQ2BI,CQ3BI,IL1,IL2,...
           DELZLK,DELSKIN,RHOICE,RHOW,TFREZ,...
           CLHVAP,SPHAIR,EMSW,SBC,TCW,DELT,HCPW,...
           CLHMLT,TCICE,HCPICE);

% fprintf(fID_82,'%f %f %f %f %f %f %f %f %f %f %f %f %f\n',T0(1),LKICEH(1),...
%     Q0(1),F0(1),QSURFL(1),FSGL(1),FLGL(1),HFSL(1),HEVL(1),QFL(1),...
%     ALVSL(1),ALIRL(1),G0(1));
dat=[NSTEP,TLAK(1),T0(1),T0old(1),QFLold(1),LKICEH(1),LKICEHold(1),QFL(1),ALVSL(1),ALIRL(1),...%10
            QSWIN(1),QLWIN(1),TA(1),QA(1),VA(1),PRES(1),RHOAIR(1),CDHL(1),...%8 (18)
            GZEROSL(1),QTRANSL(1),HTCL(1),FLS(1),FICE(1),ALBW(1),ALBI(1),...%8 (26)
            CQ1A(1),CQ1B(1),CQ2A(1),CQ2B(1),CQ3A(1),CQ3B(1),...%6 (32)
            CQ1BI(1),CQ2BI(1),CQ3BI(1),IL1(1),IL2(1),...%5 (37)
            DELZLK(1),DELSKIN(1),RHOICE(1),RHOW(1),TFREZ(1),...%5 (42)
            CLHVAP(1),SPHAIR(1),EMSW(1),SBC(1),TCW(1),DELT(1),HCPW(1),...%7 (49)
            CLHMLT(1),TCICE(1),HCPICE(1),Q0(1),F0(1),QSURFL(1),FSGL(1),FLGL(1),HFSL(1),...%9 (58)
            HEVL(1),G0(1)];%2 (60)
n=length(dat);
forStr=['%5d'];
for i=2:n
     forStr=[forStr,' %17.8e'];
end
forStr=[forStr,'\n'];
fprintf(fID_82,forStr,dat);
% C
% C-----------------------------------------------------------------------
% C --- * COMPUTE SOLAR FLUX (QFLX) FOR CURRENT TIMESTEP
% C        include attenuation through ice (including leads) if present
% C
for I=IL1:IL2%DO 200 I=IL1,IL2
    ICEBOT=RHOIW*LKICEH(I);
    ICETOP=LKICEH(I)-ICEBOT;
    for J=1:NLAK(I)%DO 180 J=1,NLAK(I)
        Z=DELSKIN + DELZLK*J; 
        if (LKICEH(I) <= 0.0) %       !NO ICE 
            ATTEN1=CQ1B(I)*Z;
            ATTEN2=CQ2B(I)*Z;
            ATTEN3=CQ3B(I)*Z;
        else 
            if (ICEBOT > Z) %        !Z inside ice
                ATTEN1=FICE(I)*(CQ1BI(I)*(Z+ICETOP)) + (1.-FICE(I))*CQ1B(I)*Z;
                ATTEN2=FICE(I)*(CQ2BI(I)*(Z+ICETOP)) + (1.-FICE(I))*CQ2B(I)*Z;
                ATTEN3=FICE(I)*(CQ3BI(I)*(Z+ICETOP)) + (1.-FICE(I))*CQ3B(I)*Z;
            else
                ATTEN1=FICE(I)*(CQ1BI(I)*LKICEH(I)+CQ1B(I)*(Z-ICEBOT)) + (1.-FICE(I))*CQ1B(I)*Z;
                ATTEN2=FICE(I)*(CQ2BI(I)*LKICEH(I)+CQ2B(I)*(Z-ICEBOT)) + (1.-FICE(I))*CQ2B(I)*Z;
                ATTEN3=FICE(I)*(CQ3BI(I)*LKICEH(I)+CQ3B(I)*(Z-ICEBOT)) + (1.-FICE(I))*CQ3B(I)*Z;
            end
        end
        %C----------- Remove fraction of shortwave flux from water column to be used ---
        %C            to heat sediment if sediment heating option activated.
        %C            Only for ice-free conditions
        %C
        if ( LKICEH(I) <= 0.0)%        !NO ICE 
            QFLX(I,J)=(1.0-QSED)*FSGL(I)*(CQ1A(I)*exp(-ATTEN1) +...
                CQ2A(I)*exp(-ATTEN2) + CQ3A(I)*exp(-ATTEN3) );
        else
            QFLX(I,J)=FSGL(I)*(CQ1A(I)*exp(-ATTEN1) + ...
                CQ2A(I)*exp(-ATTEN2) + CQ3A(I)*exp(-ATTEN3) );
        end
    end%180       CONTINUE
end%200   CONTINUE
% C
% C-----------------------------------------------------------------------
% C --- * COMPUTE CONDUCTIVE HEAT FLUX (FFLX) FOR CURRENT TIMESTEP
% C --- * THERMAL CONDUCTIVITY IS WEIGHTED AVERAGE OF ICE AND WATER
% C --- * if ICE IS PRESENT IN LAYER
% C
for I=IL1:IL2%DO 300 I=IL1,IL2
    ICEBOT=RHOIW*LKICEH(I);
    ICETOP=LKICEH(I)-ICEBOT;
    for J=1:NLAK(I)%DO 250, J=1,NLAK(I)-1;
        ZTOP=DELSKIN + DELZLK*(J -1);
        ZBOT=DELSKIN + DELZLK*J;
        if (ICEBOT >= ZBOT) 
            FFLX(I,J)=(-TCICE/DELZLK)*(TLAK(I,J+1)-TLAK(I,J));
        elseif (ICEBOT < ZBOT && ICEBOT > ZTOP)
            TC=((ICEBOT-ZTOP)*TCICE + (ZBOT-ICEBOT)*TCW)/DELZLK;
            FFLX(I,J)=(-TC/DELZLK)*(TLAK(I,J+1)-TLAK(I,J));
        else
            FFLX(I,J)=(-TCW/DELZLK)*(TLAK(I,J+1)-TLAK(I,J));
        end
    end%250     CONTINUE
    NLEV=NLAK(I);
    %C ---* COMPUTE THERMAL FLUX AT LAKE BOTTOM INTO SEDIMENT
    %C ---* ASSUMES LAKE DOES NOT FREEZE TO BOTTOM
    %C ---* (ADIABATIC BC RECOVERED WHEN TCSED=0)
    TSEDT=( (TCSED*TSED(I)/DSED)+(TCW*TLAK(I,NLEV)/DELZLK) )/( (TCSED/DSED)+(TCW/DELZLK) );
    FFLX(I,NLEV)= (-2.0*TCW/DELZLK)*(TSEDT-TLAK(I,NLEV));
    %CCC     FFLX(I,NLEV)=0.0	!adiabatic lower B.C.
end%300   CONTINUE
% C
% C-----------------------------------------------------------------------
% C ---* COMPUTE TEMPERATURE PROFILE (TLAK) BEFORE TURBULENT MIXING 
% C    * FOR NEXT TIMESTEP.  
% C    * HEAT CAPACITY OF LAYER IS WEIGHTED AVERAGE OF WATER AND ICE if
% C    * NECESSARY
% C
oldTLAK=TLAK(1,:);
oldT0=T0(1);
for I=IL1:IL2%DO 400 I=IL1,IL2
    ICEBOT=RHOIW*LKICEH(I);
    ICETOP=LKICEH(I)-ICEBOT;
    ICEBOT0=RHOIW*LKICEH0(I);
    ICETOP0=LKICEH0(I)-ICEBOT;
    % C
    % C --COLUMN  --- TEMP BEFORE MELT/FREEZE
    %TLAKOld4=TLAK(1,:);
    for J=1:NLAK(I)%DO 310, J=1,NLAK(I)
        TLAK0(J)=TLAK(I,J);
        ZTOP=DELSKIN + DELZLK*(J -1);
        ZBOT=DELSKIN + DELZLK*J;
        if (ICEBOT >= ZBOT) 
            HCAP=HCPICE;
        elseif (ICEBOT <= ZTOP) 
            HCAP=HCPW;
        else 
            Z=ICEBOT-ZTOP;
            HCAP=(Z*HCPICE + (DELZLK-Z)*HCPW)/DELZLK;
        end
        if (J == 1) 
            EFLUX=F0(I)-FFLX(I,1)+Q0(I)-QFLX(I,1);
        else
            EFLUX=FFLX(I,J-1)-FFLX(I,J)+QFLX(I,J-1)-QFLX(I,J);
        end
        TLAK(I,J)=TLAK(I,J) + (DELT/(DELZLK*HCAP))*EFLUX;
    end%310     CONTINUE
    % C --UPDATE SEDIMENT TEMPERATURE
    NLEV=NLAK(I);
    FFLXSED=(-2.0*TCSED/DSED)*(TSEDB-TSED(I));
    TSED(I)=TSED(I) + (DELT/(DSED*HCPSED))*(FFLX(I,NLEV)-FFLXSED+QSED*FSGL(I));
    % C  
    % C   ICE GROWTH OR DECAY
    % C 
    %TLAKOld5=TLAK(1,:);
    for J=1:NLAK(I)%DO 320, J=1,NLAK(I)
        ZTOP=DELSKIN + DELZLK*(J -1);
        ZBOT=DELSKIN + DELZLK*J;
        if (J == 1) 
            EFLUX=F0(I)-FFLX(I,1)+Q0(I)-QFLX(I,1);
        else
            EFLUX=FFLX(I,J-1)-FFLX(I,J)+QFLX(I,J-1)-QFLX(I,J);
        end
        % C
        % C --- FREEZING --------------------
        % C
        if (EFLUX < 0.0 && ICEBOT0 < ZBOT && ICEBOT0 >= ZTOP) 
            % C ----Net energy flux used to lower T to TFREZ
            if (TLAK0(J) > TFREZ)
                ECOOL=(ZBOT-ICEBOT0)*HCPW*(TLAK0(J)-TFREZ)/DELT;
            else
                ECOOL=0.0;
            end
            % C ----Remaining energy flux (if any) used to freeze ice
            EAVAIL=EFLUX+ECOOL;
            if (EAVAIL < 0.0) 
                NEWICE=-(DELT/(RHOICE*CLHMLT))*EAVAIL;
                LKICEH(I)=LKICEH(I)+NEWICE;
                ICEBOT=RHOIW*LKICEH(I);
                TLAK(I,J)=TFREZ;
                % C-----LIMIT ICE GROWTH TO THE CURRENT LAYER
                if (ICEBOT > ZBOT) 
                    EHEAT=(RHOICE*CLHMLT*(ICEBOT-ZBOT))/DELT;
                    TLAK(I,J)=TLAK(I,J) - (EHEAT*DELT)/(DELZLK*HCPICE);
                    LKICEH(I)=ZBOT/RHOIW;
                end
            end
        end
        % C
        % C --- MELTING --------------------
        % C
        %TLAKOld6=TLAK(1,:);
        if (EFLUX > 0.0 && ICEBOT0 > ZTOP)  
        % C ----Net energy flux used to raise T to TFREZ
            ZICE=MIN(DELZLK, ICEBOT0-ZTOP);
            EHEAT=ZICE*HCPICE*(TFREZ-TLAK0(J))/DELT;
            % C ----Remaining energy flux (if any) used to melt ice
            EAVAIL=EFLUX-EHEAT;
            if (EAVAIL > 0.0) 
                NEWICE=-(DELT/(RHOICE*CLHMLT))*EAVAIL;
                LKICEH(I)=LKICEH(I)+NEWICE;
                ICEBOT=RHOIW*LKICEH(I);
                TLAK(I,J)=TFREZ;
                % C-----LIMIT ICE MELT TO THE CURRENT LAYER
                if (ICEBOT < ZTOP) 
                    EHEAT=RHOICE*CLHMLT*(ZTOP-ICEBOT)/DELT;
                    TLAK(I,J)=TFREZ + (EHEAT*DELT)/(DELZLK*HCPW);
                    LKICEH(I)=ZTOP/RHOIW;
                end
            end
        end
    end%320     CONTINUE
end%400   CONTINUE
   
% C-----------------------------------------------------------------------
% C ---* COMPUTE MIXED LAYER DEPTH (HDPTH), TKE, SHEAR FLOW (DELU)
% C    * FOR NEXT TIMESTEP
% % C
% if NSTEP==467
%     fprintf('');
% end

oldHDPTH=HDPTH;
oldTKE=TKE;
oldDELU=DELU;
[...DTEMP,Q0,NLAK,USTAR,IL1,IL2,ILG,NLAKMAX,
              HDPTH,TKE,DELU,FQU,BFLX,DISS,...EXPW,FSGL,
              FSHEAR,FENTRA,TRAN,...HLAK,LLAK,GRED,
              ...CQ1A,CQ1B,CQ2A,CQ2B,CQ3A,CQ3B,RHOMIX,
              ] = ...FLGL,HFSL,HEVL,LKICEH
    MIXLYR(DTEMP,NLAK,USTAR,IL1,IL2,...Q0,ILG,NLAKMAX,
              HDPTH,TKE,DELU,EXPW,FSGL,...FQU,BFLX,DISS,
              HLAK,LLAK,GRED,...FSHEAR,FENTRA,TRAN,
              CQ1A,CQ1B,CQ2A,CQ2B,CQ3A,CQ3B,RHOMIX,...
              FLGL,HFSL,HEVL,LKICEH,...
              TKECN,TKECF,TKECE,TKECS,TKECL,GRAV,...%Gloabl issues from here
              DELSKIN,SPHW,DELT,DHMAX,DELZLK,HDPTHMIN,...
              TKEMIN,DUMAX);%Gloabl issues
          
% dat=[NSTEP(1),NLAK(1),oldT0(1),oldTLAK(1,1:NLAK(1)),... 
%     T0(1),TLAK(1,1:NLAK(1))];
% dat=[NSTEP(1),oldHDPTH(1),oldTKE(1),oldDELU(1),...
%       FQU(1),BFLX(1),DISS(1),FSHEAR(1),FENTRA(1),TRAN(1),...
%       DTEMP(1),NLAK(1),USTAR(1),IL1,IL2,...Q0,ILG,NLAKMAX,
%       HDPTH(1),TKE(1),DELU(1),EXPW(1),FSGL(1),...FQU,BFLX,DISS,
%       HLAK(1),LLAK(1),GRED(1),...FSHEAR,FENTRA,TRAN,
%       CQ1A(1),CQ1B(1),CQ2A(1),CQ2B(1),CQ3A(1),CQ3B(1),RHOMIX(1),...
%       FLGL(1),HFSL(1),HEVL(1),LKICEH(1),...
%       TKECN(1),TKECF(1),TKECE(1),TKECS(1),TKECL(1),GRAV(1),...%Gloabl issues from here
%       DELSKIN(1),SPHW(1),DELT(1),DHMAX(1),DELZLK(1),HDPTHMIN(1),...
%       TKEMIN(1),DUMAX(1)];
% n=length(dat);
% forStr=['%5d'];
% for i=2:n
%      forStr=[forStr,' %12.8e'];
% end
% forStr=[forStr,'\n'];
% fprintf(fID_82,forStr,dat); 
          
% C
% C-----------------------------------------------------------------------
% C ---* MIX TEMP OVER MIXED LAYER (now mass weighted)
% C
JBOT=nan;
JTOP=nan;
JMIX=nan;
TBOT=nan;
TTOP=nan;
TMIX=nan;
HBOT=nan;
HTOP=nan;


%oldTLAK=nan(1,NLAK(1));%subtle changes to TLAK past here, problem is higher up
% for i=1:NLAK(1)
%     oldTLAK(1,i)=TLAK(1,i);
% end
oldTLAK=TLAK;
oldT0=T0;
RHO1_1=single(0);
RHO2_1=single(0);
RHO1_2=single(0);
RHO1_3=single(0);
RHO1_4=single(0);
T4=T0;
for I=IL1:IL2%DO 600 I=IL1,IL2
    ICEBOT=RHOIW*LKICEH(I);
    %C
    Z2=DELZLK+DELSKIN;
    ZBOT=DELSKIN+DELZLK*NLAK(I);
    if (HDPTH(I) < Z2) %!fully stratified: no mixing
        JMIX=1;
        ZMIX=DELZLK+DELSKIN;
    elseif (HDPTH(I) >= ZBOT) %	!fully mixed
        JMIX=NLAK(I);
        ZMIX=ZBOT;
    else
        for J=2:NLAK(I)%DO 510, J=2,NLAK(I)
            Z=DELSKIN + DELZLK*J;
            if (HDPTH(I) < Z)
                break;%EXIT;
            end
        end%510      CONTINUE
        JMIX=J-1;
        ZMIX=Z-DELZLK;
    end
    % C ---------------------------------------------------------------------------
    % C -- MIXING UNDER ICE ALLOWED (ie temp mixed between bottom of mixed layer
    % C                              and bottom of ice cover)
    % C ---------------------------------------------------------------------------
    if (LKICEH(I) <= 0.0) 
        TC1=T0(I)-TFREZ;
        [~,RHO1a]=EQNST(TC1,0.05);
        TBAR=DELSKIN*RHO1a*T0(I);
        MASS=DELSKIN*RHO1a;
        RHO1_1=RHO1a;
    else
        TBAR=0.0;
        MASS=0.0;
    end
    for J=1:JMIX%DO 520, J=1,JMIX
        ZTOP=DELSKIN + (J-1)*DELZLK;
        ZBOT=ZTOP+DELZLK;
        if (ICEBOT <= ZTOP) 
            TC1=TLAK(I,J)-TFREZ;
            [~,RHO1b]=EQNST(TC1,ZBOT);
            TBAR=TBAR + RHO1b*TLAK(I,J)*DELZLK;
            MASS=MASS + RHO1b*DELZLK;
            RHO1_2=RHO1b;
        end
    end% 520     CONTINUE
    if ((JMIX >= 2) && ((ZMIX-ICEBOT)>=DELZLK) ) 
        TMIX=TBAR/MASS;
    else
        TMIX=TLAK(I,1);
    end
    % C-----------------------------------------------------------------------
    % C MIX TEMPERATURES: include skin layer if no ice
    % C-----------------------------------------------------------------------
    % C--Mix skin with first layer temperature
    % Cmdm    if ( (LKICEH(I)<=0.0) && (JMIX<2) ) THEN
    %TLAKOld7=TLAK(1,:);
    if NSTEP>=472
        %fprintf('here');
    end
    if ( (LKICEH(I)<=0.0) && (JMIX==1) )  
        TC1=T0(I)-TFREZ;
        TC2=TLAK(I,1)-TFREZ;
        [~,RHO1]=EQNST(TC1,0.05);
        [~,RHO2]=EQNST(TC2,0.5);
        T0(I) = (RHO1*DELSKIN*T0(I) + RHO2*DELZLK*TLAK(I,1))/(RHO1*DELSKIN + RHO2*DELZLK);
        TLAK(I,1)=T0(I);
        RHO1_3=RHO1;
        RHO2_1=RHO2;
    elseif ( (LKICEH(I)<=0.0) && (JMIX>=2) ) 
        T0(I)=TMIX;
        for J=1:JMIX%DO 525, J=1,JMIX
            TLAK(I,J)=TMIX;
        end%525       CONTINUE
        % C
        % C--- MIXING UNDER ICE
    elseif ((JMIX >= 2) && ((ZMIX-ICEBOT)>=DELZLK) )
        %Cmdm      print*, "MIXING UNDER ICE"
        for J=1:JMIX%DO 530, J=1,JMIX
            ZTOP=DELSKIN + (J-1)*DELZLK;
            if (ICEBOT <= ZTOP) 
                TLAK(I,J)=TMIX;
            end
        end%530       CONTINUE
    end
    % C======================================================================
    % C ---* COMPUTE TEMPERATURE, DENSITY, EXPANSIVITY OF MIXED LAYER WATER
    % C      EQN OF STATE FROM FARMER AND CARMACK, 1982
    % C
    TCMIX(I)=TMIX-TFREZ;
%     tmp=NLAK(I);
%     tmp(tmp>JMIX+1)=JMIX+1;
%     JMIXP1=tmp;
    JMIXP1=MIN(NLAK(I),JMIX+1);
    [EXPW(I),RHOMIX(I)] = EQNST(TCMIX(I),HDPTH(I));
    [~,RHO1]            = EQNST(TLAK(I,JMIXP1)-TFREZ,HDPTH(I));
    GRED(I)             = GRAV*abs(RHO1-RHOMIX(I))/RHOW;
    RHO1_4=RHO1;
    % Cmdm   GRED(I) = GRAV*ABS(RHO1-RHOMIX(I))/RHO1
    if (GRED(I)<0.) 		
        XIT('CLASSL',-2);
    end
        % C
        % C-----------------------------------------------------------------------
        % C ---* MIX TEMPERATURE IN THERMOCLINE LAYER
        % C          DISABLED FOR NOW: JAN 2013 (M.MACKAY)
        % C
        % C--Compute thermocline layer thickness (DELTHRM)
        % C--Limit allowed values
        % C
    if (LKICEH(I) <= 0.0) 
        if (GRED(I) > 0.0)
            DELTHRM=0.3*DELU(I)*DELU(I)/GRED(I);
        else
            DELTHRM=0.0;
        end
        if (DELTHRM>DELMAX)
            DELTHRM=DELMAX;
        end
        if (DELTHRM<DELMIN)
            DELTHRM=DELMIN;
        end

        HTOP=HDPTH(I)-(DELTHRM/2.);
        HBOT=HDPTH(I)+(DELTHRM/2.);
        for J=1:NLAK(I)%DO 540, J=1,NLAK(I)
            Z=DELZLK*J;
            if (HTOP < Z)
                break;%EXIT
            end
        end%540    CONTINUE
        JTOP=J-1;
        if (JTOP<1)
            JTOP=1;
        end
        ZTOP=DELZLK*JTOP;
        TTOP=TLAK(I,JTOP);

        for J=1:NLAK(I)%DO 550, J=1,NLAK(I)
            Z=DELZLK*J;
            if (HBOT < Z)
                break;%EXIT
            end
        end%550    CONTINUE
        JBOT=J-1;
        if (JBOT<1) 
            JBOT=1;
        end
        ZBOT=DELZLK*JBOT;
        TBOT=TLAK(I,JBOT);
        if (JBOT<JTOP)
             XIT('CLASSL',-3)
        end
        
        % C-DISABLED
        % C-       if (JBOT>JTOP) THEN
        % C-        DO 560, J=JTOP,JBOT
        % C-          Z=DELZLK*J
        % C-          TLAK(I,J)=TTOP+((Z-ZTOP)*(TBOT-TTOP)/(ZBOT-ZTOP))
        % C-560     CONTINUE
        % C-       end
        %
        % C--Compute temperature jump beneath mixed layer
        % C--
        if (JMIX < NLAK(I))
            DTEMP(I)=TMIX-TLAK(I,JBOT);% 	!see Spigel et al 1986
        else
            DTEMP(I)=0.;
        end
    else
        DTEMP(I)=0.;
        DELTHRM=0.;
    end
end%600   CONTINUE
% dat=[NSTEP(1),NLAK(1),oldT0(1),oldTLAK(1,1:NLAK(1)),... 
%     T0(1),TLAK(1,1:NLAK(1))];
% dat=[NSTEP,JBOT,JTOP,JMIX,JMIXP1,...
%     oldT0(1),oldTLAK(1,1),T0(1),TLAK(1,1),...
%     RHOIW,LKICEH(1),Z2,ZBOT,ZTOP,ZMIX,...
%     TTOP,TBOT,TBAR,...
%     MASS,TCMIX(1),EXPW(1),RHOMIX(1),...
%     RHO1,GRED(1),DELTHRM,DELMAX,DELMIN,...
%     DTEMP(1),ICEBOT,HTOP,HBOT,HDPTH(1)];
%dat=[NSTEP,RHO1_1,RHO2_1,RHO1_2,RHO1_3,RHO1_4];
% dat=[NSTEP,T1(1),T2(1),T3(1),T4(1),T0(1)];
% 
% n=length(dat);
% %forStr=['%6d %5d %5d %5d %5d'];
% forStr=['%6d'];
% %for i=6:n
% for i=2:n
%      forStr=[forStr,' %16.6e'];
% end
% forStr=[forStr,'\n'];
% fprintf(fID_82,forStr,dat);    

% C------------------------
% C--Compute Wedderburn number diagnostic
% C--
for I=IL1:IL2%DO 605 I=IL1,IL2
    if  (LKICEH(I) <= 0.0 && USTAR(I) > 0.0 ) 
        WEDB=GRED(I)*HDPTH(I)*HDPTH(I)/(LLAK(I)*USTAR(I)*USTAR(I));
    else
        WEDB=-999;%	!because USTAR=0
    end
    % C--Compute heat of melting/freezing for energy balance diagnostic
    % C--
    HMFL(I) = (LKICEH(I)-LKICEH0(I))*RHOICE*CLHMLT/DELT;
    OSW(I) = QSWIN(I) - FLS(I)*(QSWNS(I)-QTRANSL(I)) - FSGL(I);
    % C
    % C--RESET SNICEH if ICE HAS COMPLETELY MELTED; ENSURE CONSISTENCY BETWEEN 
    % C  FICE AND LKICEH
    if  (LKICEH(I) <= 0.0 ) 
        SNICEH(I)=0.0;
    end
    ICELIM=LLAK(I)*1.0E-5;        
    FICE(I)=MIN(ICELIM,LKICEH(I))/ICELIM;   
    
    % 
    % C--
    % C-----------------------------------------------------------------------
    % C ---* WRITE OUTPUT FOR THIS TIMESTEP
    % C
    fprintf(fID_71,format6010,IYEAR,IDAY,IHOUR,IMIN,G0(I),F0(I),FSGL(I),...
        Q0(I),FLGL(I),HFSL(I),HEVL(I),T0(I),LKICEH(I),SNO(I),...
    	RHOSNI(I)*S(I),RHOW*R(I),OSW(I),SNO(I)/RHOSNO(I),SNICEH(I));
    fprintf(fID_71,'\n');
    %WRITE(71,6010) IYEAR,IDAY,IHOUR,IMIN,G0(I),F0(I),FSGL(I),
    %1   Q0(I),FLGL(I),HFSL(I),HEVL(I),T0(I),LKICEH(I),SNO(I),
    %2   RHOSNI(I)*S(I),RHOW*R(I),OSW(I),SNO(I)/RHOSNO(I),SNICEH(I)
    % C    2   RHOSNI(I)*S(I),RHOW*R(I),OSW(I),ZSNOW(I),ROFICEH(I)
    % C    2   RHOSNI(I)*S(I),RHOW*R(I),OSW(I),ZSNOW(I),SNICEH(I)
    dat=nan(1,NLAK(I));
    for j=1:NLAK(I)
        dat(1,j)=TLAK(I,j)-TFREZ;
    end
    fprintf(fID_72,format6020,IYEAR,IDAY,IHOUR,IMIN,dat);
    fprintf(fID_72,'\n');
    %WRITE(72,6020) IYEAR,IDAY,IHOUR,IMIN,
    %1                 (TLAK(I,J)-TFREZ,J=1,NLAK(I))
    fprintf(fID_73,format6030,IYEAR,IDAY,IHOUR,IMIN,FQU(I),BFLX(I),DISS(I),...
        TKE(I),DELU(I),HDPTH(I),JMIX,ZMIX,TMIX,DTEMP(I),...
    	FSHEAR(I),FENTRA(I));
    fprintf(fID_73,'\n');
    %WRITE(73,6030) IYEAR,IDAY,IHOUR,IMIN,FQU(I),BFLX(I),DISS(I),
    %1                 TKE(I),DELU(I),HDPTH(I),JMIX,ZMIX,TMIX,DTEMP(I),
    %2                 FSHEAR(I),FENTRA(I)
    %Cmdm  WRITE(74,6040) IYEAR,IDAY,IHOUR,IMIN,USTAR(I),GRED(I),WEDB,
    fprintf(fID_74,format6040,IYEAR,IDAY,IHOUR,IMIN,USTAR(I),EXPW(I),WEDB,...
    	HDPTH(I),DELTHRM,DELU(I),ZREFM(I),VA(I),CDML(I),GRED(I));
    fprintf(fID_74,'\n');
    %WRITE(74,6040) IYEAR,IDAY,IHOUR,IMIN,USTAR(I),EXPW(I),WEDB,
    %1                 HDPTH(I),DELTHRM,DELU(I),
    %2                 ZREFM(I),VA(I),CDML(I)
    fprintf(fID_75,format6050,IYEAR,IDAY,IHOUR,IMIN,TSNOW(I)-TFREZ,WSNOW(I),...
    	RHOSNO(I), ZSNOW(I),FLS(I));
    fprintf(fID_75,'\n');
    %WRITE(75,6050) IYEAR,IDAY,IHOUR,IMIN,TSNOW(I)-TFREZ,WSNOW(I),
    %1               RHOSNO(I), ZSNOW(I),FLS(I)
    fprintf(fID_76,format6060,IYEAR,IDAY,IHOUR,IMIN,QSWNS(I),QTRANSL(I),...
    	QLWIN(I)-QLWOS(I),QSENSS(I),QEVAPS(I),QMELTS(I),GZEROSL(I),CFLUXS(I));
    fprintf(fID_76,'\n');
    %WRITE(76,6060) IYEAR,IDAY,IHOUR,IMIN,QSWNS(I),QTRANSL(I),
    %1               QLWIN(I)-QLWOS(I),QSENSS(I),QEVAPS(I),
    %2               QMELTS(I),GZEROSL(I),CFLUXS(I)
    fprintf(fID_77,format6070,IYEAR,IDAY,IHOUR,IMIN,QSWIN(I),FSGL(I),QFLX(I,2),...
    	FFLX(I,NLEV),TSED(I)-TFREZ,TLAK(I,NLEV)-TFREZ);
    fprintf(fID_77,'\n');
    %WRITE(77,6070) IYEAR,IDAY,IHOUR,IMIN,QSWIN(I),FSGL(I),QFLX(I,2),
    %1               FFLX(I,NLEV),TSED(I)-TFREZ,TLAK(I,NLEV)-TFREZ
end%605   CONTINUE
forStr='%5d %d %3d %2d %2d %f %f %f %f %f %f %f %f';
%numOfLakeTemps=10;
% for twice=1:2
%     for i=1:numOfLakeTemps
%         forStr=[forStr,' %f'];
%     end
% end
% forStr=[forStr,'\n'];
% fprintf(fID_82,forStr,NSTEP(1), IYEAR(1), IDAY(1), IHOUR(1), IMIN(1),...
%     TLAKOld1(1,1),TLAKOld2(1,1),TLAKOld3(1,1),TLAKOld4(1,1),TLAKOld5(1,1),...
%     TLAKOld6(1,1),TLAKOld7(1,1),TLAK(1,1));

% 
% 6010  FORMAT(I4,1X,3(I3,1X),8F8.2,2F7.3,5E10.3)
% 6020  FORMAT(I4,1X,3(I3,1X),200F7.2)
% 6030  FORMAT(I4,1X,3(I3,1X),4E10.2,F5.2,F6.2,I3,3F7.2,2E10.2)
% 6040  FORMAT(I4,1X,3(I3,1X),3E10.3,3F7.2,F8.3,F8.3,E10.3)
% 6050  FORMAT(I4,1X,3(I3,1X),5E10.3)
% 6060  FORMAT(I4,1X,3(I3,1X),8E10.3)
% 6070  FORMAT(I4,1X,3(I3,1X),6F7.1) 
% C                                                                       
% C     * SCREEN LEVEL DIAGNOSTICS.                                                    
% C                                                                       
for I=IL1:IL2% DO I=IL1,IL2
    CDM(I)=FLS(I)*CDMS(I)+(1.0-FLS(I))*CDML(I);
    CDH(I)=FLS(I)*CDHS(I)+(1.0-FLS(I))*CDHL(I);
    DRAG(I)=FLS(I)*DRAGS(I)+(1.0-FLS(I))*DRAGL(I);
    ZOM(I)=ZREFM(I)/exp(VKC/sqrt(DRAG(I)));
    ZOH(I)=ZOM(I)/3.0;
    TSURF(I)=FLS(I)*TZEROS(I)+(1.0-FLS(I))*T0(I);
    QSURF(I)=FLS(I)*QZEROS(I)+(1.0-FLS(I))*QSURFL(I);
    ALVS(I)=FLS(I)*ALVSSN(I)+(1.0-FLS(I))*ALVSL(I);
    ALIR(I)=FLS(I)*ALIRSN(I)+(1.0-FLS(I))*ALIRL(I);
end%ENDDO
%C
if(ISLFD==0)
    for I=IL1:IL2%DO I=IL1,IL2                                            
        FACTM=ZDIAGM(I)+ZOM(I);  
        FACTH=ZDIAGH(I)+ZOM(I);   
        RATIOM=sqrt(CDM(I))*log(FACTM/ZOM(I))/VKC;              
        RATIOM=MIN(RATIOM,1.);                                  
        RATIOH=sqrt(CDM(I))*log(FACTH/ZOH(I))/VKC;             
        RATIOH=MIN(RATIOH,1.);                                   
        if(TSURF(I)>TA(I))
            RATIOH=RATIOH*CDH(I)/CDM(I);                         
            RATIOH=MIN(RATIOH,(FACTH/ZREFH(I))^(1./3.));         
        end                                                   
        ST(I)=TSURF(I)-(MIN(RATIOH,1.))*(TSURF(I)-TA(I));       
        SQ(I)=QSURF(I)-(MIN(RATIOH,1.))*(QSURF(I)-QA(I));      
        SU(I)=RATIOM*UWIND(I);                                  
        SV(I)=RATIOM*VWIND(I);                                  
    end%ENDDO
    % C                                                                       
    %After fiddling around with this function, Murray MacKay told me it is
    %not needed.  It only changes SH, which is passed back out of CLASSL
    %and assigned as a diagnostic variabile. 
    %[SH,ST,SQ,PRES,FTOT,ILG,IL1,IL2]=SCREENRH(SH,ST,SQ,PRES,FTOT,ILG,IL1,IL2);           
    %[SH]=SCREENRH(SH,ST,SQ,PRES,FTOT,ILG,IL1,IL2);           
    %[SH]=SCREENRH(IL1,IL2);
    % C                                                                       
elseif(ISLFD==1)
    [SU,SV,ST,SQ,...    
                CDM,CDH,UWIND,VWIND,TA,QA,...                 
                TSURF,QSURF,ZOM,ZOH,FTOT,ZREFM,...                 
                ZDIAGM,ZDIAGH,ILG,IL1,IL2,JL]=...
        SLDIAG(SU,SV,ST,SQ,...    
                CDM,CDH,UWIND,VWIND,TA,QA,...                 
                TSURF,QSURF,ZOM,ZOH,FTOT,ZREFM,...                 
                ZDIAGM,ZDIAGH,ILG,IL1,IL2,JL);                     
    % C
    %See note above where this is called about this function not being needed.
    %[SH,ST,SQ,PRES,FTOT,ILG,IL1,IL2]=SCREENRH(SH,ST,SQ,PRES,FTOT,ILG,IL1,IL2);           
    %[SH]=SCREENRH(SH,ST,SQ,PRES,FTOT,ILG,IL1,IL2);           
    %[SH]=SCREENRH(IL1,IL2);
    % C                                                                       
elseif(ISLFD==2)                              
    [SU,SV,ST,SQ]=...,ILG,UWIND,VWIND,TSURF,QSURF,ZOM,ZOH,ILMO,ZREFM,HBL,UE,FTEMP,FVAP,ZDIAGM,ZDIAGH,RADJ,FTOT,IL1,IL2,JL...
        DIASURFZ(SU,SV,ST,SQ,UWIND,VWIND,TSURF,QSURF,ZOM,ZOH,ILMO,ZREFM,...ILG,
          HBL,UE,FTEMP,FVAP,ZDIAGM,ZDIAGH,RADJ,FTOT,IL1,IL2);% JL              ,
end

%C                                                                       
for I=IL1:IL2%DO I=IL1,IL2                                              
    FSGS(I) =FLS(I)*(QSWNS(I)-QTRANSL(I));           
    FLGS(I) =FLS(I)*(QLWIN(I)-QLWOS(I));            
    HFSS(I) =FLS(I)*QSENSS(I);                     
    HEVS(I) =FLS(I)*QEVAPS(I);                     
    HTCS(I) =HTCS(I)-(FLS(I)*GZEROSL(I));
    HTCL(I) =HTCL(I)+(FLS(I)*GZEROSL(I));
    QLWAVG(I)=FLS(I)*QLWOS(I)+(1.0-FLS(I))*QLWIN(I)-FLGL(I);
    QSENS(I)=HFSL(I)+HFSS(I);
    TFLUX(I)=-QSENS(I)/(RHOAIR(I)*SPHAIR);
    QEVAP(I)=HEVL(I)+HEVS(I);
    EVAP(I)=QFL(I)+QFN(I);
    QFLUX(I)=-EVAP(I)/RHOAIR(I);
    EVPPOT(I)=EVAP(I);
    EVAPB(I)=1.0;
    %C         GT(I)=(QLWAVG(I)/SBC)**0.25
    GT(I)=FLS(I)*TZEROS(I)+(1.0-FLS(I))*T0(I);
end%ENDDO
%C                                                                       
%RETURN   
%END 
% fclose(fID_71);
% fclose(fID_72);
% fclose(fID_73);
% fclose(fID_74);
% fclose(fID_75);
% fclose(fID_76);
% fclose(fID_77);
%fclose(fID_78);
%fclose(fID_79);


end%End of CLASSL function