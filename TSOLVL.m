function [T0,LKICEH,KLAK,GLAK,Q0SAT,KSTAR,LSTAR,QSENS,...
           QEVAP,EVAP,ALVS,ALIR,G0] = ...
    TSOLVL(TLAK,T0,LKICEH,EVAP,ALVS,ALIR,...
           QSWIN,QLWIN,TA,QA,VA,PRES,RHOAIR,CDH,...
           GZEROL,QTRANSL,HTCL,FLS,FICE,ALBW,ALBI,...
           CQ1A,CQ1B,CQ2A,CQ2B,CQ3A,CQ3B,...
           CQ1BI,CQ2BI,CQ3BI,IL1,IL2,...
           DELZLK,DELSKIN,RHOICE,RHOW,TFREZ,...
           CLHVAP,SPHAIR,EMSW,SBC,TCW,DELT,HCPW,...
           CLHMLT,TCICE,HCPICE)
                   %Already reduced the input/output scope
% C======================================================================
% C     * DEC 19/16 - M.LAZARE/   - NSTEP REMOVED (WASN'T USED).
% C     *             D.VERSEGHY. - BUGFIX IN ACCOUNTING FOR EFFECT
% C     *                           OF FICE ON ALBEDO AND LATENT HEAT.
% C     *                         - {ALBW,ALBI} NOW DEFINED IN CLASSL
% C     *                           AND PASSED IN, FOR CONSISTENCY.
% C     * APR 11/16 - M.MACKAY.   THERMAL EXPANSION OF ICE BUG CORRECTED
% C     *                         ICE DRAFT,FREEBOARD NOW INCLUDED
% C     *                         SW ATTENUATION THROUGH LEADS INCLUDED
% C     * SEP 02/15 - D.VERSEGHY. ADD EFFECTS OF SNOW COVERAGE ON SURFACE
% C     *                         FLUXES; ADDITIONAL DIAGNOSTIC VARIABLES; 
% C     *                         COSMETIC CHANGES TO CODE AND NAMING.
% C     * SEP  2/11 - M.MACKAY.  	ICE COMPUTED IN SKIN LAYER
% C     * SEP 28/07 - M.MACKAY.  	SFC ENERGY BALANCE NOW COMPUTED OVER
% C     *                         SKIN OF FINITE WIDTH
% C     * MAR 15/07 - M.MACKAY.   COMPUTES LAKE SURFACE ENERGY BALANCE
% C     *  
% C
%       IMPLICIT NONE
% C
% C ----* GLOBAL LAKE VARIABLES *---------------------------------------
% C
%       INTEGER NLAKMAX
%       REAL,DIMENSION(ILG) :: T0,LKICEH,KLAK,GLAK,Q0SAT
%       REAL,DIMENSION(ILG,NLAKMAX) :: TLAK
% C
% C ----* OUTPUT FIELDS *------------------------------------------------
% C
%       REAL,DIMENSION(ILG) :: KSTAR,LSTAR,QSENS,QEVAP,EVAP,ALVS,ALIR
% C
% C ----* INPUT FIELDS *------------------------------------------------
% C
%       REAL,DIMENSION(ILG) :: QSWIN,QLWIN,CSZ,TA,QA,VA,PRES,RHOAIR,CDH,
%      1                       GZEROL,QTRANSL,HTCL,FLS,FICE,ALBW,ALBI,
%      2                       CQ1A,CQ1B,CQ2A,CQ2B,CQ3A,CQ3B,
%      3                       CQ1BI,CQ2BI,CQ3BI
%       INTEGER ILG,IL1,IL2
% C
% C ----* INTERNAL ARRAYS *----------------------------------------------
% C
%       REAL,DIMENSION(ILG) :: G0,XICE
% C
% C ----* COMMON BLOCKS *------------------------------------------------
% C
%       REAL DELT,TFREZ,RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN,TCW,TCICE,
%      1     TCSAND,TCCLAY,TCOM,TCDRYS,RHOSOL,RHOOM,HCPW,HCPICE,HCPSOL,
%      2     HCPOM,HCPSND,HCPCLY,SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,
%      3     TCGLAC,CLHMLT,CLHVAP,CGRAV,CKARM,CPD,AS,ASX,CI,BS,
%      4     BETA,FACTN,HMIN
%       REAL TKECN,TKECF,TKECE,TKECS,HDPTHMIN,
%      1     TKEMIN,DELMAX,DELMIN,EMSW,DELZLK,DELSKIN,DHMAX,
%      2     TKECL,DUMAX
%       COMMON /CLASS1/ DELT,TFREZ                                       
%       COMMON /CLASS2/ RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN
%       COMMON /CLASS3/ TCW,TCICE,TCSAND,TCCLAY,TCOM,TCDRYS,
%      1                RHOSOL,RHOOM
%       COMMON /CLASS4/ HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,
%      1                SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,
%      2                TCGLAC,CLHMLT,CLHVAP
%       COMMON /LAKECON/ TKECN,TKECF,TKECE,TKECS,HDPTHMIN,          
%      1                 TKEMIN,DELMAX,DELMIN,EMSW,DELZLK,DELSKIN,DHMAX,
%      2                 TKECL,DUMAX
% C
% C ----* LOCAL VARIABLES *--------------------------------------------
% C
%       INTEGER I,J,ITMAX
%       REAL DZ,CA,CB,E0SAT,DS
%       REAL ALBTOT,T0old,NEWICE,ESKIN,ECOOL
%       REAL EAVAIL,ATTEN1,ATTEN2,ATTEN3,EHEAT,TC,CPHCH
%       REAL RHOIW,ICETOP,ICEBOT
% C
% C ----* LOCAL PARAMETER DEFINITIONS *-------------------------------
% C

DZ=DELZLK;
DS=DELSKIN;
RHOIW=RHOICE/RHOW;
XICE=single(nan(1,length(IL1:IL2)));
% C----------------------------------------------------------------------
for I=IL1:IL2%DO 100 I=IL1,IL2

    % C======================================================================
    % C COMPUTE SURFACE ENERGY BALANCE FOR CURRENT TIMESTEP
    % C======================================================================
    % C COMPUTE SW FLUXES --------------------------------------------------
    % C   KSTAR is net SW at surface of skin
    % C   KLAK is penetrating SW at base of skin
    % C
    XICE(I)=MAX(FICE(I)-FLS(I),0.0);
    if(FICE(I)>(FLS(I)+0.001))
        ALBTOT=(XICE(I)*ALBI(I)+(1.-FICE(I))*ALBW(I))/(1.0-FLS(I));
    else
        ALBTOT=ALBW(I);
    end

    KSTAR(I)=(1.-FLS(I))*(1.-ALBTOT)*QSWIN(I)+FLS(I)*QTRANSL(I);
    if ( KSTAR(I) < 0.0)  
        KSTAR(I)=0.0;
    end
    ALVS(I)=ALBTOT;
    ALIR(I)=ALBTOT;

    % C--- Attenuation through ice (now includes attenuation through leads)
    % C---
    ICEBOT=RHOIW*LKICEH(I);%                !ICE DRAFT
    ICETOP=LKICEH(I)-ICEBOT;%               !ICE FREEBOARD
    if ( LKICEH(I) <= 0.0) %                !NO ICE 
        ATTEN1=CQ1B(I)*DS;
        ATTEN2=CQ2B(I)*DS;
        ATTEN3=CQ3B(I)*DS;
    else 
        if (ICEBOT > DS)       %            !Z inside ice
            ATTEN1=FICE(I)*(CQ1BI(I)*(DS+ICETOP)) + (1.-FICE(I))*CQ1B(I)*DS;
            ATTEN2=FICE(I)*(CQ2BI(I)*(DS+ICETOP)) + (1.-FICE(I))*CQ2B(I)*DS;
            ATTEN3=FICE(I)*(CQ3BI(I)*(DS+ICETOP)) + (1.-FICE(I))*CQ3B(I)*DS;
        else
            ATTEN1=FICE(I)*(CQ1BI(I)*LKICEH(I) + CQ1B(I)*(DS-ICEBOT)) + (1.-FICE(I))*CQ1B(I)*DS;
            ATTEN2=FICE(I)*(CQ2BI(I)*LKICEH(I) + CQ2B(I)*(DS-ICEBOT)) + (1.-FICE(I))*CQ2B(I)*DS;
            ATTEN3=FICE(I)*(CQ3BI(I)*LKICEH(I) + CQ3B(I)*(DS-ICEBOT)) + (1.-FICE(I))*CQ3B(I)*DS;
        end
    end
    KLAK(I)=KSTAR(I)*(CQ1A(I)*exp(-ATTEN1) + CQ2A(I)*exp(-ATTEN2) + CQ3A(I)*exp(-ATTEN3) );
    % C
    % C COMPUTE TURBULENT FLUXES -------------------------------------------
    % C     * CALCULATION OF E0SAT CONSISTENT WITH CLASSI
    % C     * BUT CONSTANTS DifFER FROM ROGERS&YAU
    % C     * Rogers and Yau values
    % C         CA=17.67
    % C         CB=29.65
    % C
    if(T0(I)>=TFREZ)                              
        CA=17.269;                                       
        CB=35.86;                                       
    else                                              
        CA=21.874;                                    
        CB=7.66;                                   
    end                                          

    if(FICE(I)>(FLS(I)+0.001))
        CPHCH=(XICE(I)*(CLHMLT+CLHVAP)+(1.-FICE(I))*CLHVAP)/(1.0-FLS(I));
    else
        CPHCH=CLHVAP;
    end
    QSENS(I)=(1.-FLS(I))*RHOAIR(I)*SPHAIR*CDH(I)*VA(I)*(T0(I)-TA(I));
    E0SAT=611.0*exp(CA*(T0(I)-TFREZ)/(T0(I)-CB));
    Q0SAT(I)=0.622*E0SAT/(PRES(I)-0.378*E0SAT);
    EVAP(I)=(1.-FLS(I))*RHOAIR(I)*CDH(I)*VA(I)*(Q0SAT(I)-QA(I));
    % Cmdm    QEVAP(I)=CPHCH*EVAP(I)
    QEVAP(I)=(1.-FLS(I))*RHOAIR(I)*CPHCH*CDH(I)*VA(I)*(Q0SAT(I)-QA(I));
    % C
    % C COMPUTE NET LW AND NET SFC ENERGY -------------------------------
    % C GLAK IS THERMAL FLUX AT BASE OF SKIN.  THERMAL CONDUCTIVITY BASED
    % C ON WEIGHTED AVERAGE FOR WATER AND ICE if ICE PRESENT IN LAYER
    % C
    LSTAR(I)=(1.-FLS(I))*EMSW*(QLWIN(I)-SBC*T0(I)*T0(I)*T0(I)*T0(I));
    G0(I)=KSTAR(I)-KLAK(I)+LSTAR(I)-QSENS(I)-QEVAP(I)+HTCL(I)+FLS(I)*GZEROL(I);
    if (ICEBOT >= DS) 
        GLAK(I)=(-2.0*TCICE/(DZ+DS))*(TLAK(I,1)-T0(I));
    elseif (ICEBOT < DS && LKICEH(I) > 0.0)
        TC=(ICEBOT*TCICE + (DS-ICEBOT)*TCW)/DS;
        % Cmdm      TC=20.0          !mdm test
        GLAK(I)=(-2.0*TC/(DZ+DS))*(TLAK(I,1)-T0(I));
    else
        GLAK(I)=(-2.0*TCW/(DZ+DS))*(TLAK(I,1)-T0(I));
    end
    % C-----NET ENERGY FLUX INTO SKIN (W/M2)
    ESKIN= G0(I) - GLAK(I);
    % C
    % C STEP FORWARD SKIN TEMP T0
    % C
    T0old=T0(I);
    if (LKICEH(I) <= 0.0)
        T0(I) = T0(I) + (DELT/(DS*HCPW))*ESKIN;
    elseif (LKICEH(I) > 0.0 && ICEBOT <= DS) 
        T0(I) = T0(I) + (DELT/((LKICEH(I)*HCPICE) + (DS-ICEBOT)*HCPW))*ESKIN;
    else
        T0(I) = T0(I) + (DELT/((DS+ICETOP)*HCPICE))*ESKIN;
    end
    % C
    % C ICE GROWTH OR DECAY
    % C
    if (ESKIN < 0.0 && ICEBOT < DS)
        % C-----NET ENERGY FLUX USED TO LOWER T0 TO TFREZ 
        if (T0old > TFREZ)
            ECOOL=(DS-ICEBOT)*HCPW*(T0old-TFREZ)/DELT;
        else
            ECOOL=0.0;
        end
    % C-----REMAINING ENERGY FLUX (if ANY) USED TO FREEZE ICE
        EAVAIL=ESKIN+ECOOL;
        if (EAVAIL < 0.0)
            NEWICE=-(DELT/(RHOICE*CLHMLT))*EAVAIL;
            LKICEH(I)=LKICEH(I)+NEWICE;
            ICEBOT=RHOIW*LKICEH(I);
            T0(I)=TFREZ;
            % C-----LIMIT ICE GROWTH TO THE CURRENT LAYER
            if (ICEBOT > DS)
                EHEAT=(RHOICE*CLHMLT*(ICEBOT-DS))/DELT;
                T0(I)=TFREZ - (EHEAT*DELT)/(DS*HCPICE);
                LKICEH(I)=DS/RHOIW;
            end
        end
    end

    if (ESKIN > 0.0 && LKICEH(I) > 0.0)
        % C-----NET ENERGY FLUX USED FIRST TO RAISE T TO ZERO
        if (ICEBOT <= DS)
            EHEAT=LKICEH(I)*HCPICE*(TFREZ-T0old)/DELT;
        else
            EHEAT=(DS+ICETOP)*HCPICE*(TFREZ-T0old)/DELT;
        end

        % C-----NET ENERGY FLUX USED TO MELT ICE
        if (ESKIN > EHEAT)
            NEWICE=-(DELT/(RHOICE*CLHMLT))*(ESKIN-EHEAT);
            LKICEH(I)=LKICEH(I)+NEWICE;
            T0(I)=TFREZ;
            % C-----LIMIT ICE MELT TO THE CURRENT LAYER
            if (LKICEH(I) < 0.0)
                EHEAT=-(RHOICE*CLHMLT*LKICEH(I))/DELT;
                T0(I)=TFREZ + (EHEAT*DELT)/(DS*HCPW);
                LKICEH(I)=0.0;
            end
        end
    end
    % C -----------------
end%100     CONTINUE

%RETURN
%END
end
