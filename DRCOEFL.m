function [CDM,CDH,CDMN] =...,VA,T0,TA,QA,PRES,ZREFM,ZREFH,FICE,FLS,ILG,IL1,IL2
    DRCOEFL (VA,T0,TA,QA,PRES,ZREFM,ZREFH,FICE,FLS,IL1,IL2,...CDM,CDH,CDMN,ILG,
    GRAV,VKC,TFREZ) %Added these because gloabl 

% C=======================================================================
% C     * DEC 15/16 - M.LAZARE/   - NSTEP REMOVED.
% C     *             D.VERSEGHY. - BUGFIX IN ACCOUNTING FOR EFFECT
% C     *                           OF FICE ON DRAG COEFFICIENTS.
% C     * JAN 25/16 - M.MACKAY.   FRACTIONAL ICE COVER INCLUDED
% C     *                         BUG FIX FOR NEAR SFC HUMIDITY
% C     *                         THERMODYNAMIC REF HEIGHT ADDED
% C     * MAY 22/15 - D.VERSEGHY. WEIGHT DRAG COEFFICIENTS FOR PRESENCE
% C     *                         OF ICE.
% C     * NOV 27/07 - M.MACKAY.  	TURBULENT TRANSFER COEFFICIENTS
% C     *
% C=======================================================================
% C
%       IMPLICIT NONE
% C
% C ----* INPUT FIELDS *------------------------------------------------
% C
%       INTEGER ILG,IL1,IL2
%       REAL,DIMENSION(ILG) :: VA,T0,TA,QA,PRES,ZREFM,ZREFH,FICE,FLS
% C
% C ----* OUTPUT FIELDS *------------------------------------------------
% C
%       REAL,DIMENSION(ILG) ::
%      1  CDH, CDM, CDMN
% C
% C ----* COMMON BLOCK PARAMETERS *--------------------------------------
% C
%       REAL RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN,DELT,TFREZ
% C
% C ----* CLASS COMMON BLOCKS *------------------------------------------
% C
%       COMMON /CLASS1/ DELT,TFREZ                                       
%       COMMON /CLASS2/ RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN
% C
% C ----* LOCAL VARIABLES *---------------------------------------------
% C
%       INTEGER I,J,K,ITER,ITMAX
%       REAL G,TOL,CHARN,RESID,VSQ,Cold,EASAT,QASAT,CA,CB,SHF,LHF,
%      1     TVIRT,MOL,Z_L,CDHN,X,PSIM,PSIH,PI,MOLold,TOL2,
%      2     CDHRAT,CDMNW,CDHNW,CDMNI,CDHNI,DENOM,XICE
% C
% C ----* LOCAL PARAMETERS *--------------------------------------------
% C
G=GRAV;
TOL=single(1.0E-8);
TOL2=single(1.0E-2);
ITMAX=single(100);
CHARN=single(0.0175);
PI=single(3.14159);

CDM=single(0);
CDH=single(0);
CDMN=single(0);
% C
% C======================================================================
for I=IL1:IL2%DO 100 I=IL1,IL2
    % C-----------------------------------------------------------------------
    % C NEUTRAL DRAG COEFFICIENTS  
    % C  Iterative routine to solve for CDN based on Charnock relation for 
    % C  roughness
    % C
    CDHNW=single(1.35E-3);
    CDMNW=single(1.0E-3);%        !initial trial value
    VSQ=VA(I)*VA(I);
    RESID=999.0;
    while (abs(RESID) > TOL) %DO 
        %if (abs(RESID) <= TOL ) 
        %    break;%EXIT
        %end
        Cold=CDMNW;
        CDMNW=VKC*VKC/(log(ZREFM(I)*G/(CHARN*Cold*VSQ))*log(ZREFM(I)*G/(CHARN*Cold*VSQ)));
        RESID=CDMNW-Cold;
    end %END DO   
    % C
    CDMNI=(VKC/(log(ZREFM(I)/0.002)))^2;
    CDHNI=(VKC/(log(ZREFH(I)/0.00067)))^2;
    if(FICE(I)>(FLS(I)+0.001))
        XICE=MAX(FICE(I)-FLS(I),0.0);
        CDMN(I)=(XICE*CDMNI+(1.0-FICE(I))*CDMNW)/(1.0-FLS(I));
        CDHN=(XICE*CDHNI+(1.0-FICE(I))*CDHNW)/(1.0-FLS(I));
    else
        CDMN(I)=CDMNW;
        CDHN=CDHNW;
    end
    % C-----------------------------------------------------------------------
    % C INITIAL TRIAL VALUES FOR TRANSFER COEFFICIENTS: SET TO NEUTRAL VALUES
    % C
    CDH(I)=CDHN;
    CDM(I)=CDMN(I);
    % C-----------------------------------------------------------------------
    % C ITERATIVELY COMPUTE DRAG COEFFICIENTS UNTIL M.O. LENGTH CONVERGES
    % C
    RESID=999.0;
    MOL=9999.0;
    ITER=0;
    while (ITER < ITMAX)%DO  %MATLAB doesn't have DO, so while and if work here.
%         IF (ABS(RESID) .LE. TOL2 .OR. ITER .GE. ITMAX) EXIT
         if (abs(RESID) <= TOL2) 
             break;%EXIT
         end
        % C-----------------------------------------------------------------------
        % C HEAT FLUXES
        % C----------------------------------------------------
        % C     * CALCULATION OF EASAT CONSISTENT WITH CLASSI
        % C     * BUT CONSTANTS DifFER FROM ROGERS&YAU
        % C     * Rogers and Yau values
        % C         CA=17.67
        % C         CB=29.65
        % C----------------------------------------------------
        if(T0(I)>=TFREZ)
            CA=17.269;                                                       
            CB=35.86;                                                        
        else                                                                
            CA=21.874;                                                       
            CB=7.66;                                                         
        end                                                               
        SHF=CDH(I)*VA(I)*(T0(I)-TA(I));
        EASAT=611.0*exp(CA*(T0(I)-TFREZ)/(T0(I)-CB));
        QASAT=0.622*EASAT/(PRES(I)-0.378*EASAT);
        LHF=CDH(I)*VA(I)*(QASAT-QA(I));

        % C-----------------------------------------------------------------------
        % C VIRTUAL TEMPERATURE AND M.-O. LENGTH
        % C
        TVIRT=TA(I)*(1.0+0.61*QA(I));
        MOLold=MOL;
        MOL = -VA(I)*VA(I)*VA(I)*CDM(I)*sqrt(CDM(I))*TVIRT/( VKC*G*(SHF + 0.61*LHF*TA(I)) );
        Z_L = ZREFM(I)/MOL;
        RESID=MOL-MOLold;

        % C-----------------------------------------------------------------------
        % C STABILITY CORRECTIONS
        % C
        % C
        % C- UNSTABLE CASE
        % C---------------
        if (Z_L < 0.0) 
            X = (1.0 - (16.0*Z_L))^0.25;
            PSIM = 2.0*log((1.0+X)/2.0) + log((1.0+X*X)/2.0) - 2.0*atan(X) + PI/2.0;
            PSIH = 2.0*log((1.0+X*X)/2.0);
            % C
            % C- STABLE CASE
            % C-------------
        elseif (Z_L >= 0 && Z_L < 0.5)
            PSIM = -5.0*Z_L;
            PSIH=PSIM;
        elseif (Z_L >= 0.5 && Z_L < 10.0)
            PSIM = (0.5/(Z_L*Z_L)) - (4.25/Z_L) - 7.0*log(Z_L) - 0.852;
            PSIH=PSIM;
        else 
            PSIM = log(Z_L) - 0.76*Z_L - 12.093;
            PSIH=PSIM;
        end%END if

        % C-----------------------------------------------------------------------
        % C RECOMPUTE DRAG COEFFICIENTS WITH STABILTY CORRECTIONS
        % C
        DENOM = (1.0 + (CDMN(I)/(VKC*VKC))*(PSIM*PSIM - (2.0*VKC*PSIM/sqrt(CDMN(I)))) );
        if(DENOM<1.0E-6) 
            DENOM=1.0E-6;
        end
        % C        if(abs(DENOM)<1.0E-6) DENOM=SIGN(1.0E-6,DENOM)
        CDM(I) = CDMN(I)/DENOM;
        DENOM = (1.0 + (CDHN/(VKC*VKC))*(PSIM*PSIH -(VKC*PSIH/sqrt(CDMN(I)))-(VKC*PSIM*sqrt(CDMN(I))/CDHN)));
        if(DENOM<1.0E-6) 
            DENOM=1.0E-6;
        end
        % C        if(abs(DENOM)<1.0E-6) DENOM=SIGN(1.0E-6,DENOM)
        CDH(I) = CDHN/DENOM;

        ITER=ITER+1;
    end%END DO

    if(ITER>=ITMAX) 
        CDM(I)=CDMN(I);
        CDH(I)=CDHN;
    end

    % C     if (ITER >= ITMAX) print*, "^ max iters reached"

    if (CDH(I) < 0.0) 
        CDH(I)=CDHN;
    end

    if (CDM(I) < 0.0)
        CDM(I)=CDMN(I);
    end

    CDHRAT=CDH(I)/CDHN;
    if (CDHRAT >= 8.0) 
        CDH(I)=CDHN;
    end

end%100   CONTINUE

%RETURN                                                                      
%END     