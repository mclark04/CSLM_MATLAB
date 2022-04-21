function [...ISNOW,FI,
        QSWNET,QLWOUT,QTRANS,QSENS,QEVAP,EVAP,...
        TZERO,QZERO,GZERO,QMELT,CDH,CDM,RIB,CFLUX,...
        FTEMP,FVAP,ILMO,UE,H,ZRSLFH,ZRSLFM,ZOH,ZOM,FCOR,ILG,...%used in FLXSURFZ
        ...QLWIN,TPOTA,QA,VA,PADRY,RHOAIR,
        ...ALVISG,ALNIRG,CRIB,CPHCH,CEVAP,TVIRTA,
        ...ZOSCLH,ZOSCLM,
        ...GCONST,GCOEFF,TSTART,PCPR,TRSNOWG,FSSB,ALSNO,
        ...THLIQ,THLMIN,DELZW,RHOSNO,ZSNOW,ZPOND,
        ITERCT,...IWATER,IEVAP,ISAND,
        ...ISLFD,ITG,IG,IL1,IL2,JL,NBS,ISNOALB,
        TSTEP,TVIRTS,EVBETA,Q0SAT,RESID,...
        DCFLXM,CFLUXM,WZERO,TRTOP,...A,B,
        LZZ0,LZZ0T,FM,FH,ITER,NITER,JEVAP,KF] = ...
    TSOLVE(ISNOW,FI,...
        QLWOUT,QSENS,QEVAP,EVAP,...,QSWNET,QTRANS
        TZERO,QZERO,GZERO,QMELT,CDH,CDM,RIB,CFLUX,...
        FTEMP,FVAP,ILMO,UE,H,ZRSLFH,ZRSLFM,ZOH,ZOM,FCOR,...%passed through FLZSURFZ
        QLWIN,TPOTA,QA,VA,PADRY,RHOAIR,...
        ALVISG,ALNIRG,CRIB,CPHCH,CEVAP,TVIRTA,...
        ZOSCLH,ZOSCLM,...
        GCONST,GCOEFF,TSTART,TRSNOWG,FSSB,ALSNO,...PCPR,
        THLIQ,THLMIN,DELZW,RHOSNO,ZSNOW,ZPOND,...
        IWATER,IEVAP,ITERCT,ISAND,...
        ISLFD,ITG,ILG,IL1,IL2,NBS,ISNOALB,...IG,JL,
        TSTEP,TVIRTS,EVBETA,Q0SAT,RESID,...
        DCFLXM,CFLUXM,WZERO,TRTOP,...A,B,
        LZZ0,LZZ0T,FM,FH,ITER,NITER,KF,...JEVAP,
        DELT,TFREZ,RHOW,SBC,SPHAIR,GRAV,VKC) 
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
% !-------------------------------------- LICENCE END --------------------------------------
% SUBROUTINE TSOLVE(ISNOW,FI,                                       
% 1                  QSWNET,QLWOUT,QTRANS,QSENS,QEVAP,EVAP,          
% 2                  TZERO,QZERO,GZERO,QMELT,CDH,CDM,RIB,CFLUX,      
% 3                  FTEMP,FVAP,ILMO,UE,H,                           
% 4                  QLWIN,TPOTA,QA,VA,PADRY,RHOAIR,                 
% 5                  ALVISG,ALNIRG,CRIB,CPHCH,CEVAP,TVIRTA,          
% 6                  ZOSCLH,ZOSCLM,ZRSLFH,ZRSLFM,ZOH,ZOM,FCOR,       
% 7                  GCONST,GCOEFF,TSTART,PCPR,TRSNOWG,FSSB,ALSNO,   
% 8                  THLIQ,THLMIN,DELZW,RHOSNO,ZSNOW,ZPOND,
% +                  IWATER,IEVAP,ITERCT,ISAND, 
% 9                  ISLFD,ITG,ILG,IG,IL1,IL2,JL,NBS,ISNOALB,        
% A                  TSTEP,TVIRTS,EVBETA,Q0SAT,RESID,                
% B                  DCFLXM,CFLUXM,WZERO,TRTOP,A,B,                  
% C                  LZZ0,LZZ0T,FM,FH,ITER,NITER,JEVAP,KF)           
% C                                                                       
% C     * OCT 26/16 - D.VERSEGHY. ADD ZPOND TO CALCULATION OF EVPMAX.
% C     * JUL 22/15 - D.VERSEGHY. LIMIT CALCULATED EVAPORATION RATE
% C     *                         ACCORDING TO WATER AVAILABILITY.
% C     * JAN 09/15 - D.VERSEGHY. FIX TO SUPPRESS EVAPORATION FROM ROCK.
% C     * JUN 27/14 - D.VERSEGHY. CHANGE ITERATION LIMIT BACK TO 50 FOR
% C     *                         BISECTION SCHEME.
% C     * NOV 16/13 - J.COLE/     FINAL VERSION FOR GCM17:                
% C     *             M.LAZARE.   - FIX COMPUTATION OF QSWNI OVER SNOW FREE 
% C     *                           BARE SOIL for ISNOW=0 and ISNOALB=1 (NEED 
% C     *                           TO SUM OVER THE 3 NEAR-IR BANDS).     
% C     * JUN 22/13 - J.COLE/     - ADD "ISNOALB" OPTION (4-BAND SOLAR).  
% C     *             M.LAZARE.   - MODIFY ABORT CONDITION FOR TOO COLD   
% C     *                           TEMPS FROM 173 TO 123, SO WON'T       
% C     *                           BLOW UP OVER ANTARCTICA.              
% C     * OCT 14/11 - D.VERSEGHY. FOR POST-ITERATION CLEANUP WITH N-R SCHEME,
% C     *                         REMOVE CONDITION INVOLVING LAST ITERATION 
% C     *                         TEMPERATURE.                             
% C     * DEC 07/09 - D.VERSEGHY. RESTORE EVAPORATION WHEN PRECIPITATION   
% C     *                         IS OCCURRING.                            
% C     * MAR 13/09 - D.VERSEGHY. REPLACE SURFCON COMMON BLOCK WITH CLASSD2;
% C     *                         REVISED CALL TO FLXSURFZ.               
% C     * JAN 06/09 - D.VERSEGHY/M.LAZARE. SPLIT if CONDITIONS FRAMING    
% C     *                         300 LOOP.                               
% C     * FEB 25/08 - D.VERSEGHY. STREAMLINE SOME CALCULATIONS; REMOVE    
% C     *                         "ILW" SWITCH; SUPPRESS WATER VAPOUR FLUX  
% C     *                         if PRECIPITATION IS OCCURRING.            
% C     * MAY 17/06 - D.VERSEGHY. SUPPRESS EVAPORATION WHEN PONDED WATER    
% C     *                         IS FREEZING; ADD IL1 AND IL2 TO CALL TO   
% C     *                         FLXSURFZ; REMOVE JL FROM CALL TO DRCOEF.  
% C     * APR 13/05 - R.BROWN. ADD WINDLESS TRANFER COEFFICIENT TO QSENS    
% C     *                         CALCULATION FOR SNOW PACKS.               
% C     * DEC 17/04 - Y.DELAGE/D.VERSEGHY. ADD SWITCH TO USE EITHER SECANT/ 
% C     *                         BISECTION OR NEWTON-RAPHSON ITERATION     
% C     *                         SCHEME (WITH NUMBER OF ITERATIONS LIMITED 
% C     *                         TO FIVE AND CORRECTION FOR REMAINING      
% C     *                         RESIDUAL).                              
% C     * NOV 04/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.            
% C     * AUG 06/04 - Y.DELAGE/D.VERSEGHY. PROTECT SENSITIVE CALCULATIONS 
% C     *                         FROM ROUNDOFF ERRORS.                   
% C     * NOV 07/02 - Y.DELAGE/D.VERSEGHY. NEW CALL TO FLXSURFZ.          
% C     * JUL 26/02 - D.VERSEGHY. SHORTENED CLASS4 COMMON BLOCK.          
% C     * MAR 28/02 - D.VERSEGHY. STREAMLINED SUBROUTINE CALL.            
% C     *                         BYPASS EVAPORATION EFFICIENCY PARAMETER 
% C     *                         IN CASES OF CONDENSATION.               
% C     * JAN 18/02 - P.BARTLETT/D.VERSEGHY. NEW "BETA" FORMULATION FOR   
% C     *                         BARE SOIL EVAPORATION BASED ON LEE AND  
% C     *                         PIELKE.                                 
% C     * APR 11/01 - M.LAZARE.   SHORTENED "CLASS2" COMMON BLOCK.        
% C     * OCT 06/00 - D.VERSEGHY. CONDITIONAL "if" IN ITERATION SEQUENCE  
% C     *                         TO AVOID DIVIDE BY ZERO.                
% C     * DEC 07/99 - A.WU/D.VERSEGHY. NEW SOIL EVAPORATION FORMULATION.  
% C     * JUL 24/97 - D.VERSEGHY. CLASS - VERSION 2.7.                    
% C     *                         REPLACE BISECTION METHOD IN SURFACE     
% C     *                         TEMPERATURE ITERATION SCHEME WITH       
% C     *                         SECANT METHOD FOR FIRST TEN ITERATIONS.   
% C     *                         PASS QZERO,QA,ZOMS,ZOHS TO REVISED        
% C     *                         DRCOEF (ZOMS AND ZOHS ALSO NEW WORK ARRAYS
% C     *                         PASSED TO THIS ROUTINE).                 
% C     * JUN 20/97 - D.VERSEGHY. PASS IN NEW "CLASS4" COMMON BLOCK.      
% C     * JAN 02/96 - D.VERSEGHY. CLASS - VERSION 2.5.                    
% C     *                         COMPLETION OF ENERGY BALANCE            
% C     *                         DIAGNOSTICS.  ALSO, PASS SWITCH "ILW"   
% C     *                         THROUGH SUBROUTINE CALL, SPECIFYING     
% C     *                         WHETHER QLWIN REPRESENTS INCOMING       
% C     *                         (ILW=1) OR NET (ILW=2) LONGWAVE         
% C     *                         RADIATION ABOVE THE GROUND.             
% C     * NOV 30/94 - M.LAZARE.   CLASS - VERSION 2.3.                    
% C     *                         NEW DRAG COEFFICIENT AND RELATED FIELDS,
% C     *                         NOW DETERMINED IN ROUTINE "DRCOEF"      
% C     *                         "CFLUX" NOW WORK FIELD INSTEAD OF "CLIMIT". 
% C     * OCT 04/94 - D.VERSEGHY. CHANGE "CALL ABORT" TO "CALL XIT" TO       
% C     *                         ENABLE RUNNING ON PCS.                    
% C     * JAN 24/94 - M.LAZARE.   UNFORMATTED I/O COMMENTED OUT IN LOOP 200.
% C     * JUL 29/93 - D.VERSEGHY. CLASS - VERSION 2.2.                     
% C     *                         REMOVE RE-DEFINITION OF QMELT NEAR END   
% C     *                         (SINCE DONE ELSEWHERE ALREADY) AND      
% C     *                         REDEFINE QSWNET FOR DIAGNOSTIC PURPOSES 
% C     *                         TO INCLUDE TRANSMISSION THROUGH         
% C     *                         SNOWPACK.                               
% C     * OCT 15/92 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.1.           
% C     *                                  REVISED AND VECTORIZED CODE    
% C     *                                  FOR MODEL VERSION GCM7.        
% C     * AUG 12/91 - D.VERSEGHY. CODE FOR MODEL VERSION GCM7U -          
% C     *                         CLASS VERSION 2.0 (WITH CANOPY).        
% C     * APR 11/89 - D.VERSEGHY. ITERATIVE SURFACE TEMPERATURE           
% C     *                         CALCULATIONS FOR SNOW/SOIL.             
% C                                                                       
% IMPLICIT NONE                                                     
% 
% C     * INTEGER CONSTANTS.                                              
% C                                                                       
% INTEGER ISNOW,ISLFD,ITG,ILG,IG,IL1,IL2,JL,I,IB,NBS,ISNOALB        
% C                                                                       
% INTEGER NUMIT,NIT,IBAD,ITERMX                                     
% C                                                                       
% C     * OUTPUT ARRAYS.                                                  
% C                                                                       
% REAL QSWNET(ILG),    QLWOUT(ILG),    QTRANS(ILG),    QSENS (ILG), 
% 1     QEVAP (ILG),    EVAP  (ILG),    TZERO (ILG),    QZERO (ILG), 
% 2     GZERO (ILG),    QMELT (ILG),    CDH   (ILG),    CDM   (ILG), 
% 3     RIB   (ILG),    CFLUX (ILG),    FTEMP (ILG),    FVAP  (ILG), 
% 4     ILMO  (ILG),    UE    (ILG),    H     (ILG)                  
% C                                                                       
% C     * INPUT ARRAYS.                                                   
% C                                                                       
% REAL FI    (ILG),    QLWIN (ILG),                                 
% 1     TPOTA (ILG),    QA    (ILG),    VA    (ILG),    PADRY (ILG), 
% 2     RHOAIR(ILG),    ALVISG(ILG),    ALNIRG(ILG),    CRIB  (ILG), 
% 3     CPHCH (ILG),    CEVAP (ILG),    TVIRTA(ILG),                 
% 4     ZOSCLH(ILG),    ZOSCLM(ILG),    ZRSLFH(ILG),    ZRSLFM(ILG), 
% 5     ZOH   (ILG),    ZOM   (ILG),    GCONST(ILG),    GCOEFF(ILG), 
% 6     TSTART(ILG),    FCOR  (ILG),    PCPR  (ILG),
% 7     RHOSNO(ILG),    ZSNOW (ILG),    ZPOND (ILG)
% C
% REAL THLIQ (ILG,IG), THLMIN(ILG,IG), DELZW (ILG,IG)
% C                                                                       
% INTEGER          IWATER(ILG),        IEVAP (ILG)                  
% INTEGER          ITERCT(ILG,6,50),   ISAND(ILG,IG)                
% C                                                                       
% C     * BAND-DEPENDANT ARRAYS.                                          
% C                                                                       
% REAL TRSNOWG(ILG,NBS), ALSNO(ILG,NBS), FSSB(ILG,NBS),             
% 1     TRTOP  (ILG,NBS)                                             
% C                                                                       
% C     * INTERNAL WORK ARRAYS.                                           
% C                                                                       
% REAL TSTEP (ILG),    TVIRTS(ILG),    EVBETA(ILG),    Q0SAT (ILG), 
% 1     RESID (ILG),    DCFLXM(ILG),    CFLUXM(ILG),                 
% 2     A     (ILG),    B     (ILG),                                 
% 3     LZZ0  (ILG),    LZZ0T (ILG),    FM    (ILG),    FH    (ILG), 
% 4     WZERO (ILG),    EVPMAX(ILG)
% C                                                                       
% INTEGER              ITER  (ILG),    NITER (ILG),    JEVAP (ILG), 
% 1                     KF    (ILG)                                  
% C                                                                       
% C     * TEMPORARY VARIABLES.                                            
% C                                                                       
% REAL QSWNV,QSWNI,DCFLUX,DRDT0,TZEROT,QEVAPT,BOWEN,EZERO           
% C                                                                       
% C     * COMMON BLOCK PARAMETERS.                                        
% C                                                                       
% REAL DELT,TFREZ,RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN,HCPW,HCPICE,      
% 1     HCPSOL,HCPOM,HCPSND,HCPCLY,SPHW,SPHICE,SPHVEG,SPHAIR,        
% 2     RHOW,RHOICE,TCGLAC,CLHMLT,CLHVAP,DELTA,CGRAV,CKARM,CPD,      
% 3     AS,ASX,CI,BS,BETA,FACTN,HMIN,ANGMAX                          
% C                                                                       
% COMMON /CLASS1/ DELT,TFREZ                                        
% COMMON /CLASS2/ RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN                   
% COMMON /CLASS4/ HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,           
% 1                SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,            
% 2                TCGLAC,CLHMLT,CLHVAP                              
% COMMON /PHYCON/ DELTA,CGRAV,CKARM,CPD                             
% COMMON /CLASSD2/ AS,ASX,CI,BS,BETA,FACTN,HMIN,ANGMAX              
% C-----------------------------------------------------------------------
% C     * INITIALIZATION AND PRE-ITERATION SEQUENCE.                      
% C
[JEVAP,EVPMAX]=makeZeros(length(IL1:IL2),1);

if(ITG<2)
    ITERMX=50;                                                     
else                                                              
    ITERMX=5;                                                      
end                                                             
% C                                                                       
% C      if(ISNOW==0) %THEN                                              
% C          EZERO=0.0                                                    
% C      else                                                             
% C          EZERO=2.0                                                    
% C      end                                                            
EZERO=single(0.0);                                                        
% C                                                                       
for I=IL1:IL2%DO I=IL1,IL2                                                      
    QSWNET(I)=single(0.0);                                                  
    QTRANS(I)=single(0.0) ;                                                 
end%END DO                                                            
% C                                                                       
if(ISNOW == 0)%    %THEN ! Use usual snow-free bare soil formulation
    for I=IL1:IL2%DO I=IL1,IL2                                                   
        if(FI(I)>0.) %THEN                                        
            TRTOP(I,1)=0.;                                            
            QSWNV=FSSB(I,1)*(1.0-ALVISG(I));                        
            if (ISNOALB == 0) %THEN                                 
                QSWNI=FSSB(I,2)*(1.0-ALNIRG(I));                       
            elseif (ISNOALB == 1) %THEN                            
                QSWNI=0.0;                                             
                for IB = 2:NBS%DO IB = 2, NBS                                        
                    QSWNI=QSWNI+FSSB(I,IB)*(1.0-ALNIRG(I));             
                end %END DO ! IB                                           
            end                                                    
            QSWNET(I)=QSWNV+QSWNI;                                    
            QTRANS(I)=QSWNET(I)*TRTOP(I,1);                           
            QSWNET(I)=QSWNET(I)-QTRANS(I);                            
        end%END if                                                      
    end%END DO ! I                                                     
else                                                              
    if (ISNOALB == 0) %THEN ! Use the existing snow albedo and transmission 
        for I=IL1:IL2%DO I=IL1,IL2                                                
            if(FI(I)>0.) %THEN                                     
                TRTOP(I,1)=TRSNOWG(I,1);                               
                QSWNV=FSSB(I,1)*(1.0-ALSNO(I,1));                      
                QSWNI=FSSB(I,2)*(1.0-ALSNO(I,2));                      
                QSWNET(I)=QSWNV+QSWNI;                                
                QTRANS(I)=QSWNET(I)*TRTOP(I,1);                        
                QSWNET(I)=QSWNET(I)-QTRANS(I);                         
            end%END if                                                   
        end%END DO ! I                                                  
    elseif(ISNOALB == 1) %THEN ! Use the band-by-band snow albedo and transmission
        for I=IL1:IL2%DO I=IL1,IL2                                                
            QTRANS(I) = 0.0;                                          
            QSWNET(I) = 0.0;                                          
        end%END DO ! I                                                  
        for IB = 1:NBS%DO IB = 1, NBS                                              
            for I=IL1:IL2%DO I=IL1,IL2                                             
                if(FI(I)>0.) %THEN                                  
                    TRTOP(I,IB)=TRSNOWG(I,IB);                          
                    QSWNV=FSSB(I,IB)*(1.0-ALSNO(I,IB));                 
                    QSWNET(I)=QSWNET(I)+FSSB(I,IB)*(1.0-ALSNO(I,IB));   
                    QTRANS(I)=QTRANS(I)+QSWNV*TRTOP(I,IB) ;             
                end%END if                                                
            end%END DO ! I                                               
        end%END DO ! IB                                                 
        for I=IL1:IL2%DO I=IL1,IL2                                                
            if(FI(I)>0.) %THEN                                     
                QSWNET(I)=QSWNET(I)-QTRANS(I);                         
            end%END if                                                   
        end%END DO ! I                                                  
    end%END if ! ISNOALB                                               
end%END if ! ISNOW                                                    
% C                                                                       
for I=IL1:IL2%DO 50 I=IL1,IL2                                                   
    if(FI(I)>0.)                                          %THEN 
        TZERO(I)=TSTART(I);                                        
        TSTEP(I)=1.0;                                              
        ITER(I)=1;                                                 
        NITER(I)=1;                                                
        % C                                                                       
        QMELT(I)=0.0;                                              
        RESID(I)=999999.;                                          
        DCFLXM(I)=0.0;                                             
        CFLUX(I)=0.0;                                              
        if(ISNOW==1)                      %THEN                  
            KF(I)=3;                                               
            EVPMAX(I)=RHOSNO(I)*ZSNOW(I)/DELT;
        else                                                      
            KF(I)=6;                                               
            EVPMAX(I)=RHOW*(ZPOND(I)+(THLIQ(I,1)-THLMIN(I,1))*DELZW(I,1))/DELT;
        end                                                     
    end                                                         
end%50 CONTINUE                                                          
% C                                                                       
% C     * ITERATION SECTION.                                              
% C     * LOOP IS REPEATED UNTIL SOLUTIONS HAVE BEEN FOUND FOR ALL POINTS 
% C     * ON THE CURRENT LATITUDE CIRCLE(S).                              
% C  
continueLoop=true;
while continueLoop %100 CONTINUE                                                          
% C                                                                       
    NUMIT=0;                                                           
    NIT=0;                                                             
    for I=IL1:IL2%DO 150 I=IL1,IL2                                                  
        if(FI(I)>0. && ITER(I)==1)                       %THEN 
            NIT=NIT+1;                                                 
            CFLUXM(I)=CFLUX(I);                                        
            if(TZERO(I)>=TFREZ)                        %THEN         
                A(I)=17.269;                                           
                B(I)=35.86;                                            
            else                                                      
                A(I)=21.874;                                           
                B(I)=7.66;                                             
            end                                                     
            WZERO(I)=0.622*611.0*exp(A(I)*(TZERO(I)-TFREZ)/(TZERO(I)-B(I)))/PADRY(I);                           
            Q0SAT(I)=WZERO(I)/(1.0+WZERO(I));                          
            if(IWATER(I)>0)                              %THEN      
                EVBETA(I)=1.0;                                         
                QZERO(I)=Q0SAT(I);                                     
            else                                                      
                EVBETA(I)=CEVAP(I);                                    
                QZERO(I)=EVBETA(I)*Q0SAT(I)+(1.0-EVBETA(I))*QA(I);     
                if(QZERO(I)>QA(I) && IEVAP(I)==0) %THEN        
                    EVBETA(I)=0.0;                                     
                    QZERO(I)=QA(I);                                    
                end                                                 
            end                                                     
        TVIRTS(I)=TZERO(I)*(1.0+0.61*QZERO(I));                    
        end                                                         
    end%150 CONTINUE                                                          
    % C                                                                       
    if(NIT>0)                                                  %THEN
        % C                                                                       
        % C     * CALCULATE SURFACE DRAG COEFFICIENTS (STABILITY-DEPENDENT) AND   
        % C     * OTHER RELATED QUANTITIES.                                       
        % C                                                                       
        if(ISLFD<2) %THEN                                             
            % CALL DRCOEF (CDM,CDH,RIB,CFLUX,QZERO,QA,ZOSCLM,ZOSCLH,      
            % 1                   CRIB,TVIRTS,TVIRTA,VA,FI,ITER,                 
            % 2                   ILG,IL1,IL2) 
            [CDM,CDH,RIB,CFLUX] = ...,QZERO,QA,ZOSCLM,ZOSCLH,CRIB,TVIRTS,TVIRTA,VA,FI,ITER,ILG,IL1,IL2
                DRCOEF (CDM,CDH,RIB,CFLUX,ZOSCLM,ZOSCLH,CRIB,...,QZERO,QA
                    TVIRTS,TVIRTA,VA,FI,ITER,IL1,IL2,...ILG,
                    GRAV,VKC);
        else
            error('FLXSURFZ is not converted from FORTRAN, work in progress.  Set ISLFD switch to 1 to use DRCOEF instead');
            % CALL FLXSURFZ(CDM,CDH,CFLUX,RIB,FTEMP,FVAP,ILMO,            
            % 1                    UE,FCOR,TPOTA,QA,ZRSLFM,ZRSLFH,VA,            
            % 2                    TZERO,QZERO,H,ZOM,ZOH,                        
            % 3                    LZZ0,LZZ0T,FM,FH,ILG,IL1,IL2,FI,ITER,JL )     
        end                                                           
        % C                                                                       
        % C     * REMAINING CALCULATIONS.                                         
        % C                                                                       
        for I=IL1:IL2%DO 175 I=IL1,IL2                                                
            if(FI(I)>0. && ITER(I)==1)                       %THEN 
                QLWOUT(I)=SBC*TZERO(I)*TZERO(I)*TZERO(I)*TZERO(I);         
                if(TZERO(I)<TPOTA(I))                        %THEN      
                    QSENS(I)=(RHOAIR(I)*SPHAIR*CFLUX(I)+EZERO)*(TZERO(I)-TPOTA(I));                                         
                else                                                      
                    QSENS(I)=RHOAIR(I)*SPHAIR*CFLUX(I)*(TZERO(I)-TPOTA(I));                                         
                end                                                     
                EVAP(I)=RHOAIR(I)*CFLUX(I)*(QZERO(I)-QA(I));               
                if(EVAP(I)>EVPMAX(I)) 
                    EVAP(I)=EVPMAX(I);
                end
                QEVAP(I)=CPHCH(I)*EVAP(I);                                 
                GZERO(I)=GCOEFF(I)*TZERO(I)+GCONST(I);                     
                RESID(I)=QSWNET(I)+QLWIN(I)-QLWOUT(I)-QSENS(I)-QEVAP(I)-GZERO(I);                                         
                if(abs(RESID(I))<5.0)
                    ITER(I)=0;
                end
                if(abs(TSTEP(I))< 1.0E-2)                   
                    ITER(I)=0;
                end
                if(NITER(I)==ITERMX && ITER(I)==1)      
                    ITER(I)=-1;
                end
            end                                                         
        end%175     CONTINUE                                                        
    end                                                             
% C                                                                       
    if(ITG<2) %THEN                                                 
        % C                                                                       
        % C     * OPTION #1: BISECTION ITERATION METHOD.                          
        % C                                                                       
        if(NIT>0)                                                  %THEN
            for I=IL1:IL2%DO 180 I=IL1,IL2                                                
                if(FI(I)>0. && ITER(I)==1)                       %THEN 
                    if(NITER(I)==1) %THEN                                    
                        if(RESID(I)>0.0) %THEN                              
                            TZERO(I)=TZERO(I)+1.0;                             
                        else                                                  
                            TZERO(I)=TZERO(I)-1.0;                             
                        end                                                 
                    else                                                      
                        if((RESID(I)>0. && TSTEP(I)<0.) || (RESID(I)<0. && TSTEP(I)>0.))    %THEN    
                            TSTEP(I)=-TSTEP(I)/2.0;                            
                        end                                                 
                        TZERO(I)=TZERO(I)+TSTEP(I);                            
                    end                                                     
                end                                                         
% C                                                                       
                if(FI(I)>0. && ITER(I)==1)                       %THEN 
                    NITER(I)=NITER(I)+1;                                       
                    NUMIT=NUMIT+1 ;                                            
                end                                                         
            end%180   CONTINUE                                                        
        end                                                             
% C                                                                       
% c     DO 185 I=IL1,IL2                                                  
% C         if(FI(I)>0. && ITER(I)==-1)                      %THEN 
% C             WRITE(6,6250) I,JL,RESID(I),TZERO(I),RIB(I)               
% C6250         FORMAT('0GROUND ITERATION LIMIT',3X,2I3,3(F8.2,E12.4))    
% C         end                                                         
% c 185 CONTINUE                                                          
% C                                                                       
        if(NUMIT>0)
            %GO TO 100       
            continue;
        end
% C                                                                       
    else %Hard to tell if this else belongs to line 396 (ITG<2) or line 431                                                             
% C                                                                       
% C     * OPTION #2: NEWTON-RAPHSON ITERATION METHOD.                     
% C                                                                       
        if(NIT>0)                                                  %THEN
            for I=IL1:IL2%DO 190 I=IL1,IL2                                                
                if(FI(I)>0. && ITER(I)==1)                      %THEN  
                    if(NITER(I)>1)                                 %THEN    
                        DCFLUX=(CFLUX(I)-CFLUXM(I))/SIGN(MAX(.001,abs(TSTEP(I))),TSTEP(I));            
                        if(abs(TVIRTA(I)-TVIRTS(I))<0.4)                   
                            DCFLUX=MAX(DCFLUX,0.8*DCFLXM(I));
                        end
                        DCFLXM(I)=DCFLUX;                                      
                    else                                                      
                        DCFLUX=0.;                                             
                    end                                                     
                    DRDT0= -4.0*SBC*TZERO(I)^3-...                               
                        RHOAIR(I)*SPHAIR*(CFLUX(I)+MAX(0.,TZERO(I)-TPOTA(I))*...  
                        DCFLUX) -GCOEFF(I)+...                                    
                        CPHCH(I)*RHOAIR(I)*(CFLUX(I)*WZERO(I)*A(I)*...
                        EVBETA(I)*(B(I)-TFREZ)/((TZERO(I)-B(I))*...              
                        (1.0+WZERO(I)))^2-(QZERO(I)-QA(I))*DCFLUX);            
                    TSTEP(I)=-RESID(I)/DRDT0;                                  
                    TSTEP(I)=MAX(-10.,MIN(5.,TSTEP(I)));                      
                    TZERO(I)=TZERO(I)+TSTEP(I);                                
                    NITER(I)=NITER(I)+1;                                       
                    NUMIT=NUMIT+1;                                             
                end                                                         
            end%190   CONTINUE                                                        
        end                                                             
% C                                                                       
        if(NUMIT>0)
            %GO TO 100       
            continue;
        end
        % C                                                                       
        % C     * if CONVERGENCE HAS NOT BEEN REACHED, CALCULATE TEMPERATURE AND  
        % C     * FLUXES ASSUMING NEUTRAL STABILITY.                              
        % C                                                                       
        for I=IL1:IL2%DO 195 I=IL1,IL2                                                  
            NUMIT=0;                                                       
            JEVAP(I)=0;                                                    
            if(FI(I)>0. &&ITER(I)==-1)                       %THEN 
                TZEROT=TVIRTA(I)/(1.0+0.61*QZERO(I));                      
                if(abs(RESID(I))>50.) %THEN                             
                    TZERO(I)=TZEROT;                                       
                    if(TZERO(I).GE.TFREZ)                        %THEN     
                        A(I)=17.269;                                       
                        B(I)=35.86;                                        
                    else                                                  
                        A(I)=21.874;                                       
                        B(I)=7.66;                                         
                    end                                                 
                    WZERO(I)=0.622*611.0*exp(A(I)*(TZERO(I)-TFREZ)/...       
                        (TZERO(I)-B(I)))/PADRY(I);                        
                    Q0SAT(I)=WZERO(I)/(1.0+WZERO(I));                      
                    QZERO(I)=EVBETA(I)*Q0SAT(I)+(1.0-EVBETA(I))*QA(I);     
                    QLWOUT(I)=SBC*TZERO(I)*TZERO(I)*TZERO(I)*TZERO(I);     
                    GZERO(I)=GCOEFF(I)*TZERO(I)+GCONST(I);                 
                    RESID(I)=QSWNET(I)+QLWIN(I)-QLWOUT(I)-GZERO(I);        
                    if(RESID(I)>0.)                 %THEN               
                        QEVAP(I)=RESID(I);                                 
                    else                                                  
                        QEVAP(I)=RESID(I)*0.5;                             
                    end                                                 
                    if(IEVAP(I)==0) 
                        QEVAP(I)=0.0;
                    end
                    QSENS(I)=RESID(I)-QEVAP(I);                            
                    RESID(I)=0.;                                           
                    EVAP(I)=QEVAP(I)/CPHCH(I);                             
                    TVIRTS(I)=TZERO(I)*(1.0+0.61*QZERO(I));                
                    JEVAP(I)=1;                                            
                    NUMIT=NUMIT+1;                                         
                end                                                     
            end                                                         
        end%195 CONTINUE                                                          
% C                                                                       
        if(NUMIT>0)                   %THEN                             
            if(ISLFD<2) %THEN                                             
                % CALL DRCOEF (CDM,CDH,RIB,CFLUX,QZERO,QA,ZOSCLM,ZOSCLH,      
                % 1                   CRIB,TVIRTS,TVIRTA,VA,FI,JEVAP,                
                % 2                   ILG,IL1,IL2)
                [CDM,CDH,RIB,CFLUX] = ...,QZERO,QA,ZOSCLM,ZOSCLH,CRIB,TVIRTS,TVIRTA,VA,FI,ITER,ILG,IL1,IL2
                    DRCOEF (CDM,CDH,RIB,CFLUX,ZOSCLM,ZOSCLH,CRIB,...,QZERO,QA
                        TVIRTS,TVIRTA,VA,FI,ITER,IL1,IL2,...ILG,
                        GRAV,VKC);
            else     
                error('FLXSURFZ is not converted from FORTRAN, work in progress.  Set ISLFD switch to 1 to use DRCOEF instead');
                % CALL FLXSURFZ(CDM,CDH,CFLUX,RIB,FTEMP,FVAP,ILMO,            
                % 1                    UE,FCOR,TPOTA,QA,ZRSLFM,ZRSLFH,VA,            
                % 2                    TZERO,QZERO,H,ZOM,ZOH,                        
                % 3                    LZZ0,LZZ0T,FM,FH,ILG,IL1,IL2,FI,JEVAP,JL )    
            end                                                           
        end                                                             
% C                                                                       
    end
    continueLoop=false;%Break out of the while loop. Could use break.
end
% C                                                                       
% C     * CHECK FOR BAD ITERATION TEMPERATURES.                           
% C                                                                       
IBAD=0;                                                            
for I=IL1:IL2%DO 200 I=IL1,IL2                                                  
    if(FI(I)>0. && (TZERO(I)<123.16 || TZERO(I)>373.16))               %THEN
        IBAD=I;                                                    
    end                                                         
end%200  CONTINUE                                                          
% C                                                                       
if(IBAD~=0)                                                 %THEN
    % WRITE(6,6275) IBAD,JL,TZERO(IBAD),NITER(IBAD),ISNOW           
    % 6275     FORMAT('0BAD ITERATION TEMPERATURE',3X,2I3,F16.2,2I4)         
    % WRITE(6,6280) QSWNET(IBAD),QLWIN(IBAD),QSENS(IBAD),           
    % 1        QEVAP(IBAD),GZERO(IBAD),CFLUX(IBAD),RIB(IBAD)             
    % 6280     FORMAT(2X,7F12.4)                                             
    % CALL XIT('TSOLVE',-1)                                         
    XIT('TSOLVE',-1);
end                                                             
% C                                                                       
% C     * POST-ITERATION CLEAN-UP.                                        
% C                                                                       
NIT=0;                                                             
for I=IL1:IL2%DO 300 I=IL1,IL2                                                  
    if(FI(I)>0.)                                          %THEN 
        if(((IWATER(I)==1 && TZERO(I)<TFREZ) ||...         
            (IWATER(I)==2 && TZERO(I)>TFREZ)) ||...        
            (ISAND(I,1)==-4 && TZERO(I)>TFREZ))   %THEN    
            TZERO(I)=TFREZ;                                        
            WZERO(I)=0.622*611.0/PADRY(I);                         
            QZERO(I)=WZERO(I)/(1.0+WZERO(I));                      
            TVIRTS(I)=TZERO(I)*(1.0+0.61*QZERO(I));                
            ITER(I)=1;                                             
            NIT=NIT+1;                                             
        else                                                      
            ITER(I)=0;                                             
        end                                                     
    end                                                         
end%300 CONTINUE                                                          
% C                                                                       
if(NIT>0)                                                  %THEN
    % C                                                                       
    % C       * CALCULATE SURFACE DRAG COEFFICIENTS (STABILITY-DEPENDENT) AND 
    % C       * OTHER RELATED QUANTITIES.                                     
    % C                                                                       
    if(ISLFD<2) %THEN                                             
        % CALL DRCOEF (CDM,CDH,RIB,CFLUX,QZERO,QA,ZOSCLM,ZOSCLH,      
        % 1                   CRIB,TVIRTS,TVIRTA,VA,FI,ITER,                 
        % 2                   ILG,IL1,IL2)   
        [CDM,CDH,RIB,CFLUX] = ...,QZERO,QA,ZOSCLM,ZOSCLH,CRIB,TVIRTS,TVIRTA,VA,FI,ITER,ILG,IL1,IL2
            DRCOEF (CDM,CDH,RIB,CFLUX,ZOSCLM,ZOSCLH,CRIB,...,QZERO,QA
                TVIRTS,TVIRTA,VA,FI,ITER,IL1,IL2,...ILG,
                GRAV,VKC);
    else
        error('FLXSURFZ is not converted from FORTRAN, work in progress.  Set ISLFD switch to 1 to use DRCOEF instead');
    % CALL FLXSURFZ(CDM,CDH,CFLUX,RIB,FTEMP,FVAP,ILMO,            
    % 1                    UE,FCOR,TPOTA,QA,ZRSLFM,ZRSLFH,VA,            
    % 2                    TZERO,QZERO,H,ZOM,ZOH,                        
    % 3                    LZZ0,LZZ0T,FM,FH,ILG,IL1,IL2,FI,ITER,JL )     
    end                                                           
end                                                             
% C                                                                       
% C     * REMAINING CALCULATIONS.                                         
% C                                                                       
for I=IL1:IL2%DO 350 I=IL1,IL2                                                  
    if(FI(I)>0. && ITER(I)==1)                       %THEN 
        QLWOUT(I)=SBC*TZERO(I)*TZERO(I)*TZERO(I)*TZERO(I);        
        if(TZERO(I)<TPOTA(I))                        %THEN      
            QSENS(I)=(RHOAIR(I)*SPHAIR*CFLUX(I)+EZERO)*(TZERO(I)-... 
            TPOTA(I));                                         
        else                                                      
            QSENS(I)=RHOAIR(I)*SPHAIR*CFLUX(I)*(TZERO(I)-...         
                TPOTA(I));                                         
        end                                                     
        EVAP(I)=RHOAIR(I)*CFLUX(I)*(QZERO(I)-QA(I));               
        if(EVAP(I)>EVPMAX(I)) 
            EVAP(I)=EVPMAX(I);
        end
        QEVAP(I)=CPHCH(I)*EVAP(I);                                 
        GZERO(I)=GCOEFF(I)*TZERO(I)+GCONST(I);                     
        QMELT(I)=QSWNET(I)+QLWIN(I)-QLWOUT(I)-QSENS(I)-QEVAP(I)-...  
            GZERO(I);                                         
        RESID(I)=0.0;                                              
        if(QMELT(I)<0.0) %THEN                                  
            QMELT(I)=QMELT(I)+QEVAP(I);                            
            QEVAP(I)=0.0;                                          
            EVAP(I) =0.0;                                          
        end                                                     
    end                                                         
    % C                                                                       
    if(FI(I)>0.)                                 %THEN          
        if(abs(EVAP(I))<1.0E-8) %THEN                           
            RESID(I)=RESID(I)+QEVAP(I);                            
            EVAP(I)=0.0;                                           
            QEVAP(I)=0.0;                                          
        end                                                     
        if((ISNOW==1 && QMELT(I)<0.0) ||...                
                (ISNOW==0 && QMELT(I)>0.0))     %THEN          
            GZERO(I)=GZERO(I)+QMELT(I);                            
            QMELT(I)=0.0;                                          
        end                                                     
    % C              QSENS(I)=QSENS(I)+0.5*RESID(I)                           
    % C              GZERO(I)=GZERO(I)+0.5*RESID(I)                           
        QSENS(I)=QSENS(I)+RESID(I);                                
        QSWNET(I)=QSWNET(I)+QTRANS(I);                             
        EVAP(I)=EVAP(I)/RHOW;                                     
        ITERCT(I,KF(I),NITER(I))=ITERCT(I,KF(I),NITER(I))+1;       
    end                                                         
end%350 CONTINUE                                                          
% C                                                                       
%RETURN                                                            
%END                                                               
end
