function [RHOSNO,ZSNOW,HCPSNO,TSNOW,EVAP,QFN,QFG,HTCS,...
          WLOST,TRUNOF,RUNOFF,TOVRFL,OVRFLW,...
          WSNOW] = ...FI,R,S,RHOSNI,ILG,IL1,IL2,JL
    SNOVAP(RHOSNO,ZSNOW,HCPSNO,TSNOW,EVAP,QFN,QFG,HTCS,...
          WLOST,TRUNOF,RUNOFF,TOVRFL,OVRFLW,...
          FI,R,S,RHOSNI,WSNOW,IL1,IL2,...ILG,,JL
          DELT,TFREZ,RHOW,HCPICE,RHOICE,HCPW)%globals
% C
% C     * AUG 25/11 - D.VERSEGHY. CORRECT CALCULATION OF TRUNOF
% C     *                         AND TOVRFL.
% C     * FEB 22/07 - D.VERSEGHY. NEW ACCURACY LIMITS FOR R AND S.
% C     * MAR 24/06 - D.VERSEGHY. ALLOW FOR PRESENCE OF WATER IN SNOW.
% C     * SEP 23/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
% C     * JUL 26/02 - D.VERSEGHY. CHANGE RHOSNI FROM CONSTANT TO
% C     *                         VARIABLE.
% C     * APR 11/01 - M.LAZARE.   CHECK FOR EXISTENCE OF SNOW BEFORE
% C     *                         PERFORMING CALCULATIONS.
% C     * JUN 20/97 - D.VERSEGHY. CLASS - VERSION 2.7.
% C     *                         PASS IN NEW "CLASS4" COMMON BLOCK.
% C     * JAN 02/96 - D.VERSEGHY. CLASS - VERSION 2.5.
% C     *                         COMPLETION OF ENERGY BALANCE
% C     *                         DIAGNOSTICS.
% C     * AUG 16/95 - D.VERSEGHY. CLASS - VERSION 2.4.
% C     *                         INCORPORATE DIAGNOSTIC ARRAY "WLOST". 
% C     * DEC 22/94 - D.VERSEGHY. CLASS - VERSION 2.3.
% C     *                         ADDITIONAL DIAGNOSTIC CALCULATION -
% C     *                         UPDATE HTCS.
% C     * JUL 30/93 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.2.
% C     *                                  NEW DIAGNOSTIC FIELDS.
% C     * APR 24/92 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.1.
% C     *                                  REVISED AND VECTORIZED CODE
% C     *                                  FOR MODEL VERSION GCM7.
% C     * AUG 12/91 - D.VERSEGHY. CODE FOR MODEL VERSION GCM7U -
% C     *                         CLASS VERSION 2.0 (WITH CANOPY).
% C     * APR 11/89 - D.VERSEGHY. SUBLIMATION FROM SNOWPACK.
% C                                          
%       IMPLICIT NONE
% C
% C     * INTEGER CONSTANTS.
% C
%       INTEGER ILG,IL1,IL2,JL,I
% C
% C     * INPUT/OUTPUT ARRAYS.
% C
%       REAL RHOSNO(ILG),   ZSNOW (ILG),   HCPSNO(ILG),   TSNOW (ILG), 
%      1     EVAP  (ILG),   QFN   (ILG),   QFG   (ILG),   HTCS  (ILG),
%      2     WLOST (ILG),   TRUNOF(ILG),   RUNOFF(ILG),   TOVRFL(ILG),
%      3     OVRFLW(ILG)
% C
% C     * INPUT ARRAYS.
% C
%       REAL FI    (ILG),   R     (ILG),   S     (ILG),   RHOSNI(ILG),
%      1     WSNOW (ILG)   
% C
% C     * TEMPORARY VARIABLES.
% C
%       REAL ZADD,ZLOST,ZREM
% C
% C     * COMMON BLOCK PARAMETERS.
% C
%       REAL DELT,TFREZ,HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,SPHW,
%      1     SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,TCGLAC,CLHMLT,CLHVAP
% C                                       
%       COMMON /CLASS1/ DELT,TFREZ
%       COMMON /CLASS4/ HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,
%      1                SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,
%      2                TCGLAC,CLHMLT,CLHVAP
% C-----------------------------------------------------------------------

for I=IL1:IL2%DO 100 I=IL1,IL2
    if(FI(I)>0. && (S(I)<1.0E-11 || R(I)<1.0E-11)&& ZSNOW(I)>0.)
        HTCS(I)=HTCS(I)-FI(I)*HCPSNO(I)*(TSNOW(I)+TFREZ)*ZSNOW(I)/DELT;
        if(EVAP(I)<0.)
            ZADD=-EVAP(I)*DELT*RHOW/RHOSNI(I);
            RHOSNO(I)=(ZSNOW(I)*RHOSNO(I)+ZADD*RHOSNI(I))/(ZSNOW(I)+ZADD);                          
            ZSNOW (I)=ZSNOW(I)+ZADD;                                                        
            HCPSNO(I)=HCPICE*RHOSNO(I)/RHOICE+HCPW*WSNOW(I)/(RHOW*ZSNOW(I));
            EVAP  (I)=0.0 ;                                                               
        else                                                                        
            ZLOST=EVAP(I)*DELT*RHOW/RHOSNO(I);
            if(ZLOST<=ZSNOW(I))
                ZSNOW(I)=ZSNOW(I)-ZLOST;                                                   
                EVAP (I)=0.0;                                                            
                HCPSNO(I)=HCPICE*RHOSNO(I)/RHOICE+HCPW*WSNOW(I)/(RHOW*ZSNOW(I));
            else                                                                   
                ZREM=(ZLOST-ZSNOW(I))*RHOSNO(I)/RHOW;
                ZSNOW(I)=0.0;                                                           
                HCPSNO(I)=0.0;
                EVAP(I)=ZREM*(CLHMLT+CLHVAP)/(CLHVAP*DELT);
                WLOST(I)=WLOST(I)-ZREM*RHOW*CLHMLT/CLHVAP;
                if(RUNOFF(I)>0. || WSNOW(I)>0.) 
                    TRUNOF(I)=(TRUNOF(I)*RUNOFF(I)+(TSNOW(I)+TFREZ)*WSNOW(I)/RHOW)/(RUNOFF(I)+WSNOW(I)/RHOW);
                end
                RUNOFF(I)=RUNOFF(I)+WSNOW(I)/RHOW;
                if(OVRFLW(I)>0. || WSNOW(I)>0.)
                    TOVRFL(I)=(TOVRFL(I)*OVRFLW(I)+(TSNOW(I)+TFREZ)*FI(I)*WSNOW(I)/RHOW)/(OVRFLW(I)+FI(I)*WSNOW(I)/RHOW);
                end
                OVRFLW(I)=OVRFLW(I)+FI(I)*WSNOW(I)/RHOW;
                TSNOW(I)=0.0; 
                WSNOW(I)=0.0;
                QFN(I)=QFN(I)-FI(I)*ZREM*RHOW/DELT;
                QFG(I)=QFG(I)+FI(I)*EVAP(I)*RHOW;
            end
        end
        HTCS(I)=HTCS(I)+FI(I)*HCPSNO(I)*(TSNOW(I)+TFREZ)*ZSNOW(I)/DELT;
    end
end%100 CONTINUE
% C                                                                                  
%       RETURN                                                                      
%       END  
end