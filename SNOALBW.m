function [ALBSNO,RHOSNO,ZSNOW,HCPSNO,RHOMAX]= ...TSNOW,FI,S,RMELT,WSNOW,ISAND,ILG,IG,IL1,IL2,JL
    SNOALBW(ALBSNO,RHOSNO,ZSNOW,HCPSNO,TSNOW,FI,S,RMELT,WSNOW,RHOMAX,...
      IL1,IL2,...ISAND,ILG,IG,JL,
      DELT,RHOW,HCPICE,RHOICE,HCPW)%globals
% C                                                                       
% C     * APR 17/14 - D.VERSEGHY. MAKE SNOW ALBEDO REFRESHMENT VALUE
% C     *                         CONSISTENT WITH SNOADD.
% C     * MAR 07/07 - D.VERSEGHY. STREAMLINE SOME CALCULATIONS.           
% C     * MAR 24/06 - D.VERSEGHY. ALLOW FOR PRESENCE OF WATER IN SNOW.    
% C     * SEP 23/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.            
% C     * SEP 10/04 - R.HARVEY/D.VERSEGHY. INCREASE SNOW ALBEDO           
% C     *                         REFRESHMENT AND WETTING THRESHOLDS.     
% C     * AUG 04/04 - Y.DELAGE/D.VERSEGHY. PROTECT SENSITIVE              
% C     *                         CALCULATIONS AGAIST ROUNDOFF            
% C     *                         ERRORS.                                 
% C     * APR 21/04 - F.SEGLENIEKS/D.VERSEGHY. BUG FIX IN SNOW            
% C     *                         TEMPERATURE COMPARISONS.                
% C     * JUL 26/02 - D.VERSEGHY. SHORTENED CLASS4 COMMON BLOCK.          
% C     * OCT 20/00 - R.BROWN/D.VERSEGHY. MODIFIED SNOW DENSITY           
% C     *                                 CALCULATIONS, ACCOUNTING        
% C     *                                 FOR SETTLING IN WARM AND        
% C     *                                 COLD SNOW.                      
% C     * JUN 05/97 - D.VERSEGHY. CLASS - VERSION 2.7.                    
% C     *                         SPECIFY LOCATION OF ICE SHEETS          
% C     *                         BY SOIL TEXTURE ARRAY RATHER            
% C     *                         THAN BY SOIL COLOUR INDEX.              
% C     * JAN 02/96 - D.VERSEGHY. CLASS - VERSION 2.5.                    
% C     *                         COMPLETION OF ENERGY BALANCE            
% C     *                         DIAGNOSTICS.                            
% C     * MAR 13/92 - M.LAZARE.   CLASS - VERSION 2.1.                    
% C     *                         CODE FOR MODEL VERSION GCM7 -           
% C     *                         DIVIDE PREVIOUS SUBROUTINE              
% C     *                         "SNOALB" INTO "SNOALBA" AND             
% C     *                         "SNOALBW" AND VECTORIZE.                
% C     * AUG 12/91 - D.VERSEGHY. CODE FOR MODEL VERSION GCM7U -          
% C     *                         CLASS VERSION 2.0 (WITH CANOPY).        
% C     * APR 11/89 - D.VERSEGHY. CALCULATE DECREASE IN SNOW ALBEDO       
% C     *                         AND INCREASE IN DENSITY DUE TO          
% C     *                         AGING. (ASSIGN DIFFERENT LOWER          
% C     *                         SNOW ALBEDO LIMITS FOR DRY AND          
% C     *                         MELTING SNOW.)                          
% C                                                                       
%       IMPLICIT NONE                                                     
% C                                                                       
% C     * INTEGER CONSTANTS.                                              
% C                                                                       
%       INTEGER ILG,IG,IL1,IL2,JL,I,IPTBAD                                
% C                                                                       
% C     * OUTPUT ARRAYS.                                                  
% C                                                                       
%       REAL ALBSNO(ILG),   RHOSNO(ILG),   ZSNOW (ILG),   HCPSNO(ILG)     
% C                                                                       
% C     * INPUT ARRAYS.                                                   
% C                                                                       
%       REAL TSNOW (ILG),   FI    (ILG),   S     (ILG),   RMELT (ILG),    
%      1     WSNOW (ILG)                                                  
% C                                                                       
%       INTEGER             ISAND (ILG,IG)                                
% C                                                                       
% C     * WORK ARRAY.                                                     
% C                                                                       
%       REAL RHOMAX(ILG)                                                  
% C                                                                       
% C     * TEMPORARY VARIABLES.                                            
% C                                                                       
%       REAL TIMFAC,RHOOLD                                                
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
% C----------------------------------------------------------------------

IPTBAD=single(0);                                                          
for I=IL1:IL2%DO 100 I=IL1,IL2                                                  
    if(ZSNOW(I)>0. && FI  (I)>0. && S(I)*DELT<1.0E-4)
        if(ALBSNO(I)>0.5001 && (RMELT(I)>1.0E-7 || TSNOW(I)>=-0.01))
            ALBSNO(I)=(ALBSNO(I)-0.50)*exp(-0.01*DELT/3600.0)+ 0.50;                                              
        elseif(ALBSNO(I)>0.7001 && RMELT(I)<=1.0E-7)
            ALBSNO(I)=(ALBSNO(I)-0.70)*exp(-0.01*DELT/3600.0)+0.70;                                              
        end                                                     
    end                                                         
%C                                                                       
    if(FI(I)>0. && ZSNOW(I)>0.0001)
        if(TSNOW(I)<-0.01)
            RHOMAX(I)=450.0-(204.7/ZSNOW(I))*(1.0-exp(-ZSNOW(I)/0.673));                        
        else                                                      
            RHOMAX(I)=700.0-(204.7/ZSNOW(I))*(1.0-exp(-ZSNOW(I)/0.673));                        
        end                                                     
    end                                                         
%C                                                                       
    if(FI(I)>0. && ZSNOW(I)>0.0001 && RHOSNO(I)<(RHOMAX(I)-0.01))
        RHOOLD=RHOSNO(I);                                          
        RHOSNO(I)=(RHOSNO(I)-RHOMAX(I))*exp(-0.01*DELT/3600.0)+RHOMAX(I);                                             
        ZSNOW(I)=ZSNOW(I)*RHOOLD/RHOSNO(I);                        
        HCPSNO(I)=HCPICE*RHOSNO(I)/RHOICE+HCPW*WSNOW(I)/(RHOW*ZSNOW(I));                                       
    end                                                          
    if((ALBSNO(I)<0.49 || ALBSNO(I)>1.0) && ZSNOW (I)>0. && FI(I)>0.)
        IPTBAD=I;  
    end
end%100 CONTINUE                                                          
%C                                                                       
if(IPTBAD~=0)
    %WRITE(6,6100) IPTBAD,JL,ALBSNO(IPTBAD)                         
    %6100    FORMAT('0AT (I,J)= (',I3,',',I3,'), ALBSNO = ',F10.5)          
    %write error code, but for now just exit
    XIT('SNOALBW',-1)                                         
end                                                             
%C                                                                       
%RETURN                                                            
%END                                                               
end
