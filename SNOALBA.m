function [ALVSSN,ALIRSN,ALVSSC,ALIRSC,...ALBSNO,            
     TRSNOWC, ALSNO, TRSNOWG]=..., ...FSDB, FSFB, RHOSNO,REFSN,BCSN,SNO,CSZ,ZSNOW,FSNOW,ASVDAT,ASIDAT,ALVSG, ALIRG,... ILG,IG,IL1,IL2,JL,IALS,NBS,ISNOALB             
   SNOALBA(ALVSSN,ALIRSN,ALBSNO,...ALVSSC,ALIRSC,            
     ALSNO, TRSNOWG, ...TRSNOWC, FSDB, FSFB, RHOSNO,    
     ZSNOW,FSNOW,ASVDAT,ASIDAT,...REFSN,BCSN,SNO,CSZ,ALVSG, ALIRG,  
     IL1,IL2,JL,IALS,NBS,ISNOALB)% ILG,IG, 
% C                                                                       
% C     * MAR 08/16 - M.Mackay    4-BAND SOLAR RADIATION DISABLED
% C     * NOV 16/13 - J.COLE.     Final version for gcm17:                
% C     *                         - Fixes to get the proper BC mixing ratio in 
% C     *                           snow, which required passing in and using  
% C     *                           the snow density RHON.                
% C     * JUN 22/13 - J.COLE.     ADD CODE FOR "ISNOALB" OPTION,          
% C     *                         WHICH IS BASED ON 4-BAND SOLAR.         
% C     * FEB 05/07 - D.VERSEGHY. STREAMLINE CALCULATIONS OF              
% C     *                         ALVSSN AND ALIRSN.                      
% C     * APR 13/06 - D.VERSEGHY. SEPARATE ALBEDOS FOR OPEN AND           
% C     *                         CANOPY-COVERED SNOW.                    
% C     * NOV 03/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.            
% C     * MAR 18/02 - D.VERSEGHY. UPDATES TO ALLOW ASSIGNMENT OF          
% C     *                         USER-SPECIFIED VALUES TO SNOW           
% C     *                         ALBEDO.                                 
% C     * JUN 05/97 - D.VERSEGHY. CLASS - VERSION 2.7.                    
% C     *                         SPECIFY LOCATION OF ICE SHEETS          
% C     *                         BY SOIL TEXTURE ARRAY RATHER            
% C     *                         THAN BY SOIL COLOUR INDEX.              
% C     * NOV 29/94 - M.LAZARE.   CLASS - VERSION 2.3.                    
% C     *                         CALL ABORT CHANGED TO CALL XIT TO       
% C     *                         ENABLE RUNNING ON PC'S.                 
% C     * MAR 13/92 - M.LAZARE.   CODE FOR MODEL VERSION GCM7 -           
% C     *                         DIVIDE PREVIOUS SUBROUTINE              
% C     *                         "SNOALB" INTO "SNOALBA" AND             
% C     *                         "SNOALBW" AND VECTORIZE.                
% C     * AUG 12/91 - D.VERSEGHY. CODE FOR MODEL VERSION GCM7U -          
% C     *                         CLASS VERSION 2.0 (WITH CANOPY).        
% C     * APR 11/89 - D.VERSEGHY. DISAGGREGATE SNOW ALBEDO INTO           
% C     *                         VISIBLE AND NEAR-IR PORTIONS;           
% C     *                         CALCULATE TRANSMISSIVITY TO             
% C     *                         SHORTWAVE RADIATION.                    
% C                                                                       
%       IMPLICIT NONE                                                     
% C                                                                       
% C     * INTEGER CONSTANTS.                                              
% C                                                                       
%       INTEGER ILG,IG,IL1,IL2,JL,IALS,IPTBAD,I,IB,NBS,ISNOALB            
% C                                                                       
% C     * OUTPUT ARRAYS.                                                  
% C                                                                       
%       REAL   ALSNO(ILG,NBS), TRSNOWG(ILG,NBS)                           
%       REAL   ALVSSN(ILG),  ALIRSN(ILG),  ALVSSC(ILG),  ALIRSC(ILG),     
%      1       ALVSG (ILG),  ALIRG (ILG),  TRSNOWC(ILG)                   
% C                                                                       
% C     * INPUT ARRAYS.                                                   
% C                                                                       
%       REAL   FSDB(ILG,NBS), FSFB(ILG,NBS)                               
%       REAL   ALBSNO(ILG),  ZSNOW (ILG),  FSNOW (ILG),                   
%      1       ASVDAT(ILG),  ASIDAT(ILG),  REFSN (ILG),  BCSN  (ILG),     
%      2       CSZ   (ILG),  SNO   (ILG),  RHOSNO(ILG)                    
% C                                                                       
% C     * LOCAL ARRAYS                                                    
% C                                                                       
%       REAL SALBG(ILG,NBS), ALDIR(ILG,NBS), ALDIF(ILG,NBS),              
%      +                     TRDIR(ILG,NBS), TRDIF(ILG,NBS)               
%       REAL REFSNO(ILG), BCSNO(ILG)                                      
%       INTEGER C_FLAG(ILG)                                               
% C                                                                       
% C     * CONSTANTS.                                                      
% C                                                                       
%       REAL WDIRCT, WDIFF                                                
%       INTEGER SUM_C_FLAG                                                
% C------------------------------------------------------------------ 
if (ISNOALB ~= 0) %  !4-BAND SW DISABLED 
     XIT('SNOALBA',-1);%error('SNOALBA issues');                                         
end 
IPTBAD=single(0);                                                          
for I=IL1:IL2                                                  
    if(ALBSNO(I)<0.50 && ALBSNO(I)>0.499) 
        ALBSNO(I)=0.50;   
    end
    if(FSNOW(I)>0.0 && IALS==0) 
        if(ALBSNO(I)>0.70)         
            ALVSSN(I)=0.79*(ALBSNO(I)-0.70)+0.84;
            ALIRSN(I)=1.21*(ALBSNO(I)-0.70)+0.56;                   
        else                                                       
            ALVSSN(I)=0.97*(ALBSNO(I)-0.50)+0.62;
            ALIRSN(I)=1.03*(ALBSNO(I)-0.50)+0.38;                   
        end                                                      
        if(ALVSSN(I)>0.999 || ALVSSN(I)<0.001) 
            IPTBAD=I;
        end
        if(ALIRSN(I)>0.999 || ALIRSN(I)<0.001) 
            IPTBAD=I;
        end
    elseif(FSNOW(I)>0.0 && IALS==1) 
        ALVSSN(I)=ASVDAT(I);
        ALIRSN(I)=ASIDAT(I);                                        
    end                                                          
    ALVSSC(I)=ALVSSN(I);
    ALIRSC(I)=ALIRSN(I);
    TRSNOWC(I)=exp(-25.0*ZSNOW(I));
end                                                     
%C   

if(IPTBAD~=0)    
    %WRITE(6,6100) IPTBAD,JL,ALVSSN(IPTBAD),ALIRSN(IPTBAD)          
    %6100    FORMAT('0AT (I,J)= (',I3,',',I3,'), ALVSSN,ALIRSN = ',2F10.5)  
    XIT('SNOALBA',-1);%error('SNOALBA IPTBAD~=0')
end                                                             
%C                                                                       
if (ISNOALB == 0) 
    for I = IL1:IL2                                                
        ALSNO(I,1) = ALVSSN(I);
        ALSNO(I,2) = ALIRSN(I);
        ALSNO(I,3) = ALIRSN(I);
        ALSNO(I,4) = ALIRSN(I);
        TRSNOWG(I,1:NBS) = TRSNOWC(I);
    end
end
%In the fortran code, there was extensive lines beyond this that were all
%commented out.
end