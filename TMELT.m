function [ZSNOW,TSNOW,QMELT,R,TR,GZERO,RALB,...
            HMFN,HTCS,HTC,HCPSNO,WSNOW...FI,,RHOSNO,
            ] =...ISAND,IG,ILG,IL1,IL2,JL
    TMELT (ZSNOW,TSNOW,QMELT,R,TR,GZERO,RALB,...
            HMFN,HTCS,HTC,FI,HCPSNO,RHOSNO,WSNOW,...
            ISAND,IL1,IL2,...IG,ILG,JL,
            DELT,TFREZ,RHOW,HCPICE,RHOICE,HCPW,CLHMLT)
% C                                                                                  
% C     * JAN 06/09 - D.VERSEGHY/M.LAZARE. SPLIT 100 LOOP INTO TWO.
% C     * MAR 24/06 - D.VERSEGHY. ALLOW FOR PRESENCE OF WATER IN SNOW.
% C     * SEP 24/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
% C     * JUL 26/02 - D.VERSEGHY. SHORTENED CLASS4 COMMON BLOCK.
% C     * JUN 20/97 - D.VERSEGHY. CLASS - VERSION 2.7.
% C     *                         MODIFICATIONS TO ALLOW FOR VARIABLE
% C     *                         SOIL PERMEABLE DEPTH.
% C     * JAN 02/95 - D.VERSEGHY. CLASS - VERSION 2.5.
% C     *                         COMPLETION OF ENERGY BALANCE
% C     *                         DIAGNOSTICS.
% C     * AUG 18/95 - D.VERSEGHY. CLASS - VERSION 2.4.
% C     *                         REVISIONS TO ALLOW FOR INHOMOGENEITY
% C     *                         BETWEEN SOIL LAYERS.
% C     * JUL 30/93 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.2.
% C     *                                  NEW DIAGNOSTIC FIELDS.
% C     * APR 24/92 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.1.
% C     *                                  REVISED AND VECTORIZED CODE
% C     *                                  FOR MODEL VERSION GCM7.
% C     * AUG 12/91 - D.VERSEGHY. CODE FOR MODEL VERSION GCM7U -
% C     *                         CLASS VERSION 2.0 (WITH CANOPY).
% C     * APR 11/89 - D.VERSEGHY. MELTING OF SNOWPACK.
% C
%       IMPLICIT NONE
% C
% C     * INTEGER CONSTANTS.
% C
%       INTEGER IG,ILG,IL1,IL2,JL,I
% C
% C     * INPUT/OUTPUT ARRAYS.
% C
%       REAL HTC   (ILG,IG)
% 
%       REAL ZSNOW (ILG),   TSNOW (ILG),   QMELT (ILG),   R     (ILG),
%      1     TR    (ILG),   GZERO (ILG),   RALB  (ILG),   HMFN  (ILG),
%      2     HTCS  (ILG)
% C
% C     * INPUT ARRAYS.
% C
%       REAL FI    (ILG),   HCPSNO(ILG),   RHOSNO(ILG),   WSNOW (ILG)
% C
%       INTEGER             ISAND (ILG,IG)
% C                                                                                 
% C     * TEMPORARY VARIABLES.
% C
%       REAL HADD,HCONV,ZMELT,RMELT,RMELTS,TRMELT
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
    if(FI(I)>0.)
        if(QMELT(I)>0. && ZSNOW(I)>0.)
            HTCS(I)=HTCS(I)-FI(I)*HCPSNO(I)*(TSNOW(I)+TFREZ)*ZSNOW(I)/DELT;
            HADD =QMELT(I)*DELT;                                                             
            HCONV=(0.0-TSNOW(I))*HCPSNO(I)*ZSNOW(I) + CLHMLT*RHOSNO(I)*ZSNOW(I);                          
            if(HADD<=HCONV)
                ZMELT=HADD/((0.0-TSNOW(I))*HCPSNO(I)+CLHMLT*RHOSNO(I));                           
                RMELTS=ZMELT*RHOSNO(I)/(RHOW*DELT);                                          
                RMELT=RMELTS;
                TRMELT=0.0 ;                                                             
                ZSNOW(I)=ZSNOW(I)-ZMELT;                                                       
                HCPSNO(I)=HCPICE*RHOSNO(I)/RHOICE+HCPW*WSNOW(I)/(RHOW*ZSNOW(I));
                HTCS (I)=HTCS(I)-FI(I)*(QMELT(I)-CLHMLT*RMELT*RHOW);
            else                                                                        
                RMELTS=ZSNOW(I)*RHOSNO(I)/RHOW;                                                 
                RMELT=RMELTS+WSNOW(I)/RHOW;
                HADD=HADD-HCONV;                                                         
                TRMELT=HADD/(HCPW*RMELT);                                                
                RMELT=RMELT/DELT;                                                        
                RMELTS=RMELTS/DELT;
                ZSNOW (I)=0.0;                                                               
                HCPSNO(I)=0.0;
                TSNOW (I)=0.0;                                                               
                WSNOW (I)=0.0;
                HTCS (I)=HTCS(I)-FI(I)*(QMELT(I)-CLHMLT*RMELTS*RHOW-HADD/DELT);
            end                                                                       
            HMFN (I)=HMFN(I)+FI(I)*CLHMLT*RMELTS*RHOW;
            TR   (I)=(R(I)*TR(I)+RMELT*TRMELT)/(R(I)+RMELT);
            R    (I)=R(I)+RMELT;
            QMELT(I)=0.0;
            HTCS(I)=HTCS(I)+FI(I)*HCPSNO(I)*(TSNOW(I)+TFREZ)* ZSNOW(I)/DELT;
        end
        RALB(I)=R(I);
    end
end%100 CONTINUE                                                                   
%C
for I=IL1:IL2%DO 200 I=IL1,IL2
    if(FI(I)>0.)
        if(QMELT(I)>0 && ISAND(I,1)>-4)
            GZERO(I)=GZERO(I)+QMELT(I);
            HTCS (I)=HTCS(I)-FI(I)*QMELT(I);
            HTC(I,1)=HTC(I,1)+FI(I)*QMELT(I);
        end
        RALB(I)=R(I);
    end
end%200 CONTINUE                                                                   
%C                                                                                  
%RETURN                                                                      
%END        
end
