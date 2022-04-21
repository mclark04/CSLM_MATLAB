function [SRH] = ...,ST,SQ,PRESSG,FMASK,ILG,IL1,IL2
    SCREENRH(SRH,ST,SQ,PRESSG,FMASK,ILG,IL1,IL2)%

% C
% C     * DEC 16, 2014 - V.FORTIN.  REPLACE MERGE AND REWORD VARIABLE
% C     *                           DECLARATIONS FOR COMPATABILITY WITH F77.
% C     * APR 30, 2009 - M.LAZARE.
% C
% C     * CALCULATES SCREEN RELATIVE HUMIDITY BASED ON INPUT SCREEN
% C     * TEMPERATURE, SCREEN SPECIFIC HUMIDITY AND SURFACE PRESSURE.
% C     * THE FORMULAE USED HERE ARE CONSISTENT WITH THAT USED ELSEWHERE
% C     * IN THE GCM PHYSICS.
% C
%       IMPLICIT NONE
% C
% C     * OUTPUT FIELD:
% C
%       REAL   SRH(ILG)
% C
% C     * INPUT FIELDS.
% C
%       REAL   ST(ILG),SQ(ILG),PRESSG(ILG),FMASK(ILG)
% C
%       REAL FACTE,EPSLIM,FRACW,ETMP,ESTREF,ESAT,QSW
%       REAL A,B,EPS1,EPS2,T1S,T2S,AI,BI,AW,BW,SLP
%       REAL RW1,RW2,RW3,RI1,RI2,RI3                 
%       REAL ESW,ESI,ESTEFF,TTT,UUU
% C
%       INTEGER ILG,IL,IL1,IL2
%
% C
% C     * COMMON BLOCKS FOR THERMODYNAMIC CONSTANTS.
% C
%       COMMON /EPS /  A,B,EPS1,EPS2    
%       COMMON /HTCP/  T1S,T2S,AI,BI,AW,BW,SLP
% C
% C     * PARAMETERS USED IN NEW SATURATION VAPOUR PRESSURE FORMULATION.
% C
%       COMMON /ESTWI/ RW1,RW2,RW3,RI1,RI2,RI3                 
% C
% C     * COMPUTES THE SATURATION VAPOUR PRESSURE OVER WATER OR ICE.
% C

%In all runs of the FORTRAN code this function always returns 0, likely
%because the physics COMMON's are not defined elsewhere and default to
%zero, so the math defaults to zero.  See code below.
% 
% These functions all return zero so Just overwrite zero for now.
% for i=1:length(IL1:IL2)
%     %SHR(i)=0;
%     SRH(i)=single(0);
% end

%MGC - COMMENT March 24, 2022
%FORTRAN assigns 0 to new variables, but the use of globals (i.e. COMMON)
%is problematic.  I'm sure in the GCM the golabs ESTWI, EPC and HTC are
%well defined, but when I ran this code on my local PC and output the
%variabiles they were all in the default state (=0).  So I'm just going to
%assing them as such here and email MacKay for clarity.

[A,B,EPS1,EPS2,...
  T1S,T2S,AI,BI,AW,BW,SLP,...
  RW1,RW2,RW3,RI1,RI2,RI3] = makeZeros(1,1);
% TTT=single(1);
% UUU=single(1);
% C
% C     * Main function
% C
%
%These functions are called here, MATLAB convention has these defintions
%placed at the end of the subroutine.
% ESW(TTT)        = exp(RW1+RW2/TTT)*TTT^RW3;
% ESI(TTT)        = exp(RI1+RI2/TTT)*TTT^RI3;
% ESTEFF(TTT,UUU) = UUU*ESW(TTT) + (1.-UUU)*ESI(TTT);
%C========================================================================
EPSLIM=single(0.001);
FACTE=1./EPS1-1.;
for IL=IL1:IL2%DO IL=IL1,IL2
    if(FMASK(IL)>0.)    
    %   C
    %   C       * COMPUTE THE FRACTIONAL PROBABILITY OF WATER PHASE      
    %   C       * EXISTING AS A FUNCTION OF TEMPERATURE (FROM ROCKEL,     
    %   C       * RASCHKE AND WEYRES, BEITR. PHYS. ATMOSPH., 1991.)       
    %   C
        if(ST(IL) >= T1S )
            FRACW=1.0;
        else
            FRACW=0.0059+0.9941*exp(-0.003102*(T1S-ST(IL))^2);
        end
        %C
        ETMP=ESTEFF(ST(IL),FRACW);
        ESTREF=0.01*PRESSG(IL)*(1.-EPSLIM)/(1.-EPSLIM*EPS2);
        if ( ETMP < ESTREF )
            ESAT=ETMP;
        else
            ESAT=ESTREF;
        end
        %C
        QSW=EPS1*ESAT/(0.01*PRESSG(IL)-EPS2*ESAT);
        SRH(IL)=MIN(MAX((SQ(IL)*(1.+QSW*FACTE))/(QSW*(1.+SQ(IL)*FACTE)),0.),1.);
    end%ENDIF
end%ENDDO
%C
%RETURN
%END
end

function Y=ESW(TTT)
Y = exp(RW1+RW2/TTT)*TTT^RW3;
end
function Y=ESI(TTT)
Y = exp(RI1+RI2/TTT)*TTT^RI3;
end
function Y=ESTEFF(TTT,UUU) 
Y = UUU*ESW(TTT) + (1.-UUU)*ESI(TTT);
end
