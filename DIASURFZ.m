function [UZ,VZ,TZ,QZ] = ...,NI,U,V,TG,QG,Z0,Z0T,ILMO,ZA,H,UE,FTEMP,FVAP,ZU,ZT,LAT,F,IL1,IL2,JL         
     DIASURFZ (UZ,VZ,TZ,QZ,U,V,TG,QG,Z0,Z0T,ILMO,ZA,...NI,
          H,UE,FTEMP,FVAP,ZU,ZT,LAT,F,IL1,IL2)%,JL

% *Author
% *          Yves Delage  (Aug1990)
% *
% *Revision
% * 001      G. Pellerin(JUN94)
% *          Adaptation to new surface formulation
% * 002      B. Bilodeau (Nov 95) - Replace VK by KARMAN
% * 003      R. Sarrazin (Jan 96) - Prevent problems if zu < za
% * 004      G. Pellerin (Feb 96) - Rewrite stable formulation
% * 005      Y. Delage and B. Bilodeau (Jul 97) - Cleanup
% * 006      Y. Delage (Feb 98) - Addition of HMIN
% * 007      Y. Delage (Sept 00) - Change UE2 by UE
% *                              - Introduce log-linear profile for near-
% *                                 neutral cases
% * 008      D. Verseghy (Nov 02) - Remove unused constant CLM 
% *                                 from common block CLASSD2
% * 009      M. Mackay (Nov 04) - Change all occurrences of Alog
% *                               to log for greater portability.
% * 010      F. SeglenieKs (Mar 05) - Declare LAT as REAL*8 for
% *                                   consistency
% * 011      P.Bartlett (Mar 06) - Set HI to zero for unstable case
% * 012      E.Chan (Nov 06) - Bracket entire subroutine loop with 
% *                            if(F(J)>0.)
% * 013      D.Verseghy (Nov 06) - Convert LAT to regular precision
% * 014      B.Dugas (Jan 09) - "Synchronization" with diasurf2
% *
% *Object
% *          to calculate the diagnostic values of U, V, T, Q
% *          near the surface (ZU and ZT)
% *
% *Arguments
% *
% *          - Output -
% * UZ       U component of the wind at Z=ZU
% * VZ       V component of the wind at Z=ZU
% * TZ       temperature in kelvins at Z=ZT
% * QZ       specific humidity at Z=ZT
% *
% *          - Input -
% * NI       number of points to process
% * U        U component of wind at Z=ZA
% * V        V component of wind at Z=ZA
% * TG       temperature at the surface (Z=0) in Kelvins
% * QG       specific humidity
% * PS       surface pressure at the surface
% * ILMO     inverse of MONIN-OBUKHOV lenth
% * H        height of boundary layer
% * UE       friction velocity
% * Z0       roughness lenth for winds
% * Z0T      roughness lenth for temperature and moisture
% * FTEMP    temperature flux at surface
% * FVAP     vapor flux at surface
% * ZA       heights of first model level above ground
% * ZU       heights for computation of wind components
% * ZT       heights for computation of temperature and moisture
% * LAT      LATITUDE
% * F        Fraction of surface type being studied 
% 
% IMPLICIT NONE
% INTEGER NI,JL
% REAL ZT(NI),ZU(NI)
% REAL UZ(NI),VZ(NI),TZ(NI),QZ(NI),ZA(NI),U(NI),V(NI)
% REAL TG(NI),QG(NI),UE(NI),FTEMP(NI),FVAP(NI)
% REAL ILMO(NI),Z0T(NI),Z0(NI),H(NI),F(NI)
% REAL LAT(NI)
%
% REAL ANG,ANGI,VITS,LZZ0,LZZ0T
% REAL CT,DANG,CM
% REAL XX,XX0,YY,YY0,fh,fm,hi
% REAL RAC3
% INTEGER J,IL1,IL2
% *
% REAL AS,ASX,CI,BS,BETA,FACTN,HMIN,ANGMAX
% COMMON / CLASSD2 / AS,ASX,CI,BS,BETA,FACTN,HMIN,ANGMAX
% REAL DELTA,GRAV,KARMAN,CPD
% COMMON / PHYCON / DELTA,GRAV,KARMAN,CPD

RAC3=single(sqrt(3.));

warning('\nDIASURFZ has not been debugged post traslation from MATLAB');
for J=IL1:IL2%DO J=IL1,IL2
    if(F(J)>0.0)%                           THEN

        LZZ0T=log((ZT(J)+Z0(J))/Z0T(J));
        LZZ0=log(ZU(J)/Z0(J)+1);
        if(ILMO(J)<=0.)% THEN
            % *---------------------------------------------------------------------
            % *                      UNSTABLE CASE
            % *
            HI=0.;
            % *CDIR IEXPAND
            FH=FHI(ZT(J)+Z0(J),Z0T(J),LZZ0T,ILMO(J));
            % *CDIR IEXPAND
            FM=FMI(ZU(J)+Z0(J),Z0 (J),LZZ0 ,ILMO(J));
        else
        % *---------------------------------------------------------------------
        % *                        STABLE CASE
            HI=1/MAX(HMIN,H(J),(ZA(J)+10*Z0(J))*FACTN,FACTN/...
                (4*AS*BETA*ILMO(J)));
            % *CDIR IEXPAND
            FH=BETA*(LZZ0T+MIN( PSI(ZT(J)+Z0(J),HI,ILMO(J))...
                -PSI(Z0T(J),HI,ILMO(J)),...
                ASX*ILMO(J)*(ZT(J)+Z0(J)-Z0T(J))));
            % *CDIR IEXPAND
            FM=LZZ0+MIN(PSI(zu(J)+Z0(J),HI,ILMO(J))...
                -PSI(Z0(J),HI,ILMO(J)),...
                ASX*ILMO(J)*ZU(J));
        end
        % *---------------------------------------------------------------------
        CT=KARMAN/FH;
        CM=KARMAN/FM;
        TZ(J)=TZ(J)+F(J)*(TG(J)-FTEMP(J)/(CT*UE(J))-GRAV/CPD*ZT(J));
        QZ(J)=QZ(J)+F(J)*(QG(J)-FVAP(J)/(CT*UE(J)));
        VITS=UE(J)/CM;

        % * CALCULATE WIND DIRECTION CHANGE FROM TOP OF SURFACE LAYER
        DANG= (ZA(J)-ZU(J))*HI*ANGMAX*sin(LAT(J));
        ANGI=atan2(V(J),SIGN(abs(U(J))+1.e-05,U(J)));%SIGN should have 2 inputs!        
        if(ILMO(J)>0.)%    THEN
            ANG=ANGI+DANG;
        else
            ANG=ANGI;
        end

        UZ(J)=UZ(J)+F(J)*VITS*cos(ANG);
        VZ(J)=VZ(J)+F(J)*VITS*sin(ANG);

    end 
end%ENDDO

% RETURN
% CONTAINS
end

% C   The following code is taken from the RPN/CMC physics library file
% C   /usr/local/env/armnlib/modeles/PHY_shared/ops/v_4.5/RCS/stabfunc2.cdk,v
% 
% C   Internal function FMI
% C   Stability function for momentum in the unstable regime (ILMO<0)
% c   Reference: Delage Y. and Girard C. BLM 58 (19-31) Eq. 19
% c
% REAL FUNCTION FMI(Z2,Z02,LZZ02,ILMO2,X,X0)
function [FMIout]=FMI(Z2,Z02,LZZ02,ILMO2)%,X,X0
% implicit none
% C
% REAL, INTENT(IN ) :: Z2,Z02,LZZ02,ILMO2
% REAL, INTENT(OUT) :: X,X0
% c
X =(1-CI*Z2 *BETA*ILMO2)^(0.16666666);
X0=(1-CI*Z02*BETA*ILMO2)^(0.16666666);
FMIout=LZZ02+log((X0+1)^2*sqrt(X0^2-X0+1)*(X0^2+X0+1)^1.5...
    /((X+1)^2*sqrt(X^2-X+1)*(X^2+X+1)^1.5))...
    +RAC3*atan(RAC3*((X^2-1)*X0-(X0^2-1)*X)/...
    ((X0^2-1)*(X^2-1)+3*X*X0));
% c
% RETURN
end%END FUNCTION FMI

% c
% C   Internal function FHI
% C   Stability function for heat and moisture in the unstable regime (ILMO<0)
% c   Reference: Delage Y. and Girard C. BLM 58 (19-31) Eq. 17
% c
% REAL FUNCTION FHI(Z2,Z0T2,LZZ0T2,ILMO2,Y,Y0)
function [FHIout]=FHI(Z2,Z0T2,LZZ0T2,ILMO2)%,Y,Y0
% implicit none
% C
% REAL, INTENT(IN ) :: Z2,Z0T2,LZZ0T2,ILMO2
% REAL, INTENT(OUT) :: Y,Y0
% c
Y =(1-CI*Z2  *BETA*ILMO2)^(0.33333333);
Y0=(1-CI*Z0T2*BETA*ILMO2)^(0.33333333);
FHIout=BETA*(LZZ0T2+1.5*log((Y0^2+Y0+1)/(Y^2+Y+1))+RAC3*...
    atan(RAC3*2*(Y-Y0)/((2*Y0+1)*(2*Y+1)+3)));
% c
%RETURN
end%END FUNCTION FHI

% C
% C   Internal function PSI
% C   Stability function for momentum in the stable regime (unsl>0)
% c   Reference :  Y. Delage, BLM, 82 (p23-48) (Eqs.33-37)
% c
% REAL FUNCTION PSI(Z2,HI2,ILMO2)
function [PSIout]=PSI(Z2,HI2,ILMO2)
% % implicit none
% C
% REAL a,b,c,d
% REAL, INTENT(IN ) :: ILMO2,Z2,HI2
% c
D = 4*AS*BETA*ILMO2;
C = D*HI2 - HI2^2;
B = D - 2*HI2;
A = sqrt(1 + B*Z2 - C*Z2^2);
PSIout = 0.5 * (A-Z2*HI2-log(1+B*Z2*0.5+A)-...
	B/(2*sqrt(C))*asin((B-2*C*Z2)/D));
% c
% RETURN
end%END FUNCTION PSI
%END