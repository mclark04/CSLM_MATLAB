function [EXPW,RHO]=EQNST(TCEL,H,varargin)
if nargin>2
    warning('Old EQNST');
    TCEL=varargin{3};
    H=varargin{4};
end
% C======================================================================
% C     * AUG  2/16 - M.MACKAY.  	Equation of state moved to subroutine
% C
%       IMPLICIT NONE
% C
% C ----* INPUT FIELDS *------------------------------------------------
% C
%       REAL EXPW,RHO,TCEL
%       REAL GRAV,H,C0,ALPHA,P,BETA,RHO0,T0,T,Z,XI,KAPPA,S
% C
% C SALINITY AND PRESSURE EFFECTS NOT INCLUDED 
% c
Z=single(0.0);
S=single(0.0);
% Cmdm  S=0.3
% Cmdm  S=0.014
% Cmdm  S=0.03
% Cmdm  Z=MIN(10.0,H)
%C-------------------------------------------------------------------------
GRAV=single(9.80616);
C0=single(4.9388E-5);
ALPHA=single(3.3039E-7);
BETA=single(8.2545E-6);
XI=single(8.0477E-1);
KAPPA=single(0.2151);
%Cmdm  T0=3.9816
%Cmdm  T0=3.9839
T0=single(3.98275);
RHO0=999.975 + XI*S;
P=RHO0*GRAV*Z*1.0E-5;
T=TCEL-T0-KAPPA*S;%in source code kappa is lower case but fortran is not case sensitive.
% C======================================================================
% C Farmer and Carmack, 1981
% C
RHO= RHO0*(1. + P*(C0-ALPHA*T) - BETA*T*T );
EXPW = (2.*RHO0*BETA*T + RHO0*P*ALPHA )/RHO;
end