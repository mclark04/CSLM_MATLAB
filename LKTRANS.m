function [CQ1A,CQ1B,CQ2A,CQ2B,CQ3A,CQ3B,CQ1BI,CQ2BI,CQ3BI] = ...BLAK,IL1,IL2,ILG,
    LKTRANS (BLAK,IL1,IL2)%CQ1A,CQ1B,CQ2A,CQ2B,CQ3A,CQ3B,ILG,,CQ1BI,CQ2BI,CQ3BI

% C======================================================================
% C     * DEC  7/07 - M.MACKAY.  	COMPUTES LIGHT EXTINCTION COEFFICIENTS
% C     *                         THIS ROUTINE KEPT FOR FUTURE DYNAMIC
% C				CHANGES TO EXTINCTION
% C
%       IMPLICIT NONE
% C
% C ----* INPUT FIELDS *------------------------------------------------
% C
%       REAL,DIMENSION(ILG) :: BLAK, CQ1A,CQ1B,CQ2A,CQ2B,CQ3A,CQ3B,
%      >                       CQ1BI,CQ2BI,CQ3BI
%       INTEGER IL1,IL2,ILG
% C
% C ----* LOCAL VARIABLES *------------------------------------------------
% C
%       INTEGER I
%
% C======================================================================
% C                CQ1A,    CQ1B,    CQ2A,    CQ2B,    CQ3A,    CQ3B 
% C Rayner,1980    0.54,    0.561,   0.30,    6.89,    0.16,    69.0 
% C======================================================================
for I=IL1:IL2%DO 100 I=IL1,IL2
% C FIXED WATER VALUES
    CQ1A(I)=single(0.5817);
    CQ2A(I)=single(0.4183);
    CQ2B(I)=single(6.89);
    CQ3A(I)=single(0.0);
    CQ3B(I)=single(69.0);
    % C FIXED ICE VALUES (from Patterson and Hamblin, 1988, L&O)
    CQ1BI(I)=single(1.5);
    % Cmdm  CQ1BI(I)=3.75    !test - high extinction for snow ice
    CQ2BI(I)=single(20.0);
    CQ3BI(I)=single(69.0);

    % C======================================================================
    % C CQ1B NOW READ IN .INI FILE
    % C----------------------------------------------------------------------
    CQ1B(I)=BLAK(I);

end%100   CONTINUE

%       RETURN
%       END