function [T0,TLAK] = ...LKICEH,RHOIW,NLAK,NLAKMAX,ILG,IL1,IL2,IYEAR,IHOUR,IDAY,IMIN
    FREECONV(LKICEH,T0,TLAK,RHOIW,NLAK,IL1,IL2,...NLAKMAX,ILG,,IYEAR,IHOUR,IDAY,IMIN
    TFREZ,DELSKIN,DELZLK)%added because of global issues

% C=======================================================================
%       IMPLICIT NONE
% C
% C ----* LAKE MODEL VARIABLES *----------------------------------------
% C
%       INTEGER NLAKMAX
%       INTEGER,DIMENSION(ILG) :: NLAK
%       REAL,DIMENSION(ILG) :: T0, LKICEH
%       REAL,DIMENSION(ILG,NLAKMAX) :: TLAK
% C
% C ----* INPUT *-------------------------------------------
% C
%       INTEGER ILG,IL1,IL2,IYEAR,IDAY,IHOUR,IMIN
%       REAL RHOIW
% C
% C ----* CLASS COMMON BLOCKS *------------------------------------------
% C
%       REAL DELT,TFREZ
%       REAL TKECN,TKECF,TKECE,TKECS,TKECL,HDPTHMIN,
%      1     TKEMIN,DELMAX,DELMIN,EMSW,DELZLK,DELSKIN,DHMAX,DUMAX
%       COMMON /CLASS1/ DELT,TFREZ                                       
%       COMMON /LAKECON/ TKECN,TKECF,TKECE,TKECS,HDPTHMIN,        
%      2                 TKEMIN,DELMAX,DELMIN,EMSW,DELZLK,DELSKIN,
%      3                 DHMAX,TKECL,DUMAX
% C
% C ----* LOCAL VARIABLES *---------------------------------------------
% C
%       INTEGER I,J,K,NMIX
%       REAL ZTOP,ZBOT,TTEST,RHO1,RHO2,TC1,TC2,TBAR,XXX,ICEBOT
% C=======================================================================
% C

for I=IL1:IL2%DO 100 I=IL1,IL2
    if (LKICEH(I) <= 0.0)
        TC1=T0(I)-TFREZ;
        TC2=TLAK(I,1)-TFREZ;
        [~,RHO1]=EQNST(TC1,0.05);
        [~,RHO2]=EQNST(TC2,0.5);
        if (RHO1 > RHO2)
            TBAR=((DELSKIN*RHO1*T0(I))+(DELZLK*RHO2*TLAK(I,1)))/((DELSKIN*RHO1)+(DELZLK*RHO2));
            T0(I)=TBAR;
            TLAK(I,1)=TBAR;
        end
    end
    ICEBOT=RHOIW*LKICEH(I);

    NMIX=1;
    for J=1:NLAK(I)-1%DO 420, J=1,NLAK(I)-1
        ZTOP=DELSKIN + (J-1)*DELZLK;
        ZBOT=ZTOP+DELZLK;
        if (ICEBOT <= ZTOP)
            TC1=TLAK(I,J)-TFREZ;
            TC2=TLAK(I,J+1)-TFREZ;
            %Cmdm       TTEST=(TC1-3.9816)*(TC2-3.9816)
            TTEST=(TC1-3.98275)*(TC2-3.98275);
            [~,RHO1]=EQNST(TC1,ZBOT);
            [~,RHO2]=EQNST(TC2,ZBOT+DELZLK);
            % C--------- MIX LAYERS if RHO1>RHO2 OR TEMPERATURES SPAN 
            % C--------- T_MAXDENSITY=3.9816 C.
            if ((RHO1 > RHO2) || (TTEST < 0.0))
                TBAR=((NMIX*RHO1*TLAK(I,J))+(RHO2*TLAK(I,J+1)))/(NMIX*RHO1+RHO2);
                for K=J-NMIX+1:J+1%DO 430, K=J-NMIX+1,J+1
                    TLAK(I,K)=TBAR;
                end%430         CONTINUE
                NMIX=NMIX+1;
                %C           WRITE(6,6666) "static instability removed under ice:",
                %C    >                IYEAR,IDAY,IHOUR,IMIN,j*DELZLK
            else
                NMIX=1;
            end
        end
    end%420     CONTINUE
end%100   CONTINUE
%^6666  FORMAT(A37,4I5,F5.1)
%RETURN
%END
end