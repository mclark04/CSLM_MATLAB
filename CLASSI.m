function [VPD,TADP,PADRY,RHOAIR,RHOSNI,RPCP,TRPCP,SPCP,TSPCP,TA,QA,PCPR,...
    RRATE,SRATE,PRESSG,IPCP,NL,IL1,IL2] = ...
   CLASSI(VPD,TADP,PADRY,RHOAIR,RHOSNI,RPCP,TRPCP,...
    SPCP,TSPCP,TA,QA,PCPR,RRATE,SRATE,PRESSG,IPCP,...
    NL,IL1,IL2,TFREZ,RGAS,RGASV,RHOW)
%CLASSI
%[VPDROW,TADPROW,PADRROW,RHOAROW,RHSIROW,RPCPROW,TRPCROW,SPCPROW,TSPCROW,TAROW,QAROW,PREROW,RPREROW,SPREROW,PRESROW,IPCP,NLAT,~,NLTEST]=CLASSI(VPDROW,TADPROW,PADRROW,RHOAROW,RHSIROW,RPCPROW,TRPCROW,SPCPROW,TSPCROW,TAROW,QAROW,PREROW,RPREROW,SPREROW,PRESROW,IPCP,NLAT,1,NLTEST);
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
% C
% 
% C     * NOV 17/11 - M.LAZARE.   REMOVE CALCULATION OF PCPR
% C     *                         FOR IPCP=4 (REDUNDANT SINCE
% C     *                         PCPR MUST BE PASSED IN FOR
% C     *                         IF CONDITION ON LINE 100).
% C     * NOV 22/06 - P.BARTLETT. CALCULATE PCPR IF IPCP=4.
% C     * JUN 06/06 - V.FORTIN.   ADD OPTION FOR PASSING IN
% C     *                         RAINFALL AND SNOWFALL RATES
% C     *                         CALCULATED BY ATMOSPHERIC MODEL.
% C     * NOV 03/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND;
% C     *                         MOVE CLOUD FRACTION CALCULATION
% C     *                         BACK INTO DRIVER.
% C     * SEP 04/03 - D.VERSEGHY. NEW LOWER LIMIT ON PRECIPITATION
% C     *                         RATE.
% C     * AUG 09/02 - D.VERSEGHY. MOVE CALCULATION OF SOME
% C     *                         ATMOSPHERIC VARIABLES HERE
% C     *                         PRIOR TO GATHERING.
% C     * JUL 26/02 - R.BROWN/S.FASSNACHT/D.VERSEGHY. PROVIDE 
% C     *                         ALTERNATE METHODS OF ESTIMATING 
% C     *                         RAINFALL/SNOWFALL PARTITIONING.
% C     * JUN 27/02 - D.VERSEGHY. ESTIMATE FRACTIONAL CLOUD COVER
% C     *                         AND RAINFALL/SNOWFALL RATES
% C     *                         IF NECESSARY.
% C
%       IMPLICIT NONE
% C
% C     * INTEGER CONSTANTS.
% C
%       INTEGER IPCP,NL,IL1,IL2,I
% C
% C     * OUTPUT ARRAYS.
% C
%       REAL VPD   (NL),  TADP  (NL),  PADRY (NL),  RHOAIR(NL),
%      1     RHOSNI(NL),  RPCP  (NL),  TRPCP (NL),  
%      2     SPCP  (NL),  TSPCP (NL)  
% C
% C     * INPUT ARRAYS.
% C
%       REAL TA    (NL),  QA    (NL),  PRESSG(NL), 
%      1     PCPR  (NL),  RRATE (NL),  SRATE (NL)
% C
% C     * WORK ARRAYS.
% C
%       REAL PHASE (NL)
% C
% C     * TEMPORARY VARIABLES.
% C
%       REAL EA,CA,CB,EASAT,CONST
% C
% C     * COMMON BLOCK PARAMETERS.
% C
%       REAL DELT,TFREZ,RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN,HCPW,HCPICE,
%      1     HCPSOL,HCPOM,HCPSND,HCPCLY,SPHW,SPHICE,SPHVEG,SPHAIR,
%      2     RHOW,RHOICE,TCGLAC,CLHMLT,CLHVAP
%  
%       COMMON /CLASS1/ DELT,TFREZ
%       COMMON /CLASS2/ RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN
%       COMMON /CLASS4/ HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,
%      1                SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,
%      2                TCGLAC,CLHMLT,CLHVAP
% C----------------------------------------------------------------
% C
% C     * CALCULATION OF ATMOSPHERIC INPUT VARIABLES.
% C
for I=IL1:IL2
    EA=QA(I)*PRESSG(I)/(0.622+0.378*QA(I));
    if(TA(I)>=TFREZ)
        CA=17.269;
        CB=35.86;
    else
        CA=21.874;
        CB=7.66;
    end
    EASAT=611.0*exp(CA*(TA(I)-TFREZ)/(TA(I)-CB));
    VPD(I)=MAX(0.0,(EASAT-EA)/100.0);
    PADRY(I)=PRESSG(I)-EA;
    RHOAIR(I)=PADRY(I)/(RGAS*TA(I))+EA/(RGASV*TA(I));
    CONST=log(EA/611.0);
    TADP(I)=(CB*CONST-CA*TFREZ)/(CONST-CA);
    % C
    % C     * DENSITY OF FRESH SNOW.
    % C
    if(TA(I)<=TFREZ)
        RHOSNI(I)=67.92+51.25*exp((TA(I)-TFREZ)/2.59);
    else
        RHOSNI(I)=MIN((119.17+20.0*(TA(I)-TFREZ)),200.0);
    end
    % C
    % C     * PRECIPITATION PARTITIONING BETWEEN RAIN AND SNOW.
    % C
    RPCP (I)=0.0;
    TRPCP(I)=0.0;
    SPCP (I)=0.0;
    TSPCP(I)=0.0;
    if(PCPR(I)>1.0E-8)
        if(IPCP==1)
            if(TA(I)>TFREZ)
                RPCP (I)=PCPR(I)/RHOW;
                TRPCP(I)=MAX((TA(I)-TFREZ),0.0);
            else
                SPCP (I)=PCPR(I)/RHOSNI(I);
                TSPCP(I)=MIN((TA(I)-TFREZ),0.0);
            end
        elseif(IPCP==2)
            if(TA(I)<=TFREZ)
                PHASE(I)=1.0;
            elseif(TA(I)>=(TFREZ+2.0))
                PHASE(I)=0.0;
            else
                PHASE(I)=1.0-0.5*(TA(I)-TFREZ);
            end
            RPCP(I)=(1.0-PHASE(I))*PCPR(I)/RHOW;
            if(RPCP(I)>0.0)
                TRPCP(I)=MAX((TA(I)-TFREZ),0.0);
            end
            SPCP(I)=PHASE(I)*PCPR(I)/RHOSNI(I);
            if(SPCP(I)>0.0)
                TSPCP(I)=MIN((TA(I)-TFREZ),0.0);
            end
        elseif(IPCP==3)
            if(TA(I)<=TFREZ)
                PHASE(I)=1.0;
            elseif(TA(I)>=(TFREZ+6.0))
                PHASE(I)=0.0;
            else
                PHASE(I)=(0.0202*(TA(I)-TFREZ)^6-0.3660*(TA(I)-TFREZ)^5+2.0399*(TA(I)-TFREZ)^4-1.5089*(TA(I)-TFREZ)^3-15.038*(TA(I)-TFREZ)^2+4.6664*(TA(I)-TFREZ)+100.0)/100.0;
                PHASE(I)=MAX(0.0,MIN(1.0,PHASE(I)));
            end
            RPCP(I)=(1.0-PHASE(I))*PCPR(I)/RHOW;
            if(RPCP(I)>0.0)
                TRPCP(I)=MAX((TA(I)-TFREZ),0.0);
            end
            SPCP(I)=PHASE(I)*PCPR(I)/RHOSNI(I);
            if(SPCP(I)>0.0)
                TSPCP(I)=MIN((TA(I)-TFREZ),0.0);
            end
        elseif(IPCP==4)
            RPCP(I)=RRATE(I)/RHOW;
            if(RPCP(I)>0.0)
                TRPCP(I)=MAX((TA(I)-TFREZ),0.0);
                SPCP(I)=SRATE(I)/RHOSNI(I);
            end
            if(SPCP(I)>0.0)
                TSPCP(I)=MIN((TA(I)-TFREZ),0.0);
            end
        end
    end
end
end