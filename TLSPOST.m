function [GSNOW,TSNOW,WSNOW,RHOSNO,QMELTS,GZERO,TSNBOT,...
        HTCS,HMFN,QFN,EVAPS,RPCN,TRPCN,SPCN,TSPCN,...
        HCPSNO...GCONSTS,GCOEFFS,T0,ZSNOW,TCSNOW,QTRANS,
        ...RPCP,TRPCP,SPCP,TSPCP,TZEROS,RHOSNI,
        ] = ...FLS,DELSKIN,ILG,IL1,IL2,JL
    TLSPOST(GSNOW,TSNOW,WSNOW,RHOSNO,QMELTS,GZERO,TSNBOT,...
        HTCS,HMFN,EVAPS,...QFN,RPCN,TRPCN,SPCN,TSPCN,
        T0,ZSNOW,TCSNOW,HCPSNO,...GCONSTS,GCOEFFS,
        RPCP,TRPCP,SPCP,TSPCP,TZEROS,RHOSNI,...QTRANS,
        FLS,DELSKIN,IL1,IL2,...ILG,JL,
        RHOW,DELT,TFREZ,CLHMLT,HCPICE,RHOICE,HCPW)%added these (and below) to pass in on a call
        
% C
% C     * JAN 30/18 - M.LAZARE.   LAST ELEMENT IN CALL, "N", REMOVED
% C     *                         BECAUSE NOT IN CALL STATEMENT FROM 
% C     *                         ROUTINE CLASSL, AND IS NOT USED.
% C     * SEP 01/15 - D.VERSEGHY. LAKE SNOW TEMPERATURE AND HEAT FLUX
% C     *                         CALCULATIONS (BASED ON CLASS SUBROUTINE
% C     *                         TSPOST); NET SURFACE WATER FLUX TERMS
% C     *                         (BASED ON CLASS SUBROUTINE WPREP).
% C
%       IMPLICIT NONE
% C                                                                                 
% C     * INTEGER CONSTANTS.
% C
%       INTEGER ILG,IL1,IL2,JL,I,J
% C
% C     * OUTPUT ARRAYS.
% C
%       REAL GZERO (ILG),    TSNBOT(ILG),
%      1     RPCN  (ILG),    TRPCN (ILG),    SPCN  (ILG),    TSPCN (ILG)
% C
% C     * INPUT/OUTPUT ARRAYS.
% C
%       REAL GSNOW (ILG),    TSNOW (ILG),    WSNOW (ILG),    RHOSNO(ILG),
%      1     QMELTS(ILG),    HTCS  (ILG),    HMFN  (ILG),    QFN   (ILG),
%      2     EVAPS (ILG)
% C
% C     * INPUT ARRAYS.
% C
%       REAL T0    (ILG),    ZSNOW (ILG),    TCSNOW (ILG),   HCPSNO(ILG), 
%      1     QTRANS(ILG),    GCONSTS(ILG),   GCOEFFS(ILG),   
%      2     RPCP  (ILG),    TRPCP  (ILG),   SPCP   (ILG),   TSPCP (ILG),
%      3     TZEROS(ILG),    RHOSNI (ILG),   FLS    (ILG)
% C
%       REAL DELSKIN
% C
% C     * TEMPORARY VARIABLES.
% C
%       REAL RADD  (ILG),    SADD  (ILG)
%[RADD, SADD] = makeBlanks(length(IL1:IL2),1);%ILG instead of length
[RADD, SADD] = makeZeros(length(IL1:IL2),1);%ILG instead of length
% C
%       REAL HADD,HCONV,WFREZ
% C
% C     * COMMON BLOCK PARAMETERS.
% C
%       REAL DELT,TFREZ,TCW,TCICE,TCSAND,TCCLAY,TCOM,TCDRYS,
%      1     RHOSOL,RHOOM,HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,
%      2     SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,TCGLAC,CLHMLT,CLHVAP
% C
%       COMMON /CLASS1/ DELT,TFREZ                                       
%       COMMON /CLASS3/ TCW,TCICE,TCSAND,TCCLAY,TCOM,TCDRYS,
%      1                RHOSOL,RHOOM
%       COMMON /CLASS4/ HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,
%      1                SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,
%      2                TCGLAC,CLHMLT,CLHVAP
% C-----------------------------------------------------------------------
% C
for I=IL1:IL2%DO 100 I=IL1,IL2
    if(FLS(I)>0)
        TSNBOT(I)=(ZSNOW(I)*TSNOW(I)+DELSKIN*T0(I))/(ZSNOW(I)+DELSKIN);
        GZERO(I)=-2.0*TCSNOW(I)*(TSNBOT(I)-TSNOW(I))/ZSNOW(I);
%C     1                 +TCICE*(T0(I)-TSNBOT(I))/DELSKIN)
        if(QMELTS(I)<0.)
            GSNOW(I)=GSNOW(I)+QMELTS(I);                           
            QMELTS(I)=0;                                         
        end
        TSNOW(I)=TSNOW(I)+(GSNOW(I)-GZERO(I))*DELT/(HCPSNO(I)*ZSNOW(I))-TFREZ;
        if(TSNOW(I)>0)
            QMELTS(I)=QMELTS(I)+TSNOW(I)*HCPSNO(I)*ZSNOW(I)/DELT;
            GSNOW(I)=GSNOW(I)-TSNOW(I)*HCPSNO(I)*ZSNOW(I)/DELT;
            TSNOW(I)=0;
        end
%         C              GZERO(I)=GZERO(I)+QMELTS(I)
%         C              QMELTS(I)=0.0
    end
end%100 CONTINUE
%C 
for I=IL1:IL2%      DO 200 I=IL1,IL2
    if(FLS(I)>0 && TSNOW(I)<0 && WSNOW(I)>0)    
        HTCS(I)=HTCS(I)-FLS(I)*HCPSNO(I)*(TSNOW(I)+TFREZ)*ZSNOW(I)/DELT;
        HADD=-TSNOW(I)*HCPSNO(I)*ZSNOW(I);
        HCONV=CLHMLT*WSNOW(I);
        if(HADD<=HCONV)
            WFREZ=HADD/CLHMLT;
            HADD=single(0.0);
            WSNOW(I)=MAX(0.0,WSNOW(I)-WFREZ);
            TSNOW(I)=0.0;
            RHOSNO(I)=RHOSNO(I)+WFREZ/ZSNOW(I);
            HCPSNO(I)=HCPICE*RHOSNO(I)/RHOICE+HCPW*WSNOW(I)/(RHOW*ZSNOW(I));
        else
            HADD=HADD-HCONV;
            WFREZ=WSNOW(I);
            WSNOW(I)=0.0;
            RHOSNO(I)=RHOSNO(I)+WFREZ/ZSNOW(I);
            HCPSNO(I)=HCPICE*RHOSNO(I)/RHOICE;
            TSNOW(I)=-HADD/(HCPSNO(I)*ZSNOW(I));
        end
        HMFN(I)=HMFN(I)-FLS(I)*CLHMLT*WFREZ/DELT;
        HTCS(I)=HTCS(I)-FLS(I)*CLHMLT*WFREZ/DELT;
        HTCS(I)=HTCS(I)+FLS(I)*HCPSNO(I)*(TSNOW(I)+TFREZ)*ZSNOW(I)/DELT;
    end
end %200 CONTINUE
%C
for I=IL1:IL2%DO 300 I=IL1,IL2
    QFN(I)=FLS(I)*EVAPS(I)*RHOW;
    if (SPCP(I)>0 || EVAPS(I)<0)
        SADD(I)=SPCP(I)-EVAPS(I)*RHOW/RHOSNI(I);
        if(abs(SADD(I)) < 1.0E-12) 
            SADD(I)=single(0.0);
        end
        if(SADD(I) > 0.0)
            SPCN (I)=SADD(I);                            
            if(SPCP(I) > 0.0)
                TSPCN(I)=TSPCP(I);
            else                
                TSPCN(I)=MIN((TZEROS(I)-TFREZ),0.0);
            end
            EVAPS(I)=single(0.0);
        else                                            
            EVAPS(I)=-SADD(I)*RHOSNI(I)/RHOW;
            SPCN (I)=single(0.0);                               
            TSPCN(I)=single(0.0);                               
        end 
    else                                                
        SPCN (I)=single(0.0);                                   
        TSPCN(I)=single(0.0);                                  
    end
    if(RPCP(I) > 0)
        RADD(I)=RPCP(I)-EVAPS(I);                      
        if(abs(RADD(I)) < 1.0E-12) 
            RADD(I)=single(0.0);
        end
        if (RADD(I) > 0.)
            RPCN (I)=RADD(I);                         
            TRPCN(I)=TRPCP(I);
            EVAPS(I)=single(0.0);                            
        else                                          
            EVAPS(I)=-RADD(I);                       
            RPCN (I)=single(0.0);
            TRPCN(I)=single(0.0);
        end                                       
    else                                             
        RPCN (I)=single(0.0);                               
        TRPCN(I)=single(0.0);                               
    end                                             
end%300 CONTINUE
% C
%       RETURN                    
%       END
end