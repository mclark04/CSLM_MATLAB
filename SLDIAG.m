function [SUT,SVT,STT,SQT] = ...,CDM,CDH,UA,VA,TA,QA,T0,Q0,Z0M,Z0E,F,ZA,ZU,ZT,ILG,IL1,IL2,JL
       SLDIAG (SUT,SVT,STT,SQT,CDM,CDH,UA,VA,TA,QA,T0,Q0,...      
          Z0M,Z0E,F,ZA,ZU,ZT,IL1,IL2)%,ILG,JL
 
% C     * JUN 23/14 - M.LAZARE.   New version for gcm18+:                 
% C     *                         - Bugfix to calculation of              
% C     *                           screen temperature and                
% C     *                           screen specific humidity.             
% C     *                         - Accumulation removed (now             
% C     *                           done in classt/oiflux11) so           
% C     *                           that a screen relative humidity       
% C     *                           can be calculated. Therefore,        
% C     *                           "instantaneous" fields are            
% C     *                           calculated and passed out             
% C     *                           instead.                              
% C     * OCT 17/11 - D.VERSEGHY. ADD CODE TO CIRCUMVENT SPECIAL          
% C     *                         CASE WHERE TA~T0 OR QA~QO, THUS         
% C     *                         AVOIDING A DIVIDE BY ZERO.              
% C     * NOV 04/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.            
% C     * JUL 19/96 - Y. DELAGE.                                          
% C     * CALCULATES NEAR SURFACE OUTPUT VARIABLES                        
% c     * OUTPUT FIELDS ARE:                                              
% C     *   SUT  : U COMPONENT OF THE WIND AT ZU                          
% C     *   SVT  : V COMPONENT OF THE WIND AT ZU                          
% C     *   STT  : TEMPERATURE AT ZT                                      
% C     *   SQT  : SPECIFIC HUMIDITY AT ZT                                
% C     * INPUT FIELDS ARE:                                               
% C     *   CDM : DRAG COEFFICIENT                                        
% C     *   CDH : TRASFER COEFFICIENT FOR HEAT AND MOISTURE               
% C     *   UA  : U COMPONENT OF THE WIND AT ZA                           
% C     *   VA  : V COMPONENT OF THE WIND AT ZA                           
% C     *   TA  : POTENTIAL TEMPERATURE AT ZA                             
% c     *   T0  : TEMPERATURE AT BOTTOM OF SURFACE LAYER                  
% C     *   Q0  : SPECIFIC HUMIDITY AT BOTTOM OF SURFACE LAYER            
% C     *   Z0M : ROUGHNESS LENGTH FOR MOMENTUM                           
% C     *   Z0E : ROUGHNESS LENGTH FOR HEAT AND MOISTURE                  
% C     *   F   : FRACTION OF GRID POINT AFFECTED BY PROCESS              
% C     *   ZA  : TOP OF SURFACE LAYER                                    
% C     *   ZU  : HEIGHT OF OUTPUT WIND                                   
% C     *   ZT  : HEIGHT OF OUTPUT TEMPERATURE AND HUMIDITY               
% C     *   ILG : NUMBER OF POINTS TO BE TREATED                          
% C                                                                       
% IMPLICIT NONE                                                     
% c                                                                       
% INTEGER ILG,IL1,IL2,JL,I                                          
% c                                                                       
% REAL  SUT(ILG), SVT(ILG), STT(ILG), SQT(ILG), CDM(ILG), CDH(ILG), 
% 1       UA(ILG),  VA(ILG),  TA(ILG),  QA(ILG), Z0M(ILG), Z0E(ILG), 
% 2        F(ILG),  T0(ILG),  Q0(ILG),  ZA(ILG),  ZU(ILG),  ZT(ILG)  
% c                                                                       
% REAL PR,WSPD,CM,US,TS,QS,L,UVA,RATIO,UVU,TTA,CE                   
% c                                                                       
% REAL RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN                              
% c                                                                       
% REAL PSM,PSE,Y,PIM,PIE,X                                          
% c                                                                       
% COMMON /CLASS2/ RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN                   
%
%These functions need to be defined as functions in matlab, so they are placed below
%the main code of the subrute so they can only be called by this subrutine
% c     * STABILITY FUNCTIONS FOR THE STABLE CASE                         
% PSM(X)= -X -.667*(X-5/.35)*exp(-.35*X);
% PSE(X)= -(1+.667*X)^1.5 -.667*(X-5/.35)*exp(-.35*X);
% % c     * STABILITY FUNCTIONS FOR THE UNSTABLE CASE                       
% Y(X)=(1-16*X)^.25 ;                                               
% PIM(X)= log((1+X)^2*(1+X^2)) -2*atan(X);                         
% PIE(X)= 2*log(1+X^2);                                             
PR=single(1.0);                                                            
for I=IL1:IL2%DO 100 I=IL1,IL2                                                  
    if(F(I)>0.)
        % C     * CALCULATION OF SURFACE FLUXES AND MONIN-OBUKHOV LENGTH          
        WSPD=MAX(VMIN,sqrt(UA(I)^2+VA(I)^2));                          
        CM=sqrt(CDM(I));                                                 
        US=CM*WSPD ;                                                     
        if(ABS(TA(I)-T0(I))<0.01)
            TS=-0.01*CDH(I)/CM;                                          
        else                                                            
            TS=CDH(I)*(TA(I)-T0(I))/CM;                                  
        end                                                           
        if(ABS(QA(I)-Q0(I))<1.0E-7)
            QS=-1.0E-7*CDH(I)/CM;                                        
        else                                                            
            QS=CDH(I)*(QA(I)-Q0(I))/CM;                                  
        end                                                           
        L=TA(I)*US^2/(VKC*GRAV*(TS*(1+.61*QA(I))+.61*TA(I)*QS));        
        % C     * CALCULATE CORRECTION FACTORS TO TAKE INTO ACCOUNT THE APPROXIMATIONS
        % C     * IN DRCOEF                                                          
        if(L>0.) %                                                THEN
        % C     * STABLE CASE                                                     
            UVA=US/VKC*(log(ZA(I)/Z0M(I))-PSM(ZA(I)/L)+PSM(Z0M(I)/L));      
            RATIO=WSPD/UVA;                                                
            UVU=US/VKC*(log((ZU(I)+Z0M(I))/Z0M(I))-PSM((ZU(I)+Z0M(I))/L)...   
            	+PSM(Z0M(I)/L))*RATIO;                                      
            TTA=T0(I)+TS/VKC*PR*(log(ZA(I)/Z0E(I))-PSE(ZA(I)/L)+...           
            	PSE(Z0E(I)/L));                                          
            RATIO=(TA(I)-T0(I))/SIGN(MAX(abs(TTA-T0(I)),1.E-4),TTA-T0(I));  
            CE=(log((ZT(I)+Z0M(I))/Z0E(I))-PSE((ZT(I)+Z0M(I))/L)...           
            	+PSE(Z0E(I)/L))*RATIO*PR/VKC;                                
        else                                                            
        % C     * UNSTABLE CASE                                                   
            UVA=US/VKC*(log(ZA(I)/Z0M(I))-PIM(Y_fun(ZA(I)/L))+PIM(Y_fun(Z0M(I)/L)));
            RATIO=WSPD/UVA;                                                
            UVU=US/VKC*(log((ZU(I)+Z0M(I))/Z0M(I))-PIM(Y_fun((ZU(I)+Z0M(I))/L))...
            	+PIM(Y_fun(Z0M(I)/L)))*RATIO;                                  
            TTA=T0(I)+TS/VKC*PR*(log(ZA(I)/Z0E(I))-PIE(Y_fun(ZA(I)/L))+...        
            	PIE(Y_fun(Z0E(I)/L)));                                      
            RATIO=(TA(I)-T0(I))/SIGN(MAX(abs(TTA-T0(I)),1.E-4),TTA-T0(I));  
            CE=(log((ZT(I)+Z0M(I))/Z0E(I))-PIE(Y_fun((ZT(I)+Z0M(I))/L))...     
            	+PIE(Y_fun(Z0E(I)/L)))*RATIO*PR/VKC;                             
        end                                                           
        % C                                                                       
        SUT(I)=UVU*UA(I)/WSPD;                                           
        SVT(I)=UVU*VA(I)/WSPD;                                           
        STT(I)=T0(I)+TS*CE;                                              
        SQT(I)=Q0(I)+QS*CE;                                              
    end                                                             
end%100 CONTINUE                                                          
% RETURN                                                            
% END 
end

% c     * STABILITY FUNCTIONS FOR THE STABLE CASE     
function [Y]=PSM(X)
Y = -X -.667*(X-5/.35)*exp(-.35*X);
end
function [Y]= PSE(X)
Y = -(1+.667*X)^1.5 -.667*(X-5/.35)*exp(-.35*X);
end
% c     * STABILITY FUNCTIONS FOR THE UNSTABLE CASE                       
function Y=Y_fun(X)
Y=(1-16*X)^.25 ;                                               
end
function Y=PIM(X)
Y= log((1+X)^2*(1+X^2)) -2*atan(X);                         
end
function Y=PIE(X)
Y= 2*log(1+X^2);  
end