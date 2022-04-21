function RunLake(LakeIniFileName,MetFileName)
if nargin<2
    MetFileName='LAKE.MET';
end
if nargin<1
    LakeIniFileName='LAKE.ini';
end
%-------------------------------------- LICENCE BEGIN ------------------------------------
%Environment Canada - Atmospheric Science and Technology License/Disclaimer,
%                     version 3; Last Modified: May 7, 2008.
%This is free but copyrighted software; you can use/redistribute/modify it under the terms
%of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
%version 3 or (at your option) any later version that should be found at:
%http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
%
%This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
%without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%See the above mentioned License/Disclaimer for more details.
%You should have received a copy of the License/Disclaimer along with this software;
%if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
%CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
%-------------------------------------- LICENCE END --------------------------------------


% =======================================================================
%      * THE CANADIAN SMALL LAKE MODEL STANDALONE DRIVER
%      *
%      * AUTHOR: MURRAY MACKAY 		JAN 2014
%      *
% Ported to Matlab by M.G.Clark, Feb 2021
% =======================================================================

%       IMPLICIT NONE
%       INTEGER,PARAMETER :: NLAT=1,NMOS=3,ILG=NLAT*NMOS
%       INTEGER,PARAMETER :: NLAKMAX=200                              
%       INTEGER,PARAMETER :: NBS=4,IGL=1,IRSTRT=0
%       INTEGER MIDROT (NLAT,NMOS)
%       INTEGER N,NLTEST,NMTEST,I,J,K,L,M,JL1,JL2,
%      1        NCOUNT,NDAY,IHOUR,IMIN,IDAY,IYEAR,NML,NMW,
%      2        IDISP,IZREF,ISLFD,IPCP,ITG,IALS,ISNOALB
% 1-IMPLICIT NONE is a FORTRAN statement that overrides old FORTRAN code to 
% force the declaration of all variabiles. (apparently you didn't need to 
% define i,j,k,l,m as itagers in FORTRAN 77.) Don't need this in matlab
%
% 2-Matlab doesn't need to deine constants(i.e. PARAMETERs), but you could 
% override these as constants if needed for some reason.
%
% 3- the :: is only needed if additioanal sub-types are defined (i.e.
% PARAMETER)
%
% 4- Matlab allows for type specificaiton, but it is not needed.  Even
% though these are intigers in the fortran I am going to keep them as
% double.  This may cause errors if they are ever used for math functions,
% since in fortran int's are rounded down (i.e. the decimals are dropped)
% It is worth checking all refrences if the code does not run properly.
%

%
% 1-A bunch of these are static, so I'm moving the presistant call here
%
% 
% persistent DELT TFREZ
% persistent RGAS RGASV GRAV SBC VKC CT VMIN
% persistent HCPW HCPICE HCPSOL HCPOM HCPSND HCPCLY SPHW SPHICE SPHVEG ...
%     SPHAIR RHOW RHOICE TCGLAC CLHMLT CLHVAP
% persistent DELTA CGRAV CKARM VPDROW
% persistent TKECN TKECF TKECE TKECS HDPTHMIN TKEMIN DELMAX DELMIN EMSW ...
%     DELZLK DELSKIN DHMAX TKECL DUMAX
% persistent AS ASX CI BS BETA FACTN HMIN ANGMAX


NLAT=single(1);
NMOS=single(3);
ILG=single(NLAT*NMOS);
NLAKMAX=single(200);
NBS=single(4);
IGL=single(1);
IRSTRT=single(0);%This is 0 in RunLake.f, but in the CLASSL model "initial conditions are = 1 
N=single(0);%nan);
NLTEST=single(0);%nan);
NMTEST=single(0);%nan);
I=single(0);%nan);
J=single(0);%nan);
K=single(0);%nan);
L=single(0);%nan);
M=single(0);%nan);
JL1=single(0);%nan);
JL2=single(0);%nan);
NCOUNT=single(0);%nan);
NDAY=single(0);%nan);
IHOUR=single(0);%nan);
IMIN=single(0);%nan);
IDAY=single(0);%nan);
IYEAR=single(0);%nan);
NML=single(0);%nan);
NMW=single(0);%nan);
IDISP=single(0);%nan);
IZREF=single(0);%nan);
ISLFD=single(0);%nan);
IPCP=single(0);%nan);
ITG=single(0);%nan);
IALS=single(0);%nan);
ISNOALB=single(0);%nan;

% -- GATHER-SCATTER INDEX ARRAYS
%       INTEGER,DIMENSION(ILG) :: ILMOS,JLMOS,IWMOS,JWMOS
%
% 1-Same deal as above, intagers just arrays now.
%
% ILMOS=nan(ILG,1);
% JLMOS=nan(ILG,1);
% IWMOS=nan(ILG,1);
% JWMOS=nan(ILG,1);

[ILMOS,JLMOS,IWMOS,JWMOS] = makeZeros(ILG,1);

%-- LAKE TILE PARAMETERS                                            
%       INTEGER NLAKGAT(ILG), NLAKROT(NLAT,NMOS)
%       REAL,   DIMENSION(NLAT,NMOS) :: HLAKROT, LLAKROT, BLAKROT,    
%      1           HFSLROT, HEVLROT, FSGLROT, FLGLROT, HMFLROT,
%      2           ASVDROT, ASIDROT,REFROT,BCSNROT
%       REAL,   DIMENSION(ILG) :: HLAKGAT, LLAKGAT, BLAKGAT              
%       REAL,   DIMENSION(ILG) :: ASVLGAT, ASILGAT,BCSNLGAT,REFLGAT,
%      1           ZDMLGAT,ZDHLGAT
%
% 1- Arrays do not need to be defined with the sytax dimension(n), so the
% first integers defined here are a vector and an array.
%
% 2- The only change in this block is 32-bit floating point numbers (REAL), 
% so I will continue to follow the conventions to set to double (64-bit) as 
% above.
%
% NLAKROT=nan(NLAT,NMOS);
% HLAKROT=nan(NLAT,NMOS);
% LLAKROT=nan(NLAT,NMOS);
% BLAKROT=nan(NLAT,NMOS);
% HFSLROT=nan(NLAT,NMOS);
% HEVLROT=nan(NLAT,NMOS);
% FSGLROT=nan(NLAT,NMOS);
% FLGLROT=nan(NLAT,NMOS);
% HMFLROT=nan(NLAT,NMOS);
% ASVDROT=nan(NLAT,NMOS);
% ASIDROT=nan(NLAT,NMOS);
% REFROT=nan(NLAT,NMOS);
% BCSNROT=nan(NLAT,NMOS);
[NLAKROT,HLAKROT,LLAKROT,BLAKROT,HFSLROT,...
    HEVLROT,FSGLROT,FLGLROT,HMFLROT,ASVDROT,...
    ASIDROT,REFROT,BCSNROT]=makeZeros(NLAT,NMOS);
% NLAKGAT=nan(ILG,1);
% HLAKGAT=nan(ILG,1);
% LLAKGAT=nan(ILG,1);
% BLAKGAT=nan(ILG,1);
% ASVLGAT=nan(ILG,1);
% ASILGAT=nan(ILG,1);
% BCSNLGAT=nan(ILG,1);
% REFLGAT=nan(ILG,1);
% ZDMLGAT=nan(ILG,1);
% ZDHLGAT=nan(ILG,1);
[NLAKGAT,HLAKGAT,LLAKGAT,BLAKGAT,ASVLGAT,ASILGAT,...
    BCSNLGAT,REFLGAT,ZDMLGAT,ZDHLGAT] = makeZeros(ILG,1);
%-- LAKE PROGNOSTIC VARIABLES                                       
%       REAL,   DIMENSION(NLAT,NMOS,NLAKMAX) :: TLAKROT               
%       REAL,   DIMENSION(ILG,NLAKMAX) :: TLAKGAT                     
%       REAL,   DIMENSION(ILG) :: EXPW,DTEMP,HDPTH,DELU,GRED,         
%      1                    TKELAK,T0LAK,LKICEH,RHOMIX,TSED,SNICEH,ROFICEH
%
TLAKROT = single(zeros(NLAT,NMOS,NLAKMAX));
TLAKGAT = single(zeros(ILG,NLAKMAX));
% EXPW    = nan(ILG,1);
% DTEMP   = nan(ILG,1);
% HDPTH   = nan(ILG,1);
% DELU    = nan(ILG,1);
% GRED    = nan(ILG,1);
% TKELAK  = nan(ILG,1);
% T0LAK   = nan(ILG,1);
% LKICEH  = nan(ILG,1);
% RHOMIX  = nan(ILG,1);
% TSED    = nan(ILG,1);
% SNICEH  = nan(ILG,1);
% ROFICEH = nan(ILG,1);
[EXPW,DTEMP,HDPTH,DELU,GRED,TKELAK,...
    T0LAK,LKICEH,RHOMIX,TSED,SNICEH,ROFICEH] = makeZeros(ILG,1);

% -- LAKE DIAGNOSTIC VARIABLES                                       
%       REAL,   DIMENSION(ILG) :: HFSLGAT, HEVLGAT, FSGLGAT, FLGLGAT, 
%      1                       HMFLGAT,HTCLGAT,FICEGAT,FLSGAT,G0SLGAT,
%      2                       FSGSLGAT,FLGSLGAT,HFSSLGAT,HEVSLGAT,
%      3                       HMFNLGAT,HTCSLGAT
%       REAL,   DIMENSION(ILG) :: PCPLGAT,PCPNLGAT,QFLGAT,QFNLGAT,
%      1                       ROFNLGAT,SFTLGAT,SFULGAT,SFVLGAT,SFQLGAT,
%      2                       SFHLGAT,QLWOLGAT,ALVLGAT,ALILGAT,EFLGAT,
%      3                       GTLGAT,QGLGAT,DRLGAT,PETLGAT,QSENLGAT,
%      4                       TFXLGAT,QEVPLGAT,QFSLGAT,QFXLGAT,SNOLGAT,
%      5                       RHOSLGAT,TSNOLGAT,ALBSLGAT,WSNOLGAT
%       REAL,DIMENSION(NLAT,NMOS) :: PCPLROT,PCPNLROT,QFLROT,QFNLROT,
%      1                       ROFNLROT,FSGSLROT,FLGSLROT,HFSSLROT,
%      2                       HEVSLROT,HMFNLROT,HTCSLROT,HTCLROT,
%      3                       SFTLROT,SFULROT,SFVLROT,SFQLROT,SFHLROT,
%      4                       QLWOLROT,ALVLROT,ALILROT,EFLROT,GTLROT,
%      5                       QGLROT,DRLROT,PETLROT,QSENLROT,TFXLROT,
%      6                       QEVPLROT,QFSLROT,QFXLROT,SNOLROT,RHOSLROT,
%      7                       TSNOLROT,ALBSLROT,WSNOLROT
% CC    INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(15,307)
%       INTEGER, PARAMETER :: SP=kind(1.0)
%       INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(2*precision(1.0_sp))
%       REAL(dp),   DIMENSION(ILG) :: CTLSTP 
% C
%       REAL,DIMENSION(ILG,NBS) :: 
%      1        FSDBLGAT,   FSFBLGAT,   FSSBLGAT
%       REAL,DIMENSION(NLAT,NBS) ::
%      1        FSDBROL,   FSFBROL,   FSSBROL
%
% 1-Again for the first conversion I'm not going to mess with precision, so
% all the kind refrences are not needed since I'm defining REAL as doubles.
% HFSLGAT = nan(ILG,1);
% HEVLGAT = nan(ILG,1);
% FSGLGAT = nan(ILG,1);
% FLGLGAT = nan(ILG,1);
% HMFLGAT = nan(ILG,1);
% HTCLGAT = nan(ILG,1);
% FICEGAT = nan(ILG,1);
% FLSGAT = nan(ILG,1);
% G0SLGAT = nan(ILG,1);
% FSGSLGAT = nan(ILG,1);
% FLGSLGAT = nan(ILG,1);
% HFSSLGAT = nan(ILG,1);
% HEVSLGAT = nan(ILG,1);
% HMFNLGAT = nan(ILG,1);
% HTCSLGAT = nan(ILG,1);
% PCPLGAT = nan(ILG,1);
% PCPNLGAT = nan(ILG,1);
% QFLGAT = nan(ILG,1);
% QFNLGAT = nan(ILG,1);
% ROFNLGAT = nan(ILG,1);
% SFTLGAT = nan(ILG,1);
% SFULGAT = nan(ILG,1);
% SFVLGAT = nan(ILG,1);
% SFQLGAT = nan(ILG,1);
% SFHLGAT = nan(ILG,1);
% QLWOLGAT = nan(ILG,1);
% ALVLGAT = nan(ILG,1);
% ALILGAT = nan(ILG,1);
% EFLGAT = nan(ILG,1);
% GTLGAT = nan(ILG,1);
% QGLGAT = nan(ILG,1);
% DRLGAT = nan(ILG,1);
% PETLGAT = nan(ILG,1);
% QSENLGAT = nan(ILG,1);
% TFXLGAT = nan(ILG,1);
% QEVPLGAT = nan(ILG,1);
% QFSLGAT = nan(ILG,1);
% QFXLGAT = nan(ILG,1);
% SNOLGAT = nan(ILG,1);
% RHOSLGAT = nan(ILG,1);
% TSNOLGAT = nan(ILG,1);
% ALBSLGAT = nan(ILG,1);
% WSNOLGAT = nan(ILG,1);
[ HFSLGAT, HEVLGAT, FSGLGAT, FLGLGAT,...
	HMFLGAT,HTCLGAT,FICEGAT,FLSGAT,G0SLGAT,...
	FSGSLGAT,FLGSLGAT,HFSSLGAT,HEVSLGAT,...
    HMFNLGAT,HTCSLGAT,...
    PCPLGAT,PCPNLGAT,QFLGAT,QFNLGAT,....
    ROFNLGAT,SFTLGAT,SFULGAT,SFVLGAT,SFQLGAT,...
    SFHLGAT,QLWOLGAT,ALVLGAT,ALILGAT,EFLGAT,...
    GTLGAT,QGLGAT,DRLGAT,PETLGAT,QSENLGAT,...
    TFXLGAT,QEVPLGAT,QFSLGAT,QFXLGAT,SNOLGAT,...
    RHOSLGAT,TSNOLGAT,ALBSLGAT,WSNOLGAT] = makeZeros(ILG,1);

% PCPLROT = nan(NLAT,NMOS);
% PCPNLROT = nan(NLAT,NMOS);
% QFLROT = nan(NLAT,NMOS);
% QFNLROT = nan(NLAT,NMOS);
% ROFNLROT = nan(NLAT,NMOS);
% FSGSLROT = nan(NLAT,NMOS);
% FLGSLROT = nan(NLAT,NMOS);
% HFSSLROT = nan(NLAT,NMOS);
% HEVSLROT = nan(NLAT,NMOS);
% HMFNLROT = nan(NLAT,NMOS);
% HTCSLROT = nan(NLAT,NMOS);
% HTCLROT = nan(NLAT,NMOS);
% SFTLROT = nan(NLAT,NMOS);
% SFULROT = nan(NLAT,NMOS);
% SFVLROT = nan(NLAT,NMOS);
% SFQLROT = nan(NLAT,NMOS);
% SFHLROT = nan(NLAT,NMOS);
% QLWOLROT = nan(NLAT,NMOS);
% ALVLROT = nan(NLAT,NMOS);
% ALILROT = nan(NLAT,NMOS);
% EFLROT = nan(NLAT,NMOS);
% GTLROT = nan(NLAT,NMOS);
% QGLROT = nan(NLAT,NMOS);
% DRLROT = nan(NLAT,NMOS);
% PETLROT = nan(NLAT,NMOS);
% QSENLROT = nan(NLAT,NMOS);
% TFXLROT = nan(NLAT,NMOS);
% QEVPLROT = nan(NLAT,NMOS);
% QFSLROT = nan(NLAT,NMOS);
% QFXLROT = nan(NLAT,NMOS);
% SNOLROT = nan(NLAT,NMOS);
% RHOSLROT = nan(NLAT,NMOS);
% TSNOLROT = nan(NLAT,NMOS);
% ALBSLROT = nan(NLAT,NMOS);
% WSNOLROT = nan(NLAT,NMOS);
[PCPLROT,PCPNLROT,QFLROT,QFNLROT,...
    ROFNLROT,FSGSLROT,FLGSLROT,HFSSLROT,...
    HEVSLROT,HMFNLROT,HTCSLROT,HTCLROT,...
    SFTLROT,SFULROT,SFVLROT,SFQLROT,SFHLROT,...
    QLWOLROT,ALVLROT,ALILROT,EFLROT,GTLROT,...
    QGLROT,DRLROT,PETLROT,QSENLROT,TFXLROT,...
    QEVPLROT,QFSLROT,QFXLROT,SNOLROT,RHOSLROT,...
    TSNOLROT,ALBSLROT,WSNOLROT] = makeZeros(NLAT,NMOS);
CTLSTP = single(zeros(ILG,1)); %The precision of this was defined in FORTRAN so I'm leaving as double
% FSDBLGAT = nan(ILG,NBS);
% FSFBLGAT = nan(ILG,NBS);
% FSSBLGAT = nan(ILG,NBS);
[FSDBLGAT,FSFBLGAT,FSSBLGAT] = makeZeros(ILG,NBS);
% FSDBROL = nan(NLAT,NBS);
% FSFBROL = nan(NLAT,NBS);
% FSSBROL = nan(NLAT,NBS);
[FSDBROL,FSFBROL,FSSBROL] = makeZeros(NLAT,NBS);

% 
% -- ATMOSPHERIC AND GRID-CONSTANT INPUT VARIABLES
% 
%       REAL,DIMENSION(NLAT,NMOS) :: FAREROT
%       REAL,DIMENSION(NLAT) ::
%      1      ZRFMROW,   ZRFHROW,   ZDMROW ,   ZDHROW ,  
%      2      ZBLDROW,   FSVHROW,   FSIHROW,   RADJROW,
%      3      CSZROW ,   FDLROW ,   ULROW  ,   VLROW  ,   
%      4      TAROW  ,   QAROW  ,   PRESROW,   PREROW ,  
%      5      PADRROW,   VPDROW ,   TADPROW,   RHOAROW,  
%      6      UVROW  ,   GCROW  ,   RPCPROW,   TRPCROW,
%      7      SPCPROW,   TSPCROW,   RHSIROW,   RPREROW,
%      8      SPREROW   
% 
%       REAL,DIMENSION(ILG) ::
%      1      ZRFMGAT,   ZRFHGAT,   ZDMGAT ,   ZDHGAT ,  
%      2      ZBLDGAT,   FSVHGAT,   FSIHGAT,   RADJGAT,
%      3      CSZGAT ,   FDLGAT ,   ULGAT  ,   VLGAT  ,   
%      4      TAGAT  ,   QAGAT  ,   PRESGAT,   PREGAT ,  
%      5      PADRGAT,   VPDGAT ,   TADPGAT,   RHOAGAT,
%      6      RPCPLGAT,  TRPCLGAT,  SPCPLGAT,  TSPCLGAT,
%      7      RHSILGAT,  RADJLGAT,  PADRLGAT
%
FAREROT = single(zeros(NLAT,NMOS));
[ZRFMROW,   ZRFHROW,   ZDMROW ,   ZDHROW ,...  
    ZBLDROW,   FSVHROW,   FSIHROW,   RADJROW,...
    CSZROW ,   FDLROW ,   ULROW  ,   VLROW  ,...
    TAROW  ,   QAROW  ,   PRESROW,   PREROW ,...
    PADRROW,   VPDROW ,   TADPROW,   RHOAROW,...
    UVROW  ,   GCROW  ,   RPCPROW,   TRPCROW,...
    SPCPROW,   TSPCROW,   RHSIROW,   RPREROW,...
    SPREROW] = makeZeros(NLAT,1);
% ZRFMROW = nan(NLAT,1);ZRFHROW = nan(NLAT,1);ZDMROW = nan(NLAT,1);ZDHROW = nan(NLAT,1);
% ZBLDROW = nan(NLAT,1);FSVHROW = nan(NLAT,1);FSIHROW = nan(NLAT,1);RADJROW = nan(NLAT,1);
% CSZROW  = nan(NLAT,1);FDLROW  = nan(NLAT,1);ULROW   = nan(NLAT,1);VLROW = nan(NLAT,1);
% TAROW   = nan(NLAT,1);QAROW   = nan(NLAT,1);PRESROW = nan(NLAT,1);PREROW = nan(NLAT,1); 
% PADRROW = nan(NLAT,1);VPDROW  = nan(NLAT,1);TADPROW = nan(NLAT,1);RHOAROW = nan(NLAT,1);
% UVROW   = nan(NLAT,1);GCROW   = nan(NLAT,1);RPCPROW = nan(NLAT,1);TRPCROW = nan(NLAT,1);
% SPCPROW = nan(NLAT,1);TSPCROW = nan(NLAT,1);RHSIROW = nan(NLAT,1);RPREROW = nan(NLAT,1);
% SPREROW = nan(NLAT,1);  

% ZRFMGAT = nan(ILG,1);ZRFHGAT = nan(ILG,1);ZDMGAT  = nan(ILG,1);ZDHGAT = nan(ILG,1); 
% ZBLDGAT = nan(ILG,1);FSVHGAT = nan(ILG,1);FSIHGAT = nan(ILG,1);RADJGAT = nan(ILG,1);
% CSZGAT  = nan(ILG,1);FDLGAT  = nan(ILG,1);ULGAT   = nan(ILG,1);VLGAT = nan(ILG,1);
% TAGAT   = nan(ILG,1);QAGAT   = nan(ILG,1);PRESGAT = nan(ILG,1);PREGAT = nan(ILG,1);  
% PADRGAT = nan(ILG,1);VPDGAT  = nan(ILG,1);TADPGAT = nan(ILG,1);RHOAGAT = nan(ILG,1);
% RPCPLGAT = nan(ILG,1);TRPCLGAT = nan(ILG,1);SPCPLGAT = nan(ILG,1);TSPCLGAT = nan(ILG,1);
% RHSILGAT = nan(ILG,1);RADJLGAT = nan(ILG,1);PADRLGAT = nan(ILG,1);

[ZRFMGAT,   ZRFHGAT,   ZDMGAT ,   ZDHGAT ,...
    ZBLDGAT,   FSVHGAT,   FSIHGAT,   RADJGAT,...
    CSZGAT ,   FDLGAT ,   ULGAT  ,   VLGAT  ,...
    TAGAT  ,   QAGAT  ,   PRESGAT,   PREGAT ,...
    PADRGAT,   VPDGAT ,   TADPGAT,   RHOAGAT,...
    RPCPLGAT,  TRPCLGAT,  SPCPLGAT,  TSPCLGAT,...
    RHSILGAT,  RADJLGAT,  PADRLGAT] = makeZeros(ILG,1);
%
%-- LAND SURFACE DIAGNOSTIC VARIABLES
%
%       REAL,DIMENSION(NLAT,NMOS) :: CDHROT ,   CDMROT 
%       REAL,DIMENSION(ILG) ::       CDHGAT ,   CDMGAT 
%
[CDHROT,CDMROT] = makeZeros(NLAT,NMOS);
[CDHGAT,CDMGAT] = makeZeros(ILG,1);
%
%-- ARRAYS USED FOR OUTPUT AND DISPLAY PURPOSES.
%
%       CHARACTER     TITLE1*4,     TITLE2*4,     TITLE3*4,
%      1              TITLE4*4,     TITLE5*4,     TITLE6*4
%       CHARACTER     NAME1*4,      NAME2*4,      NAME3*4,
%      1              NAME4*4,      NAME5*4,      NAME6*4
%       CHARACTER     PLACE1*4,     PLACE2*4,     PLACE3*4,
%      1              PLACE4*4,     PLACE5*4,     PLACE6*4
%
% 1-In FORTRAN when defining character vectors (i.e. strings) they are
% defined in length (here it is 4 characters).  Like above you CAN do that 
% in MATLAB, or other high level languages, but don't need to.
%
TITLE1 = '    ';TITLE2 = '    ';TITLE3 = '    ';
TITLE4 = '    ';TITLE5 = '    ';TITLE6 = '    ';
NAME1 = '    ';NAME2 = '    ';NAME3 = '    ';
NAME4 = '    ';NAME5 = '    ';NAME6 = '    ';
PLACE1 = '    ';PLACE2 = '    ';PLACE3 = '    ';
PLACE4 = '    ';PLACE5 = '    ';PLACE6 = '    ';

% 
%-- LOCAL CONSTANTS AND VARIABLES.
%
%       REAL DEGLAT,DEGLON,FSDOWN,DAY,DECL,HOUR,COSZ,
%      1     EA,CA,CB,EASAT,CONST,HCAP,ZTOP,ZBOT,Z,QSUML,
%      2     RHOIW,ICETOP,ICEBOT
% %
% EGLAT = nan; DEGLON = nan; FSDOWN = nan; DAY = nan; DECL = nan; HOUR = nan; COSZ = nan; 
% EA = nan; CA = nan; CB = nan; EASAT = nan; CONST = nan; HCAP = nan; ZTOP = nan; ZBOT = nan; Z = nan; QSUML = nan; 
% RHOIW = nan; ICETOP = nan; ICEBOT = nan;

[EGLAT,DEGLON,FSDOWN,DAY,DECL,HOUR,COSZ,...
    EA,CA,CB,EASAT,CONST,HCAP,ZTOP,ZBOT,Z,QSUML,...
    RHOIW,ICETOP,ICEBOT] = makeZeros(1,1);
%
%-- COMMON BLOCK PARAMETERS.
%
%      REAL DELT,TFREZ, CLHVAP,PI,
%      1     RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN,TCW,TCICE,TCSAND,TCCLAY,
%      2     TCOM,TCDRYS,RHOSOL,RHOOM,HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,
%      3     HCPCLY,SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,TCGLAC,CLHMLT,
%      4     CKARM,CGRAV,CPD,DELTA,AS,ASX,ANGMAX,BS,BETA,CI,FACTN,HMIN


% DELT = nan; TFREZ = nan; CLHVAP = nan; PI = nan; 
% RGAS = nan; RGASV = nan; GRAV = nan; SBC = nan; VKC = nan; CT = nan; VMIN = nan; TCW = nan; TCICE = nan; TCSAND = nan; TCCLAY = nan; 
% TCOM = nan; TCDRYS = nan; RHOSOL = nan; RHOOM = nan; HCPW = nan; HCPICE = nan; HCPSOL = nan; HCPOM = nan; HCPSND = nan; 
% HCPCLY = nan; SPHW = nan; SPHICE = nan; SPHVEG = nan; SPHAIR = nan; RHOW = nan; RHOICE = nan; TCGLAC = nan; CLHMLT = nan; 
% CKARM = nan; CGRAV = nan; CPD = nan; DELTA = nan; AS = nan; ASX = nan; ANGMAX = nan; BS = nan; BETA = nan; CI = nan; FACTN = nan; HMIN = nan;
[DELT,TFREZ, CLHVAP,PI,...
    RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN,TCW,TCICE,TCSAND,TCCLAY,...
    TCOM,TCDRYS,RHOSOL,RHOOM,HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,...
    HCPCLY,SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,TCGLAC,CLHMLT,...
    CKARM,CGRAV,CPD,DELTA,AS,ASX,ANGMAX,BS,BETA,CI,FACTN,HMIN] = makeZeros(1,1);
%
%--  LAKE TILE COMMON BLOCK PARAMETERS                            
%    *                                                              
%       REAL TKECN,TKECF,TKECE,TKECS,HDPTHMIN, DUMAX, 
%      1  TKEMIN,DELMAX,DELMIN,EMSW,DELZLK,DELSKIN,DHMAX,TKECL
% TKECN = nan; TKECF = nan; TKECE = nan; TKECS = nan; HDPTHMIN = nan; DUMAX = nan;  
% TKEMIN = nan; DELMAX = nan; DELMIN = nan; EMSW = nan; DELZLK = nan; DELSKIN = nan; DHMAX = nan; TKECL = nan;
[TKECN,TKECF,TKECE,TKECS,HDPTHMIN, DUMAX,...
    TKEMIN,DELMAX,DELMIN,EMSW,DELZLK,DELSKIN,DHMAX,TKECL] = makeZeros(1,1);
%=======================================================================
%     * PHYSICAL CONSTANTS.
%     * THE FOLLOWING COMMON BLOCKS ARE DEFINED SPECIFICALLY FOR USE 
%     * IN CLASS, AND ARE RETAINED HERE TO MAINTAIN COMPATIBILITY
%
%     COMMON /CLASS1/ DELT,TFREZ
%       COMMON /CLASS2/ RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN
%       COMMON /CLASS3/ TCW,TCICE,TCSAND,TCCLAY,TCOM,TCDRYS,
%      1                RHOSOL,RHOOM
%       COMMON /CLASS4/ HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,
%      1                SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,
%      2                TCGLAC,CLHMLT,CLHVAP
%       COMMON /PHYCON/ DELTA,CGRAV,CKARM,CPD
%
% C
% C    * LAKE TILE COMMON BLOCK                                       
% C    *                                                              
%       COMMON /LAKECON/ TKECN,TKECF,TKECE,TKECS,HDPTHMIN,            
%      2             TKEMIN,DELMAX,DELMIN,EMSW,DELZLK,DELSKIN,DHMAX,  
%      3             TKECL,DUMAX                                      
% C
% C     * ADDITIONAL VALUES FOR RPN AND GCM COMMON BLOCKS.
% C
%       COMMON /CLASSD2/ AS,ASX,CI,BS,BETA,FACTN,HMIN,ANGMAX              
% C
% 1-Global variables are bad, I subscribe to the theory that they should
% never been used but the CSLM has a few common blocks defined so I'm going
% to have to find a workaround as the global function in matlab is messy.
% but for now just use the static variabile term persistent
% 2-Had to move the persistent declarations above the intitalization call
%

% 
%      * ADDITIONAL VALUES FOR RPN AND GCM COMMON BLOCKS.
% 
%       COMMON /CLASSD2/ AS,ASX,CI,BS,BETA,FACTN,HMIN,ANGMAX              
% 
%       DATA      DELTA,      AS,         ASX,        ANGMAX
%      1       /  0.608,      12.0,       4.7,        0.85   /
% 
%       DATA      CI,         BS,         BETA,       FACTN,      HMIN 
%      1       /  40.0,       1.0,        1.0,        1.2,        40.   /
%
DELTA   = single(0.608);
AS      = single(12.0);
ASX     = single(4.7);
ANGMAX  = single(0.85);
CI      = single(40.0);
BS      = single(1.0);
BETA    = single(1.0);
FACTN   = single(1.2);
HMIN    = single(40);                                           

% 
%     * PARAMETERS ORIGINALLY FROM CLASS COMMON BLOCK DATA           
%     *                                                              
%       DATA      VKC,        CT,         VMIN,     EMSW
%      1       /  0.40,       1.15E-3,    0.1,      0.97     /
%       DATA      TCW,        TCICE,      TCSAND,    TCCLAY,  TCOM
%      1       /  0.57,       2.24,       2.5,       2.5,     0.25  /
%       DATA      TCDRYS,     RHOSOL,     RHOOM
%      1       /  0.275,      2.65E3,     1.30E3  /    
%       DATA      HCPW,       HCPICE,     HCPSOL,    HCPOM
%      1       /  4.187E6,    1.9257E6,   2.25E6,    2.50E6 /
%       DATA      HCPSND,     HCPCLY,     SPHW,      SPHICE,  SPHVEG
%      1       /  2.13E6,     2.38E6,     4.186E3,   2.10E3,  2.70E3 /
%       DATA      RHOW,       RHOICE,     TCGLAC,    CLHMLT,  CLHVAP
%      1       /  1.0E3,      0.917E3,    2.24,      0.334E6, 2.501E6/
VKC    = single(0.40);
CT     = single(1.15E-3);
VMIN   = single(0.1);
EMSW   = single(0.97);
TCW    = single(0.57);
TCICE  = single(2.24);
TCSAND = single(2.5);
TCCLAY = single(2.5);
TCOM   = single(0.25);
TCDRYS = single(0.275);
RHOSOL = single(2.65E3);
RHOOM  = single(1.30E3);
HCPW   = single(4.187E6);
HCPICE = single(1.9257E6);
HCPSOL = single(2.25E6);
HCPOM  = single(2.50E6);
HCPSND = single(2.13E6);
HCPCLY = single(2.38E6);
SPHW   = single(4.186E3);
SPHICE = single(2.10E3);
SPHVEG = single(2.70E3);
RHOW   = single(1.0E3);
RHOICE = single(0.917E3);
TCGLAC = single(2.24);
CLHMLT = single(0.334E6);
CLHVAP = single(2.501E6);
% C
% C--Lake TKE process efficiencies                                    
%       DATA  TKECN,      TKECF,      TKECE,      TKECS,      TKECL   
%      1/     1.33,       0.25,       1.15,       0.20,      0.2350/ 
% CToo?1/     1.33,       0.25,       1.15,       0.20,      0.2250/
% CELA 1/     1.33,       0.25,       1.15,       0.20,      0.2350/
%
TKECN = single(1.33);
TKECF = single(0.25);
TKECE = single(1.15);
TKECS = single(0.20);
TKECL = single(0.2350);
% C
% C--Lake process parameter limits, and grid spacing                  
%       DATA  HDPTHMIN, TKEMIN,  DELMAX,  DELMIN,  DELZLK,  DELSKIN   
%      1/     0.5,      1.0E-12, 5.0,     0.5,     0.5,     0.050 /   
%       DATA  DHMAX,    DUMAX                                         
%      1/     2.0,      0.1   /        
%
HDPTHMIN = single(0.5);
TKEMIN   = single(1.0E-12);
DELMAX   = single(5.0);
DELMIN   = single(0.5);
DELZLK   = single(0.5);
DELSKIN  = single(0.050);
DHMAX    = single(2.0);
DUMAX    = single(0.1);
%
%     * ASSIGN VALUES NORMALLY SPECIFIED WITHIN THE GCM.
%
%       DATA      PI   /  3.1415926535898    /
%       DATA      GRAV,       RGAS,       SPHAIR,    RGASV
%      1/         9.80616,    287.04,     1.00464E3,  461.50   /
%       DATA      TFREZ,      SBC
%      1/         273.16,     5.66796E-8  /
%       DATA      DELT  
%      1/         600.0   /       ! 10 minutes (eg. Raksjon, ELA)
% C     1/         300.0   /       ! 5 minutes (eg. Toolik)
% C     1/         1800.0  /       ! half hrly
%
% 1- MATLAB has built in defintions for pi() but I'll redefine to keep the
% number of sig figs.
PI     = single(3.1415926535898);%Pi
GRAV   = single(9.80616);%gravity
RGAS   = single(287.04);
SPHAIR = single(1.00464E3);
RGASV  = single(461.50);
TFREZ  = single(273.16);%freezing temperature
SBC    = single(5.66796E-8);
DELT   = single(600.0); %10 minutes (eg. Raksjon, ELA), or 300.0 for 5 minutes (eg. Toolik), or 1800.0 for half hrly
CGRAV=GRAV;
CKARM=VKC;
CPD=SPHAIR;

% C
% C=======================================================================
% C     * CLASS SWITCHES select different runtime options
% C
% C     * IF IDISP=0, VEGETATION DISPLACEMENT HEIGHTS ARE IGNORED,
% C     * BECAUSE THE ATMOSPHERIC MODEL CONSIDERS THESE TO BE PART
% C     * OF THE "TERRAIN".
% C     * IF IDISP=1, VEGETATION DISPLACEMENT HEIGHTS ARE CALCULATED.
% C
% C     * IF IZREF=1, THE BOTTOM OF THE ATMOSPHERIC MODEL IS TAKEN
% C     * TO LIE AT THE GROUND SURFACE.
% C     * IF IZREF=2, THE BOTTOM OF THE ATMOSPHERIC MODEL IS TAKEN
% C     * TO LIE AT THE LOCAL ROUGHNESS HEIGHT.
% C
% C     * IF ISLFD=0, DRCOEF IS CALLED FOR SURFACE STABILITY CORRECTIONS
% C     * AND THE ORIGINAL GCM SET OF SCREEN-LEVEL DIAGNOSTIC CALCULATIONS 
% C     * IS DONE.
% C     * IF ISLFD=1, DRCOEF IS CALLED FOR SURFACE STABILITY CORRECTIONS
% C     * AND SLDIAG IS CALLED FOR SCREEN-LEVEL DIAGNOSTIC CALCULATIONS. 
% C     * IF ISLFD=2, FLXSURFZ IS CALLED FOR SURFACE STABILITY CORRECTIONS
% C     * AND DIASURF IS CALLED FOR SCREEN-LEVEL DIAGNOSTIC CALCULATIONS. 
% C
% C     * IF IPCP=1, THE RAINFALL-SNOWFALL CUTOFF IS TAKEN TO LIE AT 0 C.
% C     * IF IPCP=2, A LINEAR PARTITIONING OF PRECIPITATION BETWEEEN 
% C     * RAINFALL AND SNOWFALL IS DONE BETWEEN 0 C AND 2 C.
% C     * IF IPCP=3, RAINFALL AND SNOWFALL ARE PARTITIONED ACCORDING TO
% C     * A POLYNOMIAL CURVE BETWEEN 0 C AND 6 C.
% C     * IF IPCP=4, THE RAINFALL, SNOWFALL AND TOTAL PRECIPITATION RATES
% C     * ARE READ IN DIRECTLY.
% C
% C     * ITC, ITCG AND ITG ARE SWITCHES TO CHOOSE THE ITERATION SCHEME TO
% C     * BE USED IN CALCULATING THE CANOPY OR GROUND SURFACE TEMPERATURE
% C     * RESPECTIVELY.  IF THE SWITCH IS SET TO 1, A BISECTION METHOD IS
% C     * USED; IF TO 2, THE NEWTON-RAPHSON METHOD IS USED.
% C     
% C     * IF IPAI, IHGT, IALC, IALS AND IALG ARE ZERO, THE VALUES OF 
% C     * PLANT AREA INDEX, VEGETATION HEIGHT, CANOPY ALBEDO, SNOW ALBEDO
% C     * AND SOIL ALBEDO RESPECTIVELY CALCULATED BY CLASS ARE USED.
% C     * IF ANY OF THESE SWITCHES IS SET TO 1, THE VALUE OF THE
% C     * CORRESPONDING PARAMETER CALCULATED BY CLASS IS OVERRIDDEN BY
% C     * A USER-SUPPLIED INPUT VALUE.
% C
%       IDISP=1             !1 for stanalone with field obs.  Not used for lake model
%       IZREF=1             !1 for standalone with field obs.  Used in TLSPREP 
%       ISLFD=0
%       IPCP=1
%       ITG=1               !std
%       IALS=0
%       ISNOALB=0
% 
IDISP= single(1);            %1 for stanalone with field obs.  Not used for lake model
IZREF= single(1);            %1 for standalone with field obs.  Used in TLSPREP 
ISLFD= single(0);
IPCP= single(1);
ITG= single(1);              %std
IALS= single(0);
ISNOALB= single(0);

% <<--- The next section is all I/O formatting.  I'm going to overhaul this
% into more MATLAB friendly formatting.
%
% =======================================================================
%      * OPEN FILES FOR READING AND WRITING.
% 
%       OPEN(UNIT=50,FILE='LAKE.ini',STATUS='OLD')
%       OPEN(UNIT=51,FILE='LAKE.met',STATUS='OLD')
%       OPEN(UNIT=71,FILE='LAKE.of1')                                 
%       OPEN(UNIT=72,FILE='LAKE.of2')                                 
%       OPEN(UNIT=73,FILE='LAKE.of3')                                 
%       OPEN(UNIT=74,FILE='LAKE.of4')                                 
%       OPEN(UNIT=75,FILE='LAKE.of5')                                 
%       OPEN(UNIT=76,FILE='LAKE.of6')                                 
%       OPEN(UNIT=77,FILE='LAKE.of7')                                 
%       OPEN(UNIT=79,FILE='LAKE.of9')                                 
%
%     * READ AND PROCESS INITIALIZATION AND BACKGROUND INFORMATION.
%     * FIRST, MODEL RUN SPECIFICATIONS.
%
% 
%       READ (50,5010) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
%       READ (50,5010) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
%       READ (50,5010) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
% C                                                                   
% C   * LAKE DIAGNOSTIC OUTPUT FILES                                  
% C                                                                   
%       WRITE(71,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6      
%       WRITE(71,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6            
%       WRITE(72,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6      
%       WRITE(72,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6            
%       WRITE(73,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6      
%       WRITE(73,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6            
%       WRITE(74,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6      
%       WRITE(74,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6            
%       WRITE(75,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6      
%       WRITE(75,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6            
%       WRITE(76,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6      
%       WRITE(76,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6            
%       WRITE(77,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6      
%       WRITE(77,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6            
%       WRITE(79,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6      
%       WRITE(79,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6            
% C
%       WRITE(71,7011)                                                
% 7011  FORMAT('YEAR DAY  HR MIN     E0     F0      Q*      Q0',      
%      1       '      L*     Hs      He     T0      Lake Ice')        
%       WRITE(72,7012)                                                
% 7012  FORMAT('YEAR DAY  HR MIN     TEMP ')                          
%       WRITE(73,7013)                                                
% 7013  FORMAT('YEAR DAY  HR MIN     FQU      BFLX      DISS',        
%      1 '    TKE DELU HDPTH JMIX ZMIX  TMIX    DTEMP',               
%      2 '    FSHEAR FENTRA')                                         
%       WRITE(74,7014)                                                
% 7014  FORMAT('YEAR DAY  HR MIN     U*     GRED      WEDB      ',    
%      1       '  HDPTH DELTHRM DELU')                                
%       WRITE(75,7015)                                                
% 7015  FORMAT('YEAR DAY  HR MIN     TSNOW   WSNOW    RHOSNO',
%      1       '   ZSNOW')
%       WRITE(76,7016)                                                
% 7016  FORMAT('YEAR DAY  HR MIN     QSWNS   QTRANSL  QLNET ',
%      1       '   QSENSS  QEVAPS   QMELT   GZEROSL')
%       WRITE(77,7017)                                                
% 7017  FORMAT('YEAR DAY  HR MIN     QSWIN   FSGL    QFLX(1m)')
%       WRITE(79,7019)                                                
% 7019  FORMAT('YEAR DAY  HR MIN     EnergyBalance')
% C                                                                   
% C=======================================================================
% C   * READ GEOPHYSICAL DATA
% C
%       READ(50,5020) DEGLAT,DEGLON,ZRFMROW(1),ZRFHROW(1),ZBLDROW(1),
%      1              GCROW(1),NLTEST,NMTEST
% 
%       DO 50 I=1,NLTEST
%       DO 50 M=1,NMTEST
%           READ(50,5010) TITLE1                                      
%           READ(50,5045) MIDROT(I,M),FAREROT(I,M)                    
%           IF (MIDROT(I,M) > 0.0 ) THEN                           
%             CALL XIT('RUNLAKE',-1)
%       else
% C-- READ IN WATER TILE INFORMATION                              
%            READ(50,5095) HLAKROT(I,M),LLAKROT(I,M),BLAKROT(I,M)
%            NLAKROT(I,M)=NINT(HLAKROT(I,M)/DELZLK)
%            READ(50,5096) (TLAKROT(I,M,J),J=1,NLAKROT(I,M))          
%       end
%
% 50    CONTINUE
%       DO 100 I=1,NLTEST
%       DO 100 M=1,NMTEST
%           IF (NLAKROT(I,M) >= 1) THEN                             
%           DO 73 J=1,NLAKROT(I,M)                                    
%             IF (TLAKROT(I,M,J)<=200.0) THEN                       
%               TLAKROT(I,M,J)=TLAKROT(I,M,J)+TFREZ                   
%             end                                                   
% 73        CONTINUE                                                  
%           end                                                     
% C-- INITIALIZE LAKE SNOW CONDITIONS
%           SNOLROT(I,M)=0.0
%           ALBSLROT(I,M)=0.0
%           RHOSLROT(I,M)=0.0
%           TSNOLROT(I,M)=0.0
%           WSNOLROT(I,M)=0.0
% 100   CONTINUE
baseDir='D:\Dropbox\PostDoc\UofS_Course\Project\TestOutput';
fID_71=fopen([baseDir,'\LAKE.of1'],'W');
fID_72=fopen([baseDir,'\LAKE.of2'],'W');
fID_73=fopen([baseDir,'\LAKE.of3'],'W');
fID_74=fopen([baseDir,'\LAKE.of4'],'W');
fID_75=fopen([baseDir,'\LAKE.of5'],'W');
fID_76=fopen([baseDir,'\LAKE.of6'],'W');
fID_77=fopen([baseDir,'\LAKE.of7'],'W');
fID_82=fopen([baseDir,'\ExtraEditing.txt'],'W');

% C                                                                   
% C=======================================================================
% C   * READ GEOPHYSICAL DATA
% C
%Open Lake initialization file and read into variabiles.
fidi = fopen(LakeIniFileName);
%First three lines are the title (i.e. lake number & date)
%Then the author
%third is the place of work 
tline = fgetl(fidi);%Read in line one, should be lake # and date
TITLE1=tline;
tline = fgetl(fidi);%read in line two, should be oprator name
NAME1=tline;
tline = fgetl(fidi);%read in line three, should be institution.
PLACE1=tline;
%First line of data
tline = fgetl(fidi);%read in line 4, should be initalization data
strFromFile=split(tline,' ');
A=ParseStringIntoSingle(strFromFile);
DEGLAT=A(1);
DEGLON=A(2);
ZRFMROW(1)=A(3);
ZRFHROW(1)=A(4);
ZBLDROW(1)=A(5);
GCROW(1)=A(6);
NLTEST=A(7);
NMTEST=A(8);
for I=1:NLTEST
    for M=1:NMTEST
        tline = fgetl(fidi);%the next three lines may be repeated if there are multiple lakes
        TITLE1 = tline;%first line is the lake info
        tline = fgetl(fidi);%second line is two variabiles MIDROT and FAREROT
        strFromFile=split(tline,' ');
        A=ParseStringIntoSingle(strFromFile);
        MIDROT(I,M)=A(1);
        FAREROT(I,M)=A(2);
        if (MIDROT(I,M) > 0.0) 
            error('RUNLAKE -> Issues loading Lake data, MIDROT is > 0');%if MIDROT is >0 safely crash and throw error
        else
            %C-- READ IN WATER TILE INFORMATION                              
            tline = fgetl(fidi);%thrid line is HLAKROT, LLAKROT and BLAKROT
            strFromFile=split(tline,' ');
            A=ParseStringIntoSingle(strFromFile);
            HLAKROT(I,M)=A(1);
            LLAKROT(I,M)=A(2);
            BLAKROT(I,M)=A(3);
            NLAKROT(I,M)=round(HLAKROT(I,M)/DELZLK);%calculate NLAKROT
            tline = fgetl(fidi);
            strFromFile=split(tline,' ');
            A=ParseStringIntoSingle(strFromFile);%This function simple cuts out white space and returns only the numerical values in the string seperated by whitespace as an array
            for J=1:NLAKROT(I,M)%FORTRAN code cut this off at NLAKROT (I think, line 332 in RUNLAKE.f), but I could also use length(A) if the other number are relevent
                  TLAKROT(I,M,J)=single(A(J));%,J=1,NLAKROT(I,M))          
            end
        end
    end
end
fclose(fidi);

%Flip over lake input values, if in deg C convert to deg K by adding
%freezing temp.
for I=1:NLTEST
    for M=1:NMTEST
        if (NLAKROT(I,M)>=1)  %This should be GEn                           
            for J=1:NLAKROT(I,M)                                    
                if (TLAKROT(I,M,J)<=200.0)
                    TLAKROT(I,M,J)=TLAKROT(I,M,J)+TFREZ;                   
                end
            end
        end
    % -- INITIALIZE LAKE SNOW CONDITIONS
    SNOLROT(I,M)=0.0;
    ALBSLROT(I,M)=0.0;
    RHOSLROT(I,M)=0.0;
    TSNOLROT(I,M)=0.0;
    WSNOLROT(I,M)=0.0;
    end
end

%Next lines in fortran are formatting settings for reading the ASCII files,
%not needed with matlab.
% 5010  FORMAT(2X,6A4)
% 5020  FORMAT(5F10.2,F7.1,3I5)
% 5045  FORMAT(I8,F8.3)                                               
% 5095  FORMAT(3F10.1)                                      
% 5096  FORMAT(200F6.2)                                               
% 5300  FORMAT(1X,I2,I3,I5,I6,2F9.2,E14.4,F9.2,E12.3,F8.2,F12.2,3F9.2,
%      1       F9.4)
% 6001  FORMAT('CLASS TEST RUN:     ',6A4)
% 6002  FORMAT('RESEARCHER:         ',6A4)
% 6003  FORMAT('INSTITUTION:        ',6A4)


% C
% C------------------------------------------------------------------
% C THE CLASS GATHER-SCATTER MACHINERY IS RETAINED HERE FOR
% C COMPATIBILITY
% C------------------------------------------------------------------
%       NML=0
%       NMW=0
%       DO 110 I=1,NLTEST
%               DO 120 J=1,NMTEST
%                   if(FAREROT(I,J)>0.0)      THEN
%                       if(MIDROT(I,J)>0)    THEN
%                           NML=NML+1
%                           ILMOS(NML)=I
%                           JLMOS(NML)=J
%                       else
%                           NMW=NMW+1
%                           IWMOS(NMW)=I
%                           JWMOS(NMW)=J
%                       end
%                   end
%   120         CONTINUE
%   110 CONTINUE
%       print*, "NML=",NML, NMW
% 
NML=single(0);
NMW=single(0);
for I=1:NLTEST
    for J=1:NMTEST
        if(FAREROT(I,J)>0.0) 
            if(MIDROT(I,J)>0)
                NML=NML+1;
                ILMOS(NML)=I;
                JLMOS(NML)=J;
            else
                NMW=NMW+1;
                IWMOS(NMW)=I;
                JWMOS(NMW)=J;
            end
        end
    end
end
fprintf('\n NML= %d,%d', NML, NMW);

% C=======================================================================
% C     * LAUNCH RUN.
% 
%       N=0
%       NCOUNT=1
%       NDAY=86400/NINT(DELT)
% 
% 200   CONTINUE

N=single(0);
NCOUNT=single(1);
NDAY=86400/round(DELT);
fprintf('\n    N, Nfrac, deltaT, totalT, IYEAR, IDAY, IHOUR, IMIN, QSUML-CTLSTP(I), QSUML');%Not in FORTRAN, just used for editing output (headers for line 1212)
oldTime=now;
firstTime=oldTime;
%Fortran uses an end of file command to jump to the end of the program.
%Thats a little problematic, so instead I just get the length of the input
%met file and loop over it's length.
[LngthOfMetDat]=GetLengthOfMetData(MetFileName);
fID_51=fopen(MetFileName);
Classlin=nan(LngthOfMetDat,95);
Classlout=nan(LngthOfMetDat,108);
%wb_id=waitbar(0,'Running Lake Model');
for i=1:LngthOfMetDat    
%     if rem(i,100)==0
%         waitbar(i/L,wb_id,'Running Lake Model');
%     end
    %the continue statement implies a loop, and there is a goto:200 way below
    %so I will either need a while or a for loop here but first I need to know
    %how the code works so I'll return to this.

    %C
    % C========================================================================
    % C     * READ IN METEOROLOGICAL FORCING DATA FOR CURRENT TIME STEP;
    % C     * CALCULATE SOLAR ZENITH ANGLE AND COMPONENTS OF INCOMING SHORT-
    % C     * WAVE RADIATION FLUX; ESTIMATE FLUX PARTITIONS IF NECESSARY.
    % C
    %       N=N+1
    %       DO 250 I=1,NLTEST
    %           READ(51,5300,END=999) IHOUR,IMIN,IDAY,IYEAR,FSDOWN,FDLROW(I),
    %      1         PREROW(I),TAROW(I),QAROW(I),UVROW(I),PRESROW(I)
    %
    % 1- I moved the read command into a function, but I pull in all the data
    % in the file, may only want the one row, see how it works.
    %
    N=N+1;
    for I=1:NLTEST
        %[IHOUR,IMIN,IDAY,IYEAR,FSDOWN,FDLROW(I),PREROW(I),TAROW(I),QAROW(I),UVROW(I),PRESROW(I)]=GetMetData('LAKE.MET',N);
        %In FORTRAN the read statement pulls in the line, but when it
        %encounters an end of line statement it jumps (GOTO) 999, which is
        %end of the RunLake code.
        tline = fgetl(fID_51);%Replacing the load function with reading line by line
        %this will need to be optimized, if NLTEST > 1 since the for loop
        %will continue exicuting
        strFromFile=split(tline,' ');
        metDat=ParseStringIntoSingle(strFromFile);
        IHOUR=metDat(1);
        IMIN=metDat(2);
        IDAY=metDat(3);
        IYEAR=metDat(4);% 
        FSDOWN=metDat(5);% downwelling shortwave
        FDLROW=metDat(6);% downwelling longwave
        PREROW=metDat(7); % Not sure (2.9084E-05) ~precipitation (m I think)
        TAROW=metDat(8); %Air Temp
        QAROW=metDat(9); %not sure (1.296E-02) ~ spacific humidity
        UVROW=metDat(10); %not sure (3.45) ~ windspeed
        PRESROW=metDat(11);%pressure in Pa

    %           FSVHROW(I)=0.5*FSDOWN
    %           FSIHROW(I)=0.5*FSDOWN
    %           TAROW(I)=TAROW(I)+TFREZ
    %           ULROW(I)=UVROW(I)
    %           VLROW(I)=0.0
    %           UVROW(I)=MAX(VMIN,UVROW(I))
    %           ZDMROW(I)=10.0
    %           ZDHROW(I)=2.0
    %           FSSBROL(I,1)=FSVHROW(I)
    %           FSSBROL(I,2)=FSIHROW(I)
    %           RPREROW(I)=0.0           !only needed for IPC=4
    %           SPREROW(I)=0.0           !only needed for IPC=4
    % 250   CONTINUE
    %
    % 1 - This is in a loop in FORTRAN but I'm not sure it's needed, related to
    % loading function.
        FSVHROW=0.5*FSDOWN;
        FSIHROW=0.5*FSDOWN;
        TAROW=TAROW+TFREZ;
        ULROW=UVROW; %all wind is assumed downwind for the first timestep.
        VLROW=single(0.0); %cross-wind speed
        %UVROW(UVROW<VMIN)=VMIN;%This is more efficent than using the max function in a loop
        UVROW(I)=MAX(VMIN,UVROW(I));
        ZDMROW(I)=10.0;
        ZDHROW(I)=2.0;
        FSSBROL(I,1)=FSVHROW;
        FSSBROL(I,2)=FSIHROW;
        RPREROW(I)=0.0;           %only needed for IPC=4
        SPREROW(I)=0.0;          %only needed for IPC=4
    end

    % C
    %       DAY=REAL(IDAY)+(REAL(IHOUR)+REAL(IMIN)/60.)/24.
    %       DECL=SIN(2.*PI*(284.+DAY)/365.)*23.45*PI/180.
    %       HOUR=(REAL(IHOUR)+REAL(IMIN)/60.)*PI/12.-PI
    %       DO 300 I=1,NLTEST
    %           RADJROW(I)=DEGLAT*PI/180.
    %           COSZ=SIN(RADJROW(I))*SIN(DECL)+COS(RADJROW(I))*
    %      >                     COS(DECL)*COS(HOUR)
    %           CSZROW(I)=SIGN(MAX(ABS(COSZ),1.0E-3),COSZ)
    % 300   CONTINUE
    %
    % 1. MATLAB does not need to specify real.  I have made every numeric
    % variable a double up to this point.  If later the code is revsed to
    % seperate intagers from floating point numbres this might need to be
    % revisited but matlab likely can handle the converison. Using alternate 
    % numeric variabiles is cumbersome in MATLAB which is optomized for 64-bit 
    % calculations.  In addition the number of bits it saves is really marginal 
    % for this application.
    % 2. You do not need the '.' following floating numbers in matlab, since
    % matlab assumes all numbers are double precision. 
    DAY=IDAY+(IHOUR+IMIN/60)/24;
    DECL=sin(2.*PI*(284+DAY)/365)*23.45*PI/180;
    HOUR=(IHOUR+IMIN/60)*PI/12-PI;
    for I=1:NLTEST
        RADJROW(I)=DEGLAT*PI/180;
        COSZ=sin(RADJROW(I))*sin(DECL)+cos(RADJROW(I))*cos(DECL)*cos(HOUR);
%       OLD MATLAB WAY, now MAX function works like FORTRAN
%         tmpCOSZ=COSZ;
%         tmpCOSZ(abs(COSZ)<1.0E-3)=1.0E-3;
%         if (tmpCOSZ<0&COSZ>0)|(tmpCOSZ>0&COSZ<0)
%             tmpCOSZ=tmpCOSZ*-1;
%         end
%         CSZROW=tmpCOSZ;
        CSZROW(I)=SIGN(MAX(abs(COSZ),1.0E-3),COSZ);
%         if COSZ>0
%             CSZROW(I)=sign(MAX(abs(COSZ),1.0E-3));
%         else
%             CSZROW(I)=-sign(MAX(abs(COSZ),1.0E-3));
%         end
    end
    % C================================================================
    % C
    % C     * CALCULATION OF ATMOSPHERIC INPUT VARIABLES.
    % C
    [VPDROW,TADPROW,PADRROW,RHOAROW,RHSIROW,RPCPROW,TRPCROW,SPCPROW,TSPCROW,...
        TAROW,QAROW,PREROW,RPREROW,SPREROW,PRESROW,IPCP,NLAT,~,NLTEST] = ...
      CLASSI(VPDROW,TADPROW,PADRROW,RHOAROW,RHSIROW,RPCPROW,TRPCROW,...
        SPCPROW,TSPCROW,TAROW,QAROW,PREROW,RPREROW,SPREROW,PRESROW,IPCP,...
        NLAT,1,NLTEST,TFREZ,RGAS,RGASV,RHOW);

    % C
    % C======================================================================
    % C GATHER CODE FROM CLASSG
    % C-- GATHER LAKE RELEVANT VARIABLES INTO WATER-TILE GAT ARRAYS                                       
    % C                                                                   
    %       DO 700 K=1,NMW                                                
    %           HLAKGAT(K+NML)=HLAKROT(IWMOS(K),JWMOS(K))                 
    %           LLAKGAT(K+NML)=LLAKROT(IWMOS(K),JWMOS(K))                 
    %           BLAKGAT(K+NML)=BLAKROT(IWMOS(K),JWMOS(K))                 
    %           NLAKGAT(K+NML)=NLAKROT(IWMOS(K),JWMOS(K))                 
    %           DO 705 L=1,NLAKGAT(K+NML)                                 
    %             TLAKGAT(K+NML,L)=TLAKROT(IWMOS(K),JWMOS(K),L)           
    % 705       CONTINUE                                                  
    %           ASVLGAT(K+NML)=ASVDROT(IWMOS(K),JWMOS(K))
    %           ASILGAT(K+NML)=ASIDROT(IWMOS(K),JWMOS(K))
    %           BCSNLGAT(K+NML)=BCSNROT(IWMOS(K),JWMOS(K))
    %           REFLGAT(K+NML)=REFROT(IWMOS(K),JWMOS(K))
    %           SNOLGAT(K+NML)=SNOLROT(IWMOS(K),JWMOS(K))
    %           RHOSLGAT(K+NML)=RHOSLROT(IWMOS(K),JWMOS(K))
    %           TSNOLGAT(K+NML)=TSNOLROT(IWMOS(K),JWMOS(K))
    %           ALBSLGAT(K+NML)=ALBSLROT(IWMOS(K),JWMOS(K))
    %           WSNOLGAT(K+NML)=WSNOLROT(IWMOS(K),JWMOS(K))
    % 700   CONTINUE                                                      
    for K=1:NMW                                                
        HLAKGAT(K+NML)=HLAKROT(IWMOS(K),JWMOS(K));                
        LLAKGAT(K+NML)=LLAKROT(IWMOS(K),JWMOS(K));                 
        BLAKGAT(K+NML)=BLAKROT(IWMOS(K),JWMOS(K));                 
        NLAKGAT(K+NML)=NLAKROT(IWMOS(K),JWMOS(K));                 
        for L=1:NLAKGAT(K+NML)  
            TLAKGAT(K+NML,L)=TLAKROT(IWMOS(K),JWMOS(K),L);
        end            
        ASVLGAT(K+NML)=ASVDROT(IWMOS(K),JWMOS(K));
        ASILGAT(K+NML)=ASIDROT(IWMOS(K),JWMOS(K));
        BCSNLGAT(K+NML)=BCSNROT(IWMOS(K),JWMOS(K));
        REFLGAT(K+NML)=REFROT(IWMOS(K),JWMOS(K));
        SNOLGAT(K+NML)=SNOLROT(IWMOS(K),JWMOS(K));
        RHOSLGAT(K+NML)=RHOSLROT(IWMOS(K),JWMOS(K));
        TSNOLGAT(K+NML)=TSNOLROT(IWMOS(K),JWMOS(K));
        ALBSLGAT(K+NML)=ALBSLROT(IWMOS(K),JWMOS(K));
        WSNOLGAT(K+NML)=WSNOLROT(IWMOS(K),JWMOS(K));
    end
    % C
    % C-- ATMOSPHERIC FORCING VARIABLES NEEDED FOR LAKE TILES         
    % C   GATHERED ON TOP OF LAND TILES                               
    % C
    %       DO 800 K=1,NMW                                                
    %           FSVHGAT(K+NML)=FSVHROW(IWMOS(K))                          
    %           FSIHGAT(K+NML)=FSIHROW(IWMOS(K))                          
    %           CSZGAT (K+NML)=CSZROW (IWMOS(K))                          
    %           FDLGAT (K+NML)=FDLROW (IWMOS(K))                          
    %           ULGAT  (K+NML)=ULROW  (IWMOS(K))                          
    %           VLGAT  (K+NML)=VLROW  (IWMOS(K))                           
    %           TAGAT  (K+NML)=TAROW  (IWMOS(K))                          
    %           QAGAT  (K+NML)=QAROW  (IWMOS(K))                          
    %           PRESGAT(K+NML)=PRESROW(IWMOS(K))                          
    %           RHOAGAT(K+NML)=RHOAROW(IWMOS(K))                          
    %           ZRFMGAT(K+NML)=ZRFMROW(IWMOS(K))                          
    %           ZRFHGAT(K+NML)=ZRFHROW(IWMOS(K))                          
    %           DO L=1,NBS
    %             FSDBLGAT(K+NML,L)=FSDBROL(IWMOS(K),L) 
    %             FSFBLGAT(K+NML,L)=FSFBROL(IWMOS(K),L) 
    %             FSSBLGAT(K+NML,L)=FSSBROL(IWMOS(K),L) 
    %           ENDDO
    %           ZDMLGAT(K+NML)=ZDMROW(IWMOS(K))
    %           ZDHLGAT(K+NML)=ZDHROW(IWMOS(K))
    %           RPCPLGAT(K+NML)=RPCPROW(IWMOS(K))
    %           TRPCLGAT(K+NML)=TRPCROW(IWMOS(K))
    %           SPCPLGAT(K+NML)=SPCPROW(IWMOS(K))
    %           TSPCLGAT(K+NML)=TSPCROW(IWMOS(K))
    %           RHSILGAT(K+NML)=RHSIROW(IWMOS(K))
    %           RADJLGAT(K+NML)=RADJROW(IWMOS(K))
    %           PADRLGAT(K+NML)=PADRROW(IWMOS(K))
    % 800   CONTINUE  
    for K=1:NMW                                                
        FSVHGAT(K+NML)=FSVHROW(IWMOS(K));                       
        FSIHGAT(K+NML)=FSIHROW(IWMOS(K));                          
        CSZGAT (K+NML)=CSZROW (IWMOS(K));                          
        FDLGAT (K+NML)=FDLROW (IWMOS(K));                          
        ULGAT  (K+NML)=ULROW  (IWMOS(K));                          
        VLGAT  (K+NML)=VLROW  (IWMOS(K));                           
        TAGAT  (K+NML)=TAROW  (IWMOS(K));                          
        QAGAT  (K+NML)=QAROW  (IWMOS(K));                          
        PRESGAT(K+NML)=PRESROW(IWMOS(K));                          
        RHOAGAT(K+NML)=RHOAROW(IWMOS(K));                          
        ZRFMGAT(K+NML)=ZRFMROW(IWMOS(K));                          
        ZRFHGAT(K+NML)=ZRFHROW(IWMOS(K));                          
        for L=1:NBS
            FSDBLGAT(K+NML,L)=FSDBROL(IWMOS(K),L); 
            FSFBLGAT(K+NML,L)=FSFBROL(IWMOS(K),L); 
            FSSBLGAT(K+NML,L)=FSSBROL(IWMOS(K),L); 
        end
        ZDMLGAT(K+NML)=ZDMROW(IWMOS(K));
        ZDHLGAT(K+NML)=ZDHROW(IWMOS(K));
        RPCPLGAT(K+NML)=RPCPROW(IWMOS(K));
        TRPCLGAT(K+NML)=TRPCROW(IWMOS(K));
        SPCPLGAT(K+NML)=SPCPROW(IWMOS(K));
        TSPCLGAT(K+NML)=TSPCROW(IWMOS(K));
        RHSILGAT(K+NML)=RHSIROW(IWMOS(K));
        RADJLGAT(K+NML)=RADJROW(IWMOS(K));
        PADRLGAT(K+NML)=PADRROW(IWMOS(K));
    end
    % C========================================================================
    % C CHECK ENERGY BALANCE - start of timestep
    % C ICEBOT doesn't include weight of snow
    % C
    %       RHOIW=RHOICE/RHOW
    %       JL1=1+NML                                                     
    %       JL2=NMW+NML                                                  
    %       DO 320 I=JL1,JL2
    %         ICEBOT=RHOIW*LKICEH(I)
    %         ICETOP=LKICEH(I)-ICEBOT
    %         IF (ICEBOT >= DELSKIN) THEN 
    %           HCAP=HCPICE
    %         else IF (LKICEH(I) <= 0.0) THEN
    %           HCAP=HCPW
    %         else 
    %           HCAP=(LKICEH(I)*HCPICE + (DELSKIN-ICEBOT)*HCPW)/DELSKIN
    %         end
    %         IF (N == 1) THEN 
    %           CTLSTP(I)= -HCAP*TLAKGAT(I,1)*DELSKIN
    %         else
    %           CTLSTP(I)= -HCAP*T0LAK(I)*DELSKIN
    %         end
    % C
    %         DO 330, J=1,NLAKGAT(I)
    %           ZTOP=DELSKIN + DELZLK*(J -1)
    %           ZBOT=DELSKIN + DELZLK*J
    %           IF (ICEBOT >= ZBOT) THEN 
    %             HCAP=HCPICE
    %           else IF (ICEBOT <= ZTOP) THEN
    %             HCAP=HCPW
    %           else 
    %             Z=ICEBOT-ZTOP
    %             HCAP=(Z*HCPICE + (DELZLK-Z)*HCPW)/DELZLK
    %           end
    %           CTLSTP(I)=CTLSTP(I) - HCAP*TLAKGAT(I,J)*DELZLK
    % 330     CONTINUE
    % 320   CONTINUE
    RHOIW=RHOICE/RHOW;
    JL1=1+NML;                                            
    JL2=NMW+NML;                                                  
    for I=JL1:JL2
        ICEBOT=RHOIW*LKICEH(I);
        ICETOP=LKICEH(I)-ICEBOT;
        if (ICEBOT >= DELSKIN)
        HCAP=HCPICE;
        elseif(LKICEH(I) <= 0.0) 
            HCAP=HCPW;
        else
            HCAP=(LKICEH(I)*HCPICE + (DELSKIN-ICEBOT)*HCPW)/DELSKIN;
        end
        if (N == 1) 
            CTLSTP(I)= -HCAP*TLAKGAT(I,1)*DELSKIN;
        else
            CTLSTP(I)= -HCAP*T0LAK(I)*DELSKIN;
        end
        %White space
        for J=1:NLAKGAT(I)
            ZTOP=DELSKIN + DELZLK*(J -1);
            ZBOT=DELSKIN + DELZLK*J;
            if (ICEBOT >= ZBOT)
                HCAP=HCPICE;
            elseif(ICEBOT <= ZTOP)
                    HCAP=HCPW;
            else
                Z=ICEBOT-ZTOP;
                HCAP=(Z*HCPICE + (DELZLK-Z)*HCPW)/DELZLK;
            end
            CTLSTP(I)=CTLSTP(I) - HCAP*TLAKGAT(I,J)*DELZLK;
        end
    end
    % C
    % C========================================================================
    % C          * LAKE MODEL
    % C                                                                   
    %       CALL CLASSL (HLAKGAT, LLAKGAT, BLAKGAT, NLAKGAT, TLAKGAT,   
    %      1 T0LAK, HDPTH, LKICEH, SNICEH, ROFICEH,
    %      2 SNOLGAT,RHOSLGAT,TSNOLGAT,ALBSLGAT,WSNOLGAT,
    %      3 CDHGAT,CDMGAT,QSENLGAT,TFXLGAT,QEVPLGAT,QFSLGAT,QFXLGAT,
    %      4 PETLGAT, EFLGAT, GTLGAT, QGLGAT, DRLGAT,
    %      5 SFTLGAT,SFULGAT,SFVLGAT,SFQLGAT,SFHLGAT,QLWOLGAT,ALVLGAT,ALILGAT,
    %      6 FSGLGAT, FLGLGAT, HFSLGAT, HEVLGAT, HMFLGAT, HTCLGAT,
    %      7 FSGSLGAT, FLGSLGAT, HFSSLGAT, HEVSLGAT, HMFNLGAT, HTCSLGAT,
    %      8 PCPLGAT,PCPNLGAT,QFLGAT,QFNLGAT,ROFNLGAT,FICEGAT,FLSGAT,G0SLGAT,
    %      9 EXPW, DTEMP, TKELAK, DELU, GRED, RHOMIX,   
    %      A FSVHGAT, FSIHGAT, FDLGAT, ULGAT, VLGAT, TAGAT, QAGAT,
    %      B RHOAGAT, PADRLGAT, PRESGAT, CSZGAT, ZRFMGAT, ZRFHGAT, 
    %      C ZDMLGAT,ZDHLGAT,RPCPLGAT,TRPCLGAT,SPCPLGAT,TSPCLGAT,RHSILGAT,
    %      >   RADJLGAT,
    %      D ASVLGAT,ASILGAT,FSDBLGAT, FSFBLGAT, FSSBLGAT, REFLGAT, BCSNLGAT,
    %      E ILG, JL1, JL2, NLAT, NLAKMAX, ISLFD, IZREF, ITG,
    %      F IALS, NBS, ISNOALB, IGL, IRSTRT,
    %      G N, IYEAR, IDAY, IHOUR, IMIN, TSED  )
    Classlin(i,:)=[HLAKGAT(1), LLAKGAT(1), BLAKGAT(1), NLAKGAT(1), TLAKGAT(1),...   
            T0LAK(1), HDPTH(1), LKICEH(1), SNICEH(1), ROFICEH(1),...
            SNOLGAT(1),RHOSLGAT(1),TSNOLGAT(1),ALBSLGAT(1),WSNOLGAT(1),...
            CDHGAT(1),CDMGAT(1),QSENLGAT(1),TFXLGAT(1),QEVPLGAT(1),QFSLGAT(1),QFXLGAT(1),...
            PETLGAT(1), EFLGAT(1), GTLGAT(1), QGLGAT(1), DRLGAT(1),...
            SFTLGAT(1),SFULGAT(1),SFVLGAT(1),SFQLGAT(1),SFHLGAT(1),QLWOLGAT(1),ALVLGAT(1),ALILGAT(1),...
            FSGLGAT(1), FLGLGAT(1), HFSLGAT(1), HEVLGAT(1), HMFLGAT(1), HTCLGAT(1),...
            FSGSLGAT(1), FLGSLGAT(1), HFSSLGAT(1), HEVSLGAT(1), HMFNLGAT(1), HTCSLGAT(1),...
            PCPLGAT(1),PCPNLGAT(1),QFLGAT(1),QFNLGAT(1),ROFNLGAT(1),FICEGAT(1),FLSGAT(1),G0SLGAT(1),...
            EXPW(1), DTEMP(1), TKELAK(1), DELU(1), GRED(1), RHOMIX(1),...   %FSVHGAT(1), FSIHGAT(1), 
            FDLGAT(1), ULGAT(1), VLGAT(1), TAGAT(1), QAGAT(1),...
            RHOAGAT(1), PADRLGAT(1), PRESGAT(1), CSZGAT(1), ZRFMGAT(1), ZRFHGAT(1),... 
            ZDMLGAT(1),ZDHLGAT(1),TRPCLGAT(1),TSPCLGAT(1),RHSILGAT(1), RADJLGAT(1),...%RPCPLGAT(1),SPCPLGAT(1),
            ASVLGAT(1),ASILGAT(1),FSDBLGAT(1), FSFBLGAT(1), FSSBLGAT(1), REFLGAT(1), BCSNLGAT(1),...
            NLAT(1), NLAKMAX(1), ISLFD(1), IZREF(1), ITG(1),...ILG(1), JL1(1), JL2(1),
            IALS(1), NBS(1), ISNOALB(1), IGL(1), ...IRSTRT(1),
            TSED(1)];
        % ---------------------------
    [HLAKGAT, LLAKGAT, BLAKGAT, NLAKGAT, TLAKGAT,...   
            T0LAK, HDPTH, LKICEH, SNICEH, ROFICEH,...
            SNOLGAT,RHOSLGAT,TSNOLGAT,ALBSLGAT,WSNOLGAT,...
            CDHGAT,CDMGAT,QSENLGAT,TFXLGAT,QEVPLGAT,QFSLGAT,QFXLGAT,...
            PETLGAT, EFLGAT, GTLGAT, QGLGAT, DRLGAT,...
            SFTLGAT,SFULGAT,SFVLGAT,SFQLGAT,SFHLGAT,QLWOLGAT,ALVLGAT,ALILGAT,...
            FSGLGAT, FLGLGAT, HFSLGAT, HEVLGAT, HMFLGAT, HTCLGAT,...
            FSGSLGAT, FLGSLGAT, HFSSLGAT, HEVSLGAT, HMFNLGAT, HTCSLGAT,...
            PCPLGAT,PCPNLGAT,QFLGAT,QFNLGAT,ROFNLGAT,FICEGAT,FLSGAT,G0SLGAT,...
            EXPW, DTEMP, TKELAK, DELU, GRED, RHOMIX,...   %FSVHGAT, FSIHGAT, 
            FDLGAT, ULGAT, VLGAT, TAGAT, QAGAT,...
            RHOAGAT, PADRLGAT, PRESGAT, CSZGAT, ZRFMGAT, ZRFHGAT,... 
            ZDMLGAT,ZDHLGAT,TRPCLGAT,TSPCLGAT,RHSILGAT, RADJLGAT,...%RPCPLGAT,SPCPLGAT,
            ASVLGAT,ASILGAT,FSDBLGAT, FSFBLGAT, FSSBLGAT, REFLGAT, BCSNLGAT,...
            NLAT, NLAKMAX, ISLFD, IZREF, ITG,...ILG, JL1, JL2,
            IALS, NBS, ISNOALB, IGL, ...IRSTRT,
            TSED] = ...N, IYEAR, IDAY, IHOUR, IMIN,
        CLASSL(HLAKGAT, LLAKGAT, BLAKGAT, NLAKGAT, TLAKGAT,...   
            T0LAK, HDPTH, LKICEH, SNICEH, ROFICEH,...
            SNOLGAT,RHOSLGAT,TSNOLGAT,ALBSLGAT,WSNOLGAT,...%%CDHGAT,CDMGAT,QSENLGAT,TFXLGAT,QEVPLGAT,QFSLGAT,QFXLGAT,PETLGAT, EFLGAT, GTLGAT, QGLGAT, DRLGAT,...
            SFTLGAT,SFULGAT,SFVLGAT,SFQLGAT,SFHLGAT,...%QLWOLGAT,ALVLGAT,ALILGAT,FSGLGAT, FLGLGAT, HFSLGAT, HEVLGAT, HMFLGAT, HTCLGAT,FSGSLGAT, FLGSLGAT, HFSSLGAT, HEVSLGAT, HMFNLGAT, HTCSLGAT,PCPLGAT,PCPNLGAT,...
            QFLGAT,...%QFNLGAT,ROFNLGAT,FICEGAT,FLSGAT,G0SLGAT,...            
            EXPW, DTEMP, TKELAK, DELU, GRED, RHOMIX,...   
            FSVHGAT, FSIHGAT, FDLGAT, ULGAT, VLGAT, TAGAT, QAGAT,...
            RHOAGAT, PADRLGAT, PRESGAT, CSZGAT, ZRFMGAT, ZRFHGAT,... 
            ZDMLGAT,ZDHLGAT,RPCPLGAT,TRPCLGAT,SPCPLGAT,TSPCLGAT,RHSILGAT, RADJLGAT,...
            ASVLGAT,ASILGAT,FSDBLGAT, FSFBLGAT, FSSBLGAT, REFLGAT, BCSNLGAT,...
            ILG, JL1, JL2, NLAT, NLAKMAX, ISLFD, IZREF,...
            ITG, IALS, NBS, ISNOALB, IGL, IRSTRT, N, IYEAR, IDAY, IHOUR, IMIN,...
            TSED,TFREZ,RHOICE,RHOW,HCPW,SPHICE,CLHMLT,VMIN,...);  %these were not passed because they were global in MacKays 
            TKEMIN,DELZLK,GRAV,DELSKIN,VKC,CLHVAP,SPHAIR,EMSW,SBC,TCW,DELT,...
            TKECN,TKECF,TKECE,TKECS,TKECL,SPHW,DHMAX,HDPTHMIN,DUMAX,DELMAX,...
            DELMIN,TCICE,HCPICE,CGRAV,CPD,CKARM,...
            fID_71,fID_72,fID_73,fID_74,fID_75,fID_76,fID_77,fID_82);
    Classlout(i,:)=[HLAKGAT(1), LLAKGAT(1), BLAKGAT(1), NLAKGAT(1), TLAKGAT(1),...   
            T0LAK(1), HDPTH(1), LKICEH(1), SNICEH(1), ROFICEH(1),...
            SNOLGAT(1),RHOSLGAT(1),TSNOLGAT(1),ALBSLGAT(1),WSNOLGAT(1),...%%CDHGAT(1),CDMGAT(1),QSENLGAT(1),TFXLGAT(1),QEVPLGAT(1),QFSLGAT(1),QFXLGAT(1),PETLGAT(1), EFLGAT(1), GTLGAT(1), QGLGAT(1), DRLGAT(1),...
            SFTLGAT(1),SFULGAT(1),SFVLGAT(1),SFQLGAT(1),SFHLGAT(1),...%QLWOLGAT(1),ALVLGAT(1),ALILGAT(1),FSGLGAT(1), FLGLGAT(1), HFSLGAT(1), HEVLGAT(1), HMFLGAT(1), HTCLGAT(1),FSGSLGAT(1), FLGSLGAT(1), HFSSLGAT(1), HEVSLGAT(1), HMFNLGAT(1), HTCSLGAT(1),PCPLGAT(1),PCPNLGAT(1),...
            QFLGAT(1),...%QFNLGAT(1),ROFNLGAT(1),FICEGAT(1),FLSGAT(1),G0SLGAT(1),...            
            EXPW(1), DTEMP(1), TKELAK(1), DELU(1), GRED(1), RHOMIX(1),...   
            FSVHGAT(1), FSIHGAT(1), FDLGAT(1), ULGAT(1), VLGAT(1), TAGAT(1), QAGAT(1),...
            RHOAGAT(1), PADRLGAT(1), PRESGAT(1), CSZGAT(1), ZRFMGAT(1), ZRFHGAT(1),... 
            ZDMLGAT(1),ZDHLGAT(1),RPCPLGAT(1),TRPCLGAT(1),SPCPLGAT(1),TSPCLGAT(1),RHSILGAT(1), RADJLGAT(1),...
            ASVLGAT(1),ASILGAT(1),FSDBLGAT(1), FSFBLGAT(1), FSSBLGAT(1), REFLGAT(1), BCSNLGAT(1),...
            ILG(1), JL1(1), JL2(1), NLAT(1), NLAKMAX(1), ISLFD(1), IZREF(1), ITG(1),...
            IALS(1), NBS(1), ISNOALB(1), IGL(1), IRSTRT(1),...
            N(1), IYEAR(1), IDAY(1), IHOUR(1), IMIN(1), TSED(1),...
            TFREZ(1),RHOICE(1),RHOW(1),HCPW(1),SPHICE(1),CLHMLT(1),VMIN(1),...);  %these were not passed because they were global in MacKays 
            TKEMIN(1),DELZLK(1),GRAV(1),DELSKIN(1),VKC(1),CLHVAP(1),SPHAIR(1),EMSW(1),SBC(1),TCW(1),DELT(1),...
            TKECN(1),TKECF(1),TKECE(1),TKECS(1),TKECL(1),SPHW(1),DHMAX(1),HDPTHMIN(1),DUMAX(1),DELMAX(1),...
            DELMIN(1),TCICE(1),HCPICE(1),CGRAV(1),CPD(1),CKARM(1)];
    % C
    % C========================================================================
    % C CHECK ENERGY BALANCE - end of timestep
    % C   
    
    for I=JL1:JL2 %DO 325
        ICEBOT=RHOIW*LKICEH(I);
        ICETOP=LKICEH(I)-ICEBOT;
        if (ICEBOT >= DELSKIN)
            HCAP=HCPICE;
        elseif (LKICEH(I) <= 0.0)
            HCAP=HCPW;
        else
            HCAP=(LKICEH(I)*HCPICE + (DELSKIN-ICEBOT)*HCPW)/DELSKIN;
        end
        CTLSTP(I)= CTLSTP(I) + HCAP*T0LAK(I)*DELSKIN;
        %DO 335, J=1,NLAKGAT(I)
        for J=1:NLAKGAT(I)
            ZTOP=DELSKIN + DELZLK*(J -1);
            ZBOT=DELSKIN + DELZLK*J;
            if (ICEBOT >= ZBOT)
                HCAP=HCPICE;
            elseif (ICEBOT <= ZTOP)
                HCAP=HCPW;
            else
                Z=ICEBOT-ZTOP;
                HCAP=(Z*HCPICE + (DELZLK-Z)*HCPW)/DELZLK;
            end
            CTLSTP(I)=CTLSTP(I) + HCAP*TLAKGAT(I,J)*DELZLK;
        %335     CONTINUE
        end
        CTLSTP(I)=CTLSTP(I)/DELT;
        QSUML=FSGLGAT(I)+FLGLGAT(I)-HFSLGAT(I)-HEVLGAT(I)-HMFLGAT(I);
        %Write to file, but for now just print to console
        %This is LAKE.of9 file
        %        WRITE(79,6019) IYEAR,IDAY,IHOUR,IMIN,QSUML-CTLSTP(I),QSUML
        %6019    FORMAT(I4,1X,3(I3,1X),2F8.2)
        if rem(N,1000)==0
            newTime=now;
            deltaT=(newTime-oldTime)*60*24;
            oldTime=newTime;
            elapseTime=(newTime-firstTime)*60*24;
            fprintf('\n%#5d,  %2.2f,  %#5.3f,  %#5.1f, %4d,  %3d,    %2d,   %2d, %8.2f, %8.2f',N,N/LngthOfMetDat,deltaT,elapseTime,IYEAR,IDAY,IHOUR,IMIN,QSUML-CTLSTP(I),QSUML);
        end
        %fprintf('\nN, Nfrac, IYEAR, IDAY, IHOUR, IMIN, QSUML-CTLSTP(I), QSUML');
        %CCC     IF(ABS(CTLSTP(I)-QSUML).GE. 10.0) THEN
        %CCC           WRITE(6,6440) N,DAY,CTLSTP(I),QSUML,QSUML-CTLSTP(I)
        %CCC6440          FORMAT(2X,'LAKE ENERGY BALANCE  ',I8,F8.3,2F20.8,F8.2)
        %Cmdm          WRITE(6,6450) FSGLGAT(I),FLGLGAT(I),HFSLGAT(I),
        %Cmdm 1             HEVLGAT(I),HMFLGAT(I),LKICEH(I),T0LAK(I)
        %Cmdm6450          FORMAT(2X,7F15.6)
        %Ctest         STOP
        %CCC     ENDIF
    %325   CONTINUE
    end
    %fprintf('\n');


    % C
    % C========================================================================
    % C SCATTER CODE FROM CLASSS 
    % C-- LAKE DIAGNOSTIC VARIABLES SPLIT OUT                     
    for K=1:NMW %DO 701 K=1,NMW                                                
        HLAKROT(IWMOS(K),JWMOS(K))=HLAKGAT(K+NML);                 
        LLAKROT(IWMOS(K),JWMOS(K))=LLAKGAT(K+NML);                
        BLAKROT(IWMOS(K),JWMOS(K))=BLAKGAT(K+NML);                 
        NLAKROT(IWMOS(K),JWMOS(K))=NLAKGAT(K+NML);                 
        for L=1:NLAKGAT(K+NML)%DO 706 L=1,NLAKGAT(K+NML)                                 
            TLAKROT(IWMOS(K),JWMOS(K),L)=TLAKGAT(K+NML,L);           
        end%706       CONTINUE                                                  
        ASVDROT(IWMOS(K),JWMOS(K))=ASVLGAT(K+NML);
        ASIDROT(IWMOS(K),JWMOS(K))=ASILGAT(K+NML);
        SNOLROT(IWMOS(K),JWMOS(K))=SNOLGAT(K+NML);
        RHOSLROT(IWMOS(K),JWMOS(K))=RHOSLGAT(K+NML);
        TSNOLROT(IWMOS(K),JWMOS(K))=TSNOLGAT(K+NML);
        ALBSLROT(IWMOS(K),JWMOS(K))=ALBSLGAT(K+NML);
        WSNOLROT(IWMOS(K),JWMOS(K))=WSNOLGAT(K+NML);
    end%701   CONTINUE                                                      
    % C
    % C-- FLUX DIAGNOSTIC VARIABLES SPLIT OUT                     
    for K=1:NMW%DO 385 K=1,NMW                                                
        CDHROT (IWMOS(K),JWMOS(K))=CDHGAT (K+NML);                 
        CDMROT (IWMOS(K),JWMOS(K))=CDMGAT (K+NML);                 
        HFSLROT (IWMOS(K),JWMOS(K))=HFSLGAT (K+NML);               
        HEVLROT (IWMOS(K),JWMOS(K))=HEVLGAT(K+NML);                
        FSGLROT (IWMOS(K),JWMOS(K))=FSGLGAT(K+NML);                
        FLGLROT (IWMOS(K),JWMOS(K))=FLGLGAT(K+NML);                
        HMFLROT (IWMOS(K),JWMOS(K))=HMFLGAT(K+NML);                
        PCPLROT (IWMOS(K),JWMOS(K))=PCPLGAT(K+NML);                
        PCPNLROT (IWMOS(K),JWMOS(K))=PCPNLGAT(K+NML);                
        QFLROT (IWMOS(K),JWMOS(K))=QFLGAT(K+NML);                
        QFNLROT (IWMOS(K),JWMOS(K))=QFNLGAT(K+NML);                
        ROFNLROT (IWMOS(K),JWMOS(K))=ROFNLGAT(K+NML);                
        HFSSLROT (IWMOS(K),JWMOS(K))=HFSSLGAT (K+NML);               
        HEVSLROT (IWMOS(K),JWMOS(K))=HEVSLGAT(K+NML);                
        FSGSLROT (IWMOS(K),JWMOS(K))=FSGSLGAT(K+NML);                
        FLGSLROT (IWMOS(K),JWMOS(K))=FLGSLGAT(K+NML);                
        HMFNLROT (IWMOS(K),JWMOS(K))=HMFNLGAT(K+NML);                
        HTCSLROT (IWMOS(K),JWMOS(K))=HTCSLGAT(K+NML);                
        HTCLROT (IWMOS(K),JWMOS(K))=HTCLGAT(K+NML);                
        SFTLROT (IWMOS(K),JWMOS(K))=SFTLGAT(K+NML);                
        SFULROT (IWMOS(K),JWMOS(K))=SFULGAT(K+NML);                
        SFVLROT (IWMOS(K),JWMOS(K))=SFVLGAT(K+NML);                
        SFQLROT (IWMOS(K),JWMOS(K))=SFQLGAT(K+NML);                
        SFHLROT (IWMOS(K),JWMOS(K))=SFHLGAT(K+NML);                
        QLWOLROT (IWMOS(K),JWMOS(K))=QLWOLGAT(K+NML);                
        ALVLROT (IWMOS(K),JWMOS(K))=ALVLGAT(K+NML);                
        ALILROT (IWMOS(K),JWMOS(K))=ALILGAT(K+NML);                
        PETLROT (IWMOS(K),JWMOS(K))=PETLGAT(K+NML);                
        EFLROT (IWMOS(K),JWMOS(K))=EFLGAT(K+NML);                
        GTLROT (IWMOS(K),JWMOS(K))=GTLGAT(K+NML);                
        QGLROT (IWMOS(K),JWMOS(K))=QGLGAT(K+NML);                
        DRLROT (IWMOS(K),JWMOS(K))=DRLGAT(K+NML);                
        QSENLROT (IWMOS(K),JWMOS(K))=QSENLGAT(K+NML);                
        TFXLROT (IWMOS(K),JWMOS(K))=TFXLGAT(K+NML);                
        QEVPLROT (IWMOS(K),JWMOS(K))=QEVPLGAT(K+NML);                
        QFSLROT (IWMOS(K),JWMOS(K))=QFSLGAT(K+NML);                
        QFXLROT (IWMOS(K),JWMOS(K))=QFXLGAT(K+NML);                
    end%385   CONTINUE                                                      
    % C
    % C-- DAY COUNTER
    NCOUNT=NCOUNT+1;
    if(NCOUNT > NDAY)
        NCOUNT=1;
    end
    % C=======================================================================
    % C
end%GO TO 200

fclose(fID_51);

csvwrite('ClassIn.csv',Classlin);
csvwrite('ClassOut.csv',Classlout);

%999   CONTINUE
% C
%STOP
%END
fprintf('\nfin\n');


end%end of RunLake
