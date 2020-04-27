%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%   Copyright   Xin-Guang Zhu, Yu Wang, Donald R. ORT and Stephen P. LONG
%CAS-MPG Partner Institute for Computational Biology, Shanghai Institutes for Biological Sciences, CAS, Shanghai,200031 
%China Institute of Genomic Biology and Department of Plant Biology, Shanghai Institutes for Biological Sciences, CAS, Shanghai,200031 
%University of Illinois at Urbana Champaign
%Global Change and Photosynthesis Research Unit, USDA/ARS, 1406 Institute of Genomic Biology, Urbana, IL 61801, USA.
 
%   This file is part of e-photosynthesis.
 
%    e-photosynthesis is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; 
 
%    e-photosynthesis is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
 
%    You should have received a copy of the GNU General Public License (GPL)
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 function SUCS_Con = SUCS_Ini(begin)
global SUCRatio;
global SUCS_OLD_TIME;
global SUCS_TIME_N;
global SUCS_VEL;
global SUCS_CON;

SUCS_OLD_TIME = 0;
SUCS_TIME_N = 1;

SUCS_VEL = zeros(1,3);    % Clean memory
SUCS_CON = zeros(3,1);    % Clean memory

 

global 	KE501	;
global 	Km511	;
global 	Km512	;
global 	Km513	;
global 	Km514	;
global  KE51;
global 	Km521	;
global 	KI521	;
global 	KI522	;
global 	KI523	;
global 	KE52	;
global 	KE531	;
global 	KE541	;
global 	Km551	;
global 	Km552	;
global 	Km553	;
global 	Km554	;
global 	KE55	;
global 	Km561	;
global 	Km562	;
global 	KI561	;
global 	KI562	;
global 	KI563	;
global 	KI564	;
global 	KI565	;
global 	KE56	;
global 	Km571	;
global 	Ki572	;
global 	KE57	;
global 	Km581	;
global 	KI581	;
global 	KI582	;
global 	Km591	;
global 	Km592	;
global 	Km593	;
global 	KI591	;
global 	KI592	;
global 	KE59	;
global 	Km601	;
global 	Km602	;
global 	Km603	;
global 	Km604	;
global 	KE60	;
global 	KE61	;
global 	Km621	;


KE501	=	1/0.05*SUCRatio(16)	;	 
Km511	=	0.02*SUCRatio(17)	;	 
Km512	=	0.3*SUCRatio(18)	;	 
Km513	=	0.4*SUCRatio(19)	;	 
KE51	=	12*SUCRatio(20);      
Km514	=	0.014*SUCRatio(21)	; 
Km521	=	0.0025*SUCRatio(22)	;	 
KI521	=	0.7*SUCRatio(23)	;	 
KI522	=	12*SUCRatio(24)	;	 
KI523	=	7*10^(-5)*SUCRatio(25);	 
KE52	=	6663*SUCRatio(26)	;	 
KE531	=	2.3*SUCRatio(27)	;	 
KE541	=	0.0584*SUCRatio(28)	;	 
Km551	=	0.14*SUCRatio(29)	;	 
Km552	=	0.1*SUCRatio(30)	;	 
Km553	=	0.11*SUCRatio(31)	;	 
Km554	=	0.12*SUCRatio(32)	;	 
KE55	=	0.31*SUCRatio(33)	; 
Km561	=	0.8*SUCRatio(34)	;	 
Km562	=	2.4*SUCRatio(35)	;	 
KI561	=	0.7*SUCRatio(36);	 
KI562	=	0.8*SUCRatio(37)	;	 
KI563	=	0.4*SUCRatio(38)	;	 
KI564	=	11*SUCRatio(39)	;	 
KI565	=	50*SUCRatio(40)	;	 
KE56	=	10*SUCRatio(41)	;	 
Km571	=	0.35*SUCRatio(42);	 
Ki572	=	10*SUCRatio(43)	;	 
KE57	=	780*SUCRatio(44)	;	 
Km581	=	0.032*SUCRatio(45);	 
KI581	=	0.1*SUCRatio(46)	;	 
KI582	=	0.5*SUCRatio(47)	;	 
Km591	=	0.5*SUCRatio(48)	;	 
Km592	=	0.021*SUCRatio(49);	 
Km593	=	0.5*SUCRatio(50)	;	 
KI591	=	0.16*SUCRatio(51);	 
KI592	=	0.7*SUCRatio(52)	;	 
KE59	=	590*SUCRatio(53)	;	 
Km601	=	0.042*SUCRatio(54);	 
Km602	=	1.66*SUCRatio(55);	 
Km603	=	0.28*SUCRatio(56);	 
Km604	=	16*SUCRatio(57)	;	 
KE60	=	16*SUCRatio(58)	;	 
KE61	=	1.2*107*SUCRatio(59);	 
Km621	=	5*SUCRatio(60)	;	 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization of the initial concentration of the different component %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize the leaves of active enzyme in a dark adapted leaves;								
%	mM
    
T3Pc	=	2	;
FBPc	=	2	;
HexPc	=	5.8	;
F26BPc	=	7.8*10^(-6)	;
ATPc	=	0.4	;
ADPc	=	0.4	;
UDPGc	=	0.57	;
UTPc	=	0.75	;
SUCP	=	0;
SUC	    =	0;
PGAc    =  0.5; 


% Assign value to a variable that is transferred to the program				

SUCS_Con = zeros(3,1);

SUCS_Con	(	1		)	=	T3Pc	;
SUCS_Con	(	2		)	=	FBPc	;
SUCS_Con	(	3		)	=	HexPc	;
SUCS_Con	(	4		)	=	F26BPc	;
SUCS_Con	(	5		)	=	ATPc	;
SUCS_Con	(	6		)	=	ADPc	;
% SUCS_Con	(	7		)	=	OPOPc	;
SUCS_Con	(	8		)	=	UDPGc	;
SUCS_Con	(	9		)	=	UTPc	;
SUCS_Con	(	10		)	=	SUCP	;
SUCS_Con	(	11		)	=	SUC	;
SUCS_Con	(	12		)	=	PGAc	;


% The following calculate the total concentration of different enzymes. 

	SC = 10; 
	SC1 = 1; 

    
global GP; 
if GP==0
	global 	V51	;
	global 	V52	;
	global 	V55	;
	global 	V56	;
	global 	V57	;
	global 	V58	;
	
	% Unit: mmol l-1 s-1;
	
	V51	=	0.107376831	* SC*SUCRatio(1)	;%	DHAP+GAP --FBP          % default 0.5
	V52	=	0.063979048	* SC*SUCRatio(2)	;	%	FBP --F6P + Pi
	V55	=	0.115403205	* SC*SUCRatio(3);	%	G1P+UTP --OPOP+UDPG 
	V56	=	0.055503446	* SC*SUCRatio(4)	;	%	UDPG+F6P--SUCP + UDP
	V57	=	0.55503446	* SC1*SUCRatio(5);	%	SUCP--Pi + SUC; 0.27 DEFALT
	V58	=	0.016819226	* SC*SUCRatio(6);	%	F26BP--F6P + Pi
end	
	global 	V59	;
	global  V60;
	global  V61;
	global	V62;	
	global  Vdhap_in;
	global  Vgap_in;
	global  Vpga_in;
    
	V59	=	0.03	* SC*SUCRatio(7);	    %	F6P + ATP --ADP + F26BP % defalut 0.03  (* 0.3)
	V60	=	6.1	*SUCRatio(8);	        %	ATP+UDP --UTP + ADP
	V61	=	10000;	        %	POPO --2PO
	V62	=	2	* SC1*SUCRatio(9);	        %	SUC Sink        0.9 works.
	Vdhap_in = 1.05* SC1*SUCRatio(10);        %   DHAP export from chloroplast
	Vgap_in  = 1.05* SC1*SUCRatio(11);        %   GAP export from chloroplast
	Vpga_in  =  1.05* SC1*SUCRatio(12);       %   PGA export from chloropalst

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here is some pool values      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ATc =  1.0*SUCRatio(13);      % mM
UTc =  1.5*SUCRatio(14);          % mM
PTc = 15*SUCRatio(15);          % 


global SUCS_Pool;		
SUCS_Pool	(	1	)	=	ATc;
SUCS_Pool	(	2	)	=	UTc;
SUCS_Pool	(	3	)	=	PTc;

