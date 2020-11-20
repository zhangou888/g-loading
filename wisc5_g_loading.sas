dm 'output' clear; dm 'log' clear;
proc datasets kill;run;quit;

libname x "C:\Users\UZHANOU\Documents\WISC4\Stdz\Data";
libname k "C:\Users\UZHANOU\Documents\WISC5\stdz\data";

/* wisc4 */
data wisc4;set x.wisc4_final;if stdz = "Y";
	keep cfid stdz agegrp 
         siss vcss coss inss wcss 
         bdss pcss mrss psss dsss 
         lnsss arss cdss ssss cass;
run;

/* wisc5 */
data wisc5;set k.wisc5_stdz_final;if stdz='Y'; 
	keep teid wisc5_in_ss wisc5_si_ss wisc5_vc_ss  
     	    wisc5_co_ss 
            wisc5_bd_ss /*use bonus plan B 12/30/2013 */
      	    wisc5_vp_ss wisc5_mr_ss  wisc5_fw_ss wisc5_pc_ss 
      	    wisc5_ds_ss 
            wisc5_ps_ss /*use ps1 12/30/2013 */
            wisc5_ln_ss wisc5_ar_ss 
       	    wisc5_ss_ss   
       	    wisc5_cd_ss 
            wisc5_ca_ss 
         	wisc5_pa1_ss   
         	wisc5_pad_ss
            wisc5_pac_ss;
run;

data wisc5;set wisc5;
 /*rename variables to no more than 6 characters in length*/
 	array x wisc5_in_ss wisc5_si_ss wisc5_vc_ss  
     	    wisc5_co_ss 
            wisc5_bd_ss /*use bonus plan B 12/30/2013 */
      	    wisc5_vp_ss wisc5_mr_ss  wisc5_fw_ss wisc5_pc_ss 
      	    wisc5_ds_ss 
            wisc5_ps_ss /*use ps1 12/30/2013 */
            wisc5_ln_ss wisc5_ar_ss 
       	    wisc5_ss_ss   
       	    wisc5_cd_ss 
            wisc5_ca_ss 
         	wisc5_pa1_ss   
         	wisc5_pad_ss
            wisc5_pac_ss; 

 	array y inss siss vcss coss bdss 
        	vpss mrss fwss pcss 
            dsss psss lnss arss 
            ssss cdss cass /*rnrss rnqss*/ 
            pass padss pacss;
  	do over x; y=x;end;
run;

/* wisc4 */
ods output lineqseqstd =estimate0 fit=model0 lineqsvarexogstd=ee0;
proc calis data = wisc4 outstat=stat0 pall;
var siss vcss coss inss wcss 
         bdss pcss mrss psss dsss 
         lnsss arss cdss ssss cass;
lineqs siss=b1 f1+e1 , 
       vcss=b2 f1+e2 , 
       inss=b3 f1+e3 , 
       coss=b4 f1+e4 ,
       wcss=b5 f1+e5 , 
	   
       bdss=b6 f2+e6 , 
       mrss=b7 f2+e7 ,
       pcss=b8 f2+e8 ,
       psss=b9 f2+e9 ,
  
       dsss=b10 f3+e10 , 
       lnsss=b11 f3+e11 ,

       cdss=b12 f4+e12 ,
       ssss=b13 f4+e13 , 
       cass=b14 f4+e14 , 

	   arss=b15 f5+e15,

       f1=b16 f6+e16,
       f2=b17 f6+e17, 
       f3=b18 f6+e18, 
       f4=b19 f6+e19,
       f5=b20 f6+e20;

std e1-e20 = the1-the20, f6=1;
bounds the1-the20>0; 
run;

data estimate1; set estimate0; if variable notin ('Std Err' 't Value');
	keep variable coefficient1 parameter1 parameter2;
run;


data estimate1; set estimate1 ;
	format factor1 $8.;
    factor1=parameter1;
run;

data estimate1; set estimate1 (rename=(coefficient1=varload1 parameter2=par2));
oo=_n_;run;
proc sort; by par2;run;

data estimate1; merge estimate1 ee0; by par2;run; 

data estimate1; set estimate1; drop par2; 
proc sort; by oo;run;

data factorest1; set estimate1;keep variable varload1 ;
data factorest1; set factorest1 (rename=(varload1=cc1));
rename variable=factor1;run;

data factorest1;set factorest1;run;

proc sort; by factor1;run;

data estimate1; set estimate1; Asterisk1='*'; plus1='+';equal='=';run;
proc sort; by  factor1;run;

data est1; merge estimate1(in=a) factorest1; by factor1; if a;run;

data est1; set est1;
 gloading=round(varload1*cc1,.0001);run;

data wisc4_g; retain variable varload1 Asterisk1 factor1 plus1 e gloading equal; set est1;
	keep variable varload1 Asterisk1 factor1 variable e gloading plus1 equal;
data wisc4_g; set wisc4_g; format model $20.; model="wisc4_5factor";run;



/* wisc5 */
ods output lineqseqstd =estimate2 lineqseq =estimate22 fit=model2 lineqsvarexogstd=ee2;
proc calis data=wisc5 outstat=stat1 pall;
var inss siss vcss coss bdss vpss mrss pcss fwss arss dsss psss lnss ssss cdss cass;
lineqs siss=b1 f1+e1 , 
       vcss=b2 f1+e2 , 
       inss=b3 f1+e3 , 
       coss=b4 f1+e4 , 
	   
       bdss=b5 f2+e5 , 
       vpss=b6 f2+e6 ,
 
       mrss=b7 f3+e7 ,
       fwss=b8 f3+e8 ,
       pcss=b9 f3+e9 ,
       
       dsss=b11 f4+e11 , 
       psss=b12 f4+e12 , 
       lnss=b13 f4+e13 ,

       arss=b10 f1 + b23 f3 + b14 f4+e10,

       cdss=b15 f5+e14 ,
       ssss=b16 f5+e15 , 
       cass=b17 f5+e16 , 

       f1=b18 f6+e17,
       f2=b19 f6+e18, 
       f3=b20 f6+e19, 
       f4=b21 f6+e20, 
       f5=b22 f6+e21;

std e1-e21 = the1-the21, f6=1;
run;

/* Estimation */
data estimate3; set estimate2; if variable notin ('Std Err' 't Value');
	keep variable coefficient1 parameter1 coefficient2 parameter2 coefficient3 parameter3 parameter4;
run;

data estimate3; set estimate3;
 if variable = 'arss' then do;
	Coefficient2 = Coefficient2;
	Coefficient3 = Coefficient3;
	parameter2 = parameter2;
	parameter3 = parameter3;
	parameter4 = parameter4;
 end;
 if variable ^= 'arss' then do;
 	  Coefficient2 = .;
      Coefficient3 = .;
	  parameter4=parameter2;
	  parameter2=.;
	  parameter3=.;
 end;
run;

data estimate3; set estimate3;
	format factor1 $8.; 
	format factor2 $8.;
	format factor3 $8.;
    factor1=parameter1;
	factor2=parameter2;
	factor3=parameter3;
run;

data estimate3; set estimate3 (rename=(coefficient1=varload1 
                                       coefficient2=varload2 
                                       coefficient3=varload3
                                       parameter4=par2));
	oo=_n_;
run;


proc sort; by par2;run;

data estimate3; merge estimate3 ee2; by par2;run; 

data estimate3; set estimate3; drop par2; 
proc sort; by oo;run;

data factorest2; set estimate3;keep variable varload1 varload2 varload3 ;
data factorest2; set factorest2 (rename=(varload1=cc1 varload2=cc2 varload3=cc3));run;
data factorest2; set factorest2; rename variable=factor1;run;

data factorest2;set factorest2;run;

proc sort; by factor1;run;

data estimate4; set estimate3; Asterisk1='*';Asterisk2='*';Asterisk1='*'; plus1='+';plus2='+';plus3='+';equal='=';run;
proc sort; by  factor1;run;

data est2; merge estimate4(in=a) factorest2; by factor1; if a;run;
/* manual adjustment-for the 2nd order factor f6 on f4*/
data est2; set est2;
	if variable = 'arss' then cc2 = 1.00;
	if variable = 'arss' then cc3 = .8073;
run;

data est2; set est2;
 if variable = 'arss' then gloading=round(varload1*cc1+varload2*cc2+varload3*cc3,.0001);
 else gloading=round(varload1*cc1,.0001);run;

data wisc5_g; retain variable varload1 Asterisk1 factor1 plus1 varload2 Asterisk2 factor2 plus2 varload3 Asterisk3 factor3 plus3 e gloading equal; set est2;
	keep variable varload1 Asterisk1 factor1 variable varload2 Asterisk2 factor2 varload3 Asterisk3 factor3 e gloading plus1 plus2 plus3 equal;
data wisc5_g; set wisc5_g; format model $20.; model="wisc5_g";run;





/** Output **/
ods rtf file="C:\Users\UZHANOU\Documents\WISC5\stdz\results\factor_analysis\g_loading\WISC4 Stdz CFA Results.rtf";
proc report data=wisc4_g nowindows headline split='/' spacing=1 ;
title2 "G-Loading on wisc4 model";
column model variable equal varload1 Asterisk1 factor1 plus1 e gloading ;
define model/order=internal group '' width=20 center;
define variable/width=10 ' ' center order=data;
define equal/ width=1 ' ' center ;
define varload1/ ' ' width=10 f=4.2 center ;
define asterisk1/ width=1  ' ' center;
define factor1/width=4 ' ' center ;
define plus1/  width=1  ' ' center ;
define e/  width=5 f=4.2 ' ' center ;
define gloading/'G-Loading' format=3.2 width=10  center;
run;quit;
ods rtf close;


ods rtf file="C:\Users\UZHANOU\Documents\WISC5\stdz\results\factor_analysis\g_loading\WISC5 Stdz CFA Results.rtf";
proc report data=wisc5_g nowindows headline split='/' spacing=1 ;
title2 "G-Loading on wisc5 model";
column model variable equal varload1 Asterisk1 factor1 plus1 varload2 Asterisk2 factor2 plus2 varload3 Asterisk3 factor3 plus3 e gloading ;
define model/order=internal group '' width=20 center;
define variable/width=10 ' ' center order=data;
define equal/ width=1 ' ' center ;
define varload1/ ' ' width=10 f=4.2 center ;
define asterisk1/ width=1  ' ' center;
define factor1/width=4 ' ' center ;
define plus1/  width=1  ' ' center ;
define varload2/ ' ' width=10 f=4.2 center ;
define asterisk2/ width=1  ' ' center;
define factor2/width=4 ' ' center ;
define plus2/  width=1  ' ' center ;
define varload3/ ' ' width=10 f=4.2 center ;
define asterisk3/ width=1  ' ' center;
define factor3/width=4 ' ' center ;
define plus3/  width=1  ' ' center ;
define e/  width=5 f=4.2 ' ' center ;
define gloading/'G-Loading' format=3.2 width=10  center;
run;quit;
ods rtf close;
