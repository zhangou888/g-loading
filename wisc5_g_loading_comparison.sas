/***********************************************/
/**    WISC4/WISC5 G-loading Comparison       **/
/**    Psychometrician: Ou Zhang              **/
/**    Date: 09-30-2014                       **/
/***********************************************/
dm 'output' clear; dm 'log' clear;
proc datasets kill;run;quit;

libname x "C:\Users\UZHANOU\Documents\WISC4\Stdz\Data";
libname k "C:\Users\UZHANOU\Documents\WISC5\stdz\data";
/* wisc4 */
data wisc4;set x.wisc4_final;if stdz = "Y";
	keep cfid stdz agegrp 
         siss vcss coss inss wcss pcss
         bdss mrss psss dsss 
         lnsss arss cdss ssss cass; 
		 *psss-picture concept;
		 *pcss-picture completion;
run;

data wisc4;set wisc4;
	rename psss=pcnss; *for model use only;
	rename pcss=pcmss;
run;

/* wisc5 */
data wisc5;set k.wisc5_stdz_final;if stdz='Y'; 
	keep teid agegrp wisc5_in_ss wisc5_si_ss wisc5_vc_ss  
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
        	vpss mrss fwss pcnss 
            dsss psss lnsss arss 
            ssss cdss cass /*rnrss rnqss*/ 
            pass padss pacss;
  	do over x; y=x;end;
run;



/************** 1st comparison **********************************/
/*   same model, same construct, different data */
/* wisc4 data */
ods output lineqseqstd =estimate0 fit=model0 lineqsvarexogstd=ee0;
proc calis data = wisc4 outstat=stat0 pall;
var siss vcss coss inss bdss pcnss mrss dsss 
    lnsss arss cdss ssss cass;
lineqs siss=b1 f1+e1 , 
       vcss=b2 f1+e2 , 
       inss=b3 f1+e3 , 
       coss=b4 f1+e4 , 
	   
       bdss=b5 f2+e5 , 
       mrss=b6 f2+e6 ,
       pcnss=b7 f2+e7 ,
  
       dsss=b8 f3+e8 , 
       lnsss=b9 f3+e9 ,

       cdss=b10 f4+e10 ,
       ssss=b11 f4+e11 , 
       cass=b12 f4+e12 , 
     
	   arss=b13 f1 + b14 f2 + b15 f3+e13,

       f1=b16 f6+e14,
       f2=b17 f6+e15, 
       f3=b18 f6+e16, 
       f4=b19 f6+e17,
       f5=b20 f6+e18;

std e1-e18 = the1-the18, f6=1;
bounds the1-the18>0; 
run;

/* manipulate estimation */
data estimate1; set estimate0; if variable notin ('Std Err' 't Value');
	keep variable coefficient1 parameter1 coefficient2 parameter2 coefficient3 parameter3 parameter4;
run;

data estimate1; set estimate1;
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

data estimate1; set estimate1;
	format factor1 $8.; 
	format factor2 $8.;
	format factor3 $8.;
    factor1=parameter1;
	factor2=parameter2;
	factor3=parameter3;
run;

data estimate1; set estimate1 (rename=(coefficient1=varload1 
                                       coefficient2=varload2 
                                       coefficient3=varload3
                                       parameter4=par2));
	oo=_n_;
run;


proc sort data=estimate1; by par2;run;

data ee2;set ee0;if Parameter ^="";
	rename Variable=par2;
	drop SpecType StdErr tValue;
run;

proc sort data=ee2;by par2;run;

data estimate1; merge estimate1 ee2; by par2;run; 

data estimate1; set estimate1; drop par2; 
proc sort data=estimate1; by oo;run;

data factorest1; set estimate1;keep variable varload1 varload2 varload3 ;
data factorest1; set factorest1 (rename=(varload1=cc1 varload2=cc2 varload3=cc3));run;
data factorest1; set factorest1; run;
proc sort data=factorest1; by variable;run;

data estimate2; set estimate1; Asterisk1='*';Asterisk2='*';Asterisk3='*'; plus1='+';plus2='+';plus3='+';equal='=';run;
proc sort data=estimate2; by variable;run;

data est2; merge estimate2(in=a) factorest1; by variable;if a;run;
/* manual adjustment-for the 2nd order factor f6 on f4*/
data est2; set est2;
	if variable = 'arss' then cc2 = .9167;
	if variable = 'arss' then cc3 = .8258;
run;

data est2; set est2;
 if variable = 'arss' then gloading=round(varload1*cc1+varload2*cc2+varload3*cc3,.0001);
 else gloading=round(varload1*cc1,.0001);run;

data est2;set est2;
	rename estimate=e;
run;

data wisc4_g; retain variable varload1 Asterisk1 factor1 plus1 varload2 Asterisk2 factor2 plus2 varload3 Asterisk3 factor3 plus3 e gloading equal; set est2;
	keep variable varload1 Asterisk1 factor1 variable varload2 Asterisk2 factor2 varload3 Asterisk3 factor3 e gloading plus1 plus2 plus3 equal;
data wisc4_g; set wisc4_g; format model $20.; model="wisc4_g";run;


*******************************************************************;
/* wisc5 data */
ods output lineqseqstd =estimate3 fit=model1 lineqsvarexogstd=ee3;
proc calis data = wisc5 outstat=stat1 pall;
var siss vcss coss inss bdss pcnss mrss dsss 
    lnsss arss cdss ssss cass;
lineqs siss=b1 f1+e1 , 
       vcss=b2 f1+e2 , 
       inss=b3 f1+e3 , 
       coss=b4 f1+e4 , 
	   
       bdss=b5 f2+e5 , 
       mrss=b6 f2+e6 ,
       pcnss=b7 f2+e7 ,
  
       dsss=b8 f3+e8 , 
       lnsss=b9 f3+e9 ,

       cdss=b10 f4+e10 ,
       ssss=b11 f4+e11 , 
       cass=b12 f4+e12 , 
     
	   arss=b13 f1 + b14 f2 + b15 f3+e13,

       f1=b16 f6+e14,
       f2=b17 f6+e15, 
       f3=b18 f6+e16, 
       f4=b19 f6+e17,
       f5=b20 f6+e18;

std e1-e18 = the1-the18, f6=1;
bounds the1-the18>0; 
run;

/* Estimation */
data estimate4; set estimate3; if variable notin ('Std Err' 't Value');
	keep variable coefficient1 parameter1 coefficient2 parameter2 coefficient3 parameter3 parameter4;
run;

data estimate4; set estimate4;
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

data estimate4; set estimate4;
	format factor1 $8.; 
	format factor2 $8.;
	format factor3 $8.;
    factor1=parameter1;
	factor2=parameter2;
	factor3=parameter3;
run;

data estimate4; set estimate4 (rename=(coefficient1=varload1 
                                       coefficient2=varload2 
                                       coefficient3=varload3
                                       parameter4=par2));
	oo=_n_;
run;


proc sort data=estimate4; by par2;run;

data ee4;set ee3;if Parameter ^="";
	rename Variable=par2;
	drop SpecType StdErr tValue;
run;

proc sort data=ee4;by par2;run;

data estimate4;merge estimate4(in=a) ee4(in=b); by par2;if a;run; 

data estimate4; set estimate4; drop par2; 
proc sort; by oo;run;

data estimate4;set estimate4;
	rename estimate = e;
run;

data factorest2; set estimate4;keep variable varload1 varload2 varload3 ;
data factorest2; set factorest2 (rename=(varload1=cc1 varload2=cc2 varload3=cc3));run;

proc sort data=factorest2; by variable;run;

data estimate5; set estimate4; Asterisk1='*';Asterisk2='*';Asterisk3='*'; plus1='+';plus2='+';plus3='+';equal='=';run;
proc sort data=estimate5; by  variable;run;

data est3; merge estimate5(in=a) factorest2; by variable; if a;run;
/* manual adjustment-for the 2nd order factor f6 on f4*/
data est3; set est3;
	if variable = 'arss' then cc2 = .9833;
	if variable = 'arss' then cc3 = .8148;
run;

data est3; set est3;
 if variable = 'arss' then gloading=round(varload1*cc1+varload2*cc2+varload3*cc3,.0001);
 else gloading=round(varload1*cc1,.0001);run;

data wisc5_g; retain variable varload1 Asterisk1 factor1 plus1 varload2 Asterisk2 factor2 plus2 varload3 Asterisk3 factor3 plus3 e gloading equal; set est3;
	keep variable varload1 Asterisk1 factor1 variable varload2 Asterisk2 factor2 varload3 Asterisk3 factor3 e gloading plus1 plus2 plus3 equal;
data wisc5_g; set wisc5_g; format model $20.; model="wisc5_g";run;


/** Output **/
ods rtf file="C:\Users\UZHANOU\Documents\WISC5\stdz\results\factor_analysis\g_compare\WISC4_Stdz_CFA_same_model.rtf";
proc report data=wisc4_g nowindows headline split='/' spacing=1 ;
title2 "G-Loading on wisc4 data_same model";
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


ods rtf file="C:\Users\UZHANOU\Documents\WISC5\stdz\results\factor_analysis\g_compare\WISC5_Stdz_CFA_same_model.rtf";
proc report data=wisc5_g nowindows headline split='/' spacing=1 ;
title2 "G-Loading on wisc5 data_same_model";
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

********* end of the same model, same construct, different data comparison *************;

/***************  EFA comparison ***********************/
%let pp=K:\clinical psychometrics\ongoing projects;
%let pp=K:\ongoing projects;

%include "&pp\WISC5\stdz\sascode\macros\Macro_FactorAnalysis_v1.sas";
%varimaxfact(wisc4, 1:11, 1, 1,siss vcss coss inss bdss pcnss mrss dsss lnsss arss cdss ssss cass,1-factor 13 subtests);
%varimaxfact(wisc5, 1:11, 2, 1,siss vcss coss inss bdss pcnss mrss dsss lnsss arss cdss ssss cass,1-factor 13 subtests);

%varimaxfact(wisc4, 1:11, 3, 1,siss vcss coss inss wcss bdss mrss pcnss pcmss dsss lnsss arss cdss ssss cass,1-factor 15 subtests);
%varimaxfact(wisc5, 1:11, 4, 1,inss siss vcss coss bdss vpss mrss fwss pcnss dsss psss lnsss arss ssss cdss cass,1-factor 16 subtests);


********************************************************************************************;

/**************** Using validity data to make comparison  ******************/
data wisc5_wisc4;set k.wisc5_wisc4;run;

proc sort data=wisc5_wisc4;by agegrp;run;
proc freq data=wisc5_wisc4;
table agegrp;run;

/****  Using surveyselect to increase sample size-Monte Carlo ****/
/* I have tried this method, the data is not very good to use, so this approach is dumped. 
    Ou Zhang -10-08-2014    
    do not run this section */
%macro samp(datin,datout,set);
proc surveyselect data = &datin out = &datout&set method = srs sampsize = 10 seed = 9876;
strata agegrp;
run;
%mend;

/* Monte Carlo 20 times */
%samp(wisc5_wisc4,ws,1);
%samp(wisc5_wisc4,ws,2);
%samp(wisc5_wisc4,ws,3);
%samp(wisc5_wisc4,ws,4);
%samp(wisc5_wisc4,ws,5);
%samp(wisc5_wisc4,ws,6);
%samp(wisc5_wisc4,ws,7);
%samp(wisc5_wisc4,ws,8);
%samp(wisc5_wisc4,ws,9);
%samp(wisc5_wisc4,ws,10);
%samp(wisc5_wisc4,ws,11);
%samp(wisc5_wisc4,ws,12);
%samp(wisc5_wisc4,ws,13);
%samp(wisc5_wisc4,ws,14);
%samp(wisc5_wisc4,ws,15);
%samp(wisc5_wisc4,ws,16);
%samp(wisc5_wisc4,ws,17);
%samp(wisc5_wisc4,ws,18);
%samp(wisc5_wisc4,ws,19);
%samp(wisc5_wisc4,ws,20);

data new_w4_w5;set ws1 ws2 ws3 ws4 ws5 ws6 ws7 ws8 ws9 ws10 ws11 ws12 ws13 ws14 ws15 ws16 ws17 ws18 ws19 ws20;run;

/* wisc4 descriptive statistics */
proc means mean stddev skew kurt data=wisc5_wisc4;
var wisc4_si_ss wisc4_vc_ss wisc4_co_ss  wisc4_in_ss wisc4_wr_ss 
    wisc4_bd_ss wisc4_mr_ss wisc4_pcn_ss wisc4_ds_ss wisc4_pcm_ss
    wisc4_ln_ss wisc4_ar_ss wisc4_cd_ss  wisc4_ss_ss wisc4_ca_ss;
run;

/* new data wisc4 descriptive statistics*/
proc means mean stddev skew kurt data=new_w4_w5;
var wisc4_si_ss wisc4_vc_ss wisc4_co_ss  wisc4_in_ss wisc4_wr_ss 
    wisc4_bd_ss wisc4_mr_ss wisc4_pcn_ss wisc4_ds_ss wisc4_pcm_ss
    wisc4_ln_ss wisc4_ar_ss wisc4_cd_ss  wisc4_ss_ss wisc4_ca_ss;
run;

/* wisc5 descriptive statistics */
proc means mean stddev skew kurt data=wisc5_wisc4;
var wisc5_in_ss wisc5_si_ss wisc5_vc_ss  
     	    wisc5_co_ss 
            wisc5_bd_ss /*use bonus plan B 12/30/2013 */
      	    wisc5_vp_ss wisc5_mr_ss  wisc5_fw_ss wisc5_pc_ss 
      	    wisc5_ds_ss 
            wisc5_ps_ss /*use ps1 12/30/2013 */
            wisc5_ln_ss wisc5_ar_ss 
       	    wisc5_ss_ss   
       	    wisc5_cd_ss 
            wisc5_ca_ss ;
run;

/* new data wisc5 descriptive statistics*/
proc means mean stddev skew kurt data=new_w4_w5;
var wisc5_in_ss wisc5_si_ss wisc5_vc_ss  
     	    wisc5_co_ss 
            wisc5_bd_ss /*use bonus plan B 12/30/2013 */
      	    wisc5_vp_ss wisc5_mr_ss  wisc5_fw_ss wisc5_pc_ss 
      	    wisc5_ds_ss 
            wisc5_ps_ss /*use ps1 12/30/2013 */
            wisc5_ln_ss wisc5_ar_ss 
       	    wisc5_ss_ss   
       	    wisc5_cd_ss 
            wisc5_ca_ss ;
run;

/**************   do not run the Monte Carlo Part   ***************/

/************** Using validition data to run CFA model **********************************/
/* wisc4 5-factor model */
ods output lineqseqstd =estimate0 fit=model0 lineqsvarexogstd=ee0;
proc calis data = wisc5_wisc4 outstat=stat0 pall;
var wisc4_si_ss wisc4_vc_ss wisc4_co_ss wisc4_in_ss wisc4_wr_ss 
    wisc4_bd_ss wisc4_pcm_ss wisc4_mr_ss wisc4_pcn_ss wisc4_ds_ss 
    wisc4_ln_ss wisc4_ar_ss wisc4_cd_ss wisc4_ss_ss wisc4_ca_ss;
lineqs wisc4_si_ss=b1 f1+e1 , 
       wisc4_vc_ss=b2 f1+e2 , 
       wisc4_in_ss=b3 f1+e3 , 
       wisc4_co_ss=b4 f1+e4 ,
       wisc4_wr_ss=b5 f1+e5 , 
	   
       wisc4_bd_ss=b6 f2+e6 , 
       wisc4_mr_ss=b7 f2+e7 ,
       wisc4_pcm_ss=b8 f2+e8 ,
       wisc4_pcn_ss=b9 f2+e9 ,
  
       wisc4_ds_ss=b10 f3+e10 , 
       wisc4_ln_ss=b11 f3+e11 ,

       wisc4_cd_ss=b12 f4+e12 ,
       wisc4_ss_ss=b13 f4+e13 , 
       wisc4_ca_ss=b14 f4+e14 , 

	   wisc4_ar_ss=b15 f5+e15,

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


/******************** wisc5 Model 5e ********************/
ods output lineqseqstd =estimate2 lineqseq =estimate22 fit=model2 lineqsvarexogstd=ee2;
proc calis data=wisc5_wisc4 outstat=stat1 pall;
var wisc5_in_ss wisc5_si_ss wisc5_vc_ss wisc5_co_ss wisc5_bd_ss wisc5_vp_ss wisc5_mr_ss wisc5_pc_ss 
    wisc5_fw_ss wisc5_ar_ss wisc5_ds_ss wisc5_ps_ss wisc5_ln_ss wisc5_ss_ss wisc5_cd_ss wisc5_ca_ss;
lineqs wisc5_si_ss=b1 f1+e1 , 
       wisc5_vc_ss=b2 f1+e2 , 
       wisc5_in_ss=b3 f1+e3 , 
       wisc5_co_ss=b4 f1+e4 , 
	   
       wisc5_bd_ss=b5 f2+e5 , 
       wisc5_vp_ss=b6 f2+e6 ,
 
       wisc5_mr_ss=b7 f3+e7 ,
       wisc5_fw_ss=b8 f3+e8 ,
       wisc5_pc_ss=b9 f3+e9 ,
       
       wisc5_ds_ss=b11 f4+e11 , 
       wisc5_ps_ss=b12 f4+e12 , 
       wisc5_ln_ss=b13 f4+e13 ,

       wisc5_ar_ss=b10 f1 + b23 f3 + b14 f4+e10,

       wisc5_cd_ss=b15 f5+e14 ,
       wisc5_ss_ss=b16 f5+e15 , 
       wisc5_ca_ss=b17 f5+e16 , 

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
 if variable = 'wisc5_ar_ss' then do;
	Coefficient2 = Coefficient2;
	Coefficient3 = Coefficient3;
	parameter2 = parameter2;
	parameter3 = parameter3;
	parameter4 = parameter4;
 end;
 if variable ^= 'wisc5_ar_ss' then do;
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
	if variable = 'wisc5_ar_ss' then cc2 = 1.00;
	if variable = 'wisc5_ar_ss' then cc3 = .5060;
run;

data est2; set est2;
 if variable = 'wisc5_ar_ss' then gloading=round(varload1*cc1+varload2*cc2+varload3*cc3,.0001);
 else gloading=round(varload1*cc1,.0001);run;

data wisc5_g; retain variable varload1 Asterisk1 factor1 plus1 varload2 Asterisk2 factor2 plus2 varload3 Asterisk3 factor3 plus3 e gloading equal; set est2;
	keep variable varload1 Asterisk1 factor1 variable varload2 Asterisk2 factor2 varload3 Asterisk3 factor3 e gloading plus1 plus2 plus3 equal;
data wisc5_g; set wisc5_g; format model $20.; model="wisc5_g";run;


ods rtf file="C:\Users\UZHANOU\Documents\WISC5\stdz\results\factor_analysis\g_compare\WISC4_valid_CFA Results.rtf";
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


ods rtf file="C:\Users\UZHANOU\Documents\WISC5\stdz\results\factor_analysis\g_compare\WISC5_valid_CFA Results.rtf";
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







