; Updated source: Model I_BW.modmod. Package v2, stage 1.
; Equations, initial values, bounds and estimation settings retained.
; Data and output filenames adapted for a synthetic demonstration.
$PROB ATOMOXETINE MODEL I - BODY WEIGHT
$INPUT ID DAT2=DROP TIME AMT ADDL II DV MDV  WT SEX AGE rs1065852 rs_grp rs5030865 rs16947 AS CYP2D6 TAD
$DATA synthetic_demo.dat IGNORE=@
$SUBROUTINE ADVAN2 TRANS2

$PK
 TVCL=THETA(1)*(WT/30.7)**0.75
 CL=TVCL*EXP(ETA(1))
 TVV=THETA(2)
 V=TVV*EXP(ETA(2))*(WT/30.7)
 KA=THETA(3)
 S2=V/1000

$ERROR 
     IPRED=F
     IRES=DV-IPRED
     DEL=0
     IF(IPRED.EQ.0)DEL=1
     IWRES=IRES/(IPRED+DEL)
     Y=F*(1+EPS(1))

$THETA
     (0,14.1,20); CL
     (0,69.4,100); V
     (1,FIX); KA

$OMEGA 
     0.19; IIV IN CL
     0.359; IIV IN V

$SIGMA 
     0.0718; ERR

$ESTIMATION
METHOD=1 INTERACTION
MAXEVAL=9999 SIG=3 PRINT=5

$COV     PRINT=E
$TABLE ID TIME AMT DV PRED CWRES IPRED ADDL II MDV CL V ETA1 WT SEX AGE rs1065852 rs5030865 rs16947 AS CYP2D6 TAD NOPRINT ONEHEADER FILE=model_I_BW.fit
