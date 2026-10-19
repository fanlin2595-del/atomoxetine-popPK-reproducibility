; Updated source: 15.modmod. Package v2, stage 1.
; Equations, initial values, bounds and estimation settings retained.
; Data and output filenames adapted for a synthetic demonstration.
$PROB ATOMOXETINE MODEL II - RS1065852
$INPUT ID DAT2=DROP TIME AMT ADDL II DV MDV  WT SEX AGE rs1065852 rs_grp rs5030865 rs16947 AS CYP2D6 TAD

$DATA synthetic_demo.dat IGNORE=@
$SUBROUTINE ADVAN2 TRANS2

$PK
 TVCL=THETA(1)*(WT/30.7)**0.75
 TVV=THETA(2)*(WT/30.7)
IF (rs_grp.EQ.1) THEN
		TVCL=TVCL*EXP(THETA(4))
END IF
CL=TVCL*EXP(ETA(1))
 V=TVV*EXP(ETA(2))
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
     (0,17.1,30); CL
     (0,67.6,100); V
     (1,FIX); KA
     (-3,-0.645,3)

$OMEGA 
     0.156 ; IIV IN CL
     0.329; IIV IN V

$SIGMA 
     0.0712; ERR

$ESTIMATION
METHOD=1 INTERACTION
MAXEVAL=9999 SIG=3 PRINT=5

$COV     PRINT=E
$TABLE ID TIME TAD AMT DV PRED CWRES IPRED ADDL II MDV CL V ETA1 WT SEX AGE rs1065852 rs_grp rs5030865 rs16947 AS CYP2D6 NOPRINT ONEHEADER FILE=model_II_rs1065852.fit
