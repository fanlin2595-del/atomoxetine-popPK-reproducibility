; Source: sim_RS1.modmod. Simulation module, package v3.
; Equations, fixed parameters, seed and NSUB retained from the supplied source.
; Exploratory concentration simulation, not a clinical dosing recommendation.
$SIZES MAXIDS=250000
$PROB SIMULATION - MODEL II RS1065852

$INPUT ID TIME AMT ADDL II DV MDV EVID WT rs_grp DOSEKG TAD
$DATA synthetic_sim_RS.dat IGNORE=@
$SUBROUTINE ADVAN2 TRANS2

$PK
 TVCL = THETA(1)*(WT/30.7)**0.75
 TVV = THETA(2)*(WT/30.7)
 IF (rs_grp.EQ.1) THEN          ; Indicator for the TT group in this encoding
   TVCL = TVCL*EXP(THETA(4))
 END IF
 CL = TVCL*EXP(ETA(1))
 V = TVV*EXP(ETA(2))
 KA = THETA(3)
 S2 = V/1000

$ERROR
 IPRED = F                     ; Individual prediction without residual error
 Y = F*(1+EPS(1))               ; Simulated observation with residual error

; Fixed at the supplied, rounded Model II estimates.
$THETA
 (17.1 FIX)                     ; Typical CL/F, L/h, at reference BW
 (67.6 FIX)                     ; Typical V/F, L, at reference BW
 (1 FIX)                        ; KA, h^-1
 (-0.645 FIX)                   ; Log multiplier for CL/F when rs_grp=1
$OMEGA
 0.156 FIX                      ; Variance of ETA(CL)
 0.329 FIX                      ; Variance of ETA(V)
$SIGMA
 0.0712 FIX                     ; Proportional residual variance

; 500 simulation subproblems; no parameter estimation.
; The original input grid and downstream summary script are not yet supplied.
$SIMULATION (20260717) ONLYSIM NSUB=500
$TABLE ID TIME TAD AMT DV IPRED CL V WT rs_grp DOSEKG MDV
NOPRINT ONEHEADER NOAPPEND FILE=sim_model_II_rs1065852.tab
