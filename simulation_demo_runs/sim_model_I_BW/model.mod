; Source: sim_BW.modmod. Simulation module, package v3.
; Equations, fixed parameters, seed and NSUB retained from the supplied source.
; Exploratory concentration simulation, not a clinical dosing recommendation.
; AS is a retained input/output field but is NOT used by this BW-only model.
$SIZES MAXIDS=250000
$PROB SIMULATION - MODEL I BW

$INPUT ID TIME AMT ADDL II DV MDV EVID WT AS DOSEKG TAD
$DATA synthetic_sim_BW.dat IGNORE=@
$SUBROUTINE ADVAN2 TRANS2

$PK
 TVCL=THETA(1)*(WT/30.7)**0.75
 CL=TVCL*EXP(ETA(1))
 TVV=THETA(2)
 V=TVV*EXP(ETA(2))*(WT/30.7)
 KA=THETA(3)
 S2=V/1000

$ERROR
 IPRED = F                     ; Individual prediction without residual error
 Y = F*(1+EPS(1))               ; Simulated observation with residual error

; Fixed at the supplied, rounded Model I estimates.
$THETA
 (14.1 FIX)                     ; Typical CL/F, L/h, at reference BW
 (69.4 FIX)                     ; Typical V/F, L, at reference BW
 (1 FIX)                        ; KA, h^-1
$OMEGA
 0.19 FIX                       ; Variance of ETA(CL)
 0.359 FIX                      ; Variance of ETA(V)
$SIGMA
 0.0718 FIX                     ; Proportional residual variance

; 500 simulation subproblems; no parameter estimation.
; The original input grid and downstream summary script are not yet supplied.
$SIMULATION (20260717) ONLYSIM NSUB=500
$TABLE ID TIME TAD AMT DV IPRED CL V WT AS DOSEKG MDV
NOPRINT ONEHEADER NOAPPEND FILE=sim_model_I_BW.tab
