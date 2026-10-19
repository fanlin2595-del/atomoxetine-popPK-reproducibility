; Source: sim_AS.modmod. Simulation module, package v3.
; Equations, fixed parameters, seed and NSUB retained from the supplied source.
; Exploratory concentration simulation, not a clinical dosing recommendation.
$SIZES MAXIDS=250000
$PROB SIMULATION - MODEL III AS

$INPUT ID TIME AMT ADDL II DV MDV EVID WT AS DOSEKG TAD
$DATA synthetic_sim_AS.dat IGNORE=@
$SUBROUTINE ADVAN2 TRANS2

$PK
 TVCL=THETA(1)*(WT/30.7)**0.75
 TVCL=TVCL*(AS/1.25)**THETA(4)
 CL=TVCL*EXP(ETA(1))
 TVV=THETA(2)
 V=TVV*EXP(ETA(2))*(WT/30.7)
 KA=THETA(3)
 S2=V/1000

$ERROR
 IPRED = F                     ; Individual prediction without residual error
 Y = F*(1+EPS(1))               ; Simulated observation with residual error

; Fixed at the supplied, rounded Model III estimates.
$THETA
 (15.6 FIX)                     ; Typical CL/F, L/h, at reference BW and AS
 (67.8 FIX)                     ; Typical V/F, L, at reference BW
 (1 FIX)                        ; KA, h^-1
 (0.585 FIX)                    ; Exponent of AS/1.25 on CL/F
$OMEGA
 0.168 FIX                      ; Variance of ETA(CL)
 0.31 FIX                       ; Variance of ETA(V)
$SIGMA
 0.0712 FIX                     ; Proportional residual variance

; 500 simulation subproblems; no parameter estimation.
; The original input grid and downstream summary script are not yet supplied.
$SIMULATION (20260717) ONLYSIM NSUB=500
$TABLE ID TIME TAD AMT DV IPRED CL V WT AS DOSEKG MDV
NOPRINT ONEHEADER NOAPPEND FILE=sim_model_III_AS.tab
