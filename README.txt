#  (FFANC) A robust feedforward active noise control system withsimultaneous online secondary- and feedback-path modeling
This repository contains code relative to simulations from paper "A robust feedforward active noise control system with simultaneous online secondary- and feedback-path modeling". The simulations code (Case 1) is provided. The detail introduction can be found in the paper.
# Getting started
# Code structure
Case_1_Sys_A.m, Case_1_Sys_B.m, Case_1_Sys_C.m, Case_1_Sys_D.m.
# Notice
In this case, you should run the programs in order of "Case_1_Sys_A.m", "Case_1_Sys_B.m",  "Case_1_Sys_C.m"and "Case_1_Sys_D.m", to produce synthetical signal and paths, and save the results for finnal plot. Note that 
one hundrand independent runs were conducted in each program. Please follow the steps below.
Step 1: set FLG_v_p=0 and FLG_v_s=0 in "Case_1_Sys_A.m", and then run it.
Step 2: set FLG_awgn=0 in "Case_1_Sys_B.m", and then run it.
Step 3: run the program "Case_1_Sys_C.m".
Step 4: run the program "Case_1_Sys_D.m".

#Please note that the code is not generally perfected for performance, but is rather meant to illustrate certain results from the paper. The code is provided as-is without guarantees.