using HomotopyContinuation
using DelimitedFiles
using LinearAlgebra


@var p11 p12 p13  p14 p15 p16 p21 p22 p23 p24 p25 p26 p31 p32 p33 p34 p35 p36 

P1=[p11 p12 p13; p12 p14 p15; p13 p15 p16]
P2=[p21 p22 p23; p22 p24 p25; p23 p25 p26]
P3=[p31 p32 p33; p32 p34 p35; p33 p35 p36]


### Example 1 parameters 
#
A=[1 0 0;0 1 0;0 0 1]; B1=[1;0;0]; B2=[0;1;0]; B3=[0;0;1]
R1=1; R2=1; R3=1
Q1=[1 -1 0;-1 1 0;0 0 0]; Q2=[1 -1 0;-1 2 -1;0 -1 1]; Q3=[0 0 0;0 1 -1;0 -1 1] 
S1=B1*R1*B1'; S2=B2*R2*B2'; S3=B3*R3*B3'
A1=A-S2*P2-S3*P3; A2=A-S1*P1-S3*P3; A3=A-S1*P1-S2*P2
stblcntindx=Int16(1)
#
are1=A1'*P1+P1*A1+Q1-P1*S1*P1
are2=A2'*P2+P2*A2+Q2-P2*S2*P2
are3=A3'*P3+P3*A3+Q3-P3*S3*P3
#

### Coupled Algebriac Riccati Equations  (18 polynomial equations in 18 variables)
###
equations = [are1[1,1]; are1[1,2]; are1[1,3]; are1[2,2]; are1[2,3]; are1[3,3];
+are2[1,1]; are2[1,2]; are2[1,3]; are2[2,2]; are2[2,3]; are2[3,3];
+are3[1,1]; are3[1,2]; are3[1,3]; are3[2,2]; are3[2,3]; are3[3,3]]

###  Solve the polynomails equaitons using Homotopy Continuation method 
###
F = System(equations)
allsolutions = solve(F)

### Real solutions of CARE
###
sol=real_solutions(allsolutions)


