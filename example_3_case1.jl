using HomotopyContinuation
using DelimitedFiles
using LinearAlgebra


@var p11 p12 p13 p21 p22 p23

P1=[p11 p12;p12 p13]
P2=[p21 p22;p22 p23]

### Example 3 parameters 
#################### Case-1 #############
A=[1 0;0 1]; 
#################### Case-2 ##############
#A=[-1 0;0 -2]; 


B1=[1;0]; B2=[0;1]
Q1=[1 0;0 1]; Q2=[1 0;0 1] 
R1=1; R2=1
V1=0.9; V2=0.9
E=[1;1]

##
S1=B1*inv(R1)*B1'; S2=B2*inv(R2)*B2'; 
M1=E*inv(V1)*E'; M2=E*inv(V2)*E';

A1=A-S2*P2; A2=A-S1*P1;
#
are1=A1'*P1+P1*A1+Q1-P1*S1*P1+P1*M1*P1
are2=A2'*P2+P2*A2+Q2-P2*S2*P2+P2*M2*P2

#
###
equations = [are1[1,1]; are1[1,2]; are1[2,2]; +are2[1,1]; are2[1,2]; are2[2,2]]

F = System(equations)


allsolutions = solve(F)



sol=real_solutions(allsolutions)



