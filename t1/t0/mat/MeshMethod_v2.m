close all
clear all
clc

pkg load symbolic
R1=vpa(1004.82572563)
R2=vpa(2032.86190649)
R3=vpa(3117.85434967) 
R4=vpa(4110.16834747)
R5=vpa(3136.28987141) 
R6=vpa(2097.41692172)
R7=vpa(1031.42587047)
Va =vpa(5.1479314917)
Id =vpa(0.00104178245279)
Kb =vpa(0.0071972032444)
Kc =vpa(8178.07176177)
Z=vpa(0)
o=vpa(1)

A=[R1+R3+R4, R3, R4;
  -R4, Z, Kc-R4-R6-R7;
   Kb*R3, Kb*R3-o, Z] 

B=[Va; Z; Z]

C=A\B

Ic=C(3)
Vc=Kc*C(3)

Ib=C(2)

Vb=(C(1)+Ib)*R3
Vb=Ib/Kb

V2=Vb
V1=R1*C(1)+V2
V4=V1-Va
V3=V2+Ib*R2
V6=V4-Ic*R6
V7=V6-R7*Ic
V5=(Id-Ib)*R5

fileID=fopen('Mesh_tab.tex','w')
fprintf(fileID, "$I_a$ & %.10f \\\\ \\hline\n", double(C(1))) 
fprintf(fileID, "$I_b$ & %.10f \\\\ \\hline\n", double(C(2))) 
fprintf(fileID, "$I_c$ & %.10f \\\\ \\hline\n", double(C(3)))
fprintf(fileID, "$I_d$ & %.10f \\\\ \\hline\n", double(Id))  
fprintf(fileID, "$V_b$ & %.10f \\\\ \\hline\n", double(Vb)) 
fprintf(fileID, "$V_c$ & %.10f \\\\ \\hline\n", double(Vc)) 
fprintf(fileID, "$V_1$ & %.10f \\\\ \\hline\n", double(V1)) 
fprintf(fileID, "$V_2$ & %.10f \\\\ \\hline\n", double(V2)) 
fprintf(fileID, "$V_3$ & %.10f \\\\ \\hline\n", double(V3)) 
fprintf(fileID, "$V_4$ & %.10f \\\\ \\hline\n", double(V4)) 
fprintf(fileID, "$V_5$ & %.10f \\\\ \\hline\n", double(V5)) 
fprintf(fileID, "$V_6$ & %.10f \\\\ \\hline\n", double(V6)) 
fprintf(fileID, "$V_7$ & %.10f \\\\ \\hline\n", double(V7)) 


fclose(fileID)



