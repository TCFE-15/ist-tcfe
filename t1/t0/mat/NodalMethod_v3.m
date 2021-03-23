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

G1=vpa(1/1004.8257256262905)
G2=vpa(1/2032.8619064870104)
G3=vpa(1/3117.854349667399) 
G4=vpa(1/4110.168347468018)
G5=vpa(1/3136.2898714052836) 
G6=vpa(1/2097.4169217171297)
G7=vpa(1/1031.4258704679073)

Va =vpa(5.147931491661591)
Id =vpa(0.0010417824527867048)
Kb =vpa(0.007197203244397765)
Kc =vpa(8178.07176177153 )

Z=vpa(0)
o=vpa(1)

A=[G1, -G1-G2-G3, G2, Z, Z, Z, Z;
    Z, Kb+G2, -G2, Z, Z, Z, Z;
    Z, Kb, Z, Z, G5, Z, Z;
    Z, Z, Z, G6, Z, -G6-G7, G7;
    o, Z, Z, -o, Z, Z, Z;
    Z, Z, Z, Kc*G6, Z, -Kc*G6, o;
    G1, -G1-G2, G2, G4+G6, G5, -G6, Z]

B=[Z; Z; Id;Z; Va; Z; Id]

C=A\B


Vc = -C(7)

Ic=Vc/Kc
Ic=(C(4)-C(6))*G6

Vb=C(2)

Ia=(C(1)-C(2))*G1
Ib=Vb*Kb
Ib=(C(3)-C(2))*G2

fileID=fopen('Node_tab.tex','w')
fprintf(fileID, "$V_1$ & %.10f \\\\ \\hline\n", double(C(1))) 
fprintf(fileID, "$V_2$ & %.10f \\\\ \\hline\n", double(C(2))) 
fprintf(fileID, "$V_3$ & %.10f \\\\ \\hline\n", double(C(3))) 
fprintf(fileID, "$V_4$ & %.10f \\\\ \\hline\n", double(C(4)))
fprintf(fileID, "$V_5$ & %.10f \\\\ \\hline\n", double(C(5)))  
fprintf(fileID, "$V_6$ & %.10f \\\\ \\hline\n", double(C(6))) 
fprintf(fileID, "$V_7$ & %.10f \\\\ \\hline\n", double(C(7))) 
fprintf(fileID, "$V_b$ & %.10f \\\\ \\hline\n", double(Vb)) 
fprintf(fileID, "$V_c$ & %.10f \\\\ \\hline\n", double(Vc)) 
fprintf(fileID, "$I_a$ & %.10f \\\\ \\hline\n", double(Ia))
fprintf(fileID, "$I_b$ & %.10f \\\\ \\hline\n", double(Ib)) 
fprintf(fileID, "$I_c$ & %.10f \\\\ \\hline\n", double(Ic)) 
fprintf(fileID, "$I_d$ & %.10f \\\\ \\hline\n", double(Id)) 
fclose(fileID)





