close all
clear all
clc

pkg load symbolic
dataID= fopen('../dataGen/data.txt','r');
formatSpec='%*s %*c %f';
for K = 1 : 9;
   fgetl(dataID); 
end;
D=fscanf(dataID, formatSpec);
fclose(dataID);

R1 = 1e3*D(1); 
R2 = 1e3*D(2);
R3 = 1e3*D(3);
R4 = 1e3*D(4);
R5 = 1e3*D(5); 
R6 = 1e3*D(6);
R7 = 1e3*D(7); 
Vs = D(8);
C1 = 1e-6*D(9);
Kb = 1e-3*D(10);
Kd = 1e3*D(11);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PONTO 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ficheiro ngspice1ID
ngspice1ID=fopen('../sim/netlist1.mod','w');
fprintf(ngspice1ID, "R1 1 2 %.8f \n", R1) 
fprintf(ngspice1ID, "R2 3 2 %.8f \n", R2) 
fprintf(ngspice1ID, "R3 2 5 %.8f \n", R3) 
fprintf(ngspice1ID, "R4 5 0 %.8f \n", R4) 
fprintf(ngspice1ID, "R5 5 6 %.8f \n", R5) 
fprintf(ngspice1ID, "R6 fic 7 %.8f \n", R6) 
fprintf(ngspice1ID, "R7 7 8 %.8f \n", R7) 
fprintf(ngspice1ID, "Vs 1 0 %.8f \n", Vs)
fprintf(ngspice1ID, "C1 6 8 %.8f \n", C1)
fprintf(ngspice1ID, "Vid 0 fic 0V \n") 
fprintf(ngspice1ID, "Gib 6 3 (2,5) %.12f \n", Kb) 
fprintf(ngspice1ID, "Hvd 5 8 Vid %.18f \n", Kd)   
fclose(ngspice1ID);

%contas teóricas 1
G1=1/R1;
G2=1/R2;
G3=1/R3;
G4=1/R4;
G5=1/R5;
G6=1/R6;
G7=1/R7;

A=[-G1, G1+G2+G3, -G2, -G3, 0, 0, 0;    %nó 2
    0, -Kb-G2, G2, Kb, 0, 0, 0;         %nó 3
    0, Kb, 0, -G5-Kb, G5, 0, 0;         %nó 6
    0, 0, 0, 0, 0, G6+G7, -G7;          %nó 7
    1, 0, 0, 0, 0, 0, 0;                %Vs
    0, 0, 0,1, 0, Kd*G6, -1;            %Vd
    G1, -G1, 0, -G4, 0, G7, -G7]     %supernó 1-0-7
B=[0; 0; 0; 0; Vs; 0; 0]
C = A\B

V1=C(1)
V2=C(2)
V3=C(3)
V5=C(4)
V6=C(5)
V7=C(6)
V8=C(7)
V6_1=V6
Ia=(V1-V2)*G1
Ib=(V3-V2)*G2
Ib=Kb*(V2-V5)
Ic=0
Id=-V7*G6
Id=(V7-V8)*G7
V4=0;
Vx=V6-V8

%imprime contas 1 em latex
tabID=fopen('oct1_tab.tex','w');
fprintf(tabID, "$V_1$ & %.10f \\\\ \\hline\n", V1) 
fprintf(tabID, "$V_2$ & %.10f \\\\ \\hline\n", V2) 
fprintf(tabID, "$V_3$ & %.10f \\\\ \\hline\n", V3) 
fprintf(tabID, "$V_4$ & %.10f \\\\ \\hline\n", V4) 
fprintf(tabID, "$V_5$ & %.10f \\\\ \\hline\n", V5) 
fprintf(tabID, "$V_5$ & %.10f \\\\ \\hline\n", V5) 
fprintf(tabID, "$V_6$ & %.10f \\\\ \\hline\n", V6) 
fprintf(tabID, "$V_7$ & %.10f \\\\ \\hline\n", V7)
fprintf(tabID, "$V_8$ & %.10f \\\\ \\hline\n", V8) 
fprintf(tabID, "$I_a$ & %.10f \\\\ \\hline\n", Ia) 
fprintf(tabID, "$I_b$ & %.10f \\\\ \\hline\n", Ib) 
fprintf(tabID, "$I_c$ & %.10f \\\\ \\hline\n", Ic) 
fprintf(tabID, "$I_d$ & %.10f \\\\ \\hline\n", Id)   
fclose(tabID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PONTO 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ficheiro ngspice2ID
Vs=0
ngspice2ID=fopen('../sim/netlist2.mod','w');
fprintf(ngspice2ID, "R1 1 2 %.8f \n", R1) 
fprintf(ngspice2ID, "R2 3 2 %.8f \n", R2) 
fprintf(ngspice2ID, "R3 2 5 %.8f \n", R3) 
fprintf(ngspice2ID, "R4 5 0 %.8f \n", R4) 
fprintf(ngspice2ID, "R5 5 6 %.8f \n", R5) 
fprintf(ngspice2ID, "R6 fic 7 %.8f \n", R6) 
fprintf(ngspice2ID, "R7 7 8 %.8f \n", R7) 
fprintf(ngspice2ID, "Vs 1 0 0V \n")
fprintf(ngspice2ID, "Vx 6 8 %.8f \n", Vx)
fprintf(ngspice2ID, "Vid 0 fic 0V \n") 
fprintf(ngspice2ID, "Gib 6 3 (2,5) %.8f \n", Kb) 
fprintf(ngspice2ID, "Hvd 5 8 Vid %.18f \n", Kd)   
fclose(ngspice2ID);

%contas 2
A(3, :)=[0,0,0,0,1,0,-1]
B(3)=Vx;
B(5)=0;

C=A\B

V1=C(1)
V2=C(2)
V3=C(3)
V5=C(4)
V6=C(5)
V7=C(6)
V8=C(7)
Vx=V6-V8
Ib=(V3-V2)*G2
Ia=(V1-V2)*G1
Id=(V7-V8)*G7
Ix=(V5-V6)*G5-Ib
Ic=Ix
Req=abs(Vx/Ix)
tau=Req*C1

%imprime contas 2 em latex
tab2ID=fopen('oct2_tab.tex','w');
fprintf(tab2ID, "$V_1$ & %.10f \\\\ \\hline\n", V1) 
fprintf(tab2ID, "$V_2$ & %.10f \\\\ \\hline\n", V2) 
fprintf(tab2ID, "$V_3$ & %.10f \\\\ \\hline\n", V3) 
fprintf(tab2ID, "$V_4$ & %.10f \\\\ \\hline\n", V4) 
fprintf(tab2ID, "$V_5$ & %.10f \\\\ \\hline\n", V5) 
fprintf(tab2ID, "$V_6$ & %.10f \\\\ \\hline\n", V6) 
fprintf(tab2ID, "$V_7$ & %.10f \\\\ \\hline\n", V7)
fprintf(tab2ID, "$V_8$ & %.10f \\\\ \\hline\n", V8) 
fprintf(tab2ID, "$I_a$ & %.10f \\\\ \\hline\n", Ia) 
fprintf(tab2ID, "$I_b$ & %.10f \\\\ \\hline\n", Ib) 
fprintf(tab2ID, "$I_c$ & %.10f \\\\ \\hline\n", Ic) 
fprintf(tab2ID, "$I_x$ & %.10f \\\\ \\hline\n", Ix)
fprintf(tab2ID, "$I_d$ & %.10f \\\\ \\hline\n", Id) 
fprintf(tab2ID, "$R_{eq}$ & %.10f \\\\ \\hline\n", Req) 
fprintf(tab2ID, "$\\tau$ & %.10f \\\\ \\hline\n", tau) 
fclose(tab2ID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PONTO 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ficheiro ngspice3ID
ngspice3ID=fopen('../sim/netlist3.mod','w');
fprintf(ngspice3ID, "R1 1 2 %.8f \n", R1) 
fprintf(ngspice3ID, "R2 3 2 %.8f \n", R2) 
fprintf(ngspice3ID, "R3 2 5 %.8f \n", R3) 
fprintf(ngspice3ID, "R4 5 0 %.8f \n", R4) 
fprintf(ngspice3ID, "R5 5 6 %.8f \n", R5) 
fprintf(ngspice3ID, "R6 fic 7 %.8f \n", R6) 
fprintf(ngspice3ID, "R7 7 8 %.8f \n", R7) 
fprintf(ngspice3ID, "Vs 1 0 0V \n")
fprintf(ngspice3ID, "C1 6 8 %.12f \n", C1)
fprintf(ngspice3ID, "Vid 0 fic 0V \n") 
fprintf(ngspice3ID, "Gib 6 3 (2,5) %.8f \n", Kb) 
fprintf(ngspice3ID, "Hvd 5 8 Vid %.18f \n \n", Kd)  
fprintf(ngspice3ID, ".op \n.ic v(6)=%.8f v(8)=%.8f v(fic)=0 \n.end", V6, V8)   
fclose(ngspice3ID);

%contas 3
t=0:2e-5:20e-3; %s

v_6n= V6*exp(-t/tau);
v_8n= V8*exp(-t/tau);

%imprime figura v(6) natural
hn = figure ();
plot (t*1000, v_6n, "r");
title("Natural solution of v6 in function of time");
xlabel ("t[ms]");
ylabel ("V6_{natural}(t) [V]");
legend('Natural solution of V6');
print (hn, "ponto3teorico.eps", "-depsc");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PONTO 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ficheiro ngspice4ID
ngspice4ID=fopen('../sim/netlist4.mod','w');
fprintf(ngspice4ID, "R1 1 2 %.8f \n", R1) 
fprintf(ngspice4ID, "R2 3 2 %.8f \n", R2) 
fprintf(ngspice4ID, "R3 2 5 %.8f \n", R3) 
fprintf(ngspice4ID, "R4 5 0 %.8f \n", R4) 
fprintf(ngspice4ID, "R5 5 6 %.8f \n", R5) 
fprintf(ngspice4ID, "R6 fic 7 %.8f \n", R6) 
fprintf(ngspice4ID, "R7 7 8 %.8f \n", R7) 
fprintf(ngspice4ID, "Vs 1 0 sin(0 1 1k) ac 1.0 sin(0 1 1k) \n")
fprintf(ngspice4ID, "C1 6 8 %.12f \n", C1)
fprintf(ngspice4ID, "Vid 0 fic 0V \n") 
fprintf(ngspice4ID, "Gib 6 3 (2,5) %.8f \n", Kb) 
fprintf(ngspice4ID, "Hvd 5 8 Vid %.18f \n", Kd) 
fprintf(ngspice4ID, "\n.op \n.ic v(6)=%.8f v(8)=%.8f v(fic)=0 \n.end", V6, V8)  
fclose(ngspice4ID);

%contas 4
w=2*pi*1e3
iwc=i*w*C1
Vscomp=1*exp(-i*(pi/2))
A=[-G1, G1+G2+G3, -G2, -G3, 0, 0, 0;    %nó 2
    0, -Kb-G2, G2, Kb, 0, 0, 0;         %nó 3
    0, Kb, 0, -Kb-G5, G5+iwc, 0, -iwc;     %nó 6
    0, 0, 0, 0, 0, G6+G7, -G7;          %nó 7
    1, 0, 0, 0, 0, 0, 0;                %Vs
    0, 0, 0,1, 0, Kd*G6, -1;            %Vd
    G1, -G1, 0, -G4, 0, G7, -G7];   %supernó 1-0-7
B=[0; 0; 0; 0; Vscomp; 0; 0];
C = A\B

V1=C(1);
V1R=real(V1)
V1I=imag(V1)
v1ang=atan(V1I/V1R)
v1ang2=angle(V1)
V1abs=abs(V1)
V2=C(2);
V2R=real(V2)
V2I=imag(V2)
v2ang=atan(V2I/V2R)
V2abs=abs(V2)
V3=C(3);
V3R=real(V3)
V3I=imag(V3)
v3ang=atan(V3I/V3R)
V3abs=abs(V3)
V5=C(4);
V5R=real(V5)
V5I=imag(V5)
v5ang=atan(V5I/V5R)
V5abs=abs(V5)
V6=C(5);
V6R=real(V6)
V6I=imag(V6)
v6ang=atan(V6I/V6R)
V6R=real(V6)
V6I=imag(V6)
V6abs=abs(V6)
V7=C(6);
V7R=real(V7)
V7I=imag(V7)
v7ang=atan(V7I/V7R)
V7abs=abs(V7)
V8=C(7);
V8R=real(V8)
V8I=imag(V8)
v8ang=atan(V8I/V8R)
V8abs=abs(V8)
V4=0*exp(i*0)
%imprime valores 4 em tabela latex
tab4ID=fopen('oct4_tab.tex','w');
fprintf(tab4ID, "$\\tilde{V_1}$ & %.10f + %.10f i \\\\ \\hline\n", V1R, V1I) 
fprintf(tab4ID, "$\\tilde{V_2}$ & %.10f + %.10f i \\\\ \\hline\n", V2R, V2I) 
fprintf(tab4ID, "$\\tilde{V_3}$ & %.10f + %.10f i \\\\ \\hline\n", V3R, V3I)
fprintf(tab4ID, "$\\tilde{V_4}$ & %.10f + %.10f i \\\\ \\hline\n", V4, V4) 
fprintf(tab4ID, "$\\tilde{V_5}$ & %.10f + %.10f i \\\\ \\hline\n", V5R, V5I) 
fprintf(tab4ID, "$\\tilde{V_6}$ & %.10f + %.10f i \\\\ \\hline\n", V6R, V6I) 
fprintf(tab4ID, "$\\tilde{V_7}$ & %.10f + %.10f i \\\\ \\hline\n", V7R, V7I)
fprintf(tab4ID, "$\\tilde{V_8}$ & %.10f + %.10f i \\\\ \\hline\n", V8R, V8I)   
fclose(tab4ID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PONTO 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%contas 5
Vs=D(8)

t1=-5e-3:5e-5:0;
vs1=Vs*ones(size(t1));
v_6total1=V6_1*ones(size(t1));

t2=0:2e-5:20e-3;
v_6f=V6abs*cos(w*t2-v6ang);
v_8f=V8abs*cos(w*t2-v8ang);
v_6total2=v_6f+v_6n;
v_8total=v_8n+v_8f;
vs2=sin(w*t2);

t=[t1, t2];
vs=[vs1, vs2];
v_6total= [v_6total1, v_6total2];

ht = figure ();
plot (t*1000, vs, "b");
hold on
plot (t*1000, v_6total, "r");
title("Value of vs and v6 in function of time.");
xlabel ("t[ms]");
ylabel ("vs(t), v6(t) [V]");
legend('vs(t)','v6(t)');
print (ht, "ponto5teorico.eps", "-depsc");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PONTO 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%contas 6
Vctotal=v_6total2-v_8total;
N=7*5
f=logspace(-1, 6, N);
for n=1:N
w(n)=2*pi*f(n);
iwc(n)=i*w(n)*C1;
Vscomp=1*exp(-i*(pi/2));
A=[-G1, G1+G2+G3, -G2, -G3, 0, 0, 0;    %nó 2
    0, -Kb-G2, G2, Kb, 0, 0, 0;         %nó 3
    0, Kb, 0, -Kb-G5, G5+iwc(n), 0, -iwc(n);     %nó 6
    0, 0, 0, 0, 0, G6+G7, -G7;          %nó 7
    1, 0, 0, 0, 0, 0, 0;                %Vs
    0, 0, 0,1, 0, Kd*G6, -1;            %Vd
    G1, -G1, 0, -G4, 0, G7, -G7];   %supernó 1-0-7
B=[0; 0; 0; 0; Vscomp; 0; 0];
C = A\B;
%amplitude
vsR(n)=abs(Vscomp);
vsR_dB(n)=20*log10(vsR(n));
v6R(n)=abs(C(5));
v6R_dB(n)=20*log10(v6R(n));
vcR(n)=abs(C(5)-C(7));
vcR_dB(n)=20*log10(vcR(n));
%phase
vsReal(n)=real(Vscomp);
vsImag(n)=imag(Vscomp);
vsangle(n)=atan(vsImag(n)/vsReal(n));
vsangle_deg(n)=(180/pi)*vsangle(n);
v6Real(n)=real(C(5));
v6Imag(n)=imag(C(5));
v6angle(n)=atan(v6Imag(n)/v6Real(n));
v6angle_deg(n)=(180/pi)*v6angle(n)-180;
vcReal(n)=real(C(5)-C(7));
vcImag(n)=imag(C(5)-C(7));
vcangle(n)=atan(vcImag(n)/vcReal(n));
vcangle_deg(n)=(180/pi)*vcangle(n)-180;
endfor

hamp= figure ()
semilogx(f, v6R_dB, "r")
hold on
semilogx(f, vsR_dB, "b")
hold on
semilogx(f, vcR_dB, "y")
title("Magnitude in function of frequency");
xlabel ("f[Hz]");
ylabel ("v6_{Magnitude}(f), vs_{Magnitude}(f), vc_{Magnitude}(f) [dB]");
legend( 'v6_{Magnitude}(f)', 'vs_{Magnitude}(f)', 'vc_{Magnitude}(f)',"location", "southwest");
print (hamp, "ponto6teorico_amplitude.eps", "-depsc");
hold off

hphase= figure ()
semilogx(f, v6angle_deg, "r")
hold on
semilogx(f, vsangle_deg, "b")
hold on
semilogx(f, vcangle_deg, "y")
title("Phase in function of frequency (Angle interval = [-270, -90] degrees)");
xlabel ("f[Hz]");
ylabel ("v6_{phase}(f), vs_{phase}(f), vc_{phase}(f) [degrees]");
legend( 'v6_{phase}(f)', 'vs_{phase}(f)', 'vc_{phase}(f)',"location", "southwest");
print (hphase, "ponto6teorico_phase.eps", "-depsc");
hold off
