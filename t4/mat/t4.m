%gain stage

% Sinal de entrada.
VS = 0;       % Componente DC (que é nula, mas enfim).
vs = 10e-3;
f = 1000;
w = 2*pi*f

t = 0:1e-6:1e-2;    % Vetor de tempo.   
N=70;
f_vector = logspace(1, 8, N);

vS = VS + vs*sin(w*t);


RL=8;
VT=25e-3
BFN=178.7
VAFN=69.7
RE1=200
RC1=700
RB1=80000
RB2=20000
VBEON=0.7     %Vamos assumir isto; é o que o professor faz!
VCC=12
RS=100


% Operating Point (Gain Stage)

% Condensador = Circuito Aberto (frequências baixas).
% VS, RS e CI vão à vida!

% Equivalente de Thévenin:

% Resistência Equivalente de Thévenin (Vcc desligado => RB1 || RB2):
RB=1/(1/RB1+1/RB2)

% Tensão Equivalente de Thévenin (tensão no nó da base quando este é um circuito aberto):
% Divisor de Tensão
VEQ=RB2/(RB1+RB2)*VCC


% KVL na malha do Bias Circuit:
IB1=(VEQ-VBEON)/(RB+(1+BFN)*RE1)


% Corrente no Coletor (Definição):
IC1=BFN*IB1


% Corrente no Emissor (KCL + Definição):
IE1=(1+BFN)*IB1


% Tensão em RE1:
VE1=RE1*IE1


% Static Output Voltage (Gain Stage) = Tensão no Coletor!
% KVL na malha da direita.
VO1=VCC-RC1*IC1


% Tensão entre o Coletor e o Emissor, para confirmar que o transístor opera na F.A.R (VCE > VBEON):
VCE=VO1-VE1


% Análise OP na Gain Stage completa!!!

% Incremental Circuit: Frequências Médias => Condensadores = Curto-Circuitos!

% Definições
gm1=IC1/VT
rpi1=BFN/gm1
ro1=VAFN/IC1


% Paralelo de RS com RB (resistência equivalente do Bias), de acordo com o slide 7 da aula 17.
RSB=RB*RS/(RB+RS)

% Expressão Completa do Ganho na Gain Stage:
AV1 = RSB/RS * RC1*(RE1-gm1*rpi1*ro1)/((ro1+RC1+RE1)*(RSB+rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)

% Ganho Completo em dB:
AVI_DB = 20*log10(abs(AV1))


% Expressão do Ganho Simplificada (ro -> +inf; RB = 0). Acho que já está em módulo!
AV1simple = RB/(RB+RS) * gm1*RC1/(1+gm1*RE1)


% Ganho Simplificado em dB:
AVIsimple_DB = 20*log10(abs(AV1simple))


%Output Voltage (Gain Stage):
v_oGS = AV1 * vS;


% Implementação do Bypass Capacitor (a corrente vai toda pelo condensador em curto-circuito <=> RE1 = 0 Ohm)
% Este é que será o ganho a sério (a partir da Lower Cut-Off Frequency)!
% A expressão anterior é só para compararmos com o ganho sem o Bypass Capacitor (podemos pôr as duas no mesmo gráfico)!
RE1BC=0
AV1BC = RSB/RS * RC1*(RE1BC-gm1*rpi1*ro1)/((ro1+RC1+RE1BC)*(RSB+rpi1+RE1BC)+gm1*RE1BC*ro1*rpi1 - RE1BC^2)
AVIBC_DB = 20*log10(abs(AV1BC))
AV1BCsimple =  - RSB/RS * gm1*RC1/(1+gm1*RE1BC)
AVIBCsimple_DB = 20*log10(abs(AV1simple))


%Output Voltage (Gain Stage) com Bypass Capacitor:
v_oBypass = AV1BC * vS;

f_pequeno = logspace(1, 4, N);

ha = figure();
plot(log10(f_vector), AVI_DB, "r")
hold on
plot(log10(f_vector), AVIBC_DB, "b")
title("Voltage Gains obtained with and without the Bypass Capacitor");
xlabel ("log10(f[Hz])");
xlim ([1 8]);
ylabel ("AV1, AV1_{Bypass} [dB]");
ylim ([0 50]);
legend('AV1','AV1_{Bypass}');
print (ha, "gainGStageOctave.eps", "-depsc");


hb = figure();
plot(t*1000, v_oGS, "r")
hold on
plot(t*1000, v_oBypass, "b")
title("Output Voltages obtained with and without the Bypass Capacitor");
xlabel ("t[ms]");
ylabel ("v_o [V]");
legend('v_o (No Bypass Capacitor)','v_o (with Bypass Capacitor)');
print (hb, "voGStageOctave.eps", "-depsc");

% Determinação da Input Impedance e da Output Impedance (Gain Stage):


RE1=100;

% Input Impedance (sem o Bypass Capacitor). Mais uma vez, deve ser só para comparar os valores.
ZI1 = 1/(1/RB+1/(((ro1+RC1+RE1)*(rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)/(ro1+RC1+RE1)))
ZIII=1/(1/RB+1/RB2+1/rpi1)
% Output Impedance (sem o Bypass Capacitor).
ZX = ro1*((RSB+rpi1)*RE1/(RSB+rpi1+RE1))/(1/(1/ro1+1/(rpi1+RSB)+1/RE1+gm1*rpi1/(rpi1+RSB)))

% Esta é a expressão que está na aula 16.
ZX = ro1*(   1/RE1+1/(rpi1+RSB)+1/ro1+gm1*rpi1/(rpi1+RSB)  )/(   1/RE1+1/(rpi1+RSB) )

% Paralelo de ZX com RC1.
ZO1 = 1/(1/ZX+1/RC1)


% Output Impedance (com o Bypass Capacitor <=> RE1 = 0)
RE1=0;

ZI1 = 1/(1/RB+1/(((ro1+RC1+RE1)*(rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)/(ro1+RC1+RE1)))

% % Paralelo de ZX com RC1.Neste caso, ZX = ro1.
ZO1 = 1/(1/ro1+1/RC1)



% Análise da Gain Stage completa!!!

% Output Stage

% Transístor PNP
BFP = 227.3;
VAFP = 37.2;
RE2 = 60;
VEBON = 0.7;      %Vamos assumir isto; é o que o professor faz!


% Operating Point (Output Stage):
% O input da Output Stage é o output da Gain Stage! 
VI2 = VO1;


% Corrente no Emissor
% KVL na malha grande:
IE2 = (VCC-VEBON-VI2)/RE2


% Corrente no Coletor
% Definição (isto é preciso para calcular gm2, rpi2 e ro2):
IC2 = BFP/(BFP+1)*IE2


% Output Voltage (KVL na malha da direita):
VO2 = VCC - RE2*IE2


% Na verdade, esta fórmula é mais conveniente.
% Vou pôr aqui para ver se dá o mesmo:
VO2_Alt = VI2 + VEBON


% Análise OP (Output Stage) completa!!!


% Determinação das admitâncias/condutâncias incrementais (definições).
% As fórmulas são o inverso das das resistências respetivas!
gm2 = IC2/VT
go2 = IC2/VAFP
gpi2 = gm2/BFP
ge2 = 1/RE2


% Gain (Output Stage):
AV2 = gm2/(gm2+gpi2+go2+ge2)


% Input Impedance (Output Stage):
ZI2 = (gm2+gpi2+go2+ge2)/gpi2/(gpi2+go2+ge2)


% Output Impedance (Output Stage):
ZO2 = 1/(gm2+gpi2+go2+ge2)


% Input Voltage (Output Stage):
vi_Output = v_oBypass;


% Output Voltage (Output Stage):
vo_Output = AV2*v_oBypass;

hc = figure();
plot(t*1000, vi_Output, "r")
hold on
plot(t*1000, vo_Output, "b")
title("Input Voltage vs Output Voltage (Output Stage)");
xlabel ("t[ms]");
ylabel ("v_i v_o [V]");
legend('v_i','v_o');
print (hc, "vivoOStageOctave.eps", "-depsc");



% Circuito Total/Completo


% Condutância equivalente à entrada da Output Stage. Torna a expressão do Ganho Total mais fácil de escrever (isto é a primeira parcela do numerador).
gB = 1/(1/gpi2+ZO1)

% Ganho Total (está no PDF dos cálculos à mão).
AV = (gB+gm2/gpi2*gB)/(gB+ge2+go2+gm2/gpi2*gB)*AV1BC
AV_DB = 20*log10(abs(AV))


% Input Impedance do Circuito Total = Input Impedance da Gain Stage.
ZI=ZI1


% Output Impedance do Circuito Total (está no PDF dos cálculos à mão).
ZO=1/(go2+gm2/gpi2*gB+ge2+gB)

tabID2=fopen('oct2_tab.tex','w');
fprintf(tabID2, "$Z_{I}$ & %.10f \\\\ \\hline\n", ZI) 
fprintf(tabID2, "$Z_{O}$ & %.10f \\\\ \\hline\n", ZO) 
fclose(tabID2);


% Output Voltage (Circuito Completo):
vo_Total = AV*vS;


hd = figure();
plot(t*1000, vS, "r")
hold on
plot(t*1000, vo_Total, "b")
title("Input Signal Voltage vs Output Voltage (Output Stage)");
xlabel ("t[ms]");
ylabel ("v_s v_o [V]");
legend('v_s','v_o');
print (hd, "vsvoTotal.eps", "-depsc");

tabID=fopen('oct1_tab.tex','w');
fprintf(tabID, "$I_{B1}$ & %.10f \\\\ \\hline\n", IB1) 
fprintf(tabID, "$I_{C1}$ & %.10f \\\\ \\hline\n", IC1) 
fprintf(tabID, "$I_{E1}$ & %.10f \\\\ \\hline\n", IE1) 
fprintf(tabID, "$V_{E1}$ & %.10f \\\\ \\hline\n", VE1) 
fprintf(tabID, "$V_{O1}$ & %.10f \\\\ \\hline\n", VO1)  
fprintf(tabID, "$V_{CE}$ & %.10f \\\\ \\hline\n", VCE) 
fprintf(tabID, "$I_{E2}$ & %.10f \\\\ \\hline\n", IE2) 
fprintf(tabID, "$I_{C2}$ & %.10f \\\\ \\hline\n", IC2) 
fprintf(tabID, "$V_{O2}$ & %.10f \\\\ \\hline\n", VO2)  
fclose(tabID);

%R quivalente cb
RE1=200
VFAKE=1
Ielol=VFAKE/RE1;
I1lol=RE1*Ielol/(rpi1+1/(1/RB1+1/RS));
RBS=1/(1/RB1+1/RS);

A=[RC1, ro1; 1,1;]
B=[-VFAKE+(gm1*rpi1*I1lol-I1lol-Ielol)*ro1; I1lol+Ielol]
C=A\B;

Is=C(2);
Zeqcb=VFAKE/Is

%req co
Zeqco=ZO+RL

%req ci
Zeqci=RS+ZI

Ci=100e-6;
Cb=2500e-6;
Co=800e-6;
w_Lower=1/(Ci*Zeqci) + 1/(Cb*Zeqcb) + 1/(Co*Zeqco)
flower=w_Lower/(2*pi)

ganho=AV;
ganhodB=20*log10(ganho);
%esboço funçao
f=logspace(1,8,70);

for i=1:length(f)
  Aux(i)=20*log10(f(i)*(10^((ganhodB-3)/20))/flower);
  if(f(i)<=flower)
    T(i)=Aux(i);
  elseif(Aux(i)< ganhodB)
    T(i)=Aux(i);
  else
    T(i)=ganhodB;
  endif  
endfor

p1= figure();
plot(log10(f), T, "r")
title("Frequency response of vo(f)/vi(f) (dB)")
xlabel ("log_{10}(f[Hz])")
ylabel ("vo(f)/vi(f)")
print (p1, "octave.eps", "-depsc");

