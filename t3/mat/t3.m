% Parte 1: Envelope Detector

A = 230;
n = 4.416766044;
f = 50;
w = 2*pi*f;
R_env = 15000;
C = 60e-6;
V_ON = 0.7;

% Vetor de tempo
t=0:2e-5:0.2;

vIN = A*cos(w*t); % Tensão à entrada do transformador.

% Fonte dependente

vS = (A/n)*cos(w*t);
vS_Corrected = abs((A/n)*cos(w*t));   % Função vS corrigida pela ponte de díodos (ideais).



% --------------------------------------A RETIRAR---------------------------------------
%{
figure
plot(t*1000, vS)
title("Tensão na Fonte")
xlabel ("t[ms]")
ylabel ("vS[V]")

figure
plot(t*1000, vS_Corrected)
title("Tensão na Fonte (Corrigida)")
xlabel ("t[ms]")
ylabel ("vS Corrigida [V]")
%}


% Determinação da corrente no Díodo

Is = 1e-14;
eta = 1;

% Temperatura Nominal em Kelvin
Tnom = 300.15; 

% Constante de Boltzmann
k = 1.38064852e-23;

% Carga do eletrão
q = 1.60217662e-19;

% Thermal Voltage
Vt = (k*Tnom)/q;



% Modelação da vS Atenuada.

vS_Translacao = vS_Corrected - 2*V_ON;

% --------------------------------------A RETIRAR---------------------------------------
%{
figure
plot(t*1000, vS_Translacao)
title("Tensão na Fonte (Translacao)")
xlabel ("t[ms]")
ylabel ("vS Translacao [V]")
%}

for i=1:length(t)
  
  if vS_Translacao(i) < 0
    
    vS_Final(i) = 0;
    
  else                                    % É maior ou igual a zero.
    
    vS_Final(i) = vS_Translacao(i);
    
  endif
  
endfor

% Gráfico vS vs vS_Final.
p1= figure
plot(t*1000, vS_Corrected)
hold on
plot(t*1000, vS_Final)
title("Transformed Voltage (vS) and Diode Corrected Voltage (vS_Final)")
xlabel ("t[ms]")
ylabel ("vS vS_Final [V]")
legend('vS','vSFinal(t)');
print (p1, "p1.eps", "-depsc");

% --------------------------------------A RETIRAR---------------------------------------
%{
% Gráfico com as várias versões da tensão.

figure
plot(t*1000, vS)
hold on
plot(t*1000, vS_Corrected)
plot(t*1000, vS_Translacao)
plot(t*1000, vS_Final)
hold off
title("Tensão na Fonte")
legend('vS(t)','vSCorrected(t)', 'vSTranslacao(t)', 'vSFinal(t)');
xlabel ("t[ms]")
ylabel ("vS[V]")
%}

% Ciclo For para modelação da tensão no condensador, v_C.

T_vS = 1/(2*f);         % Período do módulo.
t_MAX = 0;              % Instante Inicial do Máximo (como é um cosseno, o máximo é logo em t = 0s!).
t_Aux = t_MAX;          % Deslocamento da função exponencial.

% Cálculo da Resistência Equivalente.

n_D = 17;
V_O = 12;        % Componente DC da tensão de saída.
V_D = V_O/n_D;   % Tensão em cada díodo.
R_reg = 5400;                          % Resistência do Regulador.
r_d = (eta*Vt)/(Is*exp(V_D/(eta*Vt)));  % Resistência Incremental do Díodo.
                               % Número de díodos.

R_eq = 1/((1/R_env)+(1/(R_reg + n_D*r_d)));

%v_C = zeros(length(t)); % Inicialização do vetor v_C.

for i=1:length(t)
  
  ramoExponencial = ((A/n)-2*V_ON)*exp(-(t - t_Aux)/(R_eq*C));
  
  
  if t(i) < t_MAX               % Antes de atingir o máximo ((A/n)-2*V_ON).
    v_C(i) = vS_Final(i);        % Mantém a sinusoide.
  elseif t(i) == t_MAX
    v_C(i) = vS_Final(i);   % Mantém a sinusoide.
    
  else 
    if ramoExponencial(i) > vS_Final(i)          % Descarga do Condensador.
      v_C(i) = ramoExponencial(i);
    elseif ramoExponencial(i) == vS_Final(i)   % Instante em que começa a carregar.
      v_C(i) = ramoExponencial(i);      
    else                                        % vCExponencial < vCSinusoidal. Volta a carregar.
      if t(i) <= t_MAX + T_vS                   % Ainda não atingiu o máximo seguinte (se não fizer isto, a função fica só a sinusoidal e cai para zero).
        v_C(i) = vS_Final(i);
      else
        t_MAX = t_MAX + T_vS;                     % O máximo avança um período.
        t_Aux = t_MAX;                            % t_Aux é o instante em que se verifica o máximo da exponencial => t_MAX - t_Aux == 0.
        
        ramoExponencial2 = ((A/n)-2*V_ON)*exp(-(t - t_Aux)/(R_env*C));
        
        v_C(i) = ramoExponencial2(i);
      endif      
    endif
    
  endif 
endfor

% --------------------------------------A RETIRAR---------------------------------------
%{
figure
plot(t*1000, v_C)
title("Tensão no Condensador")
xlabel ("t[ms]")
ylabel ("v_C [V]")
%}

h = figure();
plot(t*1000, v_C, "r");
title("Envelope Detector Output - Voltage in capacitor");
xlabel ("t[ms]");
ylabel ("v_C [V]");
print (h, "envelopeDetectorOctave.eps", "-depsc");







% Parte 2: Voltage Regulator



V_C = mean(v_C); % Componente DC da tensão no condensador.
v_c = v_C - V_C; % Componente incremental da tensão no condensador.



% Incremental Analysis
i_d = v_c/(R_reg + n_D*r_d);  % Componente incremental da corrente nos díodos.

v_o = n_D*r_d*i_d;            % Componente incremental da tensão na saída.

% Corrente no Voltage Regulator (para o caso de ser preciso).

I_D = (V_C - n_D*V_D)/R_reg;

I_d1 = I_D + i_d;

% --------------------------------------A RETIRAR---------------------------------------
%{
figure
plot(t*1000, I_d1)
title("Corrente no Voltage Regulator (Método 1)")
xlabel ("t[ms]")
ylabel ("I_d[A]")
%}



% Tensão na saída.
V_o = V_O + v_o;

% Cálculo do Ripple.

display("Ripple = ")
ripple = range(V_o)

% --------------------------------------A RETIRAR---------------------------------------
%{
figure
plot(t*1000, V_o)
title("Tensão de Saída")
xlabel ("t[ms]")
ylabel ("V_o [V]")
%}

% Gráficos da Rafaela

% Confirmado!
ha = figure();
plot(t*1000, V_o, "r");
title("Regulator Output Voltage");
xlabel ("t[ms]");
ylabel ("V_o [V]");
print (ha, "voltageRegulatorOctave.eps", "-depsc");

% Confirmado!
hb = figure();
plot(t*1000, V_o, "r");
title("DC Voltage Output"); %mesma coisa titulos diferentes
xlabel ("t[ms]");
ylabel ("V_o [V]");
print (hb, "DCvoltageOctave.eps", "-depsc");

% Confirmado!
hc = figure();
plot(t*1000, v_o, "r");
title("Output AC component + DC deviation = V_o - 12");
xlabel ("t[ms]");
ylabel ("v_o [V]");
print (hc, "V_oMenos12Octave.eps", "-depsc");

%plots opcionais que tbm tenho no ngspice para comparar ou que podem dar jeito idk

%{
% Confirmado!
hd = figure();
plot(t*1000, V_o, "r");
hold on
plot(t*1000, vS, "b");
title("Transformed Voltage (vS) vs Output Voltage (V_o)");
xlabel ("t[ms]");
ylabel ("V_o, vS [V]");
legend('V_o','vS');
print (hd, "V_ovSOctave.eps", "-depsc");

% Confirmado! Alterar o título!
he = figure();
plot(t*1000, V_o, "r");
hold on
plot(t*1000, vS, "b");
hold on
plot(t*1000, vIN, "y");
title("Voltage Comparison (vIN, vS and V_o)");
xlabel ("t[ms]");
ylabel ("V_o, vS, vIN [V]");
legend({'V_o', 'vS', 'vIN'});
print (he, "V_ovSvINOctave.eps", "-depsc");
%}
