% Sinais de entrada.

R1=1e3
R2=500
R3=100e3
R4=500
C1=220e-9
C2=220e-9

%{
wL=1/(R1*C1)
fL=wL/(2*pi)

w0=2*pi*1e3;

wH=w0*w0/wL
fH=wH/(2*pi)


%ganho vem anteriormente
AV=100
ganho=AV;
ganhodB=20*log10(ganho);
%}
%esboço funçao
f=logspace(1,8,70);

%{
for i=1:length(f)
  Auxpos(i)=20*log10(f(i)*(ganho/fL));
  Auxneg(i)=20*log10(f(i)*(-ganho/fL)+(ganho-ganho*fH/fL));
  if(f(i)<=fL)
    T(i)=Auxpos(i);
  elseif(f(i)<fH && Auxpos(i)<= ganhodB)
    T(i)=Auxpos(i);
  elseif(f(i)<fH && Auxpos(i)>ganhodB)
    T(i)=ganhodB;
  elseif(f(i)>=fH)
    T(i)=Auxneg(i);
  endif  
endfor
%}

T=(R1*C1*2*pi.*f*j)./(1+R1*C1*2*pi*f*j).*(1+R3/R4).*(1./(1+R2*C2*2*pi.*f*j));
Tdb=20*log10(abs(T));
Tp=180*angle(T)/pi;

p1= figure();
plot(log10(f), Tdb, "r")
title("Frequency response of vo(f)/vi(f) [dB]")
xlabel ("log_{10}(f[Hz])")
ylabel ("T(f) [dB]")
print (p1, "Tdb.eps", "-depsc");

p2= figure();
plot(log10(f), Tp, "r")
title("Phase response of vo(f)/vi(f) [deg]")
xlabel ("log_{10}(f[Hz])")
ylabel ("Phase(f) [deg]")
print (p2, "Tp.eps", "-depsc");


fCentral = 1000;

TCentral = (R1*C1*2*pi.*fCentral*j)./(1+R1*C1*2*pi*fCentral*j).*(1+R3/R4).*(1./(1+R2*C2*2*pi.*fCentral*j))
TCentralAbs = abs(TCentral)
TCentraldB = 20*log10(abs(TCentral))
TCentralPhase = 180*angle(TCentral)/pi

inputImp = R1 + 1/(j*2*pi*fCentral*C1)
inputImpAbs= abs(inputImp)
outputImp = 1/((1/R2) + j*2*pi*fCentral*C2)
outputImpAbs= abs(outputImp)


tabID2=fopen('oct2_tab.tex','w');
fprintf(tabID2, "$Z_{IN}$ & %f%+fj  \\\\ \\hline\n", real(inputImp), imag(inputImp))
fprintf(tabID2, "$Z_{IN}(Absolute)$ & %f  \\\\ \\hline\n", inputImpAbs)
fprintf(tabID2, "$Z_{OUT}$ & %f%+fj \\\\ \\hline\n", real(outputImp), imag(outputImp))
fprintf(tabID2, "$Z_{OUT}(Absolute)$ & %f  \\\\ \\hline\n", outputImpAbs)
fclose(tabID2);


% Análise Experimental

output1kHz = 13.5;
input1kHz = 0.111;

gain1kHz = output1kHz/input1kHz
gain1kHzdB = 20*log10(abs(gain1kHz))



output100Hz = 3.4;
input100Hz = 0.115;

gain100Hz = output100Hz/input100Hz
gain100HzdB = 20*log10(abs(gain100Hz))



output10000Hz = 0.860;
input10000Hz = 0.113;

gain10000Hz = output10000Hz/input10000Hz
gain10000HzdB = 20*log10(abs(gain10000Hz))





