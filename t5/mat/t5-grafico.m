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
title("Frequency response of vo(f)/vi(f)[dB]")
xlabel ("log_{10}(f[Hz])")
ylabel ("T(f) [dB]")
print (p1, "Tdb.eps", "-depsc");

p2= figure();
plot(log10(f), Tp, "r")
title("Phase response of vo(f)/vi(f)[deg]")
xlabel ("log_{10}(f[Hz])")
ylabel ("Phase(f) [deg]")
print (p2, "Tp.eps", "-depsc");

inputImp=1+1j;
outputImp=2+1j;

tabID2=fopen('oct2_tab.tex','w');
fprintf(tabID2, "$Z_{IN}$ & %f%+fj  \\\\ \\hline\n", real(inputImp), imag(inputImp)) 
fprintf(tabID2, "$Z_{OUT}$ & %f%+fj \\\\ \\hline\n", real(outputImp), imag(outputImp)) 
fclose(tabID2);
