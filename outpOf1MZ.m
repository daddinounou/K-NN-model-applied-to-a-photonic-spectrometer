%One MZI, diffrent wavelenghts, one Input port, tow Output ports
clear all

neff=2.5;
L1=1;
L2=1.2;
DeltaL=0.2;
%alpha=0.05; I supposed no losses, otherwise multiply the delay *exp(-alpha*L2)


sq=1/sqrt(2);
split=[[sq,-i*sq];[-i*sq, sq]];
com=split;

Beta = neff*2*pi;
Input=[1,0];

out1=[];
out2=[];
wave=[];
% out1 and out2 correpending to the outputs 1 and 2

for lambda= (1500:0.002:1510)*10^-7
    
    delay= exp(-i*Beta*L2/lambda)*[exp(-i*Beta*DeltaL/lambda),0;0,1];
    
    x=(abs(Input*split*delay*com)).^2;
    
    out1=[out1,x(1)];
    
    out2=[out2,x(2)];
    wave=[wave,lambda];
end  

figure(1)
plot(wave,out1)
figure(2)
plot(wave,out2)