% here the transfere matrix is the product of the parts of MZ 
%% losses due to index fluctuation with lambda and with width 
close all
clear all
%neff=3.49769;
L1=1e7;
%DeltaL0=365*1e3;
DeltaL0=365;
split=[1,-i;-i,1]/sqrt(2);
com=split;
Input=[0,1];

%%

w1=[];
w2=[];
MZ=[1:100];
step=0.001;
lambda= [1548:step:1552];  
%neff0=ones(length(lambda));
neff0=-1.6*1e-3*lambda+4.5197;                    %n is function of " lambda " here for width 450nm (linear dipendence lumerical)
delta_neff=0.0226;                                %fluctuation of n due to fluctuation of width(+-20nm) (deltaN to be added randomly)  

alpha= ( log(10)/20 ) * 1e-7 *0;                  %2 dB/cm

amp=0;  %the amplitude of whit noise

r1=randn(length(lambda),length(MZ));
r2=randn(length(lambda),length(MZ));
noise1=amp*r1./max(r1(:));
noise2=amp*r2./max(r2(:));


for N=1:length(MZ)
    L2=L1+N*DeltaL0;    
    DeltaL=L2-L1;
    out1=[];
    out2=[];    
    %neff=neff0-delta_neff+2*delta_neff*rand;      %fluctuation for n is added for each MZ " due to width fluctuation"

    
    for i=1:length(lambda)

        Beta=2*pi*neff0(i)/lambda(i);

        delay= exp(-alpha*L2) * exp(-1i*Beta*L2)  *[ exp(alpha*DeltaL)*exp(1i*Beta*DeltaL)  , 0;    0,1];   %2*2 matrix

        OUT=(abs(Input*split*delay*com)).^2;     %1*2 matrix
       
        out1=[out1,OUT(1)+noise1(i,N)];
        
        out2=[out2,OUT(1)+noise2(i,N)];
        
    end
    
    w1=[w1,out1'];   %matrice de transfre du bras1
    w2=[w2,out2'];   %matrice de transfre du bras2
end

%choose one of the tow matrices to work with
%% the transfer matrix for the first MZI (sinusoidal)
% (to extract Littrow condition) (the maximum periodically)

figure
plot(lambda,w1(:,1),'DisplayName','arm1')
hold on
%plot(lambda,w2(:,1),'DisplayName','arm2')
legend('show')
title('the output of the first MZ'); xlabel('wavelenghts'); ylabel('Power'); 

%% The transfer matrix 
figure
imagesc(MZ,lambda,w1); 
colormap jet
title('the transfer matrix at 0dB/cm'); xlabel('number of MZ'); ylabel('wavelengths');

%% Original input spectrum

lambda_in=1550;
sep=0.037;
Pin=[zeros(size(w1(:,1)'))];
[l0,x1] = (min(abs(lambda - lambda_in)));
[l1,x2] = (min(abs(lambda - (lambda_in+sep))));

%x=find(lambda==lambda_in);
%x2=find(lambda==lambda_in+0.0325);                   %(add small delta_lambda)

Pin(x1)=1; 
Pin(x2)=0;
Pout=Pin*w1;
%the output
figure
plot(MZ,Pout)
title('the output as function of MZ'); xlabel('number of MZ'); ylabel('Power');
legend('show')

%% spectrum retrieval
PinRet=abs(Pout*pinv(w2));

figure
plot(lambda,PinRet./max(PinRet),'DisplayName','Retrieval')
hold on 
plot(lambda,Pin,'DisplayName','Original signal')                             %plot the original input with the retrieval spectrum
hold off
title('the input spectrum retrieval') ;xlabel('wavelengths');
legend('show')
text(lambda_in,0.9,num2str((x2-x1)*step));

%% filtring noise by convolution 
sigm=0.01;

gaus=exp(-0.5*((lambda-1550)/sigm).^2)*1/(sqrt(2*pi)*sigm);
c=conv(gaus,PinRet,'same');
figure
plot(lambda,c)
title('the filttred signal')

x=rand(100,1);         %Generating a vector with 100 random elements.
x0=0.321;              %Value specified 
y=find((x-x0)<0.5);    %Finding indices of all elements whose difference with the value specified is 0.5. (Here my error tolerance is 0.5)
z=x(y);                %Creating a vector with the values that satisfy the above condition
closestval=min(z-x0);  %Finding the value from this vector that is closest to the value specified.
