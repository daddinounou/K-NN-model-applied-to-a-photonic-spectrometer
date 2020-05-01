%here we use the simplified expression for the transfere matrix at each%port

close all
%neff=3.49769;
L1=1e7;
DeltaL0=365*1e3;
%Beta = neff*2*pi;

w=[];
res1=0;
res2=0;
lambda= [1548:0.001:1552];
neff0=-1.6*1e-3*lambda+4.5197;                    %n is function of " lambda " here for width 450nm
delta_neff=0.0226;                                %fluctuation of n due to fluctuation of width(+-20nm)                                

%alpha=log(2)/(2*(L1+DeltaL0));                    
alpha= ( log(10)/20 ) * 1e-7 *  (1e-3);
tic

for N=1:32
    L2=L1+N*DeltaL0;    
    DeltaL=L2-L1;
    out1=[];
    out2=[];
    neff=neff0-delta_neff+2*delta_neff*rand;      %fluctuation for n is added for each MZ " due to width fluctuation"
    for i=1:length(lambda)                        %change in the width dosn't change the slope neff(lambda)
        Beta=2*pi*neff(i)/lambda(i);
        res1= 0.5*exp(-2*alpha*L2) * ( 0.5*(1+exp(-2*alpha*DeltaL) )-exp(-alpha*DeltaL)*cos(Beta*DeltaL));  %propagation with losses
        %res1=0.5*(1+cos(Beta*DeltaL));           %the propagation part without losses
        out1=[out1,res1];    
    end
    w=[w,out1'];
end
toc

%% the transfer matrix for the first MZI (sinusoidal)
% (to extract Littrow condition)
figure
plot(lambda,w(:,1))
title('the amplitude against the wavelength MZ1')
%% The transfer matrix 
figure
imagesc([1:32],lambda,w); 
colormap jet
title('the matrix transfer of one port config')
%% The original input spectrum

lambda_in=1550;
sep=0.045;

Pin=[zeros(size(w(:,1)'))];

[l0,x] = (min(abs(lambda - lambda_in)))
[l1,x2] = (min(abs(lambda - (lambda_in+sep))))

%x=find(lambda==lambda_in);
%x2=find(lambda==lambda_in+0.040);                   %(add small delta_lambda)

Pin(x)=1; 
Pin(x2)=1;
Pout=Pin*w;
% the output
figure
plot([1:32],Pout)
title('the output as function of MZ')

%% retrieve the spectral input from the interferogram
PinRet=abs(Pout*pinv(w));

figure
plot(lambda,PinRet./max(PinRet))
hold on 
plot(lambda,Pin)                                 %plot the original input with the retrieval spectrum
hold off
text(lambda_in,0.9,num2str(sep));
title('retrieval of monochromatic input')


%% filtring the noise by convolution 
sigm=0.03;

gaus=exp(-0.5*((lambda-1550)/sigm).^2)*1/(sqrt(2*pi)*sigm);
c=conv(gaus,PinRet,'same');
figure
plot(lambda,c)
title('the filttred signal')

% annother method to filter using 'gausswin', we can also use 'filter(PinRet,1,y)'
% y=gausswin(4001,50);             %cette fonction comme comme 'moyenne' le centre de la longueur donn�e donn�
% c=conv(y,PinRet,'same');
% figure 
% plot(lambda,y)
% 
% figure
% plot(lambda,c)

