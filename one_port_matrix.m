close all
clear all

neff=4.38;
L1=1e7;
DeltaL0=365*1e3;
Beta = neff*2*pi;

w=[];
res1=0;
res2=0;
lambda= [1548:0.001:1552];

tic
for N=1:32
    L2=L1+N*DeltaL0;    
    DeltaL=L2-L1;
    out1=[];
    out2=[];
    for i=1:length(lambda)
        res1=0.5*(1+cos(Beta*DeltaL/lambda(i)));
        out1=[out1,res1];    
    end
    w=[w,out1'];
end
toc
%% the transfer matrix for the first MZI (sinusoidal)
% to extract Littrow condition
figure(1)
plot(lambda,w(:,1))
title('the 1st colomn of the matrix MZ1')
% the Littrow condition at 1364.5 nm for 50 MZ
%% The transfer matrix 
figure(2);
imagesc([1:32],lambda,w); 
colormap jet
title('the transfer matrix'); xlabel('number of MZ'); ylabel('wavelengths');
%% The input spectrum

lambda_in=1550;

Pin=[zeros(size(w(:,1)'))];
x=find(lambda==lambda_in);
x
x2=find(lambda==lambda_in+0.046962);

Pin(x)=1; 
Pin(x2)=1;
Pout=Pin*w;
figure(7)
plot([1:32],Pout)
title('the intensity at the output'); xlabel('number of MZ1')
%% the output interferogram 


%% retrieve the spectral input from the interferogram
PinRet=Pout*pinv(w);

figure(6)
plot(lambda,abs(PinRet./max(PinRet)))
hold on 
plot(lambda,Pin)
hold off
title('retrieval spectrogram'); xlabel('wavelengths');


%% retrival for a polychromatic input

% Pout=Pin*w;
% PinRet=Pout*pinv(w);
% figure(3)
% plot(PinRet)
% title('retrieval of polychromatic input')
% %% retrieval for a lot of monochromatic input plot together
% 
% for i=1:5:1000
%     Pin=[zeros(size(w(:,1)'))];
%     Pin(i)=100;
%     Pout=Pin*w;
%     PinRet=Pout*pinv(w);
%     figure(4)
%     plot(PinRet,'DisplayName',num2str(i));
%     hold on
% end
% hold off
% %legend('show')
% title('retrieval of lot of monochromatic inputs plot together')
% 
% %% plot the output at a certain wavelength (sinusoidal output if monochromatic input)  
% 
% figure(6)
% plot(Pout)





