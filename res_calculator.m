
clear all
tic
deltaL0=365e3;                      % the path length delay increment 
MZ=0:31;                            % number of MZ interferometers
step=0.0001;
lambda=1550:step:1551;
%lambda=1549.642:step:1550.359;      % the wavelength range (FSR 0.717) 
lambda_in=1550.45;                     % an input wavelength
stdnoise=0.00;                      % the std of noise in gaussian distribution
%% Gaussian filter
sigm=0.002;                                                                % This was claibrated to get almost the
gaussFilter = exp(-(lambda-mean(lambda)).^ 2 / (2 * sigm ^ 2));            % original retireval wehn no noise is added
gaussFilter = gaussFilter / sum (gaussFilter); % normalize

%% This script gives the reolutions as function of losses for a given amount of noise

faberr=rand(1,length(MZ));                                                
resolutions=[]; ecarts=[];

for loss=0:0                                                               % swip the losses, theis loop gives the average resolution for each loss value
    vec=[];                                                                % initialize the resolutions vector
    for i=1:1                                                              % run multiple simulations (20) for each loss value in order to calculate the average resolution
        w1=F_transfer_matrix(deltaL0,MZ,lambda,stdnoise,loss,300,faberr);  % calibration matrix used to calculate the interferogram
        w2=F_transfer_matrix(deltaL0,MZ,lambda,stdnoise,loss,300,faberr);  % calibration matrix used to calculate the retrieval
        c=ones(1,length(lambda));
        [l0,x1] = (min(abs(lambda - lambda_in))); x2=x1;                   % finds the index of the input wavelength, initialize the 2nd
        resol=0;
        pos1=x1; pos2=x2;
        peak1=1; peak2=1;
        %while c(fix((pos1+pos2)/2))/max(c) > 0.5  || (c(pos1)/max(c))<0.8  || (c(pos2)/max(c))<0.8  % a very subjective way to approximate the resolution !!!
%         while c(fix((x1+x2)/2))/max(c) > 0.5  || (c(x1)/max(c))<0.8  || (c(x2)/max(c))<0.8  % more adequat when sevral peaks exist
%                                                                            %This loop is to find the resolution that fullfil the above
%                                                                            %3 criteria, by moving the input 2
%             x2=x2+1;
%             if x2>length(lambda)                                           % if the indice x2 exceeds the length of the spectrum we stop moving input2
%                 break                                                      
%             end                
            Pin=zeros(1,length(lambda));
            Pin(x1)=1;     Pin(x2)=1;                                      % the good separation input found, build the input vector
            Pout=Pin*w1;                                                   % calculate the interferogram
            PinRet=abs(Pout*(pinv(w2)));                                   % calculate the retrieval
            c=conv(gaussFilter,PinRet,'same');                             % convolution with the gaussian filter  
            [peaks,pos]=findpeaks(c,'sortStr','descend');                  % find the peaks and their positions of the retrieved signal(filtred)
            pos1=pos(1);  pos2=pos(2);                                     % get the 2 big peaks
            %resol=abs(pos2-pos1)*step;                                      % to calculate the resolution
            resolu=(x2-x1)*step;
            
        end
        vec=[vec,resol];                                                   % the vector of resolutions for a given loss
    %end
    res=mean(vec);                      ecart=std(vec);
    resolutions=[resolutions,res];      ecarts=[ecarts,ecart];             % vector of tesolution (and std) for all losses
end
toc

%% transmittance function T=f(lambda)
% figure 
% plot(lambda,w1(:,26),'DisplayName','n_{eff}(\lambda)')
% hold on 
% xlabel('wavelength(nm'); ylabel('Transmittance');
% legend('show')

%% The transfer matrix 
figure
imagesc(MZ,lambda,w1); 
colormap jet
xlabel('number of MZI'); ylabel('wavelengths(nm)');

%% Interferogram I=f(MZ) arround Littrow

% figure
% plot(MZ,w1(4000,:),'DisplayName',num2str(1549.642+step*(1000-1)))
% legend('show')
% ylim([0 1.1])    
% xlabel('number MZI'); ylabel('Transmittance');


%%
% shift1=(min(pos1,pos2)-x1)*step;  shift2=(max(pos1,pos2)-x2)*step;         % deviation of each peak from the original input (because of T..)
% figure
% plot(lambda,PinRet./max(PinRet),'DisplayName','Retrieval')
% hold on
% plot(lambda,c./max(c),'DisplayName','smoothen Ret')
% hold on
% plot(lambda,Pin, 'DisplayName','Input')
% hold on
% plot(lambda,gaussFilter./max(gaussFilter),'DisplayName','gussian');
% hold off 
% legend('show'); xlabel('wavelength(nm)')
% %title(['resolution=',num2str(res*1e3),'pm')
% %title(['resolution=',num2str(res*1e3),'pm','     loss =',num2str(loss),' dB/cm'])
% title(['noise=',num2str(stdnoise*100),'%','     resolution=',num2str(resolu*1e3),'pm','    loss dB/cm=',num2str(loss)])
 set(gca,'FontSize',14,'FontWeight','bold')


