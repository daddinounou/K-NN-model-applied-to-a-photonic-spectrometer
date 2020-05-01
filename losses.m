close all
%neff=3.49769;
L1=1e7;
DeltaL0=365*1e3;
%Beta = neff*2*pi;


res1=0;
res2=0;
lambda= [1548:0.001:1552];
neff0=-1.6*1e-3*lambda+4.5197;                    %n is function of " lambda " here for width 450nm
delta_neff=0.0226;                                 %fluctuation of n due to fluctuation of width(+-20nm)                                
% alpha=log(10)*1e-7;
alpha=[0:5:20]*(log(10)/10)*1e-7;                           % alpha varie entre 10 et 20dB


for j=1:length(alpha);
    w=[];
    for N=1:32
        L2=L1+N*DeltaL0;    
        DeltaL=L2-L1;
        out1=[];
        out2=[];
        neff=neff0-delta_neff+2*delta_neff*rand;      %fluctuation for n is added for each MZ "width"
        for i=1:length(lambda)                        %change in the width dosn't change the slope neff(lambda)
            Beta=2*pi*neff(i)/lambda(i);
            res1=0.5*exp(-2*alpha(j)*L2)*(0.5*(1+exp(-2*alpha(j)*DeltaL))+exp(-alpha(j)*DeltaL)*cos(Beta*DeltaL));
            %res1=0.5*(1+cos(Beta*DeltaL));            %the propagation part 
            out1=[out1,res1];    
        end
        w=[w,out1'];
    end
    figure
    imagesc([1:32],lambda,w); 
    colormap jet
    title('the matrix transfer of one port config')    
end   

