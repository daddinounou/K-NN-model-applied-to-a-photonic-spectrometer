clear all

neff=2.5;
L1=1;
DeltaL=0.2;

sq=1/sqrt(2);
split=[[sq,-1i*sq];[-1i*sq, sq]];
com=split;

Beta = neff*2*pi;
Input=[0,1]; %we use only the lower arm. 1 is to say that there is signal,
% no way to specify a vector of lambda at the input 'cause lambda is
% already in the delay matrix of each MZ (which all together make the transfere matrix) 

%%

w=[];
res1=[];
res2=[];

for N=1:30
    L2=L1+N*0.2;    
    DeltaL=L2-L1;
    out1=[];
    out2=[];
    wave=[];
    for lambda= (1550:0.0002:1550.5)*10^-7

        delay= exp(-1i*Beta*L2/lambda)*[exp(-1i*Beta*DeltaL/lambda),0;0,1];

        OUT=Input*(abs(split*delay*com));
        %OUT is (1*2) matrix (it is multiplied by input),ici on a la sortie
        %(2 bras) directement, le probleme c'est qu'on ne peut pas avoir
        %accès à la matrice de transfaire du système

        out1=[out1,OUT(1)];
        
        out2=[out2,OUT(2)];
        
        wave=[wave,lambda];
                
        if lambda==1.5502e-04
            %this is just to verify the quality of our matix, since we
            %injected a 1 in the first arm and 0 in the second, we have to
            %find this result when using the output with the inverse of
            %each MZ at a given wavelength.
                    
            matrix=(abs(split*delay*com));
            %matrix x is (2*2)(not multiplied by input) c'est la matrice
            %élémentaire dans la grande matrice de transfaire, c'est le problème
            % de visualiser 2 portes simultanément.
            
            Inp=(OUT/(matrix));
            %this has tow components (2*1), 
            %don't confuse this calculated Inp with the real one Input
            %it tels us that for all MZ, only the upper arm is used
            
            res1=[res1,Inp(1)];
            res2=[res2,Inp(2)];            

        end
        
    end
    
    w=[w,out1'];
end

% this corresponds only to one output, we can do the same for the other
%we represent the transfer matrix ( this is the output precisly speaking
%because we tooke the input into accont ) (wavelength, MZI, output)

figure(1);
imagesc(w); 
colormap jet

yticklabels = (1550.5:-0.05:1550)*10^-7;
yticks = linspace(1, size(w, 1), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel',flipud(yticklabels(:)),'YDir','reverse');
title('the signal at the output of arms 1')
% it is also its transfaire matrix 'cause it's just multiplied by 1 :) 

%%

% verify the input
figure(2)
plot(res1)
title('the signal at the upper arm for all the MZs')
figure(3)
plot(res2)
title('the signal at the lower arm')

%%

% representation of the output of all MZ as a function of w, the output at Littrow
% is max and is constant 
for i=1:30
    figure(4)
    plot(wave',w(:,i))
    hold on 
end
hold off

title('superposition of outputs of all MZs')

%%
%representation of the output of all MZ at Littrow wavelengt
% after multiple tries I found that the littrow condition correspends to
% the 1939th line (for this vector of lambda!)
% despite we see the outputs 1, Matlab see them slightly diffrent !

figure(5)
plot(w(1939,:))

title('the output at Littrow condition')




