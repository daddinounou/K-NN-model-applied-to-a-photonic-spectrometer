
%clear all
tic
deltaL0=365e3;                      % the path length delay increment 
MZ=0:31;                            % number of MZ interferometers
step=0.0001;
lambda=1549.3250:step:1550.0425;
stdnoise=0.00;                      % the std of noise in gaussian distribution
loss=0;
faberr=rand(1,length(MZ)); 

 

% DATA1=[];  DATA2=[];  DATA3=[];  DATA4=[];
% 
% for T=300:0.1:305
%     w=F_transfer_matrix(deltaL0,MZ,lambda,stdnoise,loss,T,faberr);
%     larg=0.01; 
%     Pin = 1-exp(-(lambda-1549.4).^ 2 / (2 * larg ^ 2));   Pin = Pin / sum (Pin);
%     Pout=Pin*w;
%     DATA1=[DATA1;Pout];
% end
% 
% 
% for T=300:0.1:305
%     w=F_transfer_matrix(deltaL0,MZ,lambda,stdnoise,loss,T,faberr);
%     larg=0.01; 
%     Pin = 1-exp(-(lambda-1549.7).^ 2 / (2 * larg ^ 2));   Pin = Pin / sum (Pin);
%     Pout=Pin*w;
%     DATA2=[DATA2;Pout];
% end
% 
% 
% for T=300:0.1:305
%     w=F_transfer_matrix(deltaL0,MZ,lambda,stdnoise,loss,T,faberr);
%     larg=0.01; 
%     Pin = 1-exp(-(lambda-1550).^ 2 / (2 * larg ^ 2));   Pin = Pin / sum (Pin);
%     Pout=Pin*w;
%     DATA3=[DATA3;Pout];
% end
% 
% 
% for T=300:0.1:305
%     w=F_transfer_matrix(deltaL0,MZ,lambda,stdnoise,loss,T,faberr);
%     larg=0.01; 
%     Pin = ones(1,length(lambda));
%     Pout=Pin*w;
%     DATA4=[DATA4;Pout];
% end
% 
% toc


%%%%%%%%%%%
w=F_transfer_matrix(deltaL0,MZ,lambda,stdnoise,loss,300,faberr);
larg=0.01; 
Pin = 1-exp(-(lambda-1549.5).^ 2 / (2 * larg ^ 2));   Pin = Pin / sum (Pin);
Pout=Pin*w;
%%%%%%%%%%%

load('DATA1');load('DATA2');load('DATA3');
DAT=[DATA1;DATA2;DATA3];

labels=1:153;      % correspendant à 300:0.1:305;
x=DAT(137,:);   %  300.9
%x=Pout;

y=zeros(1,size(DAT,1));

for i=1:size(DAT,1)
    distance=sqrt(sum((DAT(i,:)-x).^2));
    y(i)=distance;
end

[dist,label]=min(y);     
label
    
    
% 
% for i=1:length(DATA1(:,1))
%     plot(MZ,DATA1(i,:),'DisplayName',num2str(i))
%     hold on 
% end
% legend('show')

      
figure
for i=1:5:31
    plot(T,DATA1(:,i),'DisplayName',num2str(i))
    hold on
end

legend('show')


