%Transfer matrix // specifications:
%Each MZ has 2 input ports(only the upper one is used) 
%and 2 output ports(only the lower is used). 
%This matrix is specific for 450nm width strip waveguide of Si.
%for given fabrication errors: +-20 nm .
function y=F_transfer_matrix(...
    DeltaL0,...                         % increment length (nm)
    MZ,...                              % 1*N vector giving the number of interferometers
    lambda,...                          % 1*M vector of wavelengths with small step (nm) 
    stdnoise,...                        % the std of gaussian giving noise vlaues (should be smaler than the signal:1)
    loss,...                            % loss in dB/cm, 
    T,...                               % Temperature 
    faberr...                           % fabrication errors must be the same for both matrices
)          
    % Number of given paramters must exactly match the num of possible parameters
    if nargin < nargin('F_transfer_matrix')
        error('Too few arguments.');
    end
    if nargin > nargin('F_transfer_matrix')
        error('Too many arguments.');
    end
    L1=500e3;                           % the length of the upper arm
    split=[1,-1i;-1i,1]/sqrt(2);        % splitter matrix  
    com=split;                          % combiner matrix                      
    Input=[1,0];                        % before the splitter there is 2 ports, we use only the lower one
    %neff0=2.002*ones(1,length(lambda));                        %l'indice effective à 300k et 1550nm
    %neff0=-1.6694*1e-3*lambda+4.589;                            %neff at 300k 
    neff0=1.5564*1e-7*T^2+7.559*1e-5*T-1.6694*1e-3*lambda+4.55327;   % Model of neff consedring the T and lambda, 450nm width strip waveguide
    delta_neff= 0.0226;                 % fluctuation of n due to fluctuation of width(+-20nm in Lumerical) (delta_n to be added randomly) 
   
    w=zeros(length(lambda),length(MZ));            % empty matrix 
    alpha=( log(10)/20 ) * 1e-7 *loss;            % loss in dB/cm, alpha in 1/nm
    for j=1:length(MZ)
        L2=L1+MZ(j)*DeltaL0;    
        DeltaL=L2-L1;                             % the diferrential pathlength          
        
        neff=neff0-delta_neff+2*delta_neff*faberr(j);  % adding fabrication errors     
        for n=1:length(lambda)            
            Beta=2*pi*neff(n)/lambda(n);
            delay= exp(-alpha*L2) * exp(-1i*Beta*L2) *[ exp(alpha*DeltaL) * exp(1i*Beta*DeltaL)  , 0;    0,1]; 
            OUT=(abs(Input*split*delay*com)).^2;  % 1*2 vector giving the power at the output ports
            w(n,j)=OUT(2)+stdnoise.*randn;
        end
    end 
    y=w;
end    
