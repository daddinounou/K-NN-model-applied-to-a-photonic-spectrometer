%Transfer matrix // specifications:
%Each MZ has 2 input ports(only the lower is used) 
%and 2 output ports(only the upper is used). 
%This matrix is specific for 450nm width strip waveguide of Si.
%for given fabrication errors: +-20 nm .
function y=F_transfer_matrix(...
    DeltaL0,...                         % increment length (nm)
    MZ,...                              % 1*N vector giving the number of interferometers
    lambda,...                          % 1*M vector of wavelengths with small step (nm) 
    stdnoise,...                        % the std of gaussian giving noise vlaues (should be smaler than the signal:1)
    loss...                             % loss in dB/cm, 
    T...                                % Temperature 
)          
    % Number of given paramters must exactly match the num of possible parameters
    if nargin < nargin('F_transfer_matrix')
        error('Too few arguments.');
    end
    if nargin > nargin('F_transfer_matrix')
        error('Too many arguments.');
    end
    
    L1=1e7;                             % the length of the upper arm 
    split=[1,-1i;-1i,1]/sqrt(2);        % splitter matrix  
    com=split;                          % combiner matrix                      
    Input=[0,1];                        % before the splitter there is 2 ports, we use only the lower one
    neff0=1.5564*1e-7*T1^2+7.559*1e-5*T1-1.6694*1e-3*lambda+4.55327;   % Model of neff consedring the T adn lambda, 450nm width strip waveguide
    delta_neff=0;
    %delta_neff=0.0226;                  % fluctuation of n due to fluctuation of width(+-20nm in Lumerical) (delta_n to be added randomly) 
    w=[];
    alpha=( log(10)/20 ) * 1e-7 *loss;            % loss in dB/cm, alpha in 1/nm
    for j=1:length(MZ)
        L2=L1+MZ(j)*DeltaL0;    
        DeltaL=L2-L1;                             % the diferrential pathlength
        out=[];  
        neff=neff0-delta_neff+2*delta_neff*rand;  % adding fabrication errors     

        for n=1:length(lambda)            
            Beta=2*pi*neff(n)/lambda(n);
            delay= exp(-alpha*L2) * exp(-1i*Beta*L2) *[ exp(alpha*DeltaL) * exp(1i*Beta*DeltaL)  , 0;    0,1]; 
            OUT=(abs(Input*split*delay*com)).^2;  % 1*2 vector giving the power at the output ports
            out=[out,OUT(1)+stdnoise.*randn];   % the power at the output of one port one MZ for each lambda            
        end
        w=[w,out'];                            % the transfer matrix giving the output for each MZ(one port) for each lambda
    end 
    y=w;
end    



