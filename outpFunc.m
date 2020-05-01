
function y=outpFunc(lambda)
    neff=2.5;
    L1=1;
    L2=1.2;
    DeltaL=0.2;

    sq=1/sqrt(2);
    split=[[sq,-i*sq];[-i*sq, sq]];
    com=split;

    Beta = neff*2*pi;
    Input=[0,1];
        
    delay= exp(-i*Beta*L2/lambda)*[exp(-i*Beta*DeltaL/lambda),0;0,1];
    
    x=(abs(Input*split*delay*com)).^2;
    
    y=x(2);
end
        
        
  