
function [x,y]=F_Interferogram(w,lambda,lambda_in,sep)
    Pin=[zeros(size(w(:,1)'))];
    [l0,x1] = (min(abs(lambda - lambda_in)));
    [l1,x2] = (min(abs(lambda - (lambda_in+sep))));

    Pin(x1)=1; 
    Pin(x2)=1;
    Pout=Pin*w;
    x=Pin;
    y=Pout;  
end



