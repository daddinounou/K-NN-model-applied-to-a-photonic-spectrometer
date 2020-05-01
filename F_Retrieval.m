

function y=F_Retrieval(Interferogram,matrix,lambda,lambda_in,sigm)   

    PinRet=abs(Interferogram*(pinv(matrix)));
    gaus=exp(-0.5*((lambda-lambda_in)/sigm).^2)./(sqrt(2*pi)*sigm);
    c=conv(gaus,PinRet,'same');
    y=c;
end  
