I=[];
wave=[];
for lambda= (1500:0.01:1510)*10^-7
    
    e=outpFunc(lambda);
    
    I=[I,e];
    wave=[wave,lambda];
end

plot(wave,I)
