N = 100000 ;
figure
hist(rand(N,1),100)
title('rand:Uniform distribution')
figure
hist(randn(N,1),100)
title('randn:Normal distribution')