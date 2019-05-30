

a   = 0.99
X   = rand(10,1);
Y   = frft(X,a);
YY  = frft(Y,-a);

l   = length(X);
figure(1);
plot(1:l,X,1:l,abs(Y),1:l,abs(YY))
legend('X','Y','YY')