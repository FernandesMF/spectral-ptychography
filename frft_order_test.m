clear all

x = -10:0.1:10; x = x-mean(x);

f1 = abs(hermiteH(3,x)).*exp(-x.^2/(2*0.5^2));

a1 = 0.5;
s1 = 1;
s2 = 0.5;
f2 = abs(hermiteH(3,x*s2)).*exp(-(x*s2).^2/(2*0.5^2));

a2 = 2*atan( tan(pi*a1/2)*s1^2/(s2^2) )/pi;
bar = cos(a2*pi/2)/cos(a1*pi/2);

F1 = frft(f1,a1);
F2 = frft(f2,a2);
F3 = frft(f2,a1);

F1 = F1/sum(abs(F1));
F2 = F2/sum(abs(F2)/bar*s2);
F3 = F3/sum(abs(F3));

foo = 1:length(x); foo = foo-mean(foo);

figure(1)
plot(x*s1,f1,...
     x*s2,f2)

figure(2)
subplot(2,1,1)
plot(foo,abs(F1),...
     foo/bar*s2,abs(F2),...
     foo,abs(F3))
subplot(2,1,2)
plot(foo,angle(F1),...
     foo*bar,angle(F2),...
     foo,angle(F3))
 
