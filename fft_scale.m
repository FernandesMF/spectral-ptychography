clear all
close all

dt  = 0.002;
T   = 2;
t   = 0:dt:T;       % s
N   = length(t);
Fs  = 1/dt;         % Hz
df  = Fs/N;
f   = (0:N-1)*df;	% Hz
f   = f-mean(f);

w   = 60;           % Hz
x   = 1*sin(2*pi*w*t) + 0.1*rand(1,N);
X   = fft(x);
X   = fftshift(X);

figure(1)
subplot(2,1,1)
plot(t,x)
subplot(2,1,2)
plot(f,abs(X))