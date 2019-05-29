clear all
% close all



% Creating function, point-spread function and blurred signal
x = -10:0.1:10;
y = abs(hermiteH(3,x)).*exp(-x.^2); %sig = 1/sqrt(2) ~ 0.7
y = y/sum(y);

psfun = @(x) exp(-x.^2/(2*(0.4*0.7)^2));
psf = psfun(x);

yblurred = conv(y,psf); yblurred = yblurred/(sum(yblurred));
dx = abs(x(2)-x(1)); xnew = (1:length(yblurred))*dx; xnew = xnew-mean(xnew);
psfnew = psfun(xnew);



% Applying several levels of poissonian noise to the blurred signal, calculating power spectra
N = 1e03;   %number of signals to be made
avcounts = [1e02,1e03,1e04,1e05];  %average counts (total)
ynoise(1:length(avcounts)) = struct('sigs',zeros(N,length(yblurred)),...
                                    'devs',zeros(N,length(yblurred)),...
                                    'Fmean',zeros(size(yblurred)),...
                                    'Fstd',zeros(size(yblurred)),...
                                    'DFmean',zeros(size(yblurred)),...
                                    'DFstd',zeros(size(yblurred)),...
                                    'wfilt',zeros(size(yblurred)),...)
                                    'dcvsigs',zeros(N,length(yblurred)));    %reserving memory to the struct
F = abs(fft(yblurred));

for r=1:numel(ynoise)   %making N signals
    for s=1:N
        ynoise(r).sigs(s,:) = poissrnd(avcounts(r)*yblurred);
        ynoise(r).devs(s,:) = ynoise(r).sigs(s,:)-yblurred*avcounts(r);
    end
    foo = fft(ynoise(r).sigs,[],2);
    ynoise(r).Fmean = mean(abs(foo));
    ynoise(r).Fstd = std(abs(foo));
    foo = fft(ynoise(r).devs,[],2);
    ynoise(r).DFmean = mean(abs(foo));
    ynoise(r).DFstd = std(abs(foo));
end

% Calculating Wiener fiters, deconvoluting signals
fest = abs(hermiteH(3,xnew)).*exp(-xnew.^2); fest = fest/sum(fest);
Fest = abs(fft(fest));    %estimate of power spectrum of the signal
Psfnew = fft(psfnew);
for r=1:numel(ynoise)
    ynoise(r).wfilt = conj(Psfnew).*Fest./( abs(Psfnew).^2.*Fest + ynoise(r).DFmean/avcounts(r));
    foo = fft(ynoise(r).sigs,[],2);
    for s=1:N
        ynoise(r).dcvsigs(s,:) = ifft( fft(ynoise(r).sigs(s,:)).*ynoise(r).wfilt );
    end
end



%% Plots

% Original wavefunction, point-spread function and blurred wavefunction
figure(1)
subplot(3,1,1)
plot(x,y)
h = gca; h.YLabel.String = 'Original Wavefun.';
subplot(3,1,2)
plot(x,psf)
h = gca; h.YLabel.String = 'PSF';
subplot(3,1,3)
plot(xnew,yblurred)
h = gca; h.YLabel.String = 'Blurred Wavefun.';

%Poisson-noised wavefunctions (first signal)
foo = numel(ynoise);
figure(2)
subplot(foo,1,1)
title('Blurred wavefunctions with poissonian noise (example)')
for r=1:foo
    subplot(foo,1,r)
    plot(xnew,ynoise(r).sigs(1,:));
    h = gca; h.YLabel.String = [num2str(avcounts(r),'%1.1g') ' counts'];
end

% Power spectra of noised-wavefunctions (first signal)
dk = 1/dx; N = length(xnew);
k = (0:N-1)*dk; k = k-mean(k);
figure(3)
title('Power spectra of noisy wavefunctions (example)')
for r=1:foo
    subplot(foo,1,r)
    plot(k,fftshift(abs(fft(ynoise(r).sigs(1,:)))))
    h = gca; h.YLabel.String = [num2str(avcounts(r),'%1.1g') ' counts'];
%     set(gca,'YLim',[0 1.5])
end

% Power spectra of noised-wavefunctions (mean)
dk = 1/dx; N = length(xnew);
k = (0:N-1)*dk; k = k-mean(k);
figure(4)
for r=1:foo
    subplot(foo,1,r)
    plot(k,fftshift(ynoise(r).Fmean))
    h = gca; h.YLabel.String = [num2str(avcounts(r),'%1.1g') ' counts'];
%     set(gca,'YLim',[0 1.5])
end

% STD of power spectra of noised wavefunctions
figure(5)
for r=1:foo
    subplot(foo,1,r)
    plot(k,fftshift(ynoise(r).Fstd))
    h = gca; h.YLabel.String = [num2str(avcounts(r),'%1.1g') ' counts'];
%     plot(k,(fftshift(abs(fft(ynoise(r).sigs(1,:)/avcounts(r))))-F)*avcounts(r))
%     set(gca,'YLim',[0 1.5])
end

% Power spectra of deviations (mean)
dk = 1/dx; N = length(xnew);
k = (0:N-1)*dk; k = k-mean(k);
figure(6)
for r=1:foo
    subplot(foo,1,r)
    plot(k,fftshift(ynoise(r).DFmean))
%     set(gca,'YLim',[0 1.5])
end

% STD of power spectra of deviations
figure(7)
for r=1:foo
    subplot(foo,1,r)
    plot(k,fftshift(ynoise(r).DFstd))
%     plot(k,(fftshift(abs(fft(ynoise(r).sigs(1,:)/avcounts(r))))-F)*avcounts(r))
%     set(gca,'YLim',[0 1.5])
end

% Wiener filters
figure(8)
for r=1:foo
    subplot(foo,1,r)
    plot(k,fftshift(abs(ynoise(r).wfilt)));
end

% Denoised wavefunctions (first signal)
figure(9)
for r=1:foo
    subplot(foo,2,2*r-1)
    plot(xnew,ynoise(r).sigs(1,:))
    h = gca; h.YLabel.String = [num2str(avcounts(r),'%1.1g') ' counts'];
    subplot(foo,2,2*r)
    plot(xnew,fftshift(abs(ynoise(r).dcvsigs(1,:))),...
         x,sum(ynoise(r).dcvsigs(1,:))*y,'--k');
end

% seminar figures
% figure(10)
% subplot(1,2,1)
% plot(xnew,ynoise(r).sigs(1,:),'LineWidth',1.5)
% h = gca; h.YLabel.String = [num2str(avcounts(r),'%1.1g') ' counts'];
% subplot(1,2,2)
% plot(xnew,fftshift(abs(ynoise(r).dcvsigs(1,:))),...
%      x,sum(ynoise(r).dcvsigs(1,:))*y,'--k','LineWidth',1.5);
%  
% figure(11)
% subplot(3,1,1)
% plot(x,y,'LineWidth',1.5)
% h = gca; h.YLabel.String = 'Original Wavefun.';
% subplot(3,1,2)
% plot(x,psf,'LineWidth',1.5)
% h = gca; h.YLabel.String = 'PSF';
% subplot(3,1,3)
% plot(xnew,yblurred,'LineWidth',1.5)
% h = gca; h.YLabel.String = 'Blurred Wavefun.';