% 5 and 0km fibers
% filename = './measurements/20180717-5-0-km-fibers-arrival-times.csv'

% 5 and 5km fibers
filename = './measurements/20180717-5-5-km-fibers-arrival-times.csv'

% 10 and 0km fibers
% filename = './measurements/20180717-10-0-km-fibers-arrival-times.csv'

% 10 and 10km fibers
% filename = './measurements/20180717-10-10-km-fibers-arrival-times.csv'

M = csvread(filename);
dM = sort(M(1:end-1,2)) - sort(M(2:end,2));

figure(1)
subplot(2,1,1)
plot(sort(M(:,2)),'.')
subplot(2,1,2)
plot(dM,'.')