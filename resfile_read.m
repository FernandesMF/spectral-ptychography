clear all

%% Parameters
% Results file path
% res_file	= './results/noisy-data/d10_r6_a0.72_lam1.00e+07-2019-5-7_res.txt';
% res_file	= './results/noisy-data/d11_r6_a0.89_lam1.00e+07-2019-5-22_res.txt';
% res_file	= './results/noisy-data/d16_r9_a0.79_lam1.00e+07-2019-5-22_res.txt';
% res_file	= './results/noisy-data/d20_r11_a0.817_lam1.00e+07-2019-5-22_res.txt';
res_file	= './results/noisy-data/d32_r17_a0.847_lam1.00e+07-2019-5-23_res.txt';

% Figure parameters
prin	= 1;        %  print to png at the end?
font	= 16;       % fontsize for axes
wid = 800;          % width and height (in px) of the figure
hei = 600;


% Reading parameters from file name
[foo,count] = sscanf(res_file,'./results/noisy-data/d%d_r%d_a%f_lam%e-%d-%d-%d_res.txt');
par.d   = foo(1);
par.r	= foo(2);
par.a   = foo(3);
par.lam	= foo(4);
par.yy  = foo(5);
par.mm  = foo(6);
par.dd  = foo(7);



%% Opening external files

resFID  = fopen(res_file,'r');
if(resFID==-1)
    error('Could not open results file');
end



%% Reading external files
% If Data is used just like it is read:
%   Data(1,:) has the run number (1,...,10 000)
%   Data(2,:) has the Hilbert-Schmidt distance between final estimate and correct answer
%   Data(3,:) has the fidelity between final estimate and correct answer

Data	= fscanf(resFID,'%d\t\t%f\t\t%f\n',[3,inf]);
dist    = Data(2,:);
fids    = Data(3,:);
nsts    = size(Data,2);



%% Histograming fidelities
figure(1);
h   = histogram(fids,'Normalization','probability');
f   = gcf;
a   = gca;
a.XLim(2)   = 1;

% Making figure pretty
f.Position  = [f.Position(1) f.Position(2) wid hei]
a.FontSize  = font;
xlabel('Reconstruction fidelity')
ylabel('Proportion of results')

% Print to png file
if(prin)
    flag	= 1;
    % png file name
    foo	= ['./results/histograms/d' num2str(par.d) '_r' num2str(par.r) 'lam' ...
           num2str(par.lam) '.png' ]
    
    % check if file already exists, and asks to overwrite if so
    if(exist(foo,'file'))
        flag	= questdlg('PNG file already exists. Overwrite?','OVERWRITE?',0,1,0);
    end
    
    % prints to png if file doesn existe yet, or if overwriting is desired
    if(flag)
        print(foo,'-dpng')
    end
end