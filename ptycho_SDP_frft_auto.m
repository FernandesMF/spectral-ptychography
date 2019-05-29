clear all

% FIX: print num. of iterations in res. file
% FIX: implement noise

%% Setting parameters
d   = 10;      % Hilbert space dimension
r   = 6;       % rank of mask projectors; we will use the contiguous family with d projectors scheme
a   = [0.72];   % orders of frft that will be "measured" 

Nproj	= d;
Nobs	= d*Nproj*numel(a);     % ptychographic projectors * comp. basis projectors
Nstates = 10;
noise   = 0;                    % Noise mode: 0 -> none; 1-> poissonian
noise_av_counts  = 1e01;         % average counts of the poissonian noise (if desired)

foo   = clock;

if(~noise)
    file_par    = ['./results/ideal-data/d' num2str(d) '_r' num2str(r) '_a' num2str(a) '-' ...
                   num2str(foo(1)) '-' num2str(foo(2)) '-' num2str(foo(3)) '_par.txt'];
    file_res    = ['./results/ideal-data/d' num2str(d) '_r' num2str(r) '_a' num2str(a) '-' ...
                    num2str(foo(1)) '-' num2str(foo(2)) '-' num2str(foo(3)) '_res.txt']; 
else
    file_par    = ['./results/noisy-data/d' num2str(d) '_r' num2str(r) '_a' num2str(a) '_lam' ...
                   sprintf('%2.2e',noise_av_counts) '-' num2str(foo(1)) '-' num2str(foo(2)) '-' num2str(foo(3)) '_par.txt'];
    file_res    = ['./results/noisy-data/d' num2str(d) '_r' num2str(r) '_a' num2str(a) '_lam' ...
                   sprintf('%2.2e',noise_av_counts) '-' num2str(foo(1)) '-' num2str(foo(2)) '-' num2str(foo(3)) '_res.txt'];           
end



%% Opening files

% Testing for existence, asking about overwriting in case the files exist
if( exist(file_par,'file') || exist(file_res,'file') )
    foo     = questdlg('Parameter or result file already exists.  Overwrite them?', ...
                       'Overwrite existing files?','Yes','Abort!','Abort!');
    if(~strcmp(foo,'Yes'))
        error('Operation aborted');
    end
end

% Opening files
fid_par     = fopen(file_par,'w+');
fid_res     = fopen(file_res,'w+');
if( fid_par==-1 || fid_res==-1 )
    error('Could not open parameter or results file'); 
end

%% Writing parameters to file

fprintf(fid_par,'d   = %d;      %% Hilbert space dimension\n',d);
fprintf(fid_par,'r   = %d;      %% rank of mask projectors; we will use the contiguous family with d projectors scheme\n',r);
fprintf(fid_par,'a   = %f;      %% orders of frft that will be \"measured\"\n\n',a);
fprintf(fid_par,'Nproj	= %d;\n',Nproj);
fprintf(fid_par,'Nobs	= %d;   %% ptychographic projectors * comp. basis projectors\n',Nobs);
fprintf(fid_par,'Nstates = %d;\n',Nstates);
fprintf(fid_par,'noise   = %d;  %% Noise mode: 0 -> none; 1-> poissonian\n',noise);
fprintf(fid_par,'noise_av_counts  = %2.2e %% poissonian noise average counts (if desired)',noise_av_counts);


%% Ptychography

% Observables
Obs     = zeros(d,d,Nobs);
pro     = zeros(1,d);
pro(1:r)    = 1;
com     = zeros(1,d);
com(1)  = 1;

for r=1:numel(a)
    aFou    = frft_matrix(d,a(r));

    for s=1:d
        for l = 1:Nproj
            Pl  = diag(circshift(pro,[0,l-1]));      % ptychographic projector
            Com = circshift(com,[0,s-1]).'*circshift(com,[0,s-1]); %computational basis projector (intensity measurement)
            foo = (r-1)*Nproj*numel(a) + (s-1)*Nproj + l;
            Obs(:,:,foo)    = Pl*aFou*Com*aFou'*Pl;
        end
    end
end

w   = waitbar(0,'Progress:    0.0');
for q=1:Nstates
    waitbar(q/Nstates,w,['Progress:    ' num2str(q/Nstates)]);
    
    % Clearing SDPs variables
    clear F E Rho Delta
    
    % Picking a random state
    rho = cubitt_RandomDensityMatrix(d);

    % Cleaning yalmip memory
    yalmip('clear');
    F = class('double');

    % Defining the SDP variables
    Rho     = sdpvar(d,d,'hermitian','complex');
    Delta   = sdpvar(Nobs,1);

    % Standard constraints
%     F   = [Rho>=0,trace(Rho)==1]; % let's try normalizing by the end
    F   = [Rho>=0];
    F   = [F,Delta>=0];

    % Measured probabilities
    Prob        = zeros(Nobs,1);
    measured    = zeros(Nobs,1);
    
    
    % Applying poissonian noise, if desired
    if(noise)
        for m = 1:Nobs
            measured(m) = real(trace(rho*Obs(:,:,m)));
            measured(m)	= poissrnd(noise_av_counts*measured(m));
        end
%         foo	= sum(measured);
%         measured	= measured/foo;
    else
        for m =1:Nobs
            measured(m) = real(trace(rho*Obs(:,:,m)));
        end
    end
    
    for m = 1:Nobs
        Prob = trace(Rho*Obs(:,:,m));
        F = [F,Prob-measured(m)*(1-Delta(m))>=0];
        F = [F,Prob-measured(m)*(1+Delta(m))<=0];
    end

    % Cost function
    % here the cost function is the sum of the "flexibilizations" around the measured values
    E = sdpvar(1,1);
    F = [F,E>=0];
    E = sum(Delta);

    % Solving SDP
    SOLUTION = optimize(F,E,sdpsettings('solve','mosek','verbose',0));
%     SOLUTION = solvesdp(F,E,sdpsettings('solve','mosek'));
%     disp('DEBUGGING');
%     problema = double(SOLUTION.problem);
%     disp(yalmiperror(problema));

    Rho = double(Rho);
    trace(Rho)
    Rho = Rho/trace(Rho);
    Delta = double(Delta.');
    dist = hsDistance(rho,Rho);
    fid = Fidelity(rho,Rho);

    % Writing results to file
    fprintf(fid_res,'%d\t\t%1.12f\t\t%1.12f\n',q,dist,fid);
    
end
close(w)
%%


%% Closing files
if( fclose(fid_par)==-1 || fclose(fid_res)==-1 )
    error('Could not close files');
end