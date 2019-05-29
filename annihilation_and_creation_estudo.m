N = 4; n = 1:N;
a = diag(sqrt(n),1); ad = a';
n2 = [0 0 1 0 0]';
% a*n2;
% ad*n2;
% ad*a*n2;
% dbase = length(n)+1;
% 
% base_c = cell(dbase, 1); % gera base canonica
% 
% for i=1:dbase
%     
%     base_c {i} = zeros (dbase, 1);
%     base_c {i}(i,1) = 1;
% 
% end

nt = [0 0 0 1 0]';

ad*nt