function base_c = basecanon(dbase)

base_c = cell(dbase, 1); % gera base canonica

for i=1:dbase
    
    base_c {i} = zeros (dbase, 1);
    base_c {i}(i,1) = 1;

end
