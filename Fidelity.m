function fid = Fidelity(p1,p2)
fid = ( trace(sqrtm( sqrtm(p1)*p2*sqrtm(p1) )) )^2;
