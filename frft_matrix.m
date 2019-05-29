function M = frft_matrix(d,a)
    foo = zeros(d,1);
    foo(1) = 1;
    M = zeros(d,d);
    
    for r = 1:d
        frft(foo,a);
        M(:,r) = frft(foo,a);
        foo = circshift(foo,[1,0]);
    end
end