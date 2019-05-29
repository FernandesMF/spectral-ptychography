function dist = hsDistance(A,B)
    % after zyckowski, p210
    dist = 0.5*trace( (A-B)*(A'-B') );
end