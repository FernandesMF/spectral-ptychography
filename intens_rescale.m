function [Inew] = intens_rescale(x,xnew,I)
    
    % Check if scales have same beginnng and end
    if( x(1)~=xnew(1) || x(end)~=xnew(end) )
        error('Scales should have the same beginning and end')
    end
    
    % Check if both scales are crescent
    if(~issorted(x))
        error('x scale is not sorted in ascending order')
    end
    if(~issorted(xnew))
        error('xnew scale is not sorted in ascending order')
    end

    % Caculating bin edges
    xedg = (x(1:end-1)+x(2:end))/2;
    xedg = [2*x(1)-xedg(1) xedg 2*x(end)-xedg(end)];
    nxedg = length(xedg);
    
    xnewedg = (xnew(1:end-1)+xnew(2:end))/2;
    xnewedg = [2*xnew(1)-xnewedg(1) xnewedg 2*xnew(end)-xnewedg(end)];
    
    % Recalculating intensities by proportional overlap of bins
    Inew = zeros(size(xnew));
    
    if(xnewedg(1)<xedg(1))
        Inew(1) = I(1)*(xedg(1)-xnewedg(1))/(xedg(2)-xedg(1));
    end
    
    for j=2:length(xnewedg)-1
        i = find( xnewedg(j-1)<xedg & xedg<xnewedg(j));
        if(~isempty(i))
            
            if(i(1)==1)
                i(1) = [];
            end
            
            Inew(j) = I(i(1)-1)*(xedg(i(1))-xnewedg(j-1))/(xedg(i(1))-xedg(i(1)-1));
            for s=2:length(i)
                Inew(j) = Inew(j) + I(i(s));
            end
            if(i(s)<nxedg)
                Inew(j) = Inew(j) + I(i(s)+1)*(xnewedg(j)-xedg(i(s)))/(xedg(i(s)+1)-xedg(i(s)));
            elseif(~isempty(i(s)))
                Inew(j) = Inew(j) + I(end)*(xnewedg(end)-xedg(end))/(xedg(end)-xedg(end-1));
            end
            
        else
            foo = find(xedg>xnewedg(j));
            Inew(j) = I(foo(1)-1)*(xnewedg(j)-xnewedg(j-1))/(xedg(foo(1))-xedg(foo(1)-1));
        end
    end
    
end

% % test:
% xnew = 1:10
% x = 1./(1:-0.1:0.1)
% t = rand(1,10)
% [tnew] = intens_rescale(x,xnew,t)
% plot(x,t,'o-',xnew,tnew,'x-')