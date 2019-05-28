function z = zig_zag_v(x)
    ind = reshape(1:numel(x), size(x)); %# indices of elements
    ind = fliplr( spdiags( fliplr(ind) ) );     %# get the anti-diagonals
    ind(:,1:2:end) = flipud( ind(:,1:2:end) );  %# reverse order of odd columns
    ind(ind==0) = [];                           %# keep non-zero indices
    z = x(ind); 
end

