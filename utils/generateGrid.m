function grid = generateGrid(row)
% for spatial frequency.
for ii = 1:row
            if mod(ii,2) ~= 0
                grid(ii,:) = (ii-1)*row + 1:ii*row;
            else
                grid(ii,:) = ii*row:-1:(ii-1)*row+1;
           end
end