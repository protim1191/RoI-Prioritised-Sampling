 function locations = getNeighbourIndices2(m,n,steps,seed)
%get the linear indices of a neighbourhood. m is the number of rows in the
%image, steps is the radius of the neighbourhood, seed is the center pixel.
% Added boundary conditions (hopefully all)
% use this for further experiments
% created by Protim July 16 2018
vec = @(x) x(:);
if ~ismember(seed,1:m*n)
    error('seed value not inside image')
end
if mod(seed,m) == 0
    south = []; se = []; sw = [];
    east = m.*(1:steps);
    west = -m.*(1:steps);
    north = -1.*(1:steps);
    ne = vec(bsxfun(@plus, east,-1*(1:steps)'))';
    nw = vec(bsxfun(@plus, west,-1*(1:steps)'))';
    neighbour_offset = [east, west, north, south, ne, nw, se, sw];
    neighbours_index = bsxfun(@plus,seed,neighbour_offset');
    neighbours_index(neighbours_index < 0) = [];
    neighbours_index(neighbours_index > m*n) = [];
elseif mod(seed,m) == 1
    north = []; nw = []; ne = [];
    east = m.*(1:steps);
    west = -m.*(1:steps);
    south = 1:steps;
    se = vec(bsxfun(@plus, east,1*(1:steps)'))';
    sw = vec(bsxfun(@plus, west,1*(1:steps)'))';
    neighbour_offset = [east, west, north, south, ne, nw, se, sw];
    neighbours_index = bsxfun(@plus,seed,neighbour_offset');
    neighbours_index(neighbours_index < 0) = [];
    neighbours_index(neighbours_index > m*n) = [];
else
    east = m.*(1:steps);
    west = -m.*(1:steps);
    if mod(seed,m) -1 -steps < 0
        northRadius = steps - abs((mod(seed,m) - 1 -steps));
    else
        northRadius = steps;
    end
    if  m - mod(seed,m) <= steps
        southRadius = m - mod(seed,m) ;
    else
        southRadius = steps;
    end
    north = -1.*(1:northRadius);
    south = 1:southRadius;
    %     ne = vec(bsxfun(@plus, east,-1*(1:steps)'))';
    %     nw = vec(bsxfun(@plus, west,-1*(1:steps)'))';
    % se = vec(bsxfun(@plus, east,1*(1:steps)'))';
    %     sw = vec(bsxfun(@plus, west,1*(1:steps)'))';
    ne = vec(bsxfun(@plus, east,-1*(1:northRadius)'))';
    nw = vec(bsxfun(@plus, west,-1*(1:northRadius)'))';
    se = vec(bsxfun(@plus, east,1*(1:southRadius)'))';
    sw = vec(bsxfun(@plus, west,1*(1:southRadius)'))';
    neighbour_offset = [east, west, north, south, ne, nw, se, sw];
    neighbours_index = bsxfun(@plus,seed,neighbour_offset');
    neighbours_index(neighbours_index <= 0) = [];
    neighbours_index(neighbours_index > m*n) = [];
end

locations = [neighbours_index;seed];