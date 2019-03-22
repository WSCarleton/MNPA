function inductor(node_1, node_2, value)

global G b C;

% Setup new row
n = size(C,1); 
extra = n+1; 

% Expand b, G, and C
b(extra) = 0; 
G(extra, extra) = 0; 
C(extra ,extra) = 0; 

% Check if node 1 isnt 0 but node 2 is 0
if (node_1 ~= 0)
    G(node_1, extra) = G(node_1, extra) + 1;
    G(extra, node_1) = G(extra, node_1) + 1;
end

% Check if node 2 isnt 0 but node 1 is 0
if (node_2 ~= 0)
    G(node_2, extra) = G(node_2, extra)-1;
    G(extra, node_2) = G(extra, node_2)-1;
end

% Fill in C matrix
C(extra, extra) = C(extra, extra) -value;
end