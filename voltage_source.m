function voltage_source(node_1, node_2, value)

global G b C;

% make new row
n = size(G,1); 
extra = n+1; 

b(extra) = 0;  
G(extra, extra) = 0; 
C(extra, extra) = 0; 

% check if node 1 isnt 0 but node 2 is 0
if (node_1 ~= 0)
    G(node_1, extra) = 1;
    G(extra, node_1) = 1;
end

% check if node 2 isnt 0 but node 1 is 0
if (node_2 ~= 0)
    G(node_2, extra) = -1;
    G(extra, node_2) = -1;
end

% assign solution to b matrix
b(extra) = value;
end
