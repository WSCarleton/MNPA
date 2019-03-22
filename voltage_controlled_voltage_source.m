function voltage_controlled_voltage_source(node_1, node_2, ni_1, ni_2, value)

global G C b;

% add new row to G and b
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

% check if reference node 1 isnt 0 but reference node 2 is 0
if (ni_1 ~= 0)
    G(extra, ni_1) = G(extra, ni_1)-value;
end

% check if reference node 2 isnt 0 but reference node 1 is 0
if (ni_2 ~= 0)
    G(extra, ni_2) = G(extra, ni_2)+value;
end
end
