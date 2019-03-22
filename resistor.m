function resistor(node_1, node_2, value)

global G;

% Setup value of resistor
value = 1/value;

% Check if node 1 isnt 0 but node 2 is 0
if (node_1 ~= 0)
    G(node_1,node_1) = G(node_1,node_1) + value;
end

% Check if node 2 isnt 0 but node 1 is 0
if (node_2 ~= 0)
    G(node_2,node_2) = G(node_2,node_2) + value;
end

% Check if both arent 0
if (node_1 ~= 0) && (node_2 ~= 0)
    G(node_1,node_2) = G(node_1,node_2) - value;
    G(node_2,node_1) = G(node_2,node_1) - value;
end
%END

