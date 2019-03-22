function capacitor(node_1,node_2,value)

global C;

% Check if node 1 isnt 0 but node 2 is 0
if (node_1 ~= 0)
    C(node_1, node_1) = C(node_1, node_1) + value;
end

% Check if node 2 isnt 0 but node 2 is 0
if (node_2 ~= 0)
    C(node_2, node_2) = C(node_2, node_2) + value;
end

% Check if both nodes arent 0
if (node_1 ~= 0) && (node_2 ~= 0)
    C(node_1, node_2) = C(node_1, node_2) - value;
    C(node_2, node_1) = C(node_2, node_1) - value;
end
end

