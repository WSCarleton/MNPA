function current_source(node_1, node_2, value)

global b;

% Check if node 1 isnt 0 but node 2 is 0
if (node_1 ~= 0)
    b(node_1) = b(node_1) + value;
end

% Check if node 2 isnt 0 but node 1 is 0
if (node_2 ~= 0)
    b(node_2) = b(node_2) - value;
end 
end


