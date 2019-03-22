function voltage_controlled_current_source(node_1, node_2, ni_1, ni_2, value)

global G;

% check if node 1 and reference node 1 are not 0 but node 2 and reference
% node 2 are 0
if (ni_1 ~= 0) && (node_1 ~= 0)
    G(ni_1,node_1) = G(ni_1,node_1) + value;
end

% check if node 2 and reference node 2 are not 0 but node 1 and reference
% node 1 are 0
if (ni_2 ~= 0) && (node_2 ~= 0)
    G(ni_2,node_2) = G(ni_2,node_2) + value;
end

% check if reference node 1 and node 2 are not 0 but node 1 and reference
% node 2 are 0
if (ni_1 ~= 0) && (node_2 ~= 0)
    G(ni_1,node_2) = G(ni_1,node_2) - value;
end

% check if reference node 2 and node 1 are not 0 but node 2 and reference
% node 1 are 0
if (ni_2 ~= 0) && (node_1 ~= 0)
    G(ni_2,node_1) = G(ni_2,node_1) - value;
end
end


