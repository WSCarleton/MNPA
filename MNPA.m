clear
clc

% Variables
global G C b;
nodes = 5; % number of nodes in circuit

% Initialize global variables
G = zeros(nodes, nodes);
C = zeros(nodes, nodes); 
b = zeros(nodes, 1);

% Variable sizes
r1 = 1; r2 = 2; r3 = 10; r4 = 0.1; ro = 1000; 
c = 0.25; 
l = 0.2; 

alpha = 100; 
vin = -10:0.01:10;
V3_values = zeros(size(vin));
Vo_values = zeros(size(vin));
% DC sweep
for z = 1:length(vin)
    % Setup G, C and b matrices
    G = zeros(nodes, nodes);
    C = zeros(nodes, nodes); 
    b = zeros(nodes, 1);
    
    % Initialize components in circuit
    resistor(1, 2, r1); 
    capacitor(1, 2, c);
    resistor(2, 0, r2);
    inductor(2, 3, l);
    resistor(3, 0, r3); 
    resistor(4, 5, r4); 
    resistor(5, 0, ro);
    voltage_controlled_voltage_source(4, 0, 3, 0, alpha/r3);
    voltage_source(1, 0, vin(z)); 
    
    % Find solution
    V = G\b; 
    Vo_values(z) = V(5); 
    V3_values(z) = V(3);
    
    if rem(z, 50) == 0
        fprintf('DC completion %3.3f\n', z*100/length(vin));
    end
end

figure(1)
clf
plot(vin, Vo_values, 'r', vin, V3_values, 'b');
title('Vin vs Vo and V3');
xlabel('Vin (V)'); ylabel('V0, V3 (V)');
legend('V0', 'V3');
grid on; 

% AC sweep of Vo(w)
f = linspace(1, 1e6, 1e7+1);

% Setup G, C and b matrices
G = zeros(nodes, nodes);
C = zeros(nodes, nodes); 
b = zeros(nodes, 1);
vin = 1;
        
% Initialize components in circuit
% V(1-5) are nodes 1-5
resistor(1, 2, r1); 
capacitor(1, 2, c);
resistor(2, 0, r2);
inductor(2, 3, l); % V(6)
resistor(3, 0, r3); 
resistor(4, 5, r4); 
resistor(5, 0, ro);
voltage_controlled_voltage_source(4, 0, 3, 0, alpha/r3); % V(7)
voltage_source(1, 0, vin); % V(8)

gain = zeros(size(f));
for z = 1:length(f)    
    % Find solution
    omega = 2*pi*f(z);
    V = (G+j*omega.*C)\b; 
    gain(z) = 20*log(V(5));
    
    if rem(z, 1000) == 0
        fprintf('AC sweep completion %3.2f\n', z*100/length(f));
    end
end

figure(2)
clf
semilogx(f, gain);
title('Gain(w) of V0/V1');
xlabel('Frequency'); ylabel('Gain (dB)');
grid on;

% AC sweep with C perturbation
% Setup G, C and b matrices
G = zeros(nodes, nodes);
C = zeros(nodes, nodes); 
b = zeros(nodes, 1);
vin = 1;
        
% Initialize components in circuit
% V(1-5) are nodes 1-5
resistor(1, 2, r1); 
capacitor(1, 2, c);
resistor(2, 0, r2);
inductor(2, 3, l); % V(6)
resistor(3, 0, r3); 
resistor(4, 5, r4); 
resistor(5, 0, ro);
voltage_controlled_voltage_source(4, 0, 3, 0, alpha/r3); % V(7)
voltage_source(1, 0, vin); % V(8)

r = 0.05*randn(1, 1000);
omega = pi;
gain = zeros(size(r));

for z = 1:length(r)    
    % Find solution
    V = (G + 1j*2*pi*omega .* r(z) .* C) \ b; 
    gain(z) = 20*log10(abs(V(5)));
    
    if rem(r, 1) == 0
        fprintf('AC sweep with perturbations completion %3.2f\n', z*100/length(f));
    end
end

figure(3)
clf
hist(gain);
title('Gain(w) of V0/V1 with perturbations');
xlabel('Gain (db)'); ylabel('Number of Data Points');
grid on;