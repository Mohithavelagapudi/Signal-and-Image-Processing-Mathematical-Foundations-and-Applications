                                   % This code is contributed by Miss Chatgpt and Mohitha......


% 1. Define the function f(x) that you want to differentiate.
% 2.Create a vector x that represents the spatial domain values at which you want to compute the derivative.
% 3.Compute the Fourier coefficients F(k) of the function f(x) using the Fourier series formula.
% 4.Compute the Fourier coefficients G(k) of the derivative g(x) using the relationship G(k) = (i * k * F(k)).
%5. Compute the inverse Fourier coefficients g(x) using the inverse Fourier series formula.
%6. Plot the function f(x) and its derivative g(x) to visualize the results.

% Define the function to differentiate
f = @(x) sin(2*pi*x) + cos(4*pi*x);

% Define the spatial domain
L = 2; % Length of the spatial domain
N = 1000; % Number of spatial points
dx = L / N; % Spatial step size
x = linspace(0, L-dx, N);

% Compute the Fourier coefficients of f(x)
F = zeros(size(x));
for k = 1:N
    F(k) = sum(f(x) .* exp(-1i*2*pi*(k-1)*x/L)) * dx;
end

% Compute the Fourier coefficients of the derivative g(x)
G = 1i * (1:N) .* F;

% Compute the inverse Fourier coefficients of g(x)
g = zeros(size(x));
for k = 1:N
    g = g + G(k) * exp(1i*2*pi*(k-1)*x/L);
end
g = real(g) / L;

% Plot the function and its derivative
figure;
subplot(2, 1, 1);
plot(x, f(x));
xlabel('x');
ylabel('f(x)');
title('Function');
subplot(2, 1, 2);
plot(x, g);
xlabel('x');
ylabel('g(x)');
title('Derivative');
