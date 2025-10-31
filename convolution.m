                % This code is contributed by Miss Chatgpt and Mohitha......



% 1.Define the two functions f(t) and g(t) that you want to convolve.
% 2.Create a time vector t that represents the range of time values at which you want to compute the convolution.
% 3.Compute the Fourier Transform F(w) of f(t) and G(w) of g(t) using the Fourier transform formulas.
% 4.Multiply the Fourier Transforms F(w) and G(w) element-wise to obtain the Fourier Transform of the convolution result.
% 5.Compute the inverse Fourier Transform of the product to obtain the convolution result in the time domain.
% 6.Plot the original functions f(t) and g(t), as well as the convolved function.



% Define the functions to convolve
f = @(t) exp(-t.^2/2); % Gaussian function
g = @(t) sin(2*pi*t); % Sine function

% Define the time vector
T = 10; % Total time duration
N = 1000; % Number of time points
dt = T / N; % Time step
t = linspace(0, T-dt, N);

% Compute the Fourier Transform of f(t)
F = zeros(size(t));
for k = 1:N
    F(k) = sum(f(t) .* exp(-1i*2*pi*(k-1)*t/T)) * dt;
end

% Compute the Fourier Transform of g(t)
G = zeros(size(t));
for k = 1:N
    G(k) = sum(g(t) .* exp(-1i*2*pi*(k-1)*t/T)) * dt;
end

% Multiply the Fourier Transforms element-wise
H = F .* G;

% Compute the inverse Fourier Transform of the product
h = zeros(size(t));
for k = 1:N
    h = h + H(k) * exp(1i*2*pi*(k-1)*t/T);
end
h = real(h) / T;

% Plot the original functions and the convolved function
figure;
subplot(3, 1, 1);
plot(t, f(t));
xlabel('t');
ylabel('f(t)');
title('Function f(t)');
subplot(3, 1, 2);
plot(t, g(t));
xlabel('t');
ylabel('g(t)');
title('Function g(t)');
subplot(3, 1, 3);
plot(t, h);
xlabel('t');
ylabel('h(t)');
title('Convolution Result');
