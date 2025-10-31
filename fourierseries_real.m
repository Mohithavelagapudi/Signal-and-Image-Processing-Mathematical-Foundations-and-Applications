% This code is contributed by Miss Chatgpt and Mohitha......

% Signal parameters
T = 2*pi;                   % Period of the signal
w0 = 2*pi/T;                % Fundamental frequency

% Define the signal function
f = @(t) sin(2*t) + 0.5*sin(4*t);   % Example signal: sum of a fundamental and second harmonic

% Fourier series coefficients
N = 10;                     % Number of coefficients to calculate
a0 = (1/T) * integral(f, 0, T);      % DC coefficient

an = zeros(1, N);           % Initialize array for cosine coefficients
bn = zeros(1, N);           % Initialize array for sine coefficients

for n = 1:N
    % Calculate the Fourier series coefficients
    an(n) = (2/T) * integral(@(t) f(t) .* cos(n*w0*t), 0, T);
    bn(n) = (2/T) * integral(@(t) f(t) .* sin(n*w0*t), 0, T);
end

% Fourier series approximation
t = linspace(0, T, 1000);   % Time vector for approximation
approximation = a0/2;       % Initialize approximation with DC component

for n = 1:N
    approximation = approximation + an(n)*cos(n*w0*t) + bn(n)*sin(n*w0*t);
end
%print fourier series coefficients
disp('Fourier series coefficients are : ')
disp(['a0 = ', num2str(a0)]);
for n = 1:N
    disp(['a', num2str(n), ' = ', num2str(an(n))]);
    disp(['b', num2str(n), ' = ', num2str(bn(n))]);
end

% Plot the original signal and the Fourier series approximation
figure;
plot(t, f(t), 'b', 'LineWidth', 1.5);
hold on;
plot(t, approximation, 'r--', 'LineWidth', 1.5);
xlabel('Time');
ylabel('Signal');
