% This code is contributed by Miss Chatgpt and Mohitha......

% Define the continuous signal
%The Discrete Fourier Transform (DFT) is inherently designed for discrete signals, not continuous signals. 
% If you have a continuous signal, you typically need to sample it at discrete time points to apply the DFT.
t = 0:0.01:1; % Time vector
x = sin(2*pi*10*t) + sin(2*pi*20*t); % Example signal (sum of two sinusoids)

% Length of the signal
N = length(x);

% Initialize arrays to store the real and imaginary parts of the DFT
X_real = zeros(1, N);
X_imag = zeros(1, N);

% Compute the DFT
for k = 0:N-1
    for n = 0:N-1
        X_real(k+1) = X_real(k+1) + x(n+1) * cos(2*pi*k*n/N);
        X_imag(k+1) = X_imag(k+1) - x(n+1) * sin(2*pi*k*n/N);
    end
end

% Calculate the magnitude spectrum
X_mag = sqrt(X_real.^2 + X_imag.^2);

% Frequency vector
fs = 1 / (t(2) - t(1)); % Sampling frequency
f = (0:N-1) * fs / N;

% Display the results
figure;
subplot(2, 1, 1);
plot(t, x);
title('Input Signal');
xlabel('Time');
ylabel('Amplitude');
subplot(2, 1, 2);
stem(f, X_mag);
title('Magnitude Spectrum');
xlabel('Frequency');
ylabel('Magnitude');
