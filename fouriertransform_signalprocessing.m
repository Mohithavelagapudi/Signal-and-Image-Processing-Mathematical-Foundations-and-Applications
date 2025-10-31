                 % This code is contributed by Miss Chatgpt and Mohitha......

% Define the input signal
t = linspace(0, 1, 1000); % Time vector
x = sin(2*pi*10*t) + sin(2*pi*20*t); % Example signal (sum of two sinusoids)

% Length of the input signal
N = length(x);

% Compute the Fourier Transform
X_real = zeros(1, N);
X_imag = zeros(1, N);

for k = 1:N
    for n = 1:N
        X_real(k) = X_real(k) + x(n) * cos(2*pi*(k-1)*(n-1)/N);
        X_imag(k) = X_imag(k) - x(n) * sin(2*pi*(k-1)*(n-1)/N);
    end
end

% Frequency vector
fs = 1 / (t(2) - t(1)); % Sampling frequency
f = (-fs/2:fs/N:fs/2-fs/N); % Frequency range

% Apply signal processing (e.g., filtering)
% Example: Apply a low-pass filter by setting high-frequency components to zero
cutoff_frequency = 15; % Adjust the cutoff frequency as desired

for k = 1:N
    if abs(f(k)) > cutoff_frequency
        X_real(k) = 0;
        X_imag(k) = 0;
    end
end

% Compute the inverse Fourier Transform
y = zeros(1, N); % y = zeros(1, N); initializes an array y of zeros with a length of N. 
%  This array will store the samples of the filtered signal in the time domain.

for n = 1:N % iterates over each time sample 
    for k = 1:N %iterates over each frequency sample 
        y(n) = y(n) + (X_real(k) * cos(2*pi*(k-1)*(n-1)/N) - X_imag(k) * sin(2*pi*(k-1)*(n-1)/N));
    end
    y(n) = y(n) / N;
end

% Plotting
figure;
subplot(2, 1, 1);
plot(t, x);
title('Input Signal');
xlabel('Time');
ylabel('Amplitude');
subplot(2, 1, 2);
plot(t, y);
title('Filtered Signal');
xlabel('Time');
ylabel('Amplitude');
