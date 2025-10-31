% This code is contributed by Miss Chatgpt and Mohitha......

% Define the input data
t = linspace(0, 1, 1000); % Time vector
% we define the input data as a sinusoid corrupted by Gaussian noise.
x = sin(2*pi*10*t) + 0.5*randn(size(t)); % Example data (sinusoid + Gaussian noise)

% Length of the input data
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

% Apply denoising by setting low-amplitude frequency components to zero
threshold = 0.1; % Adjust the threshold as desired
% Lower values will preserve more frequency components, including lower-amplitude ones, while higher 
% values will remove more frequency components, considering them as noise.

for k = 1:N
    amplitude = sqrt(X_real(k)^2 + X_imag(k)^2);
    % This loop iterates over each frequency bin k and computes the amplitude of the complex Fourier coefficient at that frequency. 
    if amplitude < threshold 
        %If the computed amplitude is less than the defined threshold, it implies that the 
        % corresponding frequency component is of low amplitude and likely represents noise
        X_real(k) = 0;
        X_imag(k) = 0;
    end
end

% Compute the inverse Fourier Transform
y = zeros(1, N);

for n = 1:N
    for k = 1:N
        y(n) = y(n) + (X_real(k) * cos(2*pi*(k-1)*(n-1)/N) - X_imag(k) * sin(2*pi*(k-1)*(n-1)/N));
    end
    y(n) = y(n) / N;
end

% Plotting
figure;
subplot(2, 1, 1);
plot(t, x);
title('Input Data');
xlabel('Time');
ylabel('Amplitude');
subplot(2, 1, 2);
plot(t, y);
title('Denoised Data');
xlabel('Time');
ylabel('Amplitude');
