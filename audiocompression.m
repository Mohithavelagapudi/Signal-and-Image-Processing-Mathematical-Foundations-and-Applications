                                       % This code is contributed by Miss Chatgpt and Mohitha......

% Read the audio file
[input, fs] = audioread('input.wav'); % fs=sampled frequency

% Length of the input audio
N = length(input);

% Compute the Fourier Transform
X_real = zeros(1, N);
X_imag = zeros(1, N);

for k = 1:N
    for n = 1:N
        X_real(k) = X_real(k) + input(n) * cos(2*pi*(k-1)*(n-1)/N);
        X_imag(k) = X_imag(k) - input(n) * sin(2*pi*(k-1)*(n-1)/N);
    end
end

% Apply audio compression by setting lower-amplitude frequency components to zero
compression_ratio = 0.5; % Adjust the compression ratio as desired
%The compression ratio determines the amount of compression to be applied

% Sort the frequency components based on amplitude
amplitudes = sqrt(X_real.^2 + X_imag.^2);
[sorted_amplitudes, indices] = sort(amplitudes, 'descend');
%The resulting sorted amplitudes are stored in the variable sorted_amplitudes, and the 
% corresponding indices are stored in the variable indices.

% Determine the number of frequency components to keep
num_components = round(compression_ratio * N);

% Set lower-amplitude frequency components starting from index
% num_components+1 up to the end, are set to zero
X_real(indices(num_components+1:end)) = 0;
X_imag(indices(num_components+1:end)) = 0;
%The sorting step allows us to identify the most significant frequency components with higher
% amplitudes, while the subsequent zeroing step removes the less significant components with lower amplitudes. 

% Compute the inverse Fourier Transform
output = zeros(1, N);

for n = 1:N
    for k = 1:N
        output(n) = output(n) + (X_real(k) * cos(2*pi*(k-1)*(n-1)/N) - X_imag(k) * sin(2*pi*(k-1)*(n-1)/N));
    end
    output(n) = output(n) / N;
end

% Write the compressed audio to a file
audiowrite('compressed.wav', output, fs);
