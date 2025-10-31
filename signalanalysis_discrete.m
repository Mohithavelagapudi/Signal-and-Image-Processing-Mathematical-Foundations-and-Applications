% This code is contributed by Miss Chatgpt and Mohitha......

% Define the input signal
x = [1 2 3 4 5 6 7 8];

% Length of the input signal
N = length(x);

% Initialize arrays to store the real and imaginary parts of the DFT
X_real = zeros(1, N);
X_imag = zeros(1, N);
%The code iterates over each frequency bin k and calculates the real and imaginary parts of the DFT 
% coefficients by summing the product of the input signal samples x[n] and the complex exponential term e^(-j2πkn/N) for all time samples n.

% Compute the DFT
for k = 0:N-1
    for n = 0:N-1
        X_real(k+1) = X_real(k+1) + x(n+1) * cos(2*pi*k*n/N);
        X_imag(k+1) = X_imag(k+1) - x(n+1) * sin(2*pi*k*n/N);
    end
end

% Calculate the magnitude spectrum
X_mag = sqrt(X_real.^2 + X_imag.^2);

% Display the results
disp('DFT Coefficients:');
disp('k     Real       Imaginary      Magnitude');
for k = 0:N-1
    disp([num2str(k) '     ' num2str(X_real(k+1)) '      ' num2str(X_imag(k+1)) '      ' num2str(X_mag(k+1))]);
end
