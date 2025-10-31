 % This code is contributed by Miss Chatgpt and Mohitha......

% Read the input image
inputImage = imread('sp.jpg');
inputImage = rgb2gray(inputImage);

% Obtain the size of the input image
[M, N] = size(inputImage);

% Compute the 2D Fourier Transform
DFT = zeros(M, N);

for u = 0:M-1
    for v = 0:N-1
        for x = 0:M-1
            for y = 0:N-1
                DFT(u+1, v+1) = DFT(u+1, v+1) + double(inputImage(x+1, y+1)) * exp(-1i*2*pi*(u*x/M + v*y/N));
            end
        end
    end
end

% Apply image compression by setting lower-amplitude frequency components to zero
compression_ratio = 0.1; % Adjust the compression ratio as desired

% Compute the magnitudes of the frequency components
magnitudes = abs(DFT);

% Sort the magnitudes in descending order
[sorted_magnitudes, indices] = sort(magnitudes(:), 'descend');

% Determine the number of frequency components to keep
num_components = round(compression_ratio * numel(DFT));
%The numel function in MATLAB returns the total number of elements in an array

% Set lower-amplitude frequency components to zero
DFT(indices(num_components+1:end)) = 0;

% Compute the inverse Fourier Transform
compressedImage = zeros(M, N);

for x = 0:M-1
    for y = 0:N-1
        for u = 0:M-1
            for v = 0:N-1
                compressedImage(x+1, y+1) = compressedImage(x+1, y+1) + DFT(u+1, v+1) * exp(1i*2*pi*(u*x/M + v*y/N));
            end
        end
        compressedImage(x+1, y+1) = compressedImage(x+1, y+1) / (M*N);
    end
end

% Convert the compressed image back to uint8 format
compressedImage = uint8(real(compressedImage));
%To convert the complex-valued compressedImage back to a real-valued image suitable for display
% and storage, the real function is used. The real function extracts the real part of each complex 
% number, discarding the imaginary part.

% Display the original and compressed images
figure;
subplot(1, 2, 1);
imshow(inputImage);
title('Original Image');
subplot(1, 2, 2);
imshow(compressedImage);
title('Compressed Image');

