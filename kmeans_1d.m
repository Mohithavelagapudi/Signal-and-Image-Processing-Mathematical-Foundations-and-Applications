                          % This code is contributed by Miss Chatgpt and Mohitha......

% Generate sample data
rng(1); % Set random seed for reproducibility
X = randn(200, 1);

% Set the number of clusters (k)
k = 3;

% Perform k-means clustering
[idx, centers] = kmeans(X, k);

% Plot the data points and cluster centers
figure;
scatter(X, zeros(size(X)), 50, idx, 'filled');
hold on;
scatter(centers, zeros(size(centers)), 200, 'k', 'filled');
hold off;
title('K-means Clustering for 1D Data');
