                                     % This code is contributed by Miss Chatgpt and Mohitha......

% Generate sample data
% rng=random number generator

rng(1); % Set random seed for reproducibility
X = [randn(100, 2) * 0.75 + ones(100, 1);
    randn(100, 2) * 0.5 - ones(100, 1)];

% Set the number of clusters (k)
k = 2;

% Initialize cluster centroids randomly
rng(2); % Set random seed for centroid initialization
initial_indices = randperm(size(X, 1), k);% random permutation of integers
centroids = X(initial_indices, :);

% Perform k-means clustering
max_iters = 10; % Maximum number of iterations
for iter = 1:max_iters
    % Assign each data point to the nearest centroid
    distances = pdist2(X, centroids);
    [~, cluster_indices] = min(distances, [], 2);
    
    % Update the centroids based on the assigned data points
    for i = 1:k
        cluster_points = X(cluster_indices == i, :);
        centroids(i, :) = mean(cluster_points, 1);
    end
end


% Plot the data points and final cluster centroids
figure;
gscatter(X(:,1), X(:,2), cluster_indices, 'rg', 'o');
hold on;
plot(centroids(:,1), centroids(:,2), 'k*', 'MarkerSize', 10);
hold off;
title('K-means Clustering (Implemented from Scratch)');
