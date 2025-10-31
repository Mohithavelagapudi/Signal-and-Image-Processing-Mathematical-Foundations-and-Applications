                  % This code is contributed by Miss Chatgpt and Mohitha......

% Generate sample data
%rng(1); % Set random seed for reproducibility
%X = [randn(100, 2) * 0.75 + ones(100, 1);randn(100, 2) * 0.5 - ones(100, 1)];
X=[185  72 ; 170 56 ;
    168 60 ; 179 68 ;
    182 72; 188 77]

% Set the number of clusters (k)
k = 2;

% Perform k-means clustering
[idx, centers] = kmeans(X, k); 
% returns the cluster indices idx indicating which cluster each point belongs to, as well as the cluster centers centers
variances = zeros(k, 1);
for i = 1:k
    cluster_points = X(idx == i);
    variances(i) = var(cluster_points)
end
% Plot the data points and cluster centers
figure;
gscatter(X(:,1), X(:,2), idx, 'rg', 'o');
hold on;
plot(centers(:,1), centers(:,2), 'k*', 'MarkerSize', 10);
hold off;
title('K-means Clustering');

