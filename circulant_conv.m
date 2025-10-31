                                          % This code is contributed by Miss Chatgpt and Mohitha......

c=[1 2 3]; % first row of matrix c
C= toeplitz([c(1) fliplr(c(2:end))],c);% circulant matrix
% Compute the eigenvalues and eigenvectors
[V, E] = eig(C);

d=[5 0 4]; % first row of matrix d
D=toeplitz([d(1) fliplr(d(2:end))],d);% circulant matrix
% prove that CD=DC
mult=C*D;
mult1=D*C;
% proving that first row of CD is c convolution d
cd = cconv(c,d,3)

% Extract the eigenvalues and eigenvectors
eigenvalues = diag(E);
eigenvectors = V;

% Display the eigenvalues
disp('Eigenvalues:');
disp(eigenvalues);

% Display the eigenvectors
disp('Eigenvectors:');
disp(eigenvectors);

