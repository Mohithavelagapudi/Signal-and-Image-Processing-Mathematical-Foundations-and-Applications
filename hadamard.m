                                        % This code is contributed by Miss Chatgpt and Mohitha......

N = 4;      % length of the signal
A = [1 1;1 -1]; 
B = kron(A,A); 
x=[1;2;3;4];
B1=(1/sqrt(N))*B;
X=B1*x

