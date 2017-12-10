
A = [.5 -5.5 -2.5 9 1 0 0 0; .5 -1.5 -.5 1 0 1 0 0; 1 0 0 0 0 0 1 1;-10 57 9 24 0 0 0 0];

A
B = pivot_step(A,1,1)
C = pivot_step(B,2,2)
D = pivot_step(C,1,3)
E = pivot_step(D,2,4)
F = pivot_step(E,1,5)
G = pivot_step(F,2,1)
H = pivot_step(G,3,3)
