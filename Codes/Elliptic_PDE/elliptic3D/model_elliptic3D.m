function Umax = model_elliptic(S)

Umax = zeros(size(S,1),1);
for i = 1:size(S,1)
    k = 0.1;c = S(i,1);a = S(i,2);b = S(i,3);
    Umax(i,1) = semilinear_elliptic_pde(k,c,a,b);
end
end
