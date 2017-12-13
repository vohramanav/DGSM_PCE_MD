function f = red_mor(X)
    bi = [0.05;0.59;10;0.21];
    bij = [0 80 60 40;0 30 0.73 0.18;0 0 0.64 0.93;0 0 0 0.06];
    bij4 = [0 10 0.98 0.19;0 0 0.49 50;0 0 0 1;0 0 0 0];
    f = zeros(size(X,1),1);
    for p = 1:size(X,1)
    Xp = transpose(X(p,:));
    t1 = transpose(bi)*Xp;
    t2 = transpose(transpose(bij)*Xp)*Xp;
    t3 = (transpose(transpose(bij4)*Xp)*Xp)*Xp(4,1);
    f(p,1) = t1+t2+t3;
    end
end
