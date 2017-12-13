function Y = uq_borehole(X)

rw = X(:, 1);
r  = 3698.30;
Tu = X(:, 2);
Hu = X(:, 3);
Tl = X(:, 4);
Hl = X(:, 5);
L  = X(:, 6);
Kw = X(:, 7);

% Precalculate the logarithm:
Logrrw = log(r./rw);

Numerator = 2*pi*Tu.*(Hu - Hl);
Denominator = Logrrw.*(1 + (2*L.*Tu)./(Logrrw.*rw.^2.*Kw) + Tu./Tl);

Y = Numerator./Denominator;
