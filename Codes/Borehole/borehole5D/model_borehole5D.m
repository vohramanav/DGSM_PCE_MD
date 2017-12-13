function Y = model_borehole5D(X)

rw = X(:, 1);
r  = 3698.30;
Tu = 89335;
Hu = X(:, 2);
Tl = 89.55;
Hl = X(:, 3);
L  = X(:, 4);
Kw = X(:, 5);

% Precalculate the logarithm:
Logrrw = log(r./rw);

Numerator = 2*pi*Tu.*(Hu - Hl);
Denominator = Logrrw.*(1 + (2*L.*Tu)./(Logrrw.*rw.^2.*Kw) + Tu./Tl);

Y = Numerator./Denominator;
