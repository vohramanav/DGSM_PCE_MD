function g = model_oscillator(X)

m = X(:,1); % mass
k1  = X(:,2); % spring constant
k2 = X(:,3); % spring constant
r = X(:,4); % displacement
F = X(:,5); % Force
t1 = X(:,6); % time

w0 = ((k1+k2)./m).^0.5;
s = sin((w0.*t1)./2.0);
g = 3.0.*r - abs((2.0.*F.*s)./(m.*(w0.^2)));

end

