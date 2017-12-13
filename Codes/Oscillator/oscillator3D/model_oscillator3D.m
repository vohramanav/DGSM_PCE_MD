function g = model_oscillator3D(X)

m = 1.0; % mass
k1  = 1.0; % spring constant
k2 = 0.1; % spring constant
r = X(:,1); % displacement
F = X(:,2); % Force
t1 = X(:,3); % time

w0 = ((k1+k2)./m).^0.5;
s = sin((w0.*t1)./2.0);
g = 3.0.*r - abs((2.0.*F.*s)./(m.*(w0.^2)));

end
