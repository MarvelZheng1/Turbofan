% matrix and grid
A = [6 -5; 4 2];
[x1,x2] = meshgrid(linspace(-6,6,21));

% vector field F = A*[x1;x2]
U = A(1,1).*x1 + A(1,2).*x2;   % dx1/dt
V = A(2,1).*x1 + A(2,2).*x2;   % dx2/dt

% (optional) normalize arrows so they’re comparable in length
L = sqrt(U.^2 + V.^2);  U = U./L;  V = V./L;

% plot the direction field
quiver(x1, x2, U, V, 'AutoScale', 'off');
axis equal; grid on
xlabel('x_1'); ylabel('x_2'); title('Direction field for \dot{x}=Ax'); hold on

% (optional) add eigenvector lines and a sample trajectory from x(0)=[1;1]
[Vv,~] = eig(A);
for k = 1:2
    v = Vv(:,k)/norm(Vv(:,k));
    plot([-6 6]*v(1), [-6 6]*v(2), 'k-');    % eigen-directions
end

% particular solution for x(0)=[1;1]
C1 = 2/7; C2 = -1/7;
t = linspace(-0.6, 0.6, 400);
x1_sol = 4*C1*exp(6*t) + C2*exp(-t);
x2_sol = 3*C1*exp(6*t) - C2*exp(-t);
plot(x1_sol, x2_sol, 'k', 'LineWidth', 1.5);

