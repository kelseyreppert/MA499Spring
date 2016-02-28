
function [x, z, xbar, H, B, scans] = dataGeneration(x0, y0, speed, angle, radius)
dim = 3;

% KF variables
q = 1e-4; % Process noise scaling
Q = q.*eye(dim);
r = 1e-3; % Observation noise scaling
R = r.*eye(dim);

dt = 1;
t0 = 0.0;
tf = 10.0;
scans = (tf-t0)/dt;

% Model Variables
H = eye(dim);
GAMMA = eye(dim);
w = q.*randn(dim,scans+1);
v = r.*randn(dim,scans+1);

x = zeros(dim,scans+1);
x(:,1) = [x0 y0 angle];
z = zeros(dim,scans+1);
z(:,1) = x(:,1) + w(:,1);
xbar = zeros(dim,scans+1);
xbar(:,1) = z(:,1);
xhat = zeros(dim,scans+1);
xhat(:,1) = z(:,1);
Pbar = zeros(dim,dim,scans+1);
Pbar(:,:,1) = eye(dim);
Phat = zeros(dim,dim,scans+1);
Phat(:,:,1) = eye(dim);
K = zeros(dim,dim,scans+1);
K(:,:,1) = eye(dim);
B = zeros(dim,dim,scans+1);
B(:,:,1) = eye(dim);



    function out = f(x,u,speed,radius,dt)
        out = x + dt*[speed*cos(x(3)) speed*sin(x(3)) speed*(u(3) - x(3))/radius]';
    end

    function out = h(x)
        out = x;
    end

for k = 1:scans
    
    u = [0 0 angle]';
    
    % Update KF Variables
    PHI = [1 0 -dt*speed*sin(x(3,k)); 0 1 dt*speed*cos(x(3,k)); 0 0 1 - dt*speed/radius];
    x(:,k+1) = f(x(:,k),u,speed,radius,dt) + w(:,k+1);
    z(:,k+1) = h(x(:,k+1)) + v(:,k+1);
    xbar(:,k+1) = PHI*xhat(:,k);
    Pbar(:,:,k+1) = PHI*Phat(:,:,k)*PHI' + GAMMA*Q*GAMMA';
    xhat(:,k+1) = xbar(k+1) + K(:,:,k+1)*(z(:,k+1) - H*xbar(:,k+1));
    Phat(:,:,k+1) = Pbar(:,:,k+1) - Pbar(:,:,k+1)*H'*inv(H*Pbar(:,:,k+1)*H' + R)*H*Pbar(:,:,k+1);
    K(:,:,k+1) = Phat(:,:,k+1)*H'*inv(R);
    B(:,:,k+1) = H*Pbar(:,:,k+1)*H' + R;
    
end

%plot(1:scans+1,x,'k-',1:scans+1,z,'r-',1:scans+1,xhat,'b-');
%legend('Process', 'Observation', 'Estimate')

end
