% Estimation of a constant
% NOTE: k=1 is the initial scan
% Parameters
constant = 3;
tf = 2;
dim = 1;

% KF variables
H = eye(dim);
q = 1e-4; % Process noise scaling
Q = q.*eye(dim);
r = 1e-3; % Observation noise scaling
R = r.*eye(dim);

% Model Variables
PHI = eye(dim);
GAMMA = eye(dim);
w = q.*randn(dim,tf+1);
v = r.*randn(dim,tf+1);

x = zeros(dim,tf+1);
x(1) = constant;
z = zeros(dim,tf+1);
z(1) = constant + w(1);
xbar = zeros(dim,tf+1);
xbar(1) = z(1);
xhat = zeros(dim,tf+1);
xhat(1) = z(1);
Pbar = zeros(dim,dim,tf+1);
Pbar(:,:,1) = eye(dim);
Phat = zeros(dim,dim,tf+1);
Phat(:,:,1) = eye(dim);
K = zeros(dim,dim,tf+1);
K(:,:,1) = eye(dim);
B = zeros(dim,dim,tf+1);
B(:,:,1) = eye(dim);

% MHT Parameters
betaFT = 0.1;
betaNT = 0.5;
PD = 0.9; % Probability of detection

% MHT Variables
m = ones(1,tf+1); % The number of measurements at each scan k
%m(1) = 1;
Ntgt = ones(1,tf+1); % The number of detected targets at the beginning of each scan k
Ntgt(1) = 1;
Nhyp = ones(1,tf+1); % The number of hypothesis at the beginning of each scan k
Nhyp(1) = 1;

OMEGA = generateHyp(m(1),Ntgt(1)); % The set of all hypotheseses
P = ones(size(OMEGA)); % Probability initialization
P(1) = 0;
P(end) = 0;

CT = linspace(1,Ntgt(1), Ntgt(1)); % Confirmed targets
TT = linspace(Ntgt(1)+1, Ntgt(1)+1+sum(m), sum(m)+1); % Tentative targets

for k = 1:tf
    
%     [highestProb, index] = max(P(:,k)); % Changes tentative targets to confirmed targets
%     if sum(OMEGA(index,k) == TT) > 0
%        CT = [CT OMEGA(index,k)]; 
%     end
    
    % Update KF Variables
    x(k+1) = PHI*x(k) + GAMMA*w(k+1);
    z(k+1) = H*x(k+1) + v(k+1);
    xbar(k+1) = PHI*xhat(k);
    Pbar(k+1) = PHI*Phat(:,:,k)*PHI' + GAMMA*Q*GAMMA';
    xhat(k+1) = xbar(k+1) + K(:,:,k+1)*(z(k+1) - H*xbar(k+1));
    Phat(k+1) = Pbar(:,:,k+1) - Pbar(:,:,k+1)*H'*inv(H*Pbar(:,:,k+1)*H' + R)*H*Pbar(:,:,k+1);
    K(:,:,k+1) = Phat(:,:,k+1)*H'*inv(R);
    B = H*Pbar*H' + R;
    
    [OMEGA, P] = expandHyp(OMEGA, P, m(k+1), CT);
    Nft = sum(OMEGA(:,k+1)==0, 2); % Number of false targets for each hypothesis OMEGA(i,:)
    Nnt = sum(OMEGA(:,k+1)>m(k+1), 2); % Number of new targets for each hypothesis OMEGA(i,:)
    Ndt = m(k) - Nnt - Nft; % Number of detections for each hypothesis OMEGA(i,:)
    
    Nhyp = numel(OMEGA(:,k+1)); % Nhyp is the number of hypothesis at scan k
    for i = 1:Nhyp
        product = 1;
        
        DT = Ndt(i);
        for j = 1:DT
            product = product*normpdf(z(k+1) - H*xbar(k+1),B(:,:,k));
        end
        P(i,k+1) = PD^Ndt(i)*(1-PD)^(Ntgt(k) - Ndt(i))*betaFT^Nft(i)*betaNT^Nnt(i)*product*P(i,k); % (16) Reid
    end
    P(i,k+1) = P(i,k+1)./sum(P(:,k+1),1); % Normalizing the probabilities
    
    
    
    % Ntgt(k+1) = Ndt(P(:,k+1) == max(P(:,k+1))); 
end

%plot(1:tf+1,x,'k-',1:tf+1,z,'r-',1:tf+1,xhat,'b-');
%legend('Process', 'Observation', 'Estimate')
