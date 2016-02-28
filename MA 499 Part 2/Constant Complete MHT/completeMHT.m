% Estimation of a constant
% NOTE: k=1 is the initial scan

% Parameters
constant = 3;
scans = 3;
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
w = q.*randn(dim,scans+1);
v = r.*randn(dim,scans+1);

x = zeros(dim,scans+1);
x(1) = constant;
z = zeros(dim,scans+1);
z(1) = constant + w(1);
xbar = zeros(dim,scans+1);
xbar(1) = z(1);
xhat = zeros(dim,scans+1);
xhat(1) = z(1);
Pbar = zeros(dim,dim,scans+1);
Pbar(:,:,1) = eye(dim);
Phat = zeros(dim,dim,scans+1);
Phat(:,:,1) = eye(dim);
K = zeros(dim,dim,scans+1);
K(:,:,1) = eye(dim);
B = zeros(dim,dim,scans+1);
B(:,:,1) = eye(dim);




% MHT Parameters
betaFT = 0.1;
betaNT = 0.5;
PD = 0.9; % Probability of detection

% MHT Variables
m = ones(1,scans+1); % The number of measurements at each scan k
% m(1) = 1; % Number of initial measurements

iniT = 1; % Number of initial targets
CT = linspace(1, iniT, iniT); % Confirmed targets
TT = linspace(iniT+1, iniT+m(1), m(1)); % Tentative targets

OMEGA = generateHyp(m(1), CT, TT); % The set of all hypotheseses
CT = [CT TT]; % Adds the tentative targets to the list of confirmed targets

% Probability initialization 
P = zeros(size(OMEGA)); 
%P(9,m(1)) = 1;
P(2,1) = 1;
% disp(OMEGA);
% row = input('Which row of OMEGA is has probability 1? ');
% P(row,m(1)) = 1;

Ntgt = zeros(scans);
Ntgt(1) = numel(CT);


for k = 1:scans
    
    % Update KF Variables
    x(k+1) = PHI*x(k) + GAMMA*w(k+1);
    z(k+1) = H*x(k+1) + v(k+1);
    xbar(k+1) = PHI*xhat(k);
    Pbar(k+1) = PHI*Phat(:,:,k)*PHI' + GAMMA*Q*GAMMA';
    xhat(k+1) = xbar(k+1) + K(:,:,k+1)*(z(k+1) - H*xbar(k+1));
    Phat(k+1) = Pbar(:,:,k+1) - Pbar(:,:,k+1)*H'*inv(H*Pbar(:,:,k+1)*H' + R)*H*Pbar(:,:,k+1);
    K(:,:,k+1) = Phat(:,:,k+1)*H'*inv(R);
    B = H*Pbar*H' + R;
    
    indA = sum(m(1:k))+1;
    indB = sum(m(1:k+1));
    
    [OMEGA, P, newTT] = expandHyp(OMEGA, P, m(k+1), CT);
    Nft = sum(OMEGA(:,indA:indB)==0, 2); % Number of false targets for each hypothesis OMEGA(i,:)
    % Nnt = sum(OMEGA(:,indA:indB)>m(k+1), 2); % Number of new targets for each hypothesis OMEGA(i,:)
    
    Nnt = zeros(size(Nft));
    for tt = 1:numel(TT)
        Nnt = Nnt + sum(OMEGA(:,indA:indB)==TT(tt),2);
    end
    
    Ndt = m(k) - Nnt - Nft; % Number of detections for each hypothesis OMEGA(i,:)
    
    Nhyp = numel(OMEGA(:,k+1)); % Nhyp is the number of hypothesis at scan k
    for i = 1:Nhyp
        product = 1;
        
        DT = Ndt(i);
        for j = 1:DT
            product = product*normpdf(z(k+1) - H*xbar(k+1),B(:,:,k));
        end
        P(i,indB) = PD^Ndt(i)*(1-PD)^(Ntgt(k) - Ndt(i))*betaFT^Nft(i)*betaNT^Nnt(i)*product*P(i,sum(m(1:k))); % (16) Reid
    end
    P(i,indB) = P(i,indB)./sum(P(:,indB),1); % Normalizing the probabilities
    
    CT = [CT newTT]; % Adds the new tentative targets to the list of confirmed targets
    Ntgt(k+1) = numel(CT);
   
end

%plot(1:scans+1,x,'k-',1:scans+1,z,'r-',1:scans+1,xhat,'b-');
%legend('Process', 'Observation', 'Estimate')
