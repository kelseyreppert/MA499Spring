% To make things easier, this code will assume that there is a constant
% numebr of measurements at each scan for all k

% Initial conditions and data generation
x10 = 0;
y10 = 0;
speed1 = 1;
angle1 = pi/4;
radius1 = 1;
[x1, z1, Q1, R1, dt] = dataGeneration(x10, y10, speed1, angle1, radius1);

x20 = 10;
y20 = 10;
speed2 = 1;
angle2 = 3*pi/4;
radius2 = 1;
[x2, z2, Q2, R2, dt] = dataGeneration(x20, y20, speed2, angle2, radius2);

% MHT variables
[dim, scans] = size(x1);
PD = 0.9;
H1 = eye(dim);
GAMMA1 = eye(dim);
H2 = eye(dim);
GAMMA2 = eye(dim);

% Initializing the structure arrays
Xstruct = struct('x1',x1, 'x2', x2);
Zstruct = struct('z1',z1, 'z2', z2);
Xbarstruct = struct('x1bar',zeros(dim,scans),'x2bar',zeros(dim,scans));
Xhatstruct = struct('x1hat',zeros(dim,scans),'x2hat',zeros(dim,scans));
Pbarstruct = struct('Pbar1',zeros(dim,dim,scans),'Pbar2',zeros(dim,dim,scans));
Phatstruct = struct('Phat1',zeros(dim,dim,scans),'Phat2',zeros(dim,dim,scans));
Bstruct = struct('B1',zeros(dim,dim,scans),'B2',zeros(dim,dim,scans));
Kstruct = struct('K1',zeros(dim,dim,scans),'K2',zeros(dim,dim,scans));

% Initial conditions
Xbarstruct.xbar1(:,1) = Zstruct.z1(:,1); 
Xhatstruct.xhat1(:,1) = Zstruct.z1(:,1); 
Pbarstruct.Pbar1(:,:,1) = eye(dim); 
Phatstruct.Phat1(:,:,1) = eye(dim); 
Bstruct.B1(:,:,1) = eye(dim); 
Kstruct.K1(:,:,1) = eye(dim); 
Xbarstruct.xbar2(:,1) = Zstruct.z2(:,1); 
Xhatstruct.xhat2(:,1) = Zstruct.z2(:,1); 
Pbarstruct.Pbar2(:,:,1) = eye(dim); 
Phatstruct.Phat2(:,:,1) = eye(dim); 
Bstruct.B2(:,:,1) = eye(dim); 
Kstruct.K2(:,:,1) = eye(dim); 

% Organizes the measurement data
mps = 2; % measurementsPerScan
maxTargets = mps*scans;
assignmentMatrix = zeros(maxTargets+1,scans,mps);

CT = linspace(1,mps,mps);
OMEGA = generateHyp(mps, CT);

% Initializes the probability matrix
P = zeros(size(OMEGA));
disp(OMEGA);
row = input('Which row of OMEGA has probability 1? ');
P(row,mps) = 1;
toDelete = find(P(:,end)==0);
P(toDelete,:) = [];
OMEGA(toDelete,:) = [];

% Creates the cell arrays to store the matricies
OMEGAs = cell(1,scans);
OMEGAs{1} = OMEGA;
Ps = cell(1,scans);
Ps{1} = P;

NT = linspace(mps+1,maxTargets,maxTargets-mps); % All the possible new targets

Ntgt = zeros(scans);
Ntgt(1) = mps;

for k = 1:scans
    PHI1 = [1 0 -dt*speed1*sin(Xstruct.x1(3,k)); 0 1 dt*speed1*cos(Xstruct.x1(3,k)); 0 0 1 - dt*speed1/radius1];
    Xbarstruct.xbar1(:,k+1) = PHI1*Xhatstruct.xhat1(:,k);
    Pbarstruct.Pbar(:,:,k+1) = PHI1*Phatstruct.Phat1(:,:,k)*PHI1' + GAMMA1*Q1*GAMMA1';
    Xhatstruct.xhat(:,k+1) = Xbarstruct.xbar1(:,k+1) + Kstruct.K1(:,:,k+1)*(Zstruct.z1(:,k+1) - H1*Xbarstruct.xbar1(:,k+1));
    Phatstruct.Phat(:,:,k+1) = Pbarstruct.Pbar1(:,:,k+1) - Pbarstruct.Pbar1(:,:,k+1)*H1'*inv(H1*Pbarstruct.Pbar1(:,:,k+1)*H1' + R1)*H1*Pbarstruct.Pbar1(:,:,k+1);
    Kstruct.K1(:,:,k+1) = Phatstruct.Phat1(:,:,k+1)*H1'*inv(R1);
    Bstruct.B1(:,:,k+1) = H1*Pbarstruct.Pbar1(:,:,k+1)*H1' + R1;
    
    PHI2 = [1 0 -dt*speed2*sin(Xstruct.x2(3,k)); 0 1 dt*speed2*cos(Xstruct.x2(3,k)); 0 0 1 - dt*speed2/radius2];
    Xbarstruct.xbar2(:,k+1) = PHI2*Xhatstruct.xhat2(:,k);
    Pbarstruct.Pbar(:,:,k+1) = PHI2*Phatstruct.Phat2(:,:,k)*PHI2' + GAMMA2*Q2*GAMMA2';
    Xhatstruct.xhat(:,k+1) = Xbarstruct.xbar2(:,k+1) + Kstruct.K2(:,:,k+1)*(Zstruct.z2(:,k+1) - H2*Xbarstruct.xbar2(:,k+1));
    Phatstruct.Phat(:,:,k+1) = Pbarstruct.Pbar2(:,:,k+1) - Pbarstruct.Pbar2(:,:,k+1)*H2'*inv(H2*Pbarstruct.Pbar2(:,:,k+1)*H2' + R2)*H2*Pbarstruct.Pbar2(:,:,k+1);
    Kstruct.K2(:,:,k+1) = Phatstruct.Phat2(:,:,k+1)*H2'*inv(R2);
    Bstruct.B2(:,:,k+1) = H2*Pbarstruct.Pbar2(:,:,k+1)*H2' + R2;
    
    [OMEGA, P] = expandHyp(OMEGA, P, mps, CT, NT(mps*(k-1)+1:mps*k));
    [r,c] = size(OMEGA);
    
    indA = c-mps+1;
    indB = c;
    
    Nft = sum(OMEGA(:,indA:indB)==0, 2); % Number of false targets for each hypothesis OMEGA(i,:)
    Nnt = zeros(size(Nft));
    for i = 1:numel(NT(indA:indB))
        Nnt = Nnt + sum(OMEGA(:,indA:indB)==NT(i),2);
    end
    Ndt = mps - Nnt - Nft; % Number of detections for each hypothesis OMEGA(i,:)
    
    Nhyp = numel(OMEGA(:,k+1)); % Nhyp is the number of hypothesis at scan k
    for i = 1:Nhyp
        product = 1;
        
        DT = Ndt(i);
        for j = 1:DT
           % if OMEGA(i+1,indA+j-1)
                
            product = product*mvnpdf(z(:,k+1), H*Xbarcell{1}{k+1}, B(:,:,k+1)); %%% NOTE: this should be z_j?
        end
        P(i,indB) = PD^Ndt(i)*(1-PD)^(Ntgt(k) - Ndt(i))*betaFT^Nft(i)*betaNT^Nnt(i)*product*P(i,indA); % (16) Reid
    end
    P(i,indB) = P(i,indB)./sum(P(:,indB),1); % Normalizing the probabilities
    
    
    %     [val,ind] = max(P(:,end));
    
    OMEGAs(k) = {OMEGA};
    Ps(k) = {P};
end




