% To make things easier, this code will assume that there is a constant
% numebr of measurements at each scan for all k

[x1, z1, xbar1, H1, B1, scans1] = dataGeneration(0, 0, 1, pi/4, 1);
[x2, z2, xbar2, H2, B2, scans2] = dataGeneration(10, 0, 1, 3*pi/4, 1);

if isequal(H1,H2)
    H = H1;
end

if isequal(scans1,scans2)
    scans = scans1;
end

X = {num2cell(x1,1); num2cell(x2,1)}; % Measurements Cell Array
Z = {num2cell(z1,1); num2cell(z2,1)}; % Measurements Cell Array
Xbar = {num2cell(xbar1,1); num2cell(xbar2,1)}; % Measurements Cell Array
V = {num2cell(B1,1); num2cell(B2,1)}; % Measurements Cell Array

% Assuming size(x1) == size(x2)
[dim, ~] = size(x1);

% plot(x1(1,:),x1(2,:),'r',x2(1,:),x2(2,:),'b');

% Organizes the measurement data
mps = 2; % measurementsPerScan
maxTargets = mps*scans;
assignmentMatrix = zeros(maxTargets+1,scans,mps);

CT = linspace(1,mps,mps);
OMEGA = generateHyp(mps, CT);

% Initializes the probability matrix
P = zeros(size(OMEGA));
disp(OMEGA);
row = input('Which row of OMEGA is has probability 1? ');
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

for k = 1:scans
    
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
                
            product = product*mvnpdf(z(:,k+1), H*Xbar{1}{k+1}, B(:,:,k+1)); %%% NOTE: this should be z_j?
        end
        P(i,indB) = PD^Ndt(i)*(1-PD)^(Ntgt(k) - Ndt(i))*betaFT^Nft(i)*betaNT^Nnt(i)*product*P(i,indA); % (16) Reid
    end
    P(i,indB) = P(i,indB)./sum(P(:,indB),1); % Normalizing the probabilities
    
    
    %     [val,ind] = max(P(:,end));
    
    OMEGAs(k) = {OMEGA};
    Ps(k) = {P};
end




