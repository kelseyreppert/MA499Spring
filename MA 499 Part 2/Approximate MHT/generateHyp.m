function [OMEGA, NT] = generateHyp(m, CT) % NOTE: m == numel(NT)
% m is the number of measurements
% CT is the vector of confirmed targets
% TT is the vector of tentative targets

OMEGA = [0 CT]'; % Initial OMEGA

for j = 2:m
    [r, ~] = size(OMEGA);
    
    % False Target Assignments
    nextCol = zeros(r,1); 
    
    % Confirmed Target Assignments
    for i = 1:numel(CT)
        nextCol = [nextCol; CT(i)*ones(r,1)]; 
    end
    
    OMEGA = [repmat(OMEGA,numel(CT)+1,1) nextCol];
end

% Checking for duplicate target assignments 
[r, ~] = size(OMEGA);
toDelete = [];
for j = 1:r
    row = OMEGA(j,:);
    
    if numel(unique(row)) ~= numel(row) && sum(row) ~= 0
        toDelete = [toDelete j];
    end
end
OMEGA(toDelete,:) = [];

end
