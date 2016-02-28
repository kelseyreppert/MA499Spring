function [OMEGA, P] = expandHyp(OMEGA, P, m, CT, NT)
% m is the number of measurements

ct = numel(CT); % Confirmed Targets

for j = 1:m
    % Next Measurement
    [r, ~] = size(OMEGA);
    nextColOMEGA = zeros(r,1);
    
    % nCO = floor(linspace(0,ct*r-1,ct*r)./r)+1;
    
    for i = 1:ct
        nextColOMEGA = [nextColOMEGA; CT(i)*ones(r,1)];
    end
    nextColOMEGA = [nextColOMEGA; NT(j)*ones(r,1)];
    OMEGA = repmat(OMEGA,ct+2,1);
    OMEGA = [OMEGA nextColOMEGA];
    P = [repmat(P,ct+2,1) zeros((ct+2)*r,1)];
end

if m > 1
    nextPart = OMEGA(:,end-m+1:end);
    % Checking for duplicate target assignments
    [r, ~] = size(nextPart);
    toDelete = [];
    for j = 1:r
        row = nextPart(j,:);
        
        if numel(unique(row)) ~= numel(row) && sum(row) ~= 0
            toDelete = [toDelete j];
        end
    end
    OMEGA(toDelete,:) = [];
    P(toDelete,:) = [];
end

end
