function [OMEGA, P] = expandHyp(OMEGA, P, m, CT)
% m is the number of measurements

TT = linspace(1,m,m) + max(max(OMEGA)); % Tentative Targets
ct = numel(CT); % Confirmed Targets

for j = 1:m
    % Next Measurement
    [r, ~] = size(OMEGA);
    nextColOMEGA = zeros(r,1);
    
    nCO = floor(linspace(0,ct*r-1,ct*r)./r)+1;
    
    for i = 1:ct
        nextColOMEGA = [nextColOMEGA; CT(i)*ones(r,1)];
    end
    nextColOMEGA = [nextColOMEGA; TT(j)*ones(r,1)];
    OMEGA = repmat(OMEGA,ct+2,1);
    OMEGA = [OMEGA nextColOMEGA];
    P = [repmat(P,ct+2,1) zeros((ct+2)*r,1)];
end


% Checking for duplicate target assignments -- This is no longer needed, as
% assignments such as 1-1 and 2-2 are valid
% [r, ~] = size(OMEGA);
% toDelete = [];
% for j = 1:r
%     row = OMEGA(j,:);
%     row(row == 0) = [];
%     if numel(unique(row)) ~= numel(row)
%         toDelete = [toDelete j];
%     end
% end
% OMEGA(toDelete,:) = [];
% P(toDelete,:) = [];

% Sorting
%OMEGA = sortrows(OMEGA,1);
end
