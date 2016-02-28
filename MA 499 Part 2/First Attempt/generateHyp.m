function [OMEGA] = generateHyp(m, n)
% n is the number of a priori targets
% m is the number of measurements
TT = linspace(n+1, n+m,m); % Array of tentative targets
OMEGA = [0 linspace(1,n,n) n+1]'; % Initial OMEGA

for j = 2:m
    % Next Measurement
    [r, ~] = size(OMEGA);
    nextCol = zeros(r,1);
    for i = 1:n
        nextCol = [nextCol; i*ones(r,1)];
    end
    nextCol = [nextCol; TT(j)*ones(r,1)];
    OMEGA = repmat(OMEGA,n+2,1);
    OMEGA = [OMEGA nextCol];
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

% Sorting
%OMEGA = sortrows(OMEGA,1);
end
