function OMEGAsorted = plotHyp(OMEGA)
% c is a vector that contains the size of OMEGA given m measurements.
% c(i) represents the size of OMEGA at stage i-1
% c = zeros(m,1);
% c(1) = 1;
% c(2) = n+2;
% for i=2:m
% c(i+1) = c(i)*(n+1) + 2;
% end

% Shows the hypothesis tree for one measurement
% s = [ones(1,c(2))];
% t = [linspace(c(1)+1,c(2)+1,c(2)-c(1)+1)];
% names = {'Start' '0' '1' '2' '3'};
% weights = ones(size(s));
% G = graph(s,t,weights,names);
% plot(G)
[~, m] = size(OMEGA);
OMEGAsorted = sortrows(OMEGA,linspace(1,m,m));
imagesc(OMEGAsorted)
colorbar('Ticks',unique(OMEGA(:)));
axis xy
title('Chronograph of Omega: The set of all possible hyoptheses')

end