
% canonical profile

z = -4.5:0.1:0;
cax = 1.49;
zax = -1.1;
eta = 2*(z-zax)/1;

cz = cax*(1+(0.0113/2)*(exp(eta)-eta-1));

figure; plot(cz,z); %set(gca, 'YDir', 'reverse')