%% Centroids
rvals = 2*rand(246,1)-1;
elevation = asin(rvals);
azimuth = 2*pi*rand(246,1);
radii = 3*(rand(246,1).^(1/3));
[x,y,z] = sph2cart(azimuth, elevation, radii);

figure
plot3(x,y,z,'.');
axis equal

%%
p1 = find(x>0 & y>0 & z>0);
p2 = find(x>0 & y>0 & z<0);

p3 = find(x>0 & y<0 & z>0);
p4 = find(x>0 & y<0 & z<0);

p5 = find(x<0 & y>0 & z>0);
p6 = find(x<0 & y>0 & z<0);

p7 = find(x<0 & y<0 & z>0);
p8 = find(x<0 & y<0 & z<0);

x = [x(p1); x(p2); x(p3); x(p4); x(p5); x(p6); x(p7); x(p8)];
y = [y(p1); y(p2); y(p3); y(p4); y(p5); y(p6); y(p7); y(p8)];
z = [y(p1); z(p2); z(p3); z(p4); z(p5); z(p6); z(p7); z(p8)];
%% Distances
dij = sqrt((x-x').^2 + (y-y').^2 + (z-z').^2);
figure
imagesc(dij)
colorbar

%% Adjacency matrix
aij = exp(-3*dij);
aij = aij-diag(diag(aij));
aij = aij./(max(aij(:)));
figure
imagesc(aij)
colorbar

%%
dlmwrite('data/network_example.txt',aij,' ')