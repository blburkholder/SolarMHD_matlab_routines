load('mentrop2_t.mat')

%property which will paint bottom/top boundary
prop = ftmass;

%prop quantity maximum in case of extraneous values
maxi = 150;
%where to start flines
n = nz-1;

count = 500;
[row,col] = size(prop);
bs = zeros(count,2);

center = [45,80];
bs(1,:) = center;
spread = 10;
for i = 2:count
bs(i,:) = center + (2*spread*rand(1,2)-spread);
end

prop(prop > maxi) = maxi;
prop(prop < -maxi) = -maxi;
prop = prop/max(max(prop));

figure
hold on
pcolor(ypo,xpo,prop)
shading interp
colormap(jet)
colorbar
scatter(bs(:,2),bs(:,1),[2],'k.')

load('mentrop2.mat')
prop = zfin1;
prop(prop > maxi) = maxi;
prop(prop < -maxi) = -maxi;
prop = prop/max(max(prop));
figure
hold on
pcolor(ypo,xpo,prop)
shading interp
colormap(jet)
colorbar
streamslice(y,x,sy(:,:,2)./rho(:,:,2),sx(:,:,2)./rho(:,:,2))

plines = perio_stream3(y(2:end-1),x(2:end-1),z(2:end-1),...
    by(2:end-1,2:end-1,2:end-1),bx(2:end-1,2:end-1,2:end-1),...
    bz(2:end-1,2:end-1,2:end-1),bs(:,2),bs(:,1),...
    z(n)*ones(length(bs(:,1)),1));

%this will put wittle red stars at the conjugate footpoints
conj = zeros(length(plines),3);
for l = 1:length(plines)
    conj(l,1) = plines{l}(end,1);
    conj(l,2) = plines{l}(end,2);
    conj(l,3) = plines{l}(end,3);
end
scatter(conj(:,1),conj(:,2),'rp')

daspect([1 1 1])
