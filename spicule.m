figure
%hold on
colormap(jet)
clf
% for ix = 2:nx-1
%     r = rho(floor(ix),:,:);
%     r = reshape(r(:), [ny nz]);
%     pcolor(y,z,log10(r'))
%     hold on
%     bby = by(floor(ix),:,:);
%     bbz = bz(floor(ix),:,:);
%     bby = reshape(bby(:), [ny nz]);
%     bbz = reshape(bbz(:), [ny nz]);
%     h = streamslice(y,z,bby',bbz',0.5);
%     set(h,'Color','yellow')
% 
%     vy = sy(ix,:,:)./rho(ix,:,:);
%     vz = sz(ix,:,:)./rho(ix,:,:);
%     
%     vy = reshape(vy(:), [ny nz]);
%     vz = reshape(vz(:), [ny nz]);
% 
%     [cor1,cor2] = meshgrid(y,z);
% 
%     quiver(cor1(1:10:end,1:10:end),cor2(1:10:end,1:10:end),vy(1:10:end,1:10:end)',vz(1:10:end,1:10:end)')
%     
%     shading interp
%     title(strcat('slice',num2str(floor(ix))));
%     colorbar
%     saveas(gcf,strcat('spicules/spicule',num2str(ix),'.png'))
%     clf
% end

for n = 0:9
    for ix = [100,178]
        load(strcat('mhd',num2str(n),'.mat'));
        r = rho(floor(ix),:,:);
        r = reshape(r(:), [ny nz]);
        pcolor(y,z,log10(r'))
        hold on
        bby = by(floor(ix),:,:);
        bbz = bz(floor(ix),:,:);
        bby = reshape(bby(:), [ny nz]);
        bbz = reshape(bbz(:), [ny nz]);
        h = streamslice(y,z,bby',bbz',0.5);
        set(h,'Color','yellow')

        vy = sy(ix,:,:)./rho(ix,:,:);
        vz = sz(ix,:,:)./rho(ix,:,:);

        vy = reshape(vy(:), [ny nz]);
        vz = reshape(vz(:), [ny nz]);

        [cor1,cor2] = meshgrid(y,z);

        quiver(cor1(1:10:end,1:10:end),cor2(1:10:end,1:10:end),vy(1:10:end,1:10:end)',vz(1:10:end,1:10:end)')

        shading interp
        title(strcat('slice',num2str(floor(ix))));
        colorbar
        saveas(gcf,strcat('spicules/spicule',num2str(n),'.',num2str(ix),'.png'))
        clf
    end
end