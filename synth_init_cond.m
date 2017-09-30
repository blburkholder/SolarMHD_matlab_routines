figure
%sin(mx)*sin(ny)
% cmn1 = [0 0 0 0 0 0 0 0;
%         0 0 0 0 0 0 0 0;
%         0 0 0 0 0 0 0 0;
%         0 0 0 0 0 0 0 0;
%         0 0 0 0 0 0 0 0;
%         0 0 0 0 0 0 0 0;
%         0 0 0 0 0 0 0 0;
%         0 0 0 0 0 0 0 0];
% %cos(mx)*cos(ny)
% cmn2 = [0 0 0 0 0 0 0 0;
%         0 0 0 0 0 0 0 0;
%         0 0 0 0 0 0 0 0;
%         0 0 0 0 0 0 0 0;
%         0 0 0 0 0 0 0 0;
%         0 0 0 0 0 0 0 0;
%         0 0 0 0 0 0 0 0;
%         0 0 0 0 0 0 0 0];
%}

xmodes = 8;
ymodes = 8;

%cmn1 = rand(xmodes,ymodes) - 0.5;
%cmn2 = rand(xmodes,ymodes) - 0.5;

%cmn1 = randn(xmodes,ymodes);
%cmn2 = randn(xmodes,ymodes);

%cmn1(cmn1 < 0.2 & cmn1 > -0.2) = 0;
%cmn2(cmn2 < 0.2 & cmn2 > -0.2) = 0;

flux_control = 0; %this should always be zero
[a,frac,cmn1,cmn2] = synthetic_magnetogram(x,y,cmn1,cmn2,xmodes,ymodes,bz(:,:,2),flux_control);
%%flux_control%%
%%0 even-even & odd-odd modes allowed
%%1 fourier decomposition of random gaussians (doesnt work)

if frac ~= 1.0
    cmn1 = 1/frac*cmn1;
    cmn2 = 1/frac*cmn2;
    [a,~,cmn1,cmn2] = synthetic_magnetogram(x,y,cmn1,cmn2,xmodes,ymodes,bz(:,:,2),flux_control);
end

[net,tot] = flux(a,x(2)-x(1),y(2)-y(1));

f = gca;
%colormap(jet)
pcolor(x,y-y(end-1)/2,a')
title(strcat('net flux-',num2str(net),'-total flux-',num2str(tot)))
shading interp
colorbar
hold on

velocity_perturbation

cmn11 = zeros(260,256);
cmn22 = zeros(260,256);
%cmn11(257:260,:) = NaN;
cmn11(1:xmodes,1:ymodes) = cmn1;
cmn22(1:xmodes,1:ymodes) = cmn2;


% fileID = fopen('synth_fft_4','w');
% fprintf(fileID,'%15.7f %15.7f %15.7f %15.7f %15.7f\n',cmn11);
% fprintf(fileID,'%15.7f %15.7f %15.7f %15.7f %15.7f\n',cmn22);
% fclose(fileID);
