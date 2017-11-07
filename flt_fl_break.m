%this part gets all the field lines which broke from t0 to tf 
%or all parts in between if you like! It takes long b/c it has
%to compute a ton of field lines

files = 10;
%s = 1-bottom 2-middle 3-top
s = 1;

%reolution
reso = 5;
%stay away from boundaries
off = 5;

load('mhd0.mat');
[coor_x1,coor_y1] = meshgrid(x,y);
%conj = zeros(length(coor_x1(off:reso:end-off,off:reso:end-off))^2,files+1);
conj = zeros(length(coor_x1(off:reso:end-off,off:reso:end-off))^2,2);
%for j = 0:files
for j = 0:1
    %load(strcat('mhd',num2str(j),'.mat'))

    [coor_x,coor_y] = velocity_perturbation(x,y,coor_x1(off:reso:end-off,off:reso:end-off),...
            coor_y1(off:reso:end-off,off:reso:end-off),time);
    pos1 = perio_stream3(x(2:end-1),y(2:end-1),z(2:end-1),bx(2:end-1,2:end-1,2:end-1),...
        by(2:end-1,2:end-1,2:end-1),bz(2:end-1,2:end-1,2:end-1),coor_x,coor_y,0*coor_x);

    for i = 1:(length(pos1)-1)
        try
            conj(i,j+1) = pos1{i}(end,3);
        catch
            i 
        end
    end

    if j > 0
        elg = reshape(conj(:,1),[size(coor_x)]);
        gle = reshape(conj(:,j+1),[size(coor_x)]);
        figure
        pcolor(coor_x, coor_y, gle-elg);
        shading interp
        colorbar
    end
    load('mhd10.mat')
end

flt_fl_break2



