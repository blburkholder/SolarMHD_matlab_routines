function [n_hat1,n_hat2] = normal_dir_var(bxs,bys,bzs,exs,eys,ezs)

    %MVA (minimum variance B)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bx_m = mean(bxs);
    by_m = mean(bys);
    bz_m = mean(bzs);
    bx_m2 = bx_m^2;
    bx_2m = mean(bxs.^2);
    by_m2 = by_m^2;
    by_2m = mean(bys.^2);
    bz_m2 = bz_m^2;
    bz_2m = mean(bzs.^2);
    bxz_m = mean(bxs.*bzs);
    bxy_m = mean(bxs.*bys);
    byz_m = mean(bys.*bzs);

    MB = [ bx_2m - bx_m2, bxy_m-bx_m*by_m, bxz_m-bx_m*bz_m;...
        bxy_m-bx_m*by_m, by_2m-by_m2, byz_m-by_m*bz_m;...
        bxz_m-bx_m*bz_m, byz_m-by_m*bz_m, bz_2m-bz_m2];
    [v,L] = eig(MB);
    n_hat1 = v(:,1)';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %MVA (maximum variance E)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Ex_m = mean(exs);
    Ey_m = mean(eys);
    Ez_m = mean(ezs);

    Ex_m2 = Ex_m^2;
    Ex_2m = mean(exs.^2);
    Ey_m2 = Ey_m^2;
    Ey_2m = mean(eys.^2);
    Ez_m2 = Ez_m^2;
    Ez_2m = mean(ezs.^2);
    Exz_m = mean(exs.*ezs);
    Exy_m = mean(exs.*eys);
    Eyz_m = mean(eys.*ezs);

    MB = [ Ex_2m - Ex_m2, Exy_m-Ex_m*Ey_m, Exz_m-Ex_m*Ez_m;...
        Exy_m-Ex_m*Ey_m, Ey_2m-Ey_m2, Eyz_m-Ey_m*Ez_m;...
        Exz_m-Ex_m*Ez_m, Eyz_m-Ey_m*Ez_m, Ez_2m-Ez_m2];
    [v,L] = eig(MB);
    n_hat2 = v(:,3)';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
