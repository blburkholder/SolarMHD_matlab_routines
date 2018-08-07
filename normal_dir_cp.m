function [n_hat1,n_hat2] = normal_dir(binterped_x,binterped_y,binterped_z,edht_interped_x,edht_interped_y,edht_interped_z,imjdbs,u_jdb,d_jdb)

    %B cross product method
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Bu_x = mean(binterped_x(imjdbs:u_jdb));
    Bu_y = mean(binterped_y(imjdbs:u_jdb));
    Bu_z = mean(binterped_z(imjdbs:u_jdb));
    Bd_x = mean(binterped_x(d_jdb:imjdbs));
    Bd_y = mean(binterped_y(d_jdb:imjdbs));
    Bd_z = mean(binterped_z(d_jdb:imjdbs));

    BuBd =  cross([Bu_x,Bu_y,Bu_z],[Bd_x,Bd_y,Bd_z]); 
    n_hat1 = BuBd/norm(BuBd);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %E cross product method
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Eu_x = mean(edht_interped_x(imjdbs:u_jdb));
    Eu_y = mean(edht_interped_y(imjdbs:u_jdb));
    Eu_z = mean(edht_interped_z(imjdbs:u_jdb));
    Ed_x = mean(edht_interped_x(d_jdb:imjdbs));
    Ed_y = mean(edht_interped_y(d_jdb:imjdbs));
    Ed_z = mean(edht_interped_z(d_jdb:imjdbs));

    Eu_m_Ed =  [Eu_x;Eu_y;Eu_z] - [Ed_x;Ed_y;Ed_z]; 
    n_hat2 = Eu_m_Ed'/norm(Eu_m_Ed);
end
