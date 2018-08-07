function [V_DHT] = get_DHT_frame(bxs,bys,bzs,vxs,vys,vzs)
    %%%%%first find DHT frame for E
    avg_K = zeros(3,3);
    avg_KV = zeros(3,1);
    for i = 1:length(bxs)
        Ki =  [bys(i)^2+bzs(i)^2,-bxs(i)*bys(i),-bxs(i)*bzs(i);...
            -bxs(i)*bys(i),bxs(i)^2+bzs(i)^2,-bys(i)*bzs(i);...
            -bxs(i)*bzs(i),-bys(i)*bzs(i),bxs(i)^2+bys(i)^2];
        avg_K = avg_K + Ki;
        avg_KV = avg_KV+Ki*[vxs(i);vys(i);vzs(i)];
    end

    avg_K = avg_K/length(bxs);
    avg_KV = avg_KV/length(bxs);
    V_DHT = avg_KV\avg_K;
end