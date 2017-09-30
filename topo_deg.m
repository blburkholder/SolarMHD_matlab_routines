function td = topo_deg(t_8x,t_8y,t_8z)
    
    Asum = 0;

    B1 = [t_8x(1);t_8y(1);t_8z(1)];
    B2 = [t_8x(2);t_8y(2);t_8z(2)];
    B3 = [t_8x(3);t_8y(3);t_8z(3)];
    B4 = [t_8x(4);t_8y(4);t_8z(4)];
    B5 = [t_8x(5);t_8y(5);t_8z(5)];
    B6 = [t_8x(6);t_8y(6);t_8z(6)];
    B7 = [t_8x(7);t_8y(7);t_8z(7)];
    B8 = [t_8x(8);t_8y(8);t_8z(8)];

    %12 triangles with nodes ordered in a right 
    %hand manner about the cell face normal
    for i  = 1:12
        if i  == 1
            b1 = B1;
            b2 = B2;
            b3 = B6;
        elseif i == 2
            b1 = B1;
            b2 = B6;
            b3 = B4;
        elseif i == 3
            b1 = B3;
            b2 = B1;
            b3 = B4;
        elseif i == 4
            b1 = B3;
            b2 = B4;
            b3 = B7;
        elseif i == 5
            b1 = B5;
            b2 = B3;
            b3 = B8;
        elseif i == 6
            b1 = B8;
            b2 = B3;
            b3 = B7;
        elseif i == 7
            b1 = B8;
            b2 = B6;
            b3 = B5;
        elseif i == 8
            b1 = B5;
            b2 = B6; 
            b3 = B2;
        elseif i == 9
            b1 = B3;
            b2 = B2;
            b3 = B1;
        elseif i == 10
            b1 = B3;
            b2 = B5;
            b3 = B2;
        elseif i == 11
            b1 = B7;
            b2 = B4;
            b3 = B6;
        else
            b1 = B7;
            b2 = B6;
            b3 = B8;
        end

        b1_mag = sqrt(b1(1)^2 + b1(2)^2 + b1(3)^2);
        b2_mag = sqrt(b2(1)^2 + b2(2)^2 + b2(3)^2);
        b3_mag = sqrt(b3(1)^2 + b3(2)^2 + b3(3)^2);

        theta1 = acos(sum(b2.*b3)/(b2_mag*b3_mag));
        theta2 = acos(sum(b1.*b3)/(b1_mag*b3_mag));
        theta3 = acos(sum(b1.*b2)/(b1_mag*b2_mag));

        A1 = tan((theta1 + theta2 + theta3)/4);
        A2 = tan((theta1 + theta2 - theta3)/4);
        A3 = tan((theta2 + theta3 - theta1)/4);
        A4 = tan((theta3 + theta1 - theta2)/4);
        A = 4*atan(sqrt(A1*A2*A3*A4));

        if sum(b1.*cross(b2,b3)) < 0
            A = -1*A;
        end
        Asum = Asum + A;
    end

    td = Asum/(4*pi); 

    if td > 0.9999 && td < 1.0001
        td = 1;
    elseif td > -1.0001 && td < -0.9999
        td = -1;
    end



