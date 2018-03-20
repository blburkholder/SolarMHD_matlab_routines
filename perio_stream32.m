function [str,str2] = perio_stream32(x,y,z,bx,by,bz,sx,sy,sz)
    ds = 0.01;
    pos = stream3(x,y,z,bx,by,bz,sx,sy,sz);
    neg = stream3(x,y,z,-bx,-by,-bz,sx,sy,sz);
    str = cell(1,length(neg));
    str2 = cell(1,length(neg));

    kkk = 1;

    for i = 1:length(pos)
    i
        if length(pos{i}) > length(neg{i}) && (length(pos{i}) > 3)
            %recursing flag (periodic boundary)
            rec = 1;
            %maximum boundaries flag
%            count = 0;
            while (rec == 1 )
%                p2 = [];
%                count = count + 1
                rec = 0;
                hx = pos{i}(end,1);
                hy = pos{i}(end,2);
                hz = pos{i}(end,3);
                %these are the interpolated b values at the
                %point before the boundary
                bbx = interp3(x,y,z,bx,hx,hy,hz);
                bby = interp3(x,y,z,by,hx,hy,hz);
                bbz = interp3(x,y,z,bz,hx,hy,hz);
                bsx = pos{i}(end,1) + bbx*ds/sqrt(bbx^2+bby^2+bbz^2);
                bsy = pos{i}(end,2) + bby*ds/sqrt(bbx^2+bby^2+bbz^2);
                bsz = pos{i}(end,3) + bbz*ds/sqrt(bbx^2+bby^2+bbz^2);
                if bsz <= z(end) && bsz >= z(1)
                    if (bsx <= x(1))
                        rec = 1;
                        %need to do a step to bring you exactly to boundary:
                        %this is the distance along x you need to go
                        %to hit boundary
                        step_length = hx - x(1);
                        %taking an euler step
                        dt = abs(step_length/bbx);
                        pos{i}(end+1,1) = x(1);
                        pos{i}(end,2) = pos{i}(end-1,2) + bby*dt;
                        pos{i}(end,3) = pos{i}(end-1,3) + bbz*dt;

                        p2 = stream3(x,y,z,bx,by,bz,pos{i}(end,1),y(end)-pos{i}(end,2),pos{i}(end,3));
                        pos{i} = [vertcat(pos{i},p2{1})];
                        str2{kkk} = p2{1};
                        kkk = kkk + 1;
                    elseif (bsx >= x(end))
                        rec = 1;
                        step_length = x(end) - hx;
                        dt = abs(step_length/bbx);
                        pos{i}(end+1,1) = x(end);
                        pos{i}(end,2) = pos{i}(end-1,2) + bby*dt;
                        pos{i}(end,3) = pos{i}(end-1,3) + bbz*dt;

                        p2 = stream3(x,y,z,bx,by,bz,pos{i}(end,1),y(end)-pos{i}(end,2),pos{i}(end,3));
                        pos{i} = [vertcat(pos{i},p2{1})];
                        str2{kkk} = p2{1};
                        kkk = kkk + 1;
                    elseif (bsy <= y(1))
                        rec = 1;
                        step_length = hy - y(1);
                        dt = abs(step_length/bby);
                        pos{i}(end+1,1) = pos{i}(end,1) + bbx*dt;
                        pos{i}(end,2) = y(1);
                        pos{i}(end,3) = pos{i}(end-1,3) + bbz*dt;

                        p2 = stream3(x,y,z,bx,by,bz,x(end)-pos{i}(end,1),pos{i}(end,2),pos{i}(end,3));
                        pos{i} = [vertcat(pos{i},p2{1})];
                        str2{kkk} = p2{1};
                        kkk = kkk + 1;
                    elseif (bsy >= y(end))
                        rec = 1;
                        step_length = y(end) - hy;
                        dt = abs(step_length/bby);
                        pos{i}(end+1,1) = pos{i}(end,1) + bbx*dt;
                        pos{i}(end,2) = y(end);
                        pos{i}(end,3) = pos{i}(end-1,3) + bbz*dt;

                        p2 = stream3(x,y,z,bx,by,bz,x(end)-pos{i}(end,1),pos{i}(end,2),pos{i}(end,3));
                        pos{i} = [vertcat(pos{i},p2{1})];
                        str2{kkk} = p2{1};
                        kkk = kkk + 1;
                    else
                        rec = 1;
                        pos{i}(end+1,1) = pos{i}(end,1) + bbx*ds/sqrt(bbx^2+bby^2+bbz^2);
                        pos{i}(end,2) = pos{i}(end-1,2) + bby*ds/sqrt(bbx^2+bby^2+bbz^2);
                        pos{i}(end,3) = pos{i}(end-1,3) + bbz*ds/sqrt(bbx^2+bby^2+bbz^2);  
                    end
                end
            end
            str{i} = pos{i};
        elseif length(neg{i}) > 3
            rec = 1;
            while (rec == 1)
                p2 = [];

                rec = 0;
                hx = neg{i}(end,1);
                hy = neg{i}(end,2);
                hz = neg{i}(end,3);
                bbx = interp3(x,y,z,bx,hx,hy,hz);
                bby = interp3(x,y,z,by,hx,hy,hz);
                bbz = interp3(x,y,z,bz,hx,hy,hz);
                bsx = neg{i}(end,1) - bbx*ds/sqrt(bbx^2+bby^2+bbz^2);
                bsy = neg{i}(end,2) - bby*ds/sqrt(bbx^2+bby^2+bbz^2);
                bsz = neg{i}(end,3) - bbz*ds/sqrt(bbx^2+bby^2+bbz^2);
                if bsz <= z(end) && bsz >= z(1)
                    if (bsx <= x(1))
                        rec = 1;
                        step_length = hx - x(1);
                        dt = abs(step_length/bbx);
                        neg{i}(end+1,1) = x(1);
                        neg{i}(end,2) = neg{i}(end-1,2) - bby*dt;
                        neg{i}(end,3) = neg{i}(end-1,3) - bbz*dt;

                        p2 = stream3(x,y,z,-bx,-by,-bz,neg{i}(end,1),y(end)-neg{i}(end,2),neg{i}(end,3));
                        neg{i} = [vertcat(neg{i},p2{1})];
                        str2{kkk} = p2{1};
                        kkk = kkk + 1;

                    elseif (bsx >= x(end))
                        rec = 1;
                        step_length = x(end) - hx;
                        dt = abs(step_length/bbx);
                        neg{i}(end+1,1) = x(end);
                        neg{i}(end,2) = neg{i}(end-1,2) - bby*dt;
                        neg{i}(end,3) = neg{i}(end-1,3) - bbz*dt;

                        p2 = stream3(x,y,z,-bx,-by,-bz,neg{i}(end,1),y(end)-neg{i}(end,2),neg{i}(end,3));
                        neg{i} = [vertcat(neg{i},p2{1})];
                        str2{kkk} = p2{1};
                        kkk = kkk + 1;
                    elseif (bsy <= y(1))
                        rec = 1;
                        step_length = hy - y(1);
                        dt = abs(step_length/bby);
                        neg{i}(end+1,1) = neg{i}(end,1) - bbx*dt;
                        neg{i}(end,2) = y(1);
                        neg{i}(end,3) = neg{i}(end-1,3) - bbz*dt;

                        p2 = stream3(x,y,z,-bx,-by,-bz,x(end)-neg{i}(end,1),neg{i}(end,2),neg{i}(end,3));
                        neg{i} = [vertcat(neg{i},p2{1})];
                        str2{kkk} = p2{1};
                        kkk = kkk + 1;
                    elseif (bsy >= y(end))
                        rec = 1;
                        step_length = y(end) - hy;
                        dt = abs(step_length/bby);
                        neg{i}(end+1,1) = neg{i}(end,1) - bbx*dt;
                        neg{i}(end,2) = y(end);
                        neg{i}(end,3) = neg{i}(end-1,3) - bbz*dt;

                        p2 = stream3(x,y,z,-bx,-by,-bz,x(end)-neg{i}(end,1),neg{i}(end,2),neg{i}(end,3));
                        neg{i} = [vertcat(neg{i},p2{1})];
                        str2{kkk} = p2{1};
                        kkk = kkk + 1;
                    else
                        rec = 1;
                        neg{i}(end+1,1) = neg{i}(end,1) - bbx*ds/sqrt(bbx^2+bby^2+bbz^2);
                        neg{i}(end,2) = neg{i}(end-1,2) - bby*ds/sqrt(bbx^2+bby^2+bbz^2);
                        neg{i}(end,3) = neg{i}(end-1,3) - bbz*ds/sqrt(bbx^2+bby^2+bbz^2);  
                    end
                end
            end
            str{i} = neg{i};
        end
    end
