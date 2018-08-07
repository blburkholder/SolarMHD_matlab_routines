function [fx1,fy1] = fix_coords(x,y,fx1,fy1)
    [rows,cols] = size(fx1);

    for i = 1:rows*cols
        if fx1(i) > x(end)
            fx1(i) = x(end) - mod(fx1(i),x(end));
            fy1(i) = y(end) - fy1(i);
        end
        if fx1(i) < x(1)
            fx1(i) = x(1) + abs(fx1(i));
            fy1(i) = y(end) - fy1(i);
        end
        if fy1(i) > y(end)
            fy1(i) = y(end) - mod(fy1(i),y(end));
            fx1(i) = x(end) - fx1(i);
        end
        if fy1(i) < y(1)
            fy1(i) = y(1) + abs(fy1(i));
            fx1(i) = x(end) - fx1(i);
        end
    end



%     for i = 1:length(fx1)
%         for j = 1:length(fy1)
%             if fx1(i,j) > x(end)
%                 fx1(i,j) = x(end) - mod(fx1(i,j),x(end));
%                 fy1(i,j) = y(end) - fy1(i,j);
%             end
%             if fx1(i,j) < x(1)
%                 fx1(i,j) = abs(fx1(i,j));
%                 fy1(i,j) = y(end) - fy1(i,j);
%             end
%             if fy1(i,j) > y(end)
%                 fy1(i,j) = y(end) - mod(fy1(i,j),y(end));
%                 fx1(i,j) = x(end) - fx1(i,j);
%             end
%             if fy1(i,j) < y(1)
%                 fy1(i,j) = abs(fy1(i,j));
%                 fx1(i,j) = x(end) - fx1(i,j);
%             end
%         end
%     end
end