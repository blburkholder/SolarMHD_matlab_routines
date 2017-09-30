function update_plot2(z,net1,tot1,net2,tot2,net3,tot3,net4,tot4,net5,tot5,i)
    if i < 1
        net = net1;
        tot = tot1;
    elseif i >= 1 && i < 2
        net = net2;
        tot = tot2;
    elseif i >= 2 && i < 3
        net = net3;
        tot = tot3;
    elseif i >= 3 && i < 4
        net = net4;
        tot = tot4;
    elseif i >= 4 && i <= 5
        net = net5;
        tot = tot5;
    end

    subplot(3,2,5)
    plot(z(2:end-1),net - net1)
    title('net flux diff t_i - t_0')
    subplot(3,2,6)
    plot(z(2:end-1),tot - tot1)
    title('total flux diff t_i - t_0')
    subplot(3,2,3)
    plot(z(2:end-1),net)
    title('net flux t_i')
    subplot(3,2,4)
    semilogy(z(2:end-1),tot)
    title('total flux t_i')
end