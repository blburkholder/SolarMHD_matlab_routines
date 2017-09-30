h =  findobj('type','figure');
for i = 1:length(h)
    saveas(h(i),strcat('assorted_fig',num2str(i),'.png'))
end