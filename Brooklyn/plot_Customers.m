% plot the simulation results
nImages = size(CT_Number_inS, 1);
axis([0 Area_size 0 Area_size])
%figure;
for i = 1 : nImages
    Station_size = 1 + 10 * CT_Number_inS(i, :);
    Station_color = [];
    Station_color1 = [];
    for j = 1 : N
        if CT_Number_inS(i, j) >= 10
            Station_color = [Station_color; 1, 0, 0];
        else
            Station_color = [Station_color; 0, 0, 1];
        end
    end
    for j = 1 : N
        Station_color1 = [Station_color1; 0, 0, 0];
    end    

    scatter(Zones(:, 2), Zones(:, 3), Station_size, Station_color, 'filled')
    hold on
    scatter(Zones(:, 2), Zones(:, 3), 50 * ones(N, 1), Station_color1, '+')
    hold off
    %scatter(Zones(:,2),Zones(:,3))
    title(['Number of Customers in Each Station,  t = ' num2str( i) ])
    drawnow
    frame = getframe(1);
    im{i} = frame2im(frame);
end
close;

filename = 'testAnimated.gif'; % Specify the output file name
for i = 1:nImages
    [A,map] = rgb2ind(im{i},256);
    if i == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.1);
    end
end