% plot the simulation results
% nImages = size(CT_Number_inN, 1);
startTime = 2000;
endTime = 3000;
Xcoor = repmat(Zones(:, 2), MaxCL, 1);
CL = Nodes(:, 3);
%figure;
for i = 1 : (endTime - startTime)
    k = i + startTime;
    EVs_t_temp = EVs_t{k, 1};
    Rebalance = find(~isnan(EVs_t_temp(:, 5)));
    if ~isnan(Rebalance)
        scatter(EVs_t_temp(Rebalance, 2), EVs_t_temp(Rebalance, 4), 30, 'blue', 'filled', 's') % 30 is the size for rebalancing vehicles (squares)
        hold on
    end
    
    Number_of_available_ev = [];
    for j = 1 : length(Xcoor)
        Number_of_available_ev = [Number_of_available_ev, length(find(EVs_t_temp(:, 8) == 1 & EVs_t_temp(:, 2) == Xcoor(j) & ceil(EVs_t_temp(:, 4)) == CL(j)))];
    end
    Station_size = 1 + 10 * Number_of_available_ev;
    scatter(Xcoor, CL, Station_size, 'green', 'filled', 's')
    hold on
    scatter(Xcoor, CL, 50 * ones(length(Xcoor), 1), 'black', '+') 
    
    Station_size1 = 1 + 50 * CT_Number_inN(k, :);
    Station_color1 = [];
    for j = 1 : N * MaxCL
        if CT_Number_inN(k, j) >= 3
            Station_color1 = [Station_color1; 'red'];
        else
            Station_color1 = [Station_color1; 'yellow'];
        end
    end
    scatter(Xcoor, CL, Station_size1, Station_color1, 'filled')
    
    hold off
    title(['Vehicle Relocation Simulation,  t = ' num2str( k) ])
    drawnow
    frame = getframe(1);
    im{i} = frame2im(frame);
end
close;

filename = 'testAnimated.gif'; % Specify the output file name
for i = 1 : (endTime - startTime)
    [A, map] = rgb2ind(im{i}, 256);
    if i == 1
        imwrite(A, map, filename, 'gif', 'LoopCount', Inf, 'DelayTime', 0.1);
    else
        imwrite(A, map, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    end
end