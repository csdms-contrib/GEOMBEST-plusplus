function calc_marshwidth(i,j,k) %initial width, run number, RSLR
%function calc_marshwidth(i,j)

% Calculates the final marsh width from the modern runs file

if i == 1
    init_width = 0;
%     filebump = 0;
elseif i == 2
    init_width = 400;
%     filebump = 10;
elseif i == 3
    init_width = 1850;
%     filebump = 20;
end

island_width = 275;
results = zeros(1000,6);
%RSLR = j;
ii = 1;

% for OW = 1:10
%     for RIi = (1:10)
%         RI = RIi * 2;
%         for RSLR = 1:10
for a = j
            filethread=a;
            filename = ['C:/GEOMBEST+/Output' num2str(filethread) '/surface.mat'];

            if exist(filename) == 2
                load(filename);

            filename2 = ['C:/GEOMBEST+/Output' num2str(filethread) '/xcentroids.mat'];
            load(filename2);

            filename3 = ['C:/GEOMBEST+/Output' num2str(filethread) '/SL.mat'];
            load(filename3)

RSLR = k;
            total_SLR = 10*RSLR*(numel(surface(:,1))-11); % Total sea level rise, in mm
            x = xcentroids(1,:); % X location of each cell
            y_final = surface(end,:); % Surface elevation from the final time step
            numcell = numel(y_final); % Number of surface cells

            for n = 1:numcell
                if y_final(n) >= SL(end)
                    x_islandfront = n; % The first cell above sea level gives the X location of the shoreline at the island front
                    break
                else
                    n = n + 1;
                end
            end

            x_islandback = numel(x); 

            for n = x_islandfront:numcell
                if y_final(n) < SL(end) - 0.25
                    x_islandback = n - 1; % The first cell past the island front that is below sea level gives the X location of the first estuarine cell
                    break
                else
                    n = n + 1;
                end
            end
                   
            residual = ((y_final(x_islandback) - (SL(end) - 0.4)) / 0.9)*50;
            
            if x_islandback == numel(x)
                subaerial_width = 2200; % If there are no cells below sea level behind the island, the total width is at its maximum, 2200 meters
            else
                subaerial_width = x(x_islandback - 1) + residual - x(x_islandfront); % The total subaerial width is equal to the distance from the island front to the first estuarine cell
            end

            marsh_width = subaerial_width - island_width; % The marsh width is the total subaerial width minus the island width
            
            if marsh_width < 0
                final_width = 0;
            else
                final_width = marsh_width;
            end
            
            if marsh_width < 0
                marsh_width = 0; % If the islandwidth has narrowed, marsh width is zero
            end
            
            % If the estuary is either completely filled or completely
            % empty of marsh, check previous time steps to determine when
            % the marsh first enclosed the estuary / disappeared
            if marsh_width ~= init_width
                while marsh_width == 0 & numel(surface(:,1)) > 11 | marsh_width == 1850 & numel(surface(:,1)) > 11

                    %Re-calculate marsh width to see if it is still at a minimum/maximum
                    SL(end) = [];
                    surface(end,:) = [];
                    y_final = surface(end,:);
                    numcell = numel(y_final);

                    for n = 1:numcell
                        if y_final(n) >= SL(end)
                            x_islandfront = n;
                            break
                        else
                            n = n + 1;
                        end
                    end

                    x_islandback = numel(x);

                    for n = x_islandfront:numcell
                        if y_final(n) < SL(end) - 0.25
                            x_islandback = n - 1;
                            break
                        else
                            n = n + 1;
                        end
                    end

                    if x_islandback < 162
                        residual = ((y_final(x_islandback) - (SL(end) - 0.4)) / 0.9)*50;
                    end

                    if x_islandback == numel(x)
                        subaerial_width = 2200;
                    else
                        subaerial_width = x(x_islandback - 1) + residual - x(x_islandfront);
                    end

                    marsh_width = subaerial_width - island_width;

                    if marsh_width < 0
                        [total_SLR numel(surface(:,1)) RSLR]
                        marsh_width = 0;
                    end

                    total_SLR = total_SLR - RSLR*10; % Reduce the total sea level rise to the level it was at when the marsh first filled/disappeared
                end
            end

            effective_marsh_width = marsh_width;
            effective_init_width = init_width;
            
            % Width change is the difference between the marsh initial and final width
            if final_width > 1700
                effective_marsh_width = final_width / 2;
                final_width = 42570 - x(x_islandfront) - island_width
            end

            if init_width == 1850
                effective_init_width = 925;
            end
            
            width_change = effective_marsh_width - effective_init_width;

            % The rate of change is equal to the change in width divided by
            % the time it took for the sea level to rise to the final sea
            % level (1 m), which is variable depending on the rate of sea
            % level rise
            
            T = total_SLR/RSLR;
            
            if width_change == 0
                change_rate = 0;
                T = 1000/RSLR;
            else
                change_rate = width_change / T;
            end
            
            j
            final_width

            % Save the results
            results(ii,:) = [a init_width final_width width_change change_rate T];

            ii = ii + 1;
        end
end

% Write all of the results to a spreadsheet in the Modern runs folder
if i == 1
    xlswrite('C:/GEOMBEST+/Marsh widths.xls',results,'Sheet1','A2')
%elseif i == 2
%    xlswrite('C:/GEOMBEST+/Marsh widths.xls',results,'Sheet1','I2')
elseif i == 3
    xlswrite('C:/GEOMBEST+/Marsh widths.xls',results,'Sheet3','A2')
end

ri = 1;
real_results = zeros(1000,8);

% for ii = 1:1000
%     if results(ii,3)/results(ii,2) > 2.4 | results(ii,3)/results(ii,2) < 0.6
%         ii = ii + 1;
%     else
%         real_results(ri,:) = results(ii,:);
%         ri = ri + 1;
%         ii = ii + 1;
%     end
% end
% 
% real_results(ri:1000,:) = [];
% 
% length(real_results(:,1))
% 
% if i == 1
%     xlswrite('../Walters/2Modern runs/Realistic marsh widths3.xls',real_results,'Sheet1','A2')
% elseif i == 2
%     xlswrite('../Walters/2Modern runs/Realistic marsh widths3.xls',real_results,'Sheet1','I2')
% elseif i == 3
%     xlswrite('../Walters/2Modern runs/Realistic marsh widths3.xls',real_results,'Sheet1','Q2')
end