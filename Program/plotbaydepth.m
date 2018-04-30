function plotbaydepth(a,b,c) %(low wind speed, high wind speed, #time steps)

d = c + 2;
results = zeros(1000,d);
ii = 1;

time = [1:c]*10;

for U = a:b
filethread = U;
filename = ['C:/GEOMBEST+/Output' num2str(filethread) '/surface.mat'];
load(filename)

filename2 = ['C:/GEOMBEST+/Output' num2str(filethread) '/xcentroids.mat'];
load(filename2)

filename3 = ['C:/GEOMBEST+/Output' num2str(filethread) '/SL.mat'];
load(filename3)

avg_bay_depth = zeros(1,c);

for t = [1:c]
x = xcentroids(1,:); % x location of each cell
y = surface(t,:); %elevation of surface at timestep
numcell = numel(y); %number of surface cells

for n = 1:numcell
    if y(n) >= SL(t)
        x_island = n; %first cell above sea level is island edge
        break
    else
        n = n + 1;
    end
end

x_rbay = numel(x);

for n = x_island:numcell
    if y(n) < SL(t)
        x_rbay = n; %x location of first estuarine cell (below SL)
        break
    else
        n = n + 1;
    end
end

x_lbay = numel(x);

for n = x_rbay:numcell
    if y(n) >= SL(t) %other edge of bay
        x_lbay = n - 1;
        break
    else
        n = n + 1;
    end
end

if x_rbay == numel(x)
    bay_width = 0; % no cells below sea level
else
    bay_width = x(x_lbay) - x(x_rbay)
end


%%average bay depth

avg_bay_depth(t) = sum(surface(x_rbay:x_lbay))/bay_width;

surface(90)

%save the results
results(ii,:) = [U bay_width avg_bay_depth];

end

fh=figure;
plot(time,avg_bay_depth)
xlabel('Time (years)','FontSize',15);
ylabel('Average bay depth (m)','FontSize',15);

saveas(fh,['C:/GEOMBEST+/Output' num2str(filethread) '/plot bay depth' num2str(filethread) '.fig']); 

outputfilename = ['C:/GEOMBEST+/Output' num2str(filethread) '/plot bay depth num2str(filethread)']; 
print('-djpeg',outputfilename)

ii = ii + 1;

end

%write results to spreadsheet
xlswrite('C:/GEOMBEST+/Bay widths.xls', results, 'Sheet1','A2')

    
end