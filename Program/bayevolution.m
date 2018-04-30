function [tempgrid,marshvol] = bayevolution(icrest,tempgrid,t,j,TP)

global celldim;
global SL;
global zcentroids;
global L;
global maxdepth;
global fetch;
global bay;

%%% ACCRETION
% Accrete the estuary
[q_b] = getbaysedflux(t,j); % Volume of sediment added to the estuary from riverine input
mwl = SL(t,j); % Mean Water Level - the minimum depth at which the marsh will grow
tempvol = q_b * TP(t); 
baycell = 0; % Initiate baycell to count the number of bay cells to be accreted onto
a = 1;
bay = 0;

for ii = icrest : L
    %Find all topcells underneath Sea Level
    realratio = sum(tempgrid(ii,:,:),3); % the fraction of sediment in each cell of the column at the limit of the dune
    topcell = find(realratio == 1); % find the first full cell
    topcell = topcell(1) - 1; % the boundary cell that is partially filled or empty

    k=topcell;
%     if isempty(k) == 1
%         keyboard
%     end
%         if min(min(min(tempgrid)))<0
%     keyboard
% end    
    cellfill = sum(tempgrid(ii,k,:),3); % fraction of sediment in the topcell
    H =(zcentroids(k) - celldim(3,j)./ 2) + cellfill*celldim(3,j); % height to the topcell
    ii;

        
    if H <= mwl
        baycell = baycell + 1; % Count the number of bay cells below sea level
        bay(a) = ii;

        a = a+1;
    end  
    ii = ii + 1;
end

% Calculation of accretion height

targvol = tempvol / baycell; % Distribute the riverine input of sediment over each baycell
estrate = targvol / celldim(1,j); % Height to which each baycell can be accreted


%%% EROSION
% Erode the estuary (everything below MSL)
%previous erosion method commented below
%edited 6/26/15 Rebecca Lauzon

actualerosion = 0; %initialize actual erosion


for ii=icrest :L
    % Find the top horizon of the bay
    realratio = sum(tempgrid(ii,:,:),3); % The amount of sediment in a cell
    %topcell = find(realratio == 0, 1, 'last' ); 
    topcell = max(find(realratio == 0)); %old version
    topcell = topcell + 1; % The first cell with sediment
    
    k = topcell;
    realratio = sum(tempgrid(ii,k,:),3); % The amount of sediment in the topcell
    
    %Calculate depth
    H =(zcentroids(k) - celldim(3,j) ./ 2) + realratio*celldim(3,j);
    depth = SL(t,j)-H; %(m)
    
   
    %Calculate fetch
    fetch = baycell*celldim(1,j); %number of baycells below sea level * x cell dimension (m)
    
    %calculate H and T
    g = 9.81; %gravity (m/s^2)
    U = getwindspeed(t,j); %wind speed (m/s)
    Hwave = (((U^2)*0.2413)*(tanh(0.493*((g*depth)/(U^2)).^0.75).*tanh((0.00313*((g*fetch)/(U^2)).^0.57)./tanh(0.493*((g*depth)/(U^2)).^0.75))).^0.87)./g; %(m)
   
    Twave = ((7.518*U)*(tanh(0.331*((g*depth)/(U^2)).^1.01).*tanh((0.0005215*((g*fetch)/(U^2)).^0.73)./tanh(0.331*((g*depth)/(U^2)).^1.01))).^0.37)./g;  %(s)
    
    %calculate shear stress
    friction = 0.03 ;   
    density = 1000; %(kg/m^3) ranges from 997 to 1002 from fresh to salt water
    shearstress=(9.94519296*10^14).*0.5.*friction.*density.*(pi*Hwave./(Twave.*sinh((2.*pi.*depth)./(Twave.*(g*depth).^0.5)))).^2; %(Pa)(s/yr)^2 = (kg/(m*yr^2))
    
    %%%%%%if sand is eroded it is not conserved! only fine sed conserved
    %find uppermost sand layer
    ii;
    realsandratio = tempgrid(ii,:,1); %amount of sand in a cell
    %topsandcell = find(realsandratio ~= 0, 1, 'last');
    topsandcell = max(find(realsandratio ~= 0)); %old version
    topsandcell = topsandcell - 1;
    
    ksand = topsandcell;
    
if ksand ~= 0
    realsandratio = tempgrid(ii,ksand,1); %amount of sand in topcell
    
    Hsand = ((zcentroids(ksand) - celldim(3,j)) ./ 2) + realsandratio*celldim(3,j);
    depthsand = SL(t,j) - Hsand; %depth to surface of sand layer
end
    
    %calculate erosion
    A = (4.12*10^-4)/31536000; %erosion coeff (s/m)(yr/s) = (yr/m)
    %if shear stress is greater than critical shear stress, erosion occurs
    shearcrittemp = getshearcrit(t,j);
    shearcrit = shearcrittemp*(9.94519296*10^14); %(Pa) = (kg/m/s^2)(s/yr)^2  = (kg/m/yr^2)
    seddensity = 2600; %(kg/m^3)
    if shearstress > shearcrit
        toperosion = A*(shearstress - shearcrit)*TP(t)/seddensity; %(kg/m^2/yr)(yr)(m^3/kg) = m
        if isempty(toperosion)== 0 & ksand ~= 0
            if toperosion > depthsand - depth
                toperosion = depthsand - depth;
            end
        end
         
    %if shear stress is less than or equal to critical shear stress, no erosion
    elseif shearstress <= shearcrit
        toperosion = 0;
    end
    
    if depth <= 0 %
        toperosion = 0;
    end
    if exist('toperosion') == 0
        toperosion = 0
    end
    if depth + toperosion > maxdepth  % If erosion over the time step is too large, set depth equal to resuspension depth
        toperosion = maxdepth - depth;
    end
          
    temperosion = toperosion;
    ii;
    
    %Calculate height change
    dheight = estrate - toperosion; %accretion - erosion
    if dheight > 0
        tempgrid = deposition(tempgrid,t,ii,j,dheight);
%         if min(min(min(tempgrid)))<0
%     keyboard
% end
    end
    
    if depth >= 0 & shearstress > shearcrit
        [actualerosion,tempgrid] = erosion(actualerosion,tempgrid,t,ii,j,k,dheight,temperosion); 
%         if min(min(min(tempgrid)))<0
%     keyboard
% end%do i
    end %%need an if statement for this?
    
    
end

% if min(min(min(tempgrid)))<0
%     keyboard
% end

marshvol = actualerosion / 2;
if baycell <= 2
    marshvol = tempvol/2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%previous erosion calculation method%%%%%%
% [bayerosion] = getbayerosion(t,j); % Maximum rate erosion for the estuary
% tE = bayerosion*TP(t)/1000; % Maximum depth to which a surface can be eroded during a single time step
% actualerosion = 0; % Initialize actual erosion, for the actual volume of material that is eroded
% [resuspensionelevation] = getresuspensionelevation(t,j); % The depth below which the esuary can no longer be eroded
% 
% for ii = icrest :L
%     % Find the top horizon of the estuary
%     realratio = sum(tempgrid(ii,:,:),3); % The amount of sediment in a cell
%     topcell = max(find(realratio == 0)); 
%     topcell = topcell + 1; % The first cell with sediment
% 
%     k = topcell;
%     realratio = sum(tempgrid(ii,k,:),3); % The amount of sediment in the topcell
%     
%     H =(zcentroids(k) - celldim(3,j) ./ 2) + realratio*celldim(3,j);
%     depth = SL(t,j)-H; % Depth to the surface of the estuary
% 
%     % Calculation of erosion height
%     toperosion = tE * (resuspensionelevation - depth)/resuspensionelevation; % Amount of erosion weighted by depth
%     
%     % Find the uppermost layer of pure sand
%     realsandratio = tempgrid(ii,:,1); % The amount of sand in a cell
%     topsandcell = max(find(realsandratio ~= 0));
%     topsandcell = topsandcell - 1;
%       
%     ksand = topsandcell;
%     realsandratio = tempgrid(ii,ksand,1); % The amount of sand in the topcell
%     
%     Hsand =((zcentroids(ksand) - celldim(3,j)) ./ 2) + realsandratio*celldim(3,j);
%     depthsand = SL(t,j)-Hsand; % Depth to the surface of the sand horizon
% 
%     if toperosion > depthsand - depth
%         toperosion = depthsand - depth;
%     end
%     
% 
%     if depth > resuspensionelevation % If the the estuary is already deeper than resuspension depth, there is no erosion
%         toperosion = 0;
%     elseif depth + toperosion > resuspensionelevation  % If erosion over the time step is too large, set depth equal to resuspention depth
%         toperosion = resuspensionelevation - depth;
%     end
% 
%     temperosion = toperosion;
%     
%     % Calculation of effective height change
%     
%     dheight = estrate - toperosion;
%     if dheight > 0
%         tempgrid = deposition(tempgrid,t,ii,j,dheight);
%     end
%     
%     if depth >= 0 & depth < resuspensionelevation
%         [actualerosion,tempgrid] = erosion(actualerosion,tempgrid,t,ii,j,k,dheight,temperosion);
%     end
% end
% 
% marshvol = actualerosion / 2;
% 
% if baycell <= 2
%     marshvol = tempvol/2;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%