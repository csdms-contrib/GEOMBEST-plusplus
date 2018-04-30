function [tempgrid,redepvol] = erodemarsh (tempgrid,t,j,tt,TP)
% erode marsh - erodes the edges of the marsh based on wave power

%global variables go here
global celldim;
global fetch;
global zcentroids;
global SL;
global dunelimit;
global L;
% global temprmarshedge;
% global templmarshedge;
global bay;

%calc amount of sediment to be eroded
%calculate average depth of the bay

if ((t == 1 && tt > 1) || t > 1) && length(bay) >= 3 %on first time step, no marsh
    
    rbay = bay(1); %first cell in bay on right
    lbay = bay(length(bay)); %last cell in bay on left
    rmarsh = rbay - 1; %edge of marsh
    lmarsh = lbay + 1; %edge of marsh
    
    for i=rbay:lbay %for each cell in the bay
        realratio = sum(tempgrid(i,:,:),3); %amount of sed in each cell of column i
        topcell = find(realratio == 1);
        topcell = topcell(1) - 1; %this is top cell with sediment
        k = topcell;
        realratio = sum(tempgrid(i,k,:),3); %amount of sediment in topcell
        H = (zcentroids(k)-celldim(3,j)./2)+ realratio*celldim(3,j); %height
        if isempty(H) == 1
            keyboard
        end
        depth(i+1-rbay) = SL(t,j)-H; %meters
    end
    avgdepth = sum(depth)/length(depth); %average depth of the bay
    
    highwater = SL(t,j) + gethighwater(t,j);
    
    %calculate wave power
    
    U = getwindspeed(t,j); %(m/s)
    g = 9.81; %(m/s^2)
    Hwave = (((U^2)*0.2413)*(tanh(0.493*((g*avgdepth)/(U^2))^0.75)*tanh((0.00313*((g*fetch)/(U^2))^0.57)/tanh(0.493*((g*avgdepth)/(U^2))^0.75)))^0.87)/g; % wave height,(meters)
    
    gamma = 9810; %specific weight of water (N/m^3 OR kgm/s^2 OR kg/s^2m^2) - are these units right?
    Wpower = (gamma/16)*Hwave^2*(g*avgdepth)^0.5; %what are units on wave power?
    
    kerode = geterosioncoeff(t,j); %varies per site (units m3/W/yr)
    %%%%%find surface elevation of rmarsh and first bay cell and find difference?%%%%%
    %     b = bay(1);
    %     realratio = sum(tempgrid(b,:,:),3); %amount of sed in each cell of column i
    %     topcell = find(realratio == 0,1,'last');
    %     topcell = topcell+1; %this is top cell with sediment
    %     k = topcell;
    %     realratio = sum(tempgrid(b,k,:),3); %amount of sediment in topcell
    %     H = (zcentroids(k)-celldim(3,j)./2)+ realratio*celldim(3,j); %height
    %     depthbay = highwater - H;
    %
    %     m = rmarsh;
    %     realratio = sum(tempgrid(m,:,:),3); %amount of sed in each cell of column i
    %     topcell = find(realratio == 0,1,'last');
    %     topcell = topcell+1; %this is top cell with sediment
    %     k = topcell;
    %     realratio = sum(tempgrid(m,k,:),3); %amount of sediment in topcell
    %     H = (zcentroids(k)-celldim(3,j)./2)+ realratio*celldim(3,j); %height
    %     depthmarsh = highwater - H;
    
    %marshh = depthbay - depthmarsh;
    %marshh = highwater - H; %height of marsh platform - this assumes marsh is at highwater
    
    erodevol = Wpower*kerode*celldim(1,j)*TP(t); % sediment that can be eroded based on wave power; m^3
    
    %erode the right marsh
    ii = rmarsh;
    
    %find depth of cell next to marsh boundary
    b = rbay; %bay cell adjacent to marsh boundary
    realratio = sum(tempgrid(b,:,:),3);
    topcell = find(realratio == 1);
    topcell = topcell(1) - 1; %partially full boundary cell
    k = topcell;
    cellfill = sum(tempgrid(b,k,:),3); %sed in boundary cell
    H = ((zcentroids(k) - celldim(3,j) ./ 2)) + cellfill*celldim(3,j);
    baydepth = highwater - H; %depth of adjacent bay cell
    
    tempevol = erodevol;
    fillsand = 0;
    fillkeep = 0;
    
    while tempevol > 0
        %find boundary cell
        realratio = sum(tempgrid(ii,:,:),3); %fraction of sed in each cell
        topcell = find(realratio == 1); %first full cell
        topcell = topcell(1)-1; %partially full boundary cell
        
        %find volume of boundary cell
        n = topcell;
        cellfill = sum(tempgrid(ii,n,:),3); %fraction of all sed in topcell
        fillmud = tempgrid(ii,n,3);
        fillmarsh = tempgrid(ii,n,2); %marsh is half organic half mud
        
        H = ((zcentroids(n) - celldim(3,j)./2)) + cellfill*celldim(3,j); %height to surface
        
        %targete = baydepth - (highwater-H); %amount that one column can be eroded
        targete = avgdepth - (highwater - H); %does it get rid of bumpy scarp if use avgdepth?
        
        while targete > 0
            if tempgrid(ii,n,4) == 1
                tempevol = 0;
            end
            %%bottom cell of each column cannot be eroded completely. Need
            %%to change below if statement so that if tempevol >
            %%fillmarsh+fill mud but not > targete then go ahead. Otherwise
            %%need to stop at targete = 0. Maybe figure out proportion of
            %%sed to do this?
            if tempevol > (fillmarsh+fillmud)*celldim(3,j)*celldim(1,j) %if there is enough erosion potential, empty the cell
                fillsand = fillsand + tempgrid(ii,n,1);
                tempgrid(ii,n,1:3) = 0; %empties cell by setting equal to zero
                fillkeep = fillkeep + (fillmud + (fillmarsh/2))*celldim(3,j)*celldim(1,j);
                tempevol = tempevol - (fillmud+fillmarsh)*celldim(3,j)*celldim(1,j);
                n = n + 1; %proceed to cell below
                fillmud = tempgrid(ii,n,3);
                fillmarsh = tempgrid(ii,n,2);
                H = (zcentroids(n) + celldim(3,j)./2); %update height of sed surface
                targete = baydepth - (highwater - H); %update amount left to erode
                %                 if targete < 0
                %                     keyboard
                %                 end
            else
                %subtract marsh first
                %then subtract bay
                %then deposit sand - right now this includes sand from this cell, make sure to account for that
                %tempgrid(ii,n,:) = tempgrid(ii,n,:) - tempevol/(celldim(3,j)*celldim(1,j)); % if not enough erosion potential to erode whole cell, erode as much as you can
                if tempevol > fillmarsh*celldim(3,j)*celldim(1,j)
                    tempgrid(ii,n,2) = 0;
                    tempevol = tempevol - fillmarsh*celldim(3,j)*celldim(1,j);
                    fillkeep = fillkeep + (fillmarsh/2)*celldim(3,j)*celldim(1,j);
                    fillmarsh = 0;
                    
                else
                    tempgrid(ii,n,2) = tempgrid(ii,n,2) - (tempevol/(celldim(3,j)*celldim(1,j)));
                    fillkeep = fillkeep + (tempevol/2);
                    tempevol = 0;
                end
                if tempevol > fillmud*celldim(3,j)*celldim(1,j)
                    tempgrid(ii,n,3) = 0;
                    fillkeep = fillkeep + fillmud*celldim(3,j)*celldim(1,j);
                    tempevol = tempevol - fillmud*celldim(3,j)*celldim(1,j);
                    fillmud = 0;
                else
                    tempgrid(ii,n,3) = tempgrid(ii,n,3) - (tempevol/(celldim(3,j)*celldim(1,j)));
                    fillkeep = fillkeep + tempevol;
                    tempevol = 0;
                    
                end
                
                fillsand = fillsand + tempgrid(ii,n,1);
                tempgrid(ii,n,1) = 0;
                break
                
            end
        end
        %deposit sand before moving over a column
        cellfillnew = sum(tempgrid(ii,n,:),3);
        while fillsand > 0
            if fillsand > 1 - cellfillnew
                tempgrid(ii,n,1) = 1 - cellfillnew;
                %                         if sum(tempgrid(ii,n,:),3) > sum(tempgrid(ii,n+1,:),3)
                %             keyboard
                %         end
                cellfillnew;
                fillsand;
                fillsand = fillsand - (1 - cellfillnew);
                n = n - 1;
                cellfillnew = 0;
            else
                tempgrid(ii,n,1) = fillsand;
                fillsand = 0;
            end
        end
        
        if ii > dunelimit
            ii = ii - 1; %continue to next column on right
            fillsand = 0;
        else
            break
        end
    end
    
    %erode the left marsh
    
    ii = lmarsh;
    
    %find depth of cell next to marsh boundary
    b = ii-1; %bay cell adjacent to marsh boundary; what about when bay is full?
    realratio = sum(tempgrid(b,:,:),3);
    %topcell = find(realratio == 0,1,'last');
    %topcell = topcell + 1; %partially full boundary cell
    topcell = find(realratio == 1);
    topcell = topcell(1) - 1;
    k = topcell;
    cellfill = sum(tempgrid(b,k,:),3); %sed in boundary cell
    H = (zcentroids(k)-celldim(3,j)./2)+cellfill*celldim(3,j);
    lbaydepth = highwater - H; %depth of adjacent bay cell
    
    tempevol = erodevol;
    fillsand = 0;
    
    while tempevol > 0
        %find boundary cell
        realratio = sum(tempgrid(ii,:,:),3); %fraction of sed in each cell
        topcell = find(realratio == 1); %first full cell
        topcell = topcell(1)-1; %partially full boundary cell
        
        %find volume of boundary cell
        n = topcell;
        cellfill = sum(tempgrid(ii,n,:),3); %fraction of sed in topcell
        %     fillsand = fillsand + tempgrid(ii,n,1);
        fillmud = tempgrid(ii,n,3);
        fillmarsh = tempgrid(ii,n,3);
        %     fillkeep = fillkeep + fillmud+(fillmarsh/2);
        
        H = ((zcentroids(n) - celldim(3,j)./2)) + cellfill*celldim(3,j); %height to surface
        
        targete = lbaydepth - (highwater-H); %amount that one column can be eroded
        while targete > 0
            if tempgrid(ii,n,4) == 1
                tempevol = 0;
            end
            if tempevol > (fillmarsh+fillmud)*celldim(3,j)*celldim(1,j) %if there is enough erosion potential, empty the cell -- this does not include sand in erosion potential
                fillsand = fillsand + tempgrid(ii,n,1);
                tempgrid(ii,n,1:3) = 0; %empties cell by setting equal to zero
                fillkeep = fillkeep + (fillmud + (fillmarsh/2))*celldim(3,j)*celldim(1,j);
                tempevol = tempevol - (fillmud+fillmarsh)*celldim(3,j)*celldim(1,j);
                n = n + 1; %proceed to cell below
                fillmud = tempgrid(ii,n,3);
                fillmarsh = tempgrid(ii,n,2);
                H = (zcentroids(n) + celldim(3,j)./2); %update height of sed surface
                targete = lbaydepth - (highwater - H); %update amount left to erode
            else
                %subtract marsh first
                %then subtract bay
                %then deposit sand - right now this includes sand from this cell, make sure to account for that
                %tempgrid(ii,n,:) = tempgrid(ii,n,:) - tempevol/(celldim(3,j)*celldim(1,j)); % if not enough erosion potential to erode whole cell, erode as much as you can
                if tempevol > fillmarsh*celldim(3,j)*celldim(1,j)
                    tempgrid(ii,n,2) = 0;
                    tempevol = tempevol - fillmarsh*celldim(3,j)*celldim(1,j);
                    fillkeep = fillkeep + (fillmarsh/2)*celldim(3,j)*celldim(1,j);
                    fillmarsh = 0;
                    
                else
                    tempgrid(ii,n,2) = tempgrid(ii,n,2) - (tempevol/(celldim(3,j)*celldim(1,j)));
                    fillkeep = fillkeep + (tempevol/2);
                    tempevol = 0;
                end
                if tempevol > fillmud*celldim(3,j)*celldim(1,j)
                    tempgrid(ii,n,3) = 0;
                    fillkeep = fillkeep + fillmud*celldim(3,j)*celldim(1,j);
                    tempevol = tempevol - fillmud*celldim(3,j)*celldim(1,j);
                    fillmud = 0;
                else
                    tempgrid(ii,n,3) = tempgrid(ii,n,3) - (tempevol/(celldim(3,j)*celldim(1,j)));
                    fillkeep = fillkeep + tempevol;
                    tempevol = 0;
                    
                end
                
                fillsand = fillsand + tempgrid(ii,n,1);
                tempgrid(ii,n,1) = 0;
                break
                
            end
        end
        cellfillnew = sum(tempgrid(ii,n,:),3);
        while fillsand > 0
            if fillsand > 1 - cellfillnew
                tempgrid(ii,n,1) = 1 - cellfillnew;
                fillsand = fillsand - (1 - cellfillnew);
                n = n - 1;
                cellfillnew = 0;
            else
                tempgrid(ii,n,1) = fillsand;
                fillsand = 0;
            end
        end
        if ii < L
            ii = ii + 1; %continue to next column on left
        else
            break
        end
    end
    
else
    erodevol = 0; %if t=1 no marsh to erode
end

if erodevol > 0
%     if t < 11
        redepvol = (fillkeep/2);
%     else
%         redepvol = 0.90*(fillkeep/2);
%     end
else
    redepvol = 0;
end
%amount of sed to use to grow marsh


