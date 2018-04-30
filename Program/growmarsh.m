function [tempgrid] = growmarsh (tempgrid,marshvol,t,j,tt,redepvol)


% Growmarsh -- fills the left and right-most shallow cells in the bay
% with marsh. Above MSL, organic growth helps the marsh reach its
% equilibirum elevation at mean high water.

% David Walters

% Version of 26-Jun-2012
%edited 10/26/15 - RL

global marshboundary;
global lmarshboundary;
global marshedge;
global celldim;
global SL;
global zcentroids;
global dunelimit;
global L;
global limits;
global maxdepth;
global temprmarshedge;
global templmarshedge;

%Fill the right marsh

[highwater] = SL(t,j) + gethighwater(t,j);
mwl = SL(t,j); % Mean Water Level - the minimum depth at which the marsh will grow
potmarshht = gethighwater(t,j) + maxdepth; % The depth below which the estuary can no longer be eroded

k = 0;

ii = dunelimit;
limits(t,tt) = dunelimit;
dunelimit;
L;

pc = 0.97; %percent sed conserved

marshvol;
redepvol;

% if t > 10
% redepvol = 0.60*redepvol;
% end
redepvol;

tempvol = marshvol+redepvol; % the volume of fine grain sediments eroded from the bay available to grow the marsh
if t > 10
    tempvol = pc*tempvol;
end

while tempvol > 0
    realratio = sum(tempgrid(ii,:,:),3); % the fraction of sediment in each cell of the column at the limit of the dune
    topcell = find(realratio == 1); % find the first full cell
    topcell = topcell(1) - 1; % the boundary cell that is partially filled or empty

    k=topcell;
    cellfill = sum(tempgrid(ii,k,:),3); % fraction of sediment in the topcell
    H =((zcentroids(k) - celldim(3,j) ./ 2)) + cellfill*celldim(3,j); % height to the topcell
    
    targetf = highwater-H; % the height of empty cells that can be filled by marsh/bay
    
    while targetf > 0

        subtarg1 = mwl - H; % the height of empty cells that can be filled with bay
        
        % First fill the lower part with %100 fine grain sediments
        
        while subtarg1 > 0

            s = 3; % fill using bay sediment
            if tempvol > (1 - cellfill) * celldim(3,j) * celldim(1,j) % if there is enough available sediment, fill the topcell
                cellfill = sum(tempgrid(ii,k,:),3);
                tempgrid(ii,k,s) = tempgrid(ii,k,s) + (1 - cellfill);
%                                 if max(max(max(tempgrid))) > 1.01
%                     keyboard
%                 end
%                 if min(min(min(tempgrid)))<0
%     keyboard
% end
                k = k - 1; % Proceed to the cell above
                H =((zcentroids(k) - celldim(3,j) ./ 2)); % Update the height of the sediment surface
                subtarg1 = mwl-H; % Update the subtarget
                targetf = highwater-H; % update the height left to fill with marsh

                tempvol = tempvol - (1 - cellfill) * celldim(3,j) * celldim(1,j); % subtract the volume added to the topcell from the available sediment
%                 if tempvol < 0
%                 keyboard
%             end
                cellfill = 0;
                
            else
               
                tempgrid(ii,k,s) = tempgrid(ii,k,s) + tempvol /(celldim(3,j) * celldim(1,j)); % if not enough sediment available, fill the topcell with as much as you have
%                                if max(max(max(tempgrid))) > 1.01
%                     keyboard
%                 end
%                 if min(min(min(tempgrid)))<0
%     keyboard
% end
                tempvol = 0; % there is no available sediment left
                break;
                
            end
        end        
        
        % Fill the marsh above Mean Sea Level using 50% organic and 50% fine grained sediments
                
         s = 2; % fill using marsh sediment
            if tempvol > 0.5*(1 - cellfill) * celldim(3,j) * celldim(1,j) % if there is enough inorganic sediment, fill the topcell 
                cellfill = sum(tempgrid(ii,k,:),3);
                tempgrid(ii,k,s) = tempgrid(ii,k,s) + (1 - cellfill);  
%                                 if max(max(max(tempgrid))) > 1.01
%                     keyboard
%                 end
%                 if min(min(min(tempgrid)))<0
%     keyboard
% end
                k = k - 1;
                H =((zcentroids(k) - celldim(3,j) ./ 2)); 
                targetf = highwater-H; % Update the height of empty cells that can be filled with marsh
         
                tempvol = tempvol - 0.5 * (1 - cellfill) * celldim(3,j) * celldim(1,j); % Update the volume of available sediment
%                 if tempvol < 0
%                 keyboard
%             end
                cellfill = 0; % the first marsh cell is filled
                
            else
               
                tempgrid(ii,k,s) = tempgrid(ii,k,s) + 2 * tempvol /(celldim(3,j) * celldim(1,j)); % if there isn't enough available sediment, fill with as much as you have
% if min(min(min(tempgrid)))<0
%     keyboard
% end
%                 if max(max(max(tempgrid))) > 1.01
%                     keyboard
%                 end
tempvol = 0;    % exit big loop
                break;          % exit column loop
                
            end

        
    end

    if ii < L
        ii = ii + 1; % Continue to the next column to the left
    else
       ii=ii;
        break
    end
end

if marshvol > 0
    temprmarshedge = ii;
else
    temprmarshedge = temprmarshedge;
end

% Save the position of the marsh boundary
if k > 0
    rightmarshedge = ii*celldim(1,j) + ((potmarshht - k*celldim(3,j))/potmarshht)*celldim(1,j);
else
    if t > 1
        rightmarshedge = marshboundary(t-1);
    else
        rightmarshedge = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Fill the left marsh

ii = L;
limits(t,tt) = L;
tempvol = marshvol+redepvol; % the volume of fine grain sediments eroded from the bay are available to grow the marsh

if t > 10
    tempvol = pc*tempvol;
end

while tempvol > 0

    realratio = sum(tempgrid(ii,:,:),3); % the fraction of sediment in each cell of the column at the limit of the dune
    topcell = find(realratio == 1); % find the first full cell
    topcell = topcell(1) - 1; % the boundary cell that is partially filled or empty

    k=topcell;
    cellfill = sum(tempgrid(ii,k,:),3); % fraction of sediment in the topcell
    H =((zcentroids(k) - celldim(3,j) ./ 2)) + cellfill*celldim(3,j);  % height to the topcell
    
    targetf = highwater-H; % the height of empty cells that can be filled by marsh/bay
    
    while targetf > 0

        subtarg1 = mwl - H; % the height of empty cells that can be filled with bay
        
        % First fill the lower part with %100 fine grain sediments
        
        while subtarg1 > 0

            s = 3; % fill using bay sediment
            if tempvol > (1 - cellfill) * celldim(3,j) * celldim(1,j) % if there is enough available sediment, fill the topcell
                
                tempgrid(ii,k,s) = tempgrid(ii,k,s) + (1 - cellfill); 
%                             if min(min(min(tempgrid)))< 0
%                 keyboard
%                             end
%                             if max(max(max(tempgrid))) > 1.01
%                     keyboard
%                 end
                k = k - 1;
                H =((zcentroids(k) - celldim(3,j) ./ 2));
                subtarg1 = mwl-H;
                targetf = highwater-H; % update the height left to fill with marsh

                tempvol = tempvol - (1 - cellfill) * celldim(3,j) * celldim(1,j); % subtract the volume added to the topcell from the available sediment
%                 if tempvol < 0
%                 keyboard
%             end
                cellfill = 0; % the topcell is filled to the top
                
            else
               
                tempgrid(ii,k,s) = tempgrid(ii,k,s) + tempvol /(celldim(3,j) * celldim(1,j)); % if not enough sediment available, fill the topcell with as much as you have
%                             if min(min(min(tempgrid)))< 0
%                 keyboard
%                             end
%                             if max(max(max(tempgrid))) > 1.01
%                     keyboard
%                 end
                tempvol = 0; % the topcell is filled to the top
                break;
                
            end
        end

        % Fill the marsh above Mean Sea Level using 50% organic and 50% fine grained sediments
                
         s = 2; % fill using marsh sediment
            %if tempvol > (1 - cellfill) * 0.5*(1 - cellfill) * celldim(3,j) * celldim(1,j)
            if tempvol > 0.5*(1-cellfill)*celldim(3,j)*celldim(1,j)
                
                tempgrid(ii,k,s) = tempgrid(ii,k,s) + (1 - cellfill) ;              
%                             if min(min(min(tempgrid)))< 0
%                 keyboard
%                             end
%                             if max(max(max(tempgrid))) > 1.01
%                     keyboard
%                 end
                k = k - 1;
                H =((zcentroids(k) - celldim(3,j) ./ 2));
                targetf = highwater-H; % Update the height of empty cells that can be filled with marsh

                tempvol = tempvol - 0.5 * (1 - cellfill) * celldim(3,j) * celldim(1,j); % Update the volume of available sediment
%                 if tempvol < 0
%                 keyboard
%             end
                cellfill = 0;  % the first marsh cell is filled
                
            else
               
                tempgrid(ii,k,s) = tempgrid(ii,k,s) + 2 * tempvol /(celldim(3,j) * celldim(1,j));
%                             if min(min(min(tempgrid)))< 0
%                 keyboard
%                             end
%                             if max(max(max(tempgrid))) > 1.01
%                     keyboard
%                 end
                tempvol = 0;    % exit big loop
                break;          % exit column loop
                
            end

        
    end
    
    if ii > dunelimit

    ii = ii - 1; % Continue with the next column to the right
    
    else
        ii = ii;
        break
    end

end

if marshvol > 0
    templmarshedge = ii;
else
    templmarshedge = templmarshedge;
end

% Save the position of the marsh boundary
if k > 0
    leftmarshedge = ii*celldim(1,j) + ((potmarshht - k*celldim(3,j))/potmarshht)*celldim(1,j);
else
    if t > 1
        leftmarshedge = lmarshboundary(t-1);
    else
        leftmarshedge = 0;
    end
end
marshedge = [leftmarshedge rightmarshedge];