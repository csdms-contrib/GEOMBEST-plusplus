function [tempgrid,rightmarshedge,leftmarshedge] = erodemarsh (tempgrid,t,j,tt,rightmarshedge,leftmarshedge)
%Erodemarsh - erodes the left and right marsh boundaries due to wind waves

%2/23/2015 - RL

global dunelimit;
global L;
global celldim;
global fetch;
global zcentroids;
global SL;
global marshboundary;
global lmarshboundary;

U = getwindspeed(t,j); %(m/s)
g = 9.81; %(m/s^2)

%calculate volume of sediment to be eroded

%calc Hwave using avg bay depth - make depth vector to fill, then find
%sum/baycell? make baycell global?

for i=marshboundary:lmarshboundary %for each cell in the bay - should these be -1?
realratio = sum(tempgrid(i,:,:),3); %amount of sediment in each cell of column ii
topcell = find(realratio == 0,1,'last');
topcell = topcell+1; %this is the first cell with sediment
k = topcell;
realratio = sum(tempgrid(i,k,:),3); %the amount of sediment in the topcell
H = (zcentroids(k)-celldim(3,j)./2)+ realratio*celldim(3,j);
depth=SL(t,j)-H; %m
end
length(depth) %check if depth is making a vector
depth = sum(depth)/baycell; %average bay depth

Hwave = (((U^2)*0.2413)*(tanh(0.493*((g*depth)/(U^2)).^0.75).*tanh((0.00313*((g*fetch)/(U^2)).^0.57)./tanh(0.493*((g*depth)/(U^2)).^0.75))).^0.87)./g; %(m)

gamma = 9810; %specific weight of water - N/m^3 or (kgm/s^2)/m^3 or kg/s^2m^2
Wpower = (gamma/16)*Hwave.^2*(g*depth)^0.5; %units?


kerode = 0.14; %varies per site - make input parameter - think this is how units will need to be fixed
marshh = 1; %not a real value - assume MHW? average depth? what about half-full columns?
erodevol = Wpower*kerode/marshh; %sediment that can be eroded based on wave power - units = m/yr - fix - how did I do this for bay evolution?

%%erode the right marsh

ii = marshboundary;
limits(t,tt) = marshboundary;

tempevol = erodevol;

while tempevol > 0
%find boundary cell
realratio = sum(tempgrid(ii,:,:),3); %fraction of sediment in each cell of column at limit of marsh
topcell = find(realratio == 1); %find first full cell - difference between this and realratio == 0,1,'last')
topcell = topcell(1) - 1; %partially full or empty boundary cell

%find volume of boundary cell
%if = 0, then use full column??
n = topcell;
cellfill = sum(tempgrid(ii,n,:),3); %fraction of sed in topcell
H =((zcentroids(n) - celldim(3,j) ./ 2)) + cellfill*celldim(3,j); % height to the topcell (from where to where?)

targete = H - ? %where to stop eroding?; %is this included in where height calculation stops? %height of full cells that can be eroded
    
while targete > 0 %erode that volume
            if tempevol > (1 - cellfill) * celldim(3,j) * celldim(1,j) % if there is enough available sediment, fill the topcell - need to change something here to make it full part instead of empty 
               
                tempgrid(ii,k,3) = tempgrid(ii,k,3) + (1 - cellfill); %this is how the cell is filled - make it into erosion
               

                k = k + 1; % Proceed to the cell below
                H =((zcentroids(k) - celldim(3,j) ./ 2)); % Update the height of the sediment surface - sign need to be changed?
                targete =H-?; % update the height left to erode
                tempevol = tempevol - (1 - cellfill) * celldim(3,j) * celldim(1,j); % subtract the volume added to the topcell from the available sediment %how to subtract volume
                %eroded?
                cellfill = 0; % the topcell is filled to the top - since we are eroding down...?
                
            else %if not enough to erode completely, erode as much as you can
               
                %tempgrid(ii,k,3) = tempgrid(ii,k,3) + tempevol /(celldim(3,j) * celldim(1,j)); % if not enough sediment available, fill the topcell with as much as you have
                %tempevol = 0; % there is no available sediment left
                break;
                
            end
end

%if more move on to next cell - tell it when to stop eroding down, and
%start eroding back (use adjacent bay cell - when it has eroded down to
%that depth, stop and move on.)
if ii > dunelimit
    ii = ii - 1; %continue to next cell to the right unless dune limit reached
else
    ii = ii;
   break


end
end
%save marsh boundary location
if n > 0
    rightmarshedge = ii*celldim(1,j) + ((potmarshht - n*celldim(3,j))/potmarshht)*celldim(1,j); %fix this....potmarshht needs to be defined
else
    if t > 1
        rightmarshedge = marshboundary(t-1);
    else
        rightmarshedge = 0;
    end
end
%see grow marsh - k = 0 is bay surface?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Erode the left marsh

%see grow marsh - need this stuff?

ii = lmarshboundary;

tempevol = erodevol; 

while tempevol > 0
%find boundary cell
realratio = sum(tempgrid(II,:,:),3); %fraction of sediment in each cell of column at limit of basin
topcell = find(realratio == 1); %first full cell
topcell = topcell(1) - 1; %partially full boundary cell

%find volume of boundary cell
n = topcell;
cellfill = sum(tempgrid(ii,n,:),3); %fraction of sediment in topcell
%here find difference between height and bay surface in next column over?

targete = 1;%height of full cells that can be eroded

while targete > 0 %erode that volume
%if not enough to erode completely?
%erode topcell

end
%if more move on to next cell down - tell it when to stop eroding down, and
%start eroding back (depth of next cell on bay side)
if ii < L
    ii = ii + 1; %continue to next cell to the left unless mainland reached
else
    ii = L; 
    break
end
end

%save marsh boundary location
if n > 0
    leftmarshedge = ii*celldim(1,j) + ((potmarshht - n*celldim(3,j))/potmarshht)*celldim(1,j);
else
    if t > 1
        leftmarshedge = lmarshboundary(t-1);
    else
        leftmarshedge = 0;
    end
end


