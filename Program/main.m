function  main(filethread) 
%[equil surface] = main(filethread) 

% main -- main function of the Quicksand program
% Dr David Stolper dstolper@usgs.gov
% filethread specifies the input and output files used
% i.e. inpout1 to input4, output1 to  output4 

% Version of 26-Dec-2002 09:41
% Updated    30-Jun-2003 10:42

clear global;

global gridstrat; % 5-D raster representation of morphology/stratigraphy - grid(T,L,M,N,S)
global SL; % sea level relative to the original elevation SL(t)
global T; % # timestep
global M; % # grid cells along y axis 
global erosionresponse; % depth-dependant response rate (m/yr)
global accretionresponse; % depth-dependant response rate (m/yr)
global zcentroids; % vector of cell centroids along the z axis
global residual; % residual sediment volumes from the crossshore search algorithm
global zresponseerosion; % maximum shoreface response rate at the elevations specified by zcentroids
global zresponseaccretion; % maximum shoreface response rate at the elevations specified by zcentroids
global stormcount;

loadrunfiles(filethread); 
loadstrat(filethread);   
loadrate(filethread);
buildgrid;
zresponseerosion = interp1(erosionresponse(:,1),erosionresponse(:,2),zcentroids(:)); % calculate shoreface response rate at 
% elevations corresponding to the z centroids of each cell within the grid 
zresponseaccretion = interp1(accretionresponse(:,1),accretionresponse(:,2),zcentroids(:)); % calculate shoreface response rate at 
% elevations corresponding to the z centroids of each cell within the grid

residual = zeros(M); % residual sediment volume from each time step 

for j = 1:M
    findfirstshoreline(0,j);
%     findmarshboundary(0,j);
end

tic
savetimestep(1,1,filethread);

storms = zeros(T,5);
for t = 1:T % loop through time
   
    t % Optional - print timestep to screen

    for j = 1:M % loop through tracts  
        gridstrat(2,:,j,:,:) = gridstrat(1,:,j,:,:); 
        solvecross(t,j); % balance crossshore sediment budget
    end
    savetimestep(2,t + 1,filethread); 
    gridstrat(1,:,:,:,:) = gridstrat(2,:,:,:,:);
    zresponseerosion = interp1(erosionresponse(:,1),erosionresponse(:,2),zcentroids(:) - SL(t)); % update shoreface response rate 
    zresponseaccretion = interp1(accretionresponse(:,1),accretionresponse(:,2),zcentroids(:) - SL(t)); % update shoreface response rate
    
    storms(t,:) = stormcount;    
%     if t == 13 
%         
%         block = 0
%         savedata(filethread) %added 5_20_08 LJM  If using debugger to stop 
%         %the program at a desired time step, step forward to execute
%         %the savedata command and then plotting functions can still be
%         %used.
% 
%     end    
end 

savedata(filethread)

toc