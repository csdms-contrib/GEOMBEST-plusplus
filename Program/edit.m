function edit(a,b)

% Folder = 'C:\GEOMBEST+';
% filePattern = fullfile(Folder, 'Input*');
% files = dir(filePattern);
for i = a:b
filename = ['C:\GEOMBEST+\Input' num2str(i) '\run1.xls'];

run1=xlsread(filename);

% for k=a:b
%     f=sprintf('C:\GEOMBEST+\Input' num2str(k) '\run1.xls',k,k)  
%    run1 = xlsread(f)
% end

c = size(run1);
n = c(1);

run1(:,4) = [];
run1(:,6) = [];
run1(3:n,9) = 8;
run1(3:n,10) = 0.2;
run1(3:n,11) = 0.14;
A = run1(3:n,7)./10;
run1(3:n,7) = A;
run1(1:2,:) = [];

col_header={'Time','Sea-level change (m)','Riverine Input(m^3)','Backbarrier Width','Exogenous sand volume (m3/m tract width/yr)','Overwash rate (mm/yr)','Overwash volume','High water','wind speed','shear crit','erosion coeff'};
xlswrite(['C:\GEOMBEST+\Input' num2str(i) '\run1.xls'],run1,'Sheet1','A3');
xlswrite(['C:\GEOMBEST+\Input' num2str(i) '\run1.xls'],col_header,'Sheet1','A2');


end
end

