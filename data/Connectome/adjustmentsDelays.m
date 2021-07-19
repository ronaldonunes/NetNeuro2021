
% Load Table
filename = 'mouse_sbri_distance_matrix_1.0.csv';

% Table with distances
A = readtable(filename);

% Areas sorted according to Kennedy's paper
AreaList={'GU','SSp-un','VISC','MOp','PL','SSs','SSp-bfd','ACAd',...
    'RL','AL','DP','AUDpo','AM','V1','MM','LM','PM','RSPd','P'};

% List of index 
index=zeros(1,19);

% Index of areas in table
AreasInTable=table2cell(A(:,1));

for i=1:size(AreaList,2) % To avoid problem when I have just one letter (ex: P, PM) 
    temp1=find(strcmp(AreasInTable,AreaList(i)));
    index(i)=temp1;
end

% Convert table to matrix
distanceMatrix=table2array(A(:,2:end));
% Select index 
distanceMatrix=distanceMatrix(index,index);



