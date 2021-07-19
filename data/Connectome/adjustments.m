% Salvar os arrays dos dados do Kennedy e do Allen institute, com o mesmo nome
% Porém o nome do arquivo deve ser diferente

% Verificar se esta igual ao artigo do Kennedy 


% Load Table
filename = '1-s2.0-S0896627317311856-mmc2.xlsx';
A = readtable(filename);

% SourceTargetFLn
SourceTargetFln=table2cell(A(:,[7,8,11]));

% Areas sorted according to Kennedy's paper
AreaList={'GU','SSp-un','VISC','MOp','PL','SSs','SSp-bfd','ACAd',...
    'RL','AL','DP','AUDpo','AM','V1','MM','LM','PM','RSPd','P'};

% Change Area Name to Area Number (1:19)
for i=1:size(AreaList,2) % To avoid problem when I have just one letter (ex: P, PM) 
    % Source
    temp1=find(strcmp(SourceTargetFln(:,1),AreaList(i)));
    SourceTargetFln(temp1,1)={num2str(i)};
    % Target
    temp2=find(strcmp(SourceTargetFln(:,2),AreaList(i)));
    SourceTargetFln(temp2,2)={num2str(i)};
    
end

% Convert String Cell to Double Array
SourceTargetFln=str2double(SourceTargetFln);

% Remove NaN (Areas that not belong to AreaList)
SourceTargetFln(isnan(SourceTargetFln(:,1)), :)=[];

% Create matrix of Weights 
Fln=zeros(size(AreaList,1),size(AreaList,1));

% Source
for i=1:size(AreaList,2)%1:size(AreaList,1)
    % Target
    for j=1:size(AreaList,2)
        if(i~=j)
            temp=SourceTargetFln(SourceTargetFln(:,1)==i,:);
            w=mean(temp(temp(:,2)==j,3));
            Fln(i,j)=w;
        end
    end
end
