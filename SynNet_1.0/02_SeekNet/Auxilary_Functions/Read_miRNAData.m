% This reads in the input file, here's how an example file looks like:


% # The format is simple, here how you make comments!
% # The first line includes Unique IDs. This can be anything, text, number,
% whatever you feel like today, but it needs to be unique for each column. This unique ID is used
% as the Key to relate this data to the metadata provided in the Metadata file.

% Pejman 2014 June
%---------------------


function [miRNA_Data, MetaData] = Read_miRNAData(DataFile, MetaDataFile)
miRNA_Data = Pej_Read_Expression_Table(DataFile);
MetaData   = Pej_Read_MetaData(MetaDataFile);

%% Match IDs and check if they're unique
iID   = find(strcmpi(MetaData.Labels, 'UniqueID'));

if isempty(iID)
    error('Could not identify any column labeled as "UniqueID" in %s !', MetaDataFile)
end 

if length(iID)>1
    error('There are multiple rows columns labeled as "UniqueID" in %s !', MetaDataFile)
end

ID1 = miRNA_Data.SampleLabels;
ID2 = MetaData.Values(:,iID);
if length(unique(ID1)) < length(ID1)
        error('The IDs provided in in %s are not unique!', DataFile)
end

if length(unique(ID2)) < length(ID2)
        error('The IDs provided in in %s are not unique!', MetaDataFile)
end

[tmp, I1, I2] = intersect(ID1, ID2);

if length(tmp) < length(ID1)
        error('Some of the IDs provided in %s do not exist in the metadata file: %s', DataFile, MetaDataFile)
end

clear ID2 ID1
[~, tmpI] = sort(I1, 'ascend'); % This mimics the intersect with the option 'stable' 
MetaData.Values         = MetaData.Values(I2(tmpI),:);
miRNA_Data.SampleID     = miRNA_Data.SampleLabels';

%% Look for Annotations
Annotation_col = find(strcmpi(MetaData.Labels, 'AnnotationID'));
if isempty(Annotation_col)
    error('Could not identify any column labled as "AnnotationID" in %s !', MetaDataFile)
end

if length(Annotation_col)>1
    error('There are multiple columns labled as "AnnotationID" in %s !', MetaDataFile)
end
miRNA_Data.Annotation = str2double(MetaData.Values(:,Annotation_col));
if any(isnan(miRNA_Data.Annotation))
        error('The "AnnotationID" column in %s should only contain numeric values!', MetaDataFile)
end

%% Look for SampleNames
SampleName_col = find(strcmpi(MetaData.Labels, 'SampleName'));
if isempty(SampleName_col)
    warning('Could not identify any column labled as "SampleName" in %s, Annotation will be used as SampleName!', MetaDataFile)
    SampleName_col = Annotation_col;
end

if length(SampleName_col)>1
    warning('There are multiple columns labled as "SampleName" in %s, I use the fist one!', MetaDataFile)
    SampleName_col(2:end) = [];
end

miRNA_Data.SampleLabels = MetaData.Values(:,SampleName_col);

%%
MetaData.Values(:,[Annotation_col SampleName_col])=[];
MetaData.Labels(  [Annotation_col SampleName_col])=[];

end