function DataAggregation(directory,times,Xpos,Ypos,cluster_IDs,px_size)


NumberOfClusters = max(cluster_IDs);

BurstDurations = zeros(1,NumberOfClusters);
BurstCounts = zeros(1,NumberOfClusters);
ClusterStartTime = zeros(1,NumberOfClusters);
ClusterCentroid = zeros(2,NumberOfClusters);  
ClusterRadius = zeros(1,NumberOfClusters);

XLSArray = cell(NumberOfClusters+1,6);
XLSArray{1,1} = 'Burst Durations';
XLSArray{1,2} = 'Burst Counts';
XLSArray{1,3} = 'Cluster Start Time';
XLSArray{1,4} = 'Cluster Centroid X (nm)';
XLSArray{1,5} = 'Cluster Centroid Y (nm)';
XLSArray{1,6} = 'Cluster RMS Radius (nm)';

for i = 1:NumberOfClusters
    Start = min(times(cluster_IDs==i));
    End = max(times(cluster_IDs==i));
    
    %I am evaluating the centroid in the most naive way.
    ClusterCentroid(:,i) = 1/sum(cluster_IDs==i)*[sum(Xpos(cluster_IDs==i))*px_size;sum(Ypos(cluster_IDs==i))*px_size]; 
    
    %RMS distance
    ClusterRadius(i) = sqrt(sum(sum(([(Xpos(cluster_IDs==i)*px_size-ClusterCentroid(1,i));Ypos(cluster_IDs==i)*px_size-ClusterCentroid(2,i)]).^2))/sum(cluster_IDs==i));
    
    ClusterStartTime(i) = Start;
    BurstDurations(i) = End-Start+1;
    BurstCounts(i) = sum(cluster_IDs==i);
    
    XLSArray{i+1,1} = BurstDurations(i);
    XLSArray{i+1,2} = BurstCounts(i);
    XLSArray{i+1,3} = ClusterStartTime(i);
    XLSArray{i+1,4} = ClusterCentroid(1,i);
    XLSArray{i+1,5} = ClusterCentroid(2,i);
    XLSArray{i+1,6} = ClusterRadius(i);
end

%% Write to an excel Spreadsheet

current_directory=cd;
cd(directory)

test_name = 'AggregatedData.csv';
count=1;
while exist(test_name,'file')==2
    count = count+1;
    test_name = ['AggregatedData_',num2str(count),'.csv'];
end

filehandle = fopen([directory,test_name],'w');
[a,b] = size(XLSArray);
for i = 1:a
    for j =1:b
        if isnumeric(XLSArray{i,j})
            fprintf(filehandle,num2str(XLSArray{i,j}));
        else
            fprintf(filehandle,XLSArray{i,j});
        end
        if j~=b
            fprintf(filehandle,',');
        end
    end
    fprintf(filehandle,'\n');
end

cd(current_directory)

%% Plot
figure
scatter(ClusterCentroid(1,:),ClusterCentroid(2,:),3,'filled')
figure
hist(ClusterRadius,100)
title('ClusterRadius')
figure
hist(BurstDurations,20)
title('BurstDurations')
figure
hist(BurstCounts,20)
title('BurstCounts')

end