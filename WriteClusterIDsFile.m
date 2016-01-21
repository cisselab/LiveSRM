function WriteClusterIDsFile(FileName,ClusterIDs)
    filehandle = fopen(FileName,'w');
    fprintf(filehandle,'ClusterIDs');
    fprintf(filehandle,'\n');
    for i = 1:length(ClusterIDs)
        fprintf(filehandle,num2str(ClusterIDs(i)));
        fprintf(filehandle,'\n');
    end