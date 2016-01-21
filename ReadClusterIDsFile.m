function ClusterIDs=ReadClusterIDsFile(FileName)
    ClusterIDs = csvread(FileName,1,0);
    ClusterIDs=uint16(ClusterIDs');