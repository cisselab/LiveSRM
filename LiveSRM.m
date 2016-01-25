function varargout = LiveSRM(varargin)
% LIVESRM MATLAB code for LiveSRM.fig
%      LIVESRM, by itself, creates a new LIVESRM or raises the existing
%      singleton*.
%
%      H = LIVESRM returns the handle to a new LIVESRM or the handle to
%      the existing singleton*.
%
%      LIVESRM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LIVESRM.M with the given input arguments.
%
%      LIVESRM('Property','Value',...) creates a new LIVESRM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LiveSRM_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LiveSRM_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LiveSRM

% Last Modified by GUIDE v2.5 25-Jan-2016 17:11:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LiveSRM_OpeningFcn, ...
                   'gui_OutputFcn',  @LiveSRM_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT



%%%%    GUI Initialization / Termination Functions    %%%%

% --- Executes just before LiveSRM is made visible.
function LiveSRM_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LiveSRM (see VARARGIN)

% Choose default command line output for LiveSRM
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes LiveSRM wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = LiveSRM_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

 
% --- Executes during object creation, after setting all properties.
function Isolation_Filter_Threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Isolation_Filter_Threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%    User Interactable Function    %%%%

function [Times,Xpos,Ypos]=ReadDataFile(filename)
    %Opens a .csv file containing a list of localizations and creates Matlab
    %variables labeling the Time, Xposition and Yposition of each
    %localization. The .csv file is assumed to contain one row per
    %localization with the first column representing the time, the second
    %column representing the X position, and the third column representing
    %the Y position. The first row of the file is assumed to be labels for
    %the columns. 
    
    if ~isempty(strfind(filename((end-31):end),'MTT_without_tracking_results.mat')) %Output from MTT Software
        load(filename)
        Times = matrice_results(1,:);
        Xpos = matrice_results(2,:);
        Ypos = matrice_results(3,:);
    elseif ~isempty(strfind(filename((end-6):end),'SRL.csv'))
        FileContents=csvread(filename,1,0); 
        Times = FileContents(:,1)';
        Xpos=FileContents(:,2)';
        Ypos=FileContents(:,3)';
    else
        error('File Type Not Supported')
    end
               
% --- Executes on button press in Load_Data.
function Load_Data_Callback(hObject, eventdata, handles)
% hObject    handle to Load_Data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Load Data
[filename,dirName]=uigetfile('*.*','Select the Cell Data');
handles.filename = filename;
handles.directory = dirName;
[Frames,Xpos,Ypos]=ReadDataFile([dirName,filename]);

%% Initialize Variables and Update Handles
handles.filter = true(1,length(Frames)); %Initialize for later use.
handles.fitParams = [0,0,0]; %Initialize Variables
handles.NuclearArea = 1;
handles.FreehandROI = [];
handles.Xpos=Xpos;
handles.Ypos=Ypos;
handles.Frames=Frames;
handles.InNucleus=true(1,length(Frames));
handles.MinClusterSize=2;
handles.ClusterIDs = zeros(1,length(Frames),'uint16');
handles.HaveLoadedAnnotation = false;
guidata(hObject, handles);

%% Initialize Plot
um_px_size = str2num(get(handles.px_size,'String'))/1000; %Px size in microns
axes(handles.axes1)
plot(handles.Xpos*um_px_size,handles.Ypos*um_px_size,'.k','Markersize',2)
xlabel('um')
ylabel('um')

% --- Executes on button press in Isolation_Filter_On_Off.
function Isolation_Filter_On_Off_Callback(hObject, eventdata, handles)
% hObject    handle to Isolation_Filter_On_Off (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Isolation_Filter_On_Off

Filter_On_Off= get(handles.Isolation_Filter_On_Off,'Value');

if Filter_On_Off
    threshold = str2num(get(handles.Isolation_Filter_Threshold,'String'));
    IsolatedPointIndices = IsolatedDetectionFilter(handles.Frames,handles.Xpos,handles.Ypos,threshold);
    handles.filter = ~IsolatedPointIndices;
    guidata(hObject, handles);
else
    handles.filter = true(1,length(handles.Frames));
    guidata(hObject, handles);
end

GraphUpdateCode(hObject,eventdata,handles)

% --- Executes on button press in NucleusSelection.
function NucleusSelection_Callback(hObject, eventdata, handles)
% hObject    handle to NucleusSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

FreehandROIhandle = imfreehand; %Allows the user to draw a boundary for the nucleus
handles.FreehandROI = getPosition(FreehandROIhandle); %Returns an ordered list of the x and y coordinates that defines the boundary.

um_px_size = str2num(get(handles.px_size,'String'))/1000; %px size in microns

InNucleus = inpolygon(handles.Xpos*um_px_size,handles.Ypos*um_px_size,handles.FreehandROI(:,1),handles.FreehandROI(:,2)); %Returns a logical defining which points lie within the ROI. 
handles.InNucleus = InNucleus;

NuclearArea = polyarea(handles.FreehandROI(:,1),handles.FreehandROI(:,2));

GraphUpdateCode(hObject,eventdata,handles)

Times=1:max(handles.Frames);
Counts = zeros(1,length(Times));
for i = Times
    Counts(i)=sum(handles.Frames(InNucleus)==i);
end
CumulativeCounts = cumsum(Counts);

modelfnhandle = @(params,x)(params(1)*(1-exp(-x/params(2)))+params(3)*x); %Assumes an exonential decay of detection rate, and a constant false positive rate
fitParams = nlinfit(Times,CumulativeCounts,modelfnhandle,[length(handles.Frames),0.5,0]);

figure
plot(Times,modelfnhandle(fitParams,Times),'c')
hold on
plot(Times,CumulativeCounts,'r')
xlabel('Time (Frames)')
ylabel('Cumulative Localizations')
legend('Fit Data','Raw Data')
drawnow

handles.fitParams = fitParams;
handles.NuclearArea = NuclearArea;
guidata(hObject, handles);
    
% --- Executes on button press in Render.
function Render_Callback(hObject, eventdata, handles)
% hObject    handle to Render (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


    Xmin = min(handles.Xpos);Xmax=max(handles.Xpos);Ymin=min(handles.Ypos);Ymax=max(handles.Ypos);
    px_size = str2num(get(handles.px_size,'String'));
    resolution = str2num(get(handles.render_resolution,'String'));
    render_precision = str2num(get(handles.render_precision,'String'));
    sigma_render = render_precision/px_size;
    
    dx=(Xmax-Xmin)/resolution;
    dy=(Ymax-Ymin)/resolution;
    Edges{1}=Xmin:dx:Xmax;
    Edges{2}=Ymin:dy:Ymax;

    Im = hist3([handles.Xpos',handles.Ypos'],'Edges',Edges);
    TempX=-1:dx:1;
    TempY=-1:dy:1;
    ConVecX = exp(-0.5*(TempX/sigma_render).^2); 
    ConVecY = exp(-0.5*(TempY/sigma_render).^2);
    Im2 = conv2(ConVecX,ConVecY,Im);
    Im2=Im2/max(max(Im2));
    Im2=Im2(:,end:-1:1)';

    figure
    imshow(Im2);
    colormap(hot)
    imcontrast(gca)

% --- Executes on button press in AnnotationLoader.
function AnnotationLoader_Callback(hObject, eventdata, handles)
% hObject    handle to AnnotationLoader (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[ClusterIDsFile,dirName]=uigetfile('*.*','Select the Data Annotation File');
handles.ClusterIDsFile=[dirName,ClusterIDsFile];

ClusterIDs=ReadClusterIDsFile(handles.ClusterIDsFile);

handles.ClusterIDs = ClusterIDs;
handles.HaveLoadedAnnotation = true;

guidata(hObject, handles);

GraphUpdateCode(hObject,eventdata,handles)

% --- Executes on button press in DataReviewSelection.
function DataReviewSelection_Callback(hObject, eventdata, handles)
% hObject    handle to DataReviewSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

um_px_size = str2num(get(handles.px_size,'String'))/1000;

rectangle = imrect;
rectangleCorners = getPosition(rectangle);

InXRange = and((handles.Xpos*um_px_size>rectangleCorners(1)),(handles.Xpos*um_px_size<(rectangleCorners(1)+rectangleCorners(3))));
InYRange = and((handles.Ypos*um_px_size>rectangleCorners(2)),(handles.Ypos*um_px_size<(rectangleCorners(2)+rectangleCorners(4))));
temporaryROIindices = and(InXRange,InYRange);

if handles.HaveLoadedAnnotation
    Update = ReviewDataGraphs(handles.Frames,handles.Xpos,handles.Ypos,temporaryROIindices,handles.ClusterIDs,handles.ClusterIDsFile);
else
    Update = ReviewDataGraphs(handles.Frames,handles.Xpos,handles.Ypos,temporaryROIindices,handles.ClusterIDs);
end
if Update
    handles.ClusterIDs=ReadClusterIDsFile(handles.ClusterIDsFile);
end

guidata(hObject, handles);

GraphUpdateCode(hObject,eventdata,handles)

% --- Executes on button press in ROI_Selector.
function ROI_Selector_Callback(hObject, eventdata, handles)
% hObject    handle to ROI_Selector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

um_px_size = str2num(get(handles.px_size,'String'))/1000; %px size in microns

rectangle = imrect;
rectangleCorners_um = getPosition(rectangle);
rectangleCorners=rectangleCorners_um/um_px_size;


AlreadyAnalyzed = (handles.ClusterIDs~=0);

temporaryROIindices = and(~AlreadyAnalyzed,((handles.Xpos>rectangleCorners(1))&(handles.Xpos<(rectangleCorners(1)+rectangleCorners(3))))&((handles.Ypos>rectangleCorners(2))&(handles.Ypos<(rectangleCorners(2)+rectangleCorners(4)))));
[ClustersMatrix,MinClusterSize] = ROIAnalysis(handles.Frames,handles.Xpos*um_px_size,handles.Ypos*um_px_size,temporaryROIindices,handles.fitParams,handles.NuclearArea,handles.MinClusterSize);
handles.MinClusterSize=MinClusterSize;

[a,~] = size(ClustersMatrix);
for i = 1:a
    LargestClusterID = max(handles.ClusterIDs);
    handles.ClusterIDs = (LargestClusterID+1)*uint16(ClustersMatrix(i,:))+handles.ClusterIDs; %Each ROI will be indexed by a unique integer
end

if handles.HaveLoadedAnnotation
    WriteClusterIDsFile(handles.ClusterIDsFile,handles.ClusterIDs)
else
    testname='ClusterIDs.csv';
    count=1;
    while exist(testname,'file') == 2
        count = count+1;
        testname = ['ClusterIDs_',num2str(count),'.csv'];
    end
    handles.ClusterIDsFile = testname;
    handles.HaveLoadedAnnotation=true;
    WriteClusterIDsFile(handles.ClusterIDsFile,handles.ClusterIDs)
end

guidata(hObject, handles);

GraphUpdateCode(hObject,eventdata,handles)
 
% --- Executes on button press in runpcPALM.
function runpcPALM_Callback(hObject, eventdata, handles)
% hObject    handle to runpcPALM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

px_size = str2num(get(handles.px_size,'String'));
pc_bin_size = str2num(get(handles.pc_bin_size,'String'));
pc_max_length = str2num(get(handles.pc_max_length,'String'));
pc_loc_precision=str2num(get(handles.loc_precision,'String'));

[r,g,included_points]=FreehandPairCorrelationRegion(handles.Frames,handles.Xpos,handles.Ypos,pc_bin_size/px_size,pc_max_length/px_size,pc_loc_precision/px_size,px_size);

%% Save PC Output

% I should also save an image of which points were selected, and a list of
% which points were selected.

pc_data_array=cell(length(r)+1,2);
pc_data_array{1,1} = 'Radius (nm)';
pc_data_array{1,2} = 'g(r)';
for i = 1:length(r)
    pc_data_array{i+1,1}=r(i)*px_size;
    pc_data_array{i+1,2}=g(i);
end

included_data_array=cell(length(included_points)+1,1);
included_data_array{1} = 'True/False This Data point was used in the evaluation of the pair correlation function';
for i = 1:length(included_points)
    included_data_array{i+1}=included_points(i);
end

current_directory=cd;
cd(handles.directory)

test_name = 'pair_correlation.csv';
included_data_filename = 'pc_included_points.csv';
count=1;
while exist(test_name,'file')==2
    count = count+1;
    test_name = ['pair_correlation_',num2str(count),'.csv'];
    included_data_filename = ['pc_included_points_',num2str(count),'.csv'];
end

filehandle = fopen([handles.directory,test_name],'w');
[a,b] = size(pc_data_array);
for i = 1:a
    for j =1:b
        if isnumeric(pc_data_array{i,j})
            fprintf(filehandle,num2str(pc_data_array{i,j}));
        else
            fprintf(filehandle,pc_data_array{i,j});
        end
        if j~=b
            fprintf(filehandle,',');
        end
    end
    fprintf(filehandle,'\n');
end

filehandle = fopen([handles.directory,included_data_filename],'w');
[a,b] = size(included_data_array);
for i = 1:a
    for j =1:b
        if islogical(included_data_array{i,j})
            fprintf(filehandle,num2str(double(included_data_array{i,j})));
        else
            fprintf(filehandle,included_data_array{i,j});
        end
        if j~=b
            fprintf(filehandle,',');
        end
    end
    fprintf(filehandle,'\n');
end

%% Save Image of PC Region

figure('Visible','Off')
plot(handles.Xpos*px_size/1000,handles.Ypos*px_size/1000,'.','Markersize',2,'Color',[0;0;0]/255)
xlabel('um')
ylabel('um')
title('Data Points used for pair correlation function')
hold on
plot(handles.Xpos(included_points)*px_size/1000,handles.Ypos(included_points)*px_size/1000,'.','Markersize',6,'Color',[213;94;0]/255)
plot_file_name=[included_data_filename(1:end-4),'.jpg'];
saveas(gcf,plot_file_name,'jpg')
close(gcf)

cd(current_directory)

GraphUpdateCode(hObject,eventdata,handles)

function Isolation_Filter_Threshold_Callback(hObject, eventdata, handles)
% hObject    handle to Isolation_Filter_Threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Isolation_Filter_Threshold as text
%        str2double(get(hObject,'String')) returns contents of Isolation_Filter_Threshold as a double
    
 



%%%%    Helper Functions    %%%%    
    
function GraphUpdateCode(hObject,eventdata,handles)

AlreadyAnalyzed = (handles.ClusterIDs~=0);
axes(handles.axes1)
CurrentAxisValues = axis;
um_px_size = str2num(get(handles.px_size,'String'))/1000; %px size in microns
 
if ~isempty(handles.FreehandROI)
    hold off
    minx = min(handles.Xpos*um_px_size);miny=min(handles.Ypos*um_px_size);maxx=max(handles.Xpos*um_px_size);maxy=max(handles.Ypos*um_px_size);
    fill([minx,maxx,maxx,minx,minx],[maxy,maxy,miny,miny,maxy],[1,0.9,0.9])
    hold on
    fill(handles.FreehandROI(:,1),handles.FreehandROI(:,2),[1,1,1])
    plot(handles.Xpos((handles.filter&(~AlreadyAnalyzed)))*um_px_size,handles.Ypos((handles.filter&(~AlreadyAnalyzed)))*um_px_size,'.','Markersize',2,'Color',[0;0;0]/255)
else
    hold off
    plot(handles.Xpos((handles.filter&(~AlreadyAnalyzed)))*um_px_size,handles.Ypos((handles.filter&(~AlreadyAnalyzed)))*um_px_size,'.','Markersize',2,'Color',[0;0;0]/255)
    hold on
end

plot(handles.Xpos((handles.filter&(AlreadyAnalyzed)))*um_px_size,handles.Ypos((handles.filter&(AlreadyAnalyzed)))*um_px_size,'.','Markersize',6,'Color',[213;94;0]/255)
xlabel('um')
ylabel('um')
hold off
axis(CurrentAxisValues)




%%%%    Vestigial Code    %%%%


% --- Executes on button press in AutoDetect.
% function AutoDetect_Callback(hObject, eventdata, handles)
% % hObject    handle to AutoDetect (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% if false  %%% Don't want this to run until code is finished. 
%     FastJetPath = uigetdir('/','Choose FastJet Directory');
% 
%     suffixindex=strfind(handles.filename,'.mat');
%     WriteFileName= [handles.directory,handles.filename(1:suffixindex-1),'_FastJet_Input.txt'];
%     
%     if ~exist(handles.InNucleus,'var')
%         [handles.NuclearArea,handles.InNucleus]=FreeHandNucleus(matrice_results); % Should add a pop-up telling the user that their input is wanted at this point in time. 
%     end
%     
% 
%     FastJetNuclearClusters=AutoDetectROIs(FastJetPath,FastJetInputFileName,FastJetOutputFileName,NuclearData,LengthScale,MinSize)
%     
%     FastJetClusters=zeros(1,length(handles.Frames));
%     FastJetClusters(handles.InNucleus)=FastJetNuclearClusters;
%     
%     guidata(hObject, handles);
% end

function pc_bin_size_Callback(hObject, eventdata, handles)
% hObject    handle to pc_bin_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pc_bin_size as text
%        str2double(get(hObject,'String')) returns contents of pc_bin_size as a double


% --- Executes during object creation, after setting all properties.
function pc_bin_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pc_bin_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pc_max_length_Callback(hObject, eventdata, handles)
% hObject    handle to pc_max_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pc_max_length as text
%        str2double(get(hObject,'String')) returns contents of pc_max_length as a double


% --- Executes during object creation, after setting all properties.
function pc_max_length_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pc_max_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function loc_precision_Callback(hObject, eventdata, handles)
% hObject    handle to loc_precision (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of loc_precision as text
%        str2double(get(hObject,'String')) returns contents of loc_precision as a double


% --- Executes during object creation, after setting all properties.
function loc_precision_CreateFcn(hObject, eventdata, handles)
% hObject    handle to loc_precision (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function render_resolution_Callback(hObject, eventdata, handles)
% hObject    handle to render_resolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of render_resolution as text
%        str2double(get(hObject,'String')) returns contents of render_resolution as a double


% --- Executes during object creation, after setting all properties.
function render_resolution_CreateFcn(hObject, eventdata, handles)
% hObject    handle to render_resolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function render_precision_Callback(hObject, eventdata, handles)
% hObject    handle to render_precision (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of render_precision as text
%        str2double(get(hObject,'String')) returns contents of render_precision as a double


% --- Executes during object creation, after setting all properties.
function render_precision_CreateFcn(hObject, eventdata, handles)
% hObject    handle to render_precision (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px_size_Callback(hObject, eventdata, handles)
% hObject    handle to px_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px_size as text
%        str2double(get(hObject,'String')) returns contents of px_size as a double


% --- Executes during object creation, after setting all properties.
function px_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in aggregate_data.
function aggregate_data_Callback(hObject, eventdata, handles)
% hObject    handle to aggregate_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

px_size = str2num(get(handles.px_size,'String'));
DataAggregation(handles.directory,handles.Frames,handles.Xpos,handles.Ypos,handles.ClusterIDs,px_size)



function render_pixel_size_Callback(hObject, eventdata, handles)
% hObject    handle to render_pixel_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of render_pixel_size as text
%        str2double(get(hObject,'String')) returns contents of render_pixel_size as a double


% --- Executes during object creation, after setting all properties.
function render_pixel_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to render_pixel_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
