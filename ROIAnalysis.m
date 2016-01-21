function varargout = ROIAnalysis(varargin)
% ROIANALYSIS MATLAB code for ROIAnalysis.fig
%      ROIANALYSIS, by itself, creates a new ROIANALYSIS or raises the existing
%      singleton*.
%
%      H = ROIANALYSIS returns the handle to a new ROIANALYSIS or the handle to
%      the existing singleton*.
%
%      ROIANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ROIANALYSIS.M with the given input arguments.
%
%      ROIANALYSIS('Property','Value',...) creates a new ROIANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ROIAnalysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ROIAnalysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ROIAnalysis

% Last Modified by GUIDE v2.5 28-Jul-2015 15:55:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ROIAnalysis_OpeningFcn, ...
                   'gui_OutputFcn',  @ROIAnalysis_OutputFcn, ...
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

% --- Outputs from this function are returned to the command line.
function varargout = ROIAnalysis_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if handles.savedata
    varargout{1} = handles.Clusters; 
    varargout{2}=get(handles.Cluster_Cutoff_Input,'String');
    
    %% Poster Graph Output
     ClusterIDs = handles.Clusters;
     ROIIDs = handles.ROIindices;
    %%

    CurrentDirectory = cd;
    DirectoryName = strcat(CurrentDirectory,'ROIs');
    if ~isdir(DirectoryName)
        mkdir(DirectoryName);
    end
    
    cd(DirectoryName)
    if exist([DirectoryName,'/1.mat'],'file')
        DataSets = dir('*.mat');
        DataSetName = [num2str(length(DataSets)+1),'.mat'];
    else
        DataSetName = '1.mat';
    end
    X = handles.X;
    Y = handles.Y;
    save (DataSetName,'X','Y');
    cd(CurrentDirectory)
    
else
    varargout{1}=false(1,length(handles.ROIindices));
    varargout{2}=get(handles.Cluster_Cutoff_Input,'String');
end

delete(handles.figure1)

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end

% --- Executes during object creation, after setting all properties.
function Cluster_Number_Selector_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Cluster_Number_Selector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes during object creation, after setting all properties.
function Cluster_Cutoff_Input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Cluster_Cutoff_Input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in DisplayCumDetTrace.
function DisplayCumDetTrace_Callback(hObject, eventdata, handles)
% hObject    handle to DisplayCumDetTrace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DisplayCumDetTrace

GraphUpdateCode(hObject,eventdata,handles)

function Cluster_Cutoff_Input_Callback(hObject, eventdata, handles)
% hObject    handle to Cluster_Cutoff_Input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Cluster_Cutoff_Input as text
%        str2double(get(hObject,'String')) returns contents of Cluster_Cutoff_Input as a double

GraphUpdateCode(hObject,eventdata,handles)

% --- Executes on button press in FitTrace.
function FitTrace_Callback(hObject, eventdata, handles)
% hObject    handle to FitTrace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of FitTrace

GraphUpdateCode(hObject,eventdata,handles)

% --- Executes on button press in Save_Data.
function Save_Data_Callback(hObject, eventdata, handles)
% hObject    handle to Save_Data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.savedata = true;
guidata(hObject, handles);
uiresume

% --- Executes just before ROIAnalysis is made visible.
function ROIAnalysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ROIAnalysis (see VARARGIN)

% Choose default command line output for ROIAnalysis
handles.output = hObject;

handles.Frames = varargin{1};
handles.Xpos = varargin{2};
handles.Ypos = varargin{3};
ROIindices = varargin{4}; % After Analysis, this variable will be passed out as a matrix with each separate cluster expressed as a different line.
DetectionRateModelParameters = varargin{5}; %[Limit Value, TimeConstant, False Detection Rate]
NuclearArea = varargin{6};
InitialMinClusterSize=varargin{7};

set(handles.Cluster_Cutoff_Input,'String',InitialMinClusterSize)

handles.LimitValue=DetectionRateModelParameters(1);
handles.TimeConstant=DetectionRateModelParameters(2);
handles.FalseDetRate=DetectionRateModelParameters(3);
handles.NucArea=NuclearArea;
handles.WinArea=(max(handles.Xpos(ROIindices))-min(handles.Xpos(ROIindices)))*(max(handles.Ypos(ROIindices))-min(handles.Ypos(ROIindices)));
handles.ROIindices = ROIindices; % Will not be divided into a matrix of separate clusters.
handles.Clusters = ROIindices;
axes(handles.Spatial_Axes);

plot(handles.Xpos(ROIindices),handles.Ypos(ROIindices),'ok','MarkerSize',2)

X = 1:max(handles.Frames);
Y = zeros(1,length(X));
for i = 1:length(X)
    NumberOfDetections = sum(handles.Frames(ROIindices)==i);
    Y(i) = NumberOfDetections;
end

axes(handles.Detection_Axes);
plot(X,Y,'k')
axes(handles.Cumulative_Axes);
plot(X,cumsum(Y),'k')

handles.X = X;
handles.Y = Y;

handles.savedata = false;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ROIAnalysis wait for user response (see UIRESUME)
uiwait(handles.figure1);

% --- Executes on button press in Subsection_Selector.
function Subsection_Selector_Callback(hObject, eventdata, handles)
% hObject    handle to Subsection_Selector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.Spatial_Axes);
rectangle = imrect;
rectangleCorners = getPosition(rectangle);

ROIindices = and(handles.ROIindices,((handles.Xpos>rectangleCorners(1))&(handles.Xpos<(rectangleCorners(1)+rectangleCorners(3))))&((handles.Ypos>rectangleCorners(2))&(handles.Ypos<(rectangleCorners(2)+rectangleCorners(4)))));

handles.ROIindices = ROIindices;

handles.WinArea=(max(handles.Xpos(ROIindices))-min(handles.Xpos(ROIindices)))*(max(handles.Ypos(ROIindices))-min(handles.Ypos(ROIindices)));

guidata(hObject, handles);

GraphUpdateCode(hObject,eventdata,handles)

uiwait(handles.figure1);

% --- Executes on slider movement.
function Cluster_Number_Selector_Callback(hObject, eventdata, handles)
% hObject    handle to Cluster_Number_Selector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

Number_Slider_Value = get(handles.Cluster_Number_Selector,'Value'); %Reads the value of Number slider, which was just adjusted.  
set(handles.Cluster_Number_Display,'string',['Detection Sensitivity:',num2str(Number_Slider_Value)])
GraphUpdateCode(hObject,eventdata,handles)

function GraphUpdateCode(hObject,eventdata,handles)

display(handles.FalseDetRate)

displayFit = get(handles.FitTrace,'Value');
displayBackground = get(handles.DisplayCumDetTrace,'Value');

ClusterSizeCutoff = str2num(get(handles.Cluster_Cutoff_Input,'String'));
Number_Slider_Value = get(handles.Cluster_Number_Selector,'Value'); %Reads the value of Number slider, which was just adjusted.  
[Start_Times,End_Times] = HierarchichalClusterIdentification(handles.Frames,handles.ROIindices,Number_Slider_Value,ClusterSizeCutoff);

%%% Graph Update Code %%%
NumberOfClusters = length(Start_Times);
if NumberOfClusters>0
    Clusters = false(NumberOfClusters,length(handles.Frames));
    for i = 1:length(Start_Times)
        Clusters(i,:)=(((handles.Frames>=Start_Times(i))&(handles.Frames<=End_Times(i)))&handles.ROIindices);
    end
else
    Clusters = false(1,length(handles.Frames));
end
handles.Clusters = Clusters;
guidata(hObject, handles);

axes(handles.Spatial_Axes);
hold off
plot(handles.Xpos(handles.ROIindices),handles.Ypos(handles.ROIindices),'ok','MarkerSize',2)

X = 1:max(handles.Frames);
Y = zeros(1,length(X));
for i = 1:length(X)
    NumberOfDetections = sum(handles.Frames(handles.ROIindices)==i);
    Y(i) = NumberOfDetections;
end
Z = cumsum(Y);

axes(handles.Detection_Axes);
hold off
plot(X,Y,'k')
axes(handles.Cumulative_Axes);
hold off
plot(X,Z,'k')

% Color Code Graphs %

ColorScheme = [213,94,0;... Vermillion
    86,180,233;... Sky Blue
    240,228,66;... Yellow
    204,121,167;... Reddish Purple
    0,158,115;... Bluish Green
    230,159,0;... Orange
    0,114,178;... Blue
    0,0,0]; %Black

Markers = '^osd';

% Spatial
axes(handles.Spatial_Axes)
if NumberOfClusters <= 1000
    hold on
    for i = 1:NumberOfClusters
        plot(handles.Xpos(Clusters(i,:)),handles.Ypos(Clusters(i,:)),Markers(mod(ceil(i/4),4)+1),'MarkerFaceColor',ColorScheme(mod(i,7)+1,:)/255,'Color',ColorScheme(mod(i,7)+1,:)/255,'MarkerSize',4)
    end
else
    display('Need a Larger Color Scheme!!!!!!')
end 

axes(handles.Cumulative_Axes)
if NumberOfClusters <= 1000
    hold on
    for i = 1:NumberOfClusters
        plotIndices = and(X>=Start_Times(i),X<=End_Times(i));
        plot(X(plotIndices),Z(plotIndices),'LineWidth',3,'Color',ColorScheme(mod(i,7)+1,:)/255)
    end
end 

if displayFit
    functionHandle = @(params,x)(params*(1-exp(-x{1}/x{2}))+x{1}*x{3});
    handles.A=nlinfit({X,handles.TimeConstant,handles.WinArea/handles.NucArea*handles.FalseDetRate},Z,functionHandle,10);
    AreaFraction = handles.WinArea/handles.NucArea;
    ExponentialTerm = (1-exp(-X/handles.TimeConstant));
    Zfit = handles.A*ExponentialTerm+AreaFraction*(handles.FalseDetRate*X);
    %Zfit = A*ExponentialTerm;
    hold on
    plot(X,Zfit,'g')
else
    handles.A=0;
end

if displayBackground
    AreaFraction = handles.WinArea/handles.NucArea;
    ExponentialTerm = (1-exp(-X/handles.TimeConstant));
    AverageCumTrace = AreaFraction*(handles.LimitValue*ExponentialTerm+handles.FalseDetRate*X);
    plot(X,AverageCumTrace,'r')
end
