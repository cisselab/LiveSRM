function varargout = ReviewDataGraphs(varargin)
% REVIEWDATAGRAPHS MATLAB code for ReviewDataGraphs.fig
%      REVIEWDATAGRAPHS, by itself, creates a new REVIEWDATAGRAPHS or raises the existing
%      singleton*.
%
%      H = REVIEWDATAGRAPHS returns the handle to a new REVIEWDATAGRAPHS or the handle to
%      the existing singleton*.
%
%      REVIEWDATAGRAPHS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in REVIEWDATAGRAPHS.M with the given input arguments.
%
%      REVIEWDATAGRAPHS('Property','Value',...) creates a new REVIEWDATAGRAPHS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ReviewDataGraphs_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ReviewDataGraphs_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ReviewDataGraphs

% Last Modified by GUIDE v2.5 28-Jul-2015 14:51:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ReviewDataGraphs_OpeningFcn, ...
                   'gui_OutputFcn',  @ReviewDataGraphs_OutputFcn, ...
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

% --- Executes just before ReviewDataGraphs is made visible.
function ReviewDataGraphs_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ReviewDataGraphs (see VARARGIN)

% Choose default command line output for ReviewDataGraphs
handles.output = hObject;
handles.update=false;

handles.Frames = varargin{1};
handles.Xpos = varargin{2};
handles.Ypos = varargin{3};
handles.ROIindices = varargin{4};  
ClusterIdentities = varargin{5};
if numel(varargin)>5
    handles.ClusterIDsFile = varargin{6};
end

handles.ClusterID = ClusterIdentities;

axes(handles.SpatialAxes);
plot(handles.Xpos(handles.ROIindices),handles.Ypos(handles.ROIindices),'ok','MarkerSize',2)

X = 1:max(handles.Frames);
Y = zeros(1,length(X));
for i = 1:length(X)
    NumberOfDetections = sum(handles.Frames(handles.ROIindices)==i);
    Y(i) = NumberOfDetections;
end
Z = cumsum(Y);

handles.X=X;
handles.Y=Y;
handles.Z=Z;

axes(handles.TemporalAxes);
plot(X,Y,'k')
axes(handles.CumulativeAxes);
plot(X,cumsum(Y),'k')

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
axes(handles.SpatialAxes)
hold on
for i = 1:max(ClusterIdentities)
    Filter = and(handles.ROIindices,ClusterIdentities==i);
    plot(handles.Xpos(Filter),handles.Ypos(Filter),Markers(mod(ceil(i/4),4)+1),'MarkerFaceColor',ColorScheme(mod(i,7)+1,:)/255,'Color',ColorScheme(mod(i,7)+1,:)/255,'MarkerSize',4)
end

axes(handles.CumulativeAxes)
hold on
for i = 1:max(ClusterIdentities)
    Filter = and(handles.ROIindices,ClusterIdentities==i);
    if sum(Filter)>0
        StartTime = min(handles.Frames(Filter));
        EndTime = max(handles.Frames(Filter));
        plotIndices = and(X>=StartTime,X<=EndTime);
        plot(X(plotIndices),Z(plotIndices),'LineWidth',3,'Color',ColorScheme(mod(i,7)+1,:)/255)
    end
end

%%% End Graph Update Code %%%

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ReviewDataGraphs wait for user response (see UIRESUME)
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = ReviewDataGraphs_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure


varargout{1}= handles.update; %The output tells the user whether or not to update their graph.

delete(hObject);

% --- Executes on button press in DeleteCluster.
function DeleteCluster_Callback(hObject, eventdata, handles)
% hObject    handle to DeleteCluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.SpatialAxes)

rectangle = imrect;
rectangleCorners = getPosition(rectangle);

InROI = and(and(handles.Xpos>rectangleCorners(1),handles.Xpos<(rectangleCorners(1)+rectangleCorners(3))),and(handles.Ypos>rectangleCorners(2),handles.Ypos<(rectangleCorners(2)+rectangleCorners(4))));

maxID = max(handles.ClusterID(InROI));

while maxID>0
    handles.ClusterID(handles.ClusterID==maxID)=0;
    maxID = max(handles.ClusterID(InROI));
end



axes(handles.TemporalAxes);
hold off
plot(handles.X,handles.Y,'k')
axes(handles.CumulativeAxes);
hold off
plot(handles.X,handles.Z,'k')

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
axes(handles.SpatialAxes)
hold off

plot(handles.Xpos(handles.ROIindices),handles.Ypos(handles.ROIindices),'.k','Markersize',2)
hold on
for i = 1:max(handles.ClusterID)
    Filter = and(handles.ROIindices,handles.ClusterID==i);
    plot(handles.Xpos(Filter),handles.Ypos(Filter),Markers(mod(ceil(i/4),4)+1),'MarkerFaceColor',ColorScheme(mod(i,7)+1,:)/255,'Color',ColorScheme(mod(i,7)+1,:)/255,'MarkerSize',4)
end

axes(handles.CumulativeAxes)
hold on
for i = 1:max(handles.ClusterID)
    Filter = and(handles.ROIindices,handles.ClusterID==i);
    if sum(Filter)>0
        StartTime = min(handles.Frames(Filter));
        EndTime = max(handles.Frames(Filter));
        plotIndices = and(handles.X>=StartTime,handles.X<=EndTime);
        plot(handles.X(plotIndices),handles.Z(plotIndices),'LineWidth',3,'Color',ColorScheme(mod(i,7)+1,:)/255)
    end
end

guidata(hObject, handles);

uiwait(handles.figure1);

% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'ClusterIDsFile')

    handles.update=true;
    WriteClusterIDsFile(handles.ClusterIDsFile,handles.ClusterID)
    guidata(hObject,handles);
    
end

uiwait(handles.figure1);

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end
