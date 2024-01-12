function varargout = craterdating(varargin)
% CRATERDATING MATLAB code for craterdating.fig
%      CRATERDATING, by itself, creates a new CRATERDATING or raises the existing
%      singleton*.
%
%      H = CRATERDATING returns the handle to a new CRATERDATING or the handle to
%      the existing singleton*.
%
%      CRATERDATING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CRATERDATING.M with the given input arguments.
%
%      CRATERDATING('Property','Value',...) creates a new CRATERDATING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before craterdating_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to craterdating_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help craterdating

% Last Modified by GUIDE v2.5 11-Jan-2024 23:42:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @craterdating_OpeningFcn, ...
    'gui_OutputFcn',  @craterdating_OutputFcn, ...
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
% pathname={};
% End initialization code - DO NOT EDIT


% --- Executes just before craterdating is made visible.
function craterdating_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to craterdating (see VARARGIN)
% Choose default command line output for craterdating

if 1
    %Width and height of displays, returned as an n-by-4 matrix, where n is the number of displays. Each row corresponds to one display and is a four-element vector of the form [x y width height]. For example, if there are two displays, then the matrix has this form:
    %[x1 y1 width1 height1
    % x2 y2 width2 height2]
    scnsize = get(0,'MonitorPosition'); %[1 1  1920 1080]
    [~,index]=min(scnsize(:,1));
    scnsize=scnsize(index,:);
    
    
    f_width=1200/scnsize(3);
    f_height=500/scnsize(4);
else
    ScreenPixelsPerInch = java.awt.Toolkit.getDefaultToolkit().getScreenResolution();
    ScreenDevices = java.awt.GraphicsEnvironment.getLocalGraphicsEnvironment().getScreenDevices();
    MainScreen = java.awt.GraphicsEnvironment.getLocalGraphicsEnvironment().getDefaultScreenDevice().getScreen()+1;
    MainBounds = ScreenDevices(MainScreen).getDefaultConfiguration().getBounds();
    f_width=MainBounds.width;
    f_height=MainBounds.height;
end
% handles.figure1.Position=[0.05,0.15,0.9,0.8];
handles.axes1.FontUnits='normalized';
handles.axes1.Box='on';
handles.axes1.FontSize=0.025;
handles.axes1.TitleFontSizeMultiplier=1.7;%set(handles.axes1.TitleFontSizeMultiplier,1.5)
handles.axes1.LabelFontSizeMultiplier,1.6;%set(handles.axes1.LabelFontSizeMultiplier,1.4)
% handles.axes1.XLabel='Diameter (km)';%,'Fontsize',16
% handles.axes1.YLabel='N(1) (km^{-2})';

% Choose default command line output for craterdating
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
% p=[f_width*0.05 f_height*0.25 f_width*0.8 f_height*0.7]*0.8;
xlabel('Diameter (km)','FontName','Times New Roman','Fontsize',20)%,'Fontsize',16
ylabel('Density of craters larger than 1 km in diameter (km^{-2})','FontName','Times New Roman','Fontsize',20)

% --- Outputs from this function are returned to the command line.
function varargout = craterdating_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1
handles.data.index_selected= get(handles.listbox1,'value'); %index的值代表我们选的是第几个选项
isjustPlot=1;
handles=main(handles,isjustPlot);

% Choose default command line output for craterdating
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in open.
function open_Callback(hObject, eventdata, handles)
% hObject    handle to open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'data')
    handles.data=[];
end

[filename, pathname2] = uigetfile( ...
    {'*.scc;*.shp;*.dbf;*.csfd','All Image Files'},...
    '请选择要修改的图片（可多选）', '.\data', ...
    'MultiSelect', 'on');
if ~(ischar(filename)||iscell(filename))
    return;
end
if ischar(filename)
    filename={filename};
elseif length(filename)>1
    filename=filename';
end

handles=editdata(hObject, eventdata, handles,pathname2,filename);

isjustPlot=1;
handles=main(handles,isjustPlot);

% Choose default command line output for craterdating
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);




% --- Executes on button press in togglebutton1.
function togglebutton1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton1


% --- Executes on button press in Fitting_method.
function handles=Fitting_method_Callback(hObject, eventdata, handles)
% hObject    handle to Fitting_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Fitting_method
type_string=get(handles.Fitting_method,'String');
Fitting_method=get(handles.Fitting_method,'value');

index= get(handles.listbox1,'value'); %index的值代表我们选的是第几个选项
for i=1:length(index)
    handles.data.Fitting_method{index(i)}=type_string{Fitting_method};
end
% Choose default command line output for craterdating
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

function handles=Dmin_Callback(hObject, eventdata, handles)
% hObject    handle to Dmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Dmin as text
%        str2double(get(hObject,'String')) returns contents of Dmin as a double
index= get(handles.listbox1,'value'); %index的值代表我们选的是第几个选项
for i=1:length(index)
    handles.data.Dmin{index(i)}=str2double(get(handles.Dmin,'string'));
end
% Choose default command line output for craterdating
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Dmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Dmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function handles=Dmax_Callback(hObject, eventdata, handles)
% hObject    handle to Dmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Dmax as text
%        str2double(get(hObject,'String')) returns contents of Dmax as a double
index= get(handles.listbox1,'value'); %index的值代表我们选的是第几个选项
for i=1:length(index)
    handles.data.Dmax{index(i)}=str2double(get(handles.Dmax,'string'));
end
% Choose default command line output for craterdating
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Dmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Dmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in SFD_method.
function handles=SFD_method_Callback(hObject, eventdata, handles)
% hObject    handle to SFD_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns SFD_method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SFD_method
type_string=get(handles.SFD_method,'String');
SFD_method=get(handles.SFD_method,'value');

index_selected= get(handles.listbox1,'value'); %index的值代表我们选的是第几个选项
for i=1:length(index_selected)
    handles.data.SFD_method{index_selected(i)}=type_string{SFD_method};
end
% Choose default command line output for craterdating
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function SFD_method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SFD_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in LineColor.
function handles=LineColor_Callback(hObject, eventdata, handles)
% hObject    handle to LineColor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns LineColor contents as cell array
%        contents{get(hObject,'Value')} returns selected item from LineColor
String=get(handles.LineColor,'String');
index_temp=get(handles.LineColor,'value');
String=String{index_temp};

index_selected=get(handles.listbox1,'value');
for i=1:length(index_selected)
    handles.data.LineColor{index_selected(i)}=String(1);
    if handles.data.getsetting==0
        set(handles.data.bestfit_handle{index_selected(i)},'color',String(1));
    end
end
% Choose default command line output for craterdating
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function LineColor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LineColor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in line.
function handles=line_Callback(hObject, eventdata, handles)
% hObject    handle to line (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns line contents as cell array
%        contents{get(hObject,'Value')} returns selected item from line
String=get(handles.line,'String');
index_temp=get(handles.line,'value');
String=String{index_temp};

index_selected=get(handles.listbox1,'value');
for i=1:length(index_selected)
    handles.data.LineType{index_selected(i)}=String(1:2);
    if handles.data.getsetting==0
        set(handles.data.bestfit_handle{index_selected(i)},'LineStyle',String(1:2));
    end
end
% Choose default command line output for craterdating
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function line_CreateFcn(hObject, eventdata, handles)
% hObject    handle to line (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Marker.
function handles=Marker_Callback(hObject, eventdata, handles)
% hObject    handle to Marker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Marker contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Marker
String=get(handles.Marker,'String');
index_temp=get(handles.Marker,'value');
String=String{index_temp};

index_selected=get(handles.listbox1,'value');
for i=1:length(index_selected)
    handles.data.MarkerType{index_selected(i)}=String(1);
    if handles.data.getsetting==0
        set(handles.data.datasymbol_handle{index_selected(i)},'Marker',String(1));
    end
end
% Choose default command line output for craterdating
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Marker_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Marker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in LineWidth.
function handles=LineWidth_Callback(hObject, eventdata, handles)
% hObject    handle to LineWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns LineWidth contents as cell array
%        contents{get(hObject,'Value')} returns selected item from LineWidth
String=get(handles.LineWidth,'String');
LineWidth=str2double(String);

index_selected=get(handles.listbox1,'value');
for i=1:length(index_selected)
    handles.data.LineWidth{index_selected(i)}=LineWidth;
    if handles.data.getsetting==0
        set(handles.data.bestfit_handle{index_selected(i)},'LineWidth',LineWidth);
    end
end
% Choose default command line output for craterdating
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function LineWidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LineWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in MarkerSize.
function handles=MarkerSize_Callback(hObject, eventdata, handles)
% hObject    handle to MarkerSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MarkerSize contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MarkerSize

String=get(handles.MarkerSize,'String');
MarkerSize=str2double(String);

index_selected=get(handles.listbox1,'value');
for i=1:length(index_selected)
    handles.data.MarkerSize{index_selected(i)}=MarkerSize;
    if handles.data.getsetting==0
        set(handles.data.datasymbol_handle{index_selected(i)},'MarkerSize',MarkerSize);
    end
end
% Choose default command line output for craterdating
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function MarkerSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MarkerSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in MarkerColor.
function handles=MarkerColor_Callback(hObject, eventdata, handles)
% hObject    handle to MarkerColor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MarkerColor contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MarkerColor
String=get(handles.MarkerColor,'String');
index_temp=get(handles.MarkerColor,'value');
String=String{index_temp};

index_selected=get(handles.listbox1,'value');
for i=1:length(index_selected)
    handles.data.MarkerColor{index_selected(i)}=String(1);
    if handles.data.getsetting==0
        set(handles.data.datasymbol_handle{index_selected(i)},'color',String(1));
    end
end
% Choose default command line output for craterdating
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function MarkerColor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MarkerColor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in PF.
function handles=PF_Callback(hObject, eventdata, handles)
% hObject    handle to PF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns PF contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PF
type_string=get(handles.PF,'String');
index_type=get(handles.PF,'value');

index_selected= get(handles.listbox1,'value'); %index的值代表我们选的是第几个选项
for i=1:length(index_selected)
    handles.data.PF{index_selected(i)}=type_string{index_type};
end
% Choose default command line output for craterdating
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function PF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ChronologyFunction.
function handles=ChronologyFunction_Callback(hObject, eventdata, handles)
% hObject    handle to ChronologyFunction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ChronologyFunction contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ChronologyFunction
type_string=get(handles.ChronologyFunction,'String');
index_type=get(handles.ChronologyFunction,'value');

index_selected= get(handles.listbox1,'value'); %index的值代表我们选的是第几个选项
for i=1:length(index_selected)
    handles.data.ChronologyFunction{index_selected(i)}=type_string{index_type};
end
% Choose default command line output for craterdating
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function ChronologyFunction_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ChronologyFunction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function uitoggletool6_OffCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function uitoggletool6_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in Axis_limits.
function Axis_limits_Callback(hObject, eventdata, handles)
% hObject    handle to Axis_limits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Axis_limits contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Axis_limits
axes(handles.axes1);
index=get(handles.Axis_limits,'value');
if index==2
    Xlim=get(gca,'xlim');
    set(handles.Axis_lowerlimits,'string',num2str(Xlim(1)));
    set(handles.Axis_upperlimits,'string',num2str(Xlim(2)));
else
    Ylim=get(gca,'ylim');
    set(handles.Axis_lowerlimits,'string',num2str(Ylim(1)));
    set(handles.Axis_upperlimits,'string',num2str(Ylim(2)));
end


function Axis_lowerlimits_Callback(hObject, eventdata, handles)
% hObject    handle to Axis_lowerlimits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Axis_lowerlimits as text
%        str2double(get(hObject,'String')) returns contents of Axis_lowerlimits as a double
axes(handles.axes1);
index=get(handles.Axis_limits,'value');
if index==2
    Xlim2=get(gca,'xlim');
    Xlim(1)=str2double(get(handles.Axis_lowerlimits,'string'));
    if ~isinf(Xlim(1))
        Xlim2(1)=Xlim(1);
    end
    Xlim(2)=str2double(get(handles.Axis_upperlimits,'string'));
    if ~isinf(Xlim(2))
        Xlim2(2)=Xlim(2);
    end
    xlim(Xlim2);
else
    Ylim2=get(gca,'ylim');
    Ylim(1)=str2double(get(handles.Axis_lowerlimits,'string'));
    if ~isinf(Ylim(1))
        Ylim2(1)=Ylim(1);
    end
    Ylim(2)=str2double(get(handles.Axis_upperlimits,'string'));
    if ~isinf(Ylim(2))
        Ylim2(2)=Ylim(2);
    end
    ylim(Ylim2);
end

% --- Executes during object creation, after setting all properties.
function Axis_lowerlimits_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Axis_lowerlimits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Axis_upperlimits_Callback(hObject, eventdata, handles)
% hObject    handle to Axis_upperlimits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Axis_upperlimits as text
%        str2double(get(hObject,'String')) returns contents of Axis_upperlimits as a double
Axis_lowerlimits_Callback(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function Axis_upperlimits_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Axis_upperlimits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Delete.
function Delete_Callback(hObject, eventdata, handles)
% hObject    handle to Delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if 1
    index=get(handles.listbox1,'value'); %选中数据的索引
    handles=editdata(hObject, eventdata, handles,index,'Delete');
else
    str=get(handles.listbox1,'string'); %获取列表中的所有数据
    if ischar(str)
        str={str};
    end
    index=get(handles.listbox1,'value'); %选中数据的索引
    data_table=get(handles.table,'data');
    for i=length(index):-1:1
        str(index(i))=[];
        data_table(index(i),:)=[];
    end
    set(handles.listbox1,'string',str); %
    set(handles.table,'data',data_table);
end

% --- Executes on button press in popupmenu13.
function popupmenu13_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function popupmenu13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in Legend.
function Legend_Callback(hObject, eventdata, handles)
% hObject    handle to Legend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Legend contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Legend


% --- Executes during object creation, after setting all properties.
function Legend_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Legend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in LegendFontSize.
function LegendFontSize_Callback(hObject, eventdata, handles)
% hObject    handle to LegendFontSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns LegendFontSize contents as cell array
%        contents{get(hObject,'Value')} returns selected item from LegendFontSize


% --- Executes during object creation, after setting all properties.
function LegendFontSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LegendFontSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Legend_name.
function Legend_name_Callback(hObject, eventdata, handles)
% hObject    handle to Legend_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(handles.table,'visible'),'off')
    set(handles.table,'visible','on');
    return;
    data_table=get(handles.table,'data');
    str=get(handles.listbox1,'String');
    if isempty(str)
        msgbox('listbox is empty')
        set(handles.table,'visible','off');
        set(handles.table,'data',[]);
        return;
    end
    if isempty(data_table{1})
        if isempty(str)
            msgbox('listbox is empty')
            set(handles.table,'visible','off');
            set(handles.table,'data',[]);
        end
        if ischar(str)
            str={str};
        end
        for i=1:length(str)
            temp=strrep(str{i,1},'_',' ');
            temp=temp(1:end-4);
            if temp(end)=='.'
                temp(end)=[];
            end
            str{i,2}=temp;
        end
        
        columnformat = {'char','char'};%这句是关键
        set(handles.table,'columnformat',columnformat,'data',str);
    end
    
else
    data_table=get(handles.table,'data');
    if isempty(data_table{1,1})
        set(handles.table,'visible','off');
    else
        for i=1:length(handles.data.legendname)
            handles.data.legendname{i}=data_table{i,2};
            %         data.data_table{i,1}=data_table{i,1};
            handles.data.data_table{i,2}=data_table{i,2};
        end
        set(handles.table,'visible','off');

        % Choose default command line output for craterdating
        handles.output = hObject;

        % Update handles structure
        guidata(hObject, handles);

        listbox1_Callback(hObject, eventdata, handles);
    end
end



% --- Executes when entered data in editable cell(s) in table.
function table_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to table (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in File.
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns File contents as cell array
%        contents{get(hObject,'Value')} returns selected item from File
data=handles.data;
index=get(handles.File,'value');
set(handles.File,'value',1);
if index==2
    [filename, pathname] = uiputfile('Project.mat', 'Save Project as');
    if ischar(filename)||iscell(filename)
        save([pathname filename],'-struct','data');
    end
elseif index==3
    [filename, pathname] = uigetfile('Project.mat', 'Open Project','MultiSelect', 'off');
    data=load([pathname filename]);
    set(handles.listbox1,'string',data.filename);
    set(handles.listbox1,'value',1:length(data.filename));
    set(handles.table,'data',data.data_table);
    listbox1_Callback(hObject, eventdata, handles);
elseif index==4 %print figure
    if 0
        handles2=handles;
        new_f_handle=figure('visible','off');
        set(new_f_handle,'position',[400,80 600,350]);
        axes_new=axes;
        handles2.axes1=axes_new;
        listbox1_Callback(hObject, eventdata, handles2)
        
    else
        
        new_f_handle=figure('visible','on');
        set(new_f_handle,'position',[400,80 600,400]);
        new_axes=copyobj(handles.axes1,new_f_handle); %picture是GUI界面绘图的坐标系句柄
        [h,object_h]=legend(new_axes,data.legend.plothandle,data.legend.text);
        set(h,'FontSize',11,'box','off')
        set(new_axes,'units','default','position','default');
        p=get(gca,'position');
        p(1)=p(1)-0.02;
        p(2)=p(2)+0.02;
        p(3)=p(3)+0.05;
        p(4)=p(4)-0.02;
        set(gca,'position',p);
        
        index_selected=get(handles.listbox1,'value');
        objhdl= findobj(object_h, 'type', 'line');  %Find objects of legend of type line.
        for i=1:length(objhdl)/2
            set(objhdl((i-1)*2+1),'Marker',data.MarkerType{index_selected(i)},'MarkerEdgeColor',data.MarkerColor{index_selected(i)});
        end
    end
    [filename,pathname]=uiputfile({'*.jpg';'*.eps';'*.pdf';'*.emf'},'save picture as');
    
    set(gcf,'paperpositionmode','auto');%print命令保存eps文件与figure文件显示不同之解决方法
    print(new_f_handle,'-dmeta','-r300',[pathname filename(1:end-4) '.emf'])%保存图片
    print(new_f_handle,'-djpeg100','-r300',[pathname filename(1:end-4) '.jpg'])%保存图片
    
    delete(new_f_handle);
    h=msgbox('Done...');
    pause(0.5);
    close(h);
end

handles.data=data;
% Choose default command line output for craterdating
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function File_CreateFcn(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Batch.
function Batch_Callback(hObject, eventdata, handles)
% hObject    handle to Batch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Reset.
function Reset_Callback(hObject, eventdata, handles)
% hObject    handle to Reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=handles.data;
index=1:length(data.pathname);
if ~isempty(index)
    editdata(hObject, eventdata, handles,index,'Delete');
end
cla


% --- Executes on button press in N1.
function N1_Callback(hObject, eventdata, handles)
% hObject    handle to N1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of N1
if 0
    data=handles.data;
    h=legend(data.legend.plothandle,data.legend.text);
    set(h,'FontSize',11);
else
    listbox1_Callback(hObject, eventdata, handles);
end


% --- Executes on button press in Equilibrium_checkbox.
function Equilibrium_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to Equilibrium_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Equilibrium_checkbox
listbox1_Callback(hObject, eventdata, handles);

% --- Executes on selection change in EquilibriumFunction.
function EquilibriumFunction_Callback(hObject, eventdata, handles)
% hObject    handle to EquilibriumFunction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns EquilibriumFunction contents as cell array
%        contents{get(hObject,'Value')} returns selected item from EquilibriumFunction


% --- Executes during object creation, after setting all properties.
function EquilibriumFunction_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EquilibriumFunction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when plotting is resized.
function plotting_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to plotting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes on selection change in target_type.
function handles=target_type_Callback(hObject, eventdata, handles)
% hObject    handle to target_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns target_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from target_type
type_string=get(handles.target_type,'String');
target_type=get(handles.target_type,'value');

index_selected= get(handles.listbox1,'value'); %index的值代表我们选的是第几个选项
for i=1:length(index_selected)
    handles.data.target_type{index_selected(i)}=type_string{target_type};
end
% Choose default command line output for craterdating
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function target_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to target_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function handles=editdata(varargin)%(isadd,pathname,filename,hObject, eventdata, handles)
%nargin=4;delete,editdata(hObject, eventdata, handles,index_delete)
%nargin=5;add,editdata(hObject, eventdata, handles,pathname,filename)
hObject=varargin{1};
eventdata=varargin{2};
handles=varargin{3};

if ~ischar(varargin{4})
    index_delete=varargin{4};
    type=varargin{5};
elseif nargin==5
    type='Add';
    pathname=varargin{4};
    filename=varargin{5};
end

str=get(handles.listbox1,'string'); %获取列表中已经存在的所有数据
handles.data.getsetting=0;
switch type
    case 'Add'
        filename_added=[str;filename];
        set(handles.listbox1,'string',filename_added);
        
        index= length(str)+1:length(filename_added); %index的值代表我们选的是第几个选项
        set(handles.listbox1,'value',index);
        handles.data.index_selected=length(filename_added);
        
        handles.data.getsetting=1;
    case 'Delete'
        %     index=get(handles.listbox1,'value');
        if length(index_delete)==length(handles.data.pathname)%reset
            set(handles.listbox1,'string',[]);
        else
            set(handles.listbox1,'value',max(1,index_delete-1));
            str(index_delete)=[];
            set(handles.listbox1,'string',str);
        end
end


axes(handles.axes1);
% xlabel('Diameter (km)','Fontsize',15)
% ylabel('N(1) (km^{-2})','Fontsize',15)
switch type
    case {'Add','Update'}
        for i=(1:length(filename))+length(str)
            if strcmp(type,'Add')
                handles.data.pathname{i}=pathname;
                handles.data.filename{i}=filename_added{i,1};
            end
            %-----------------fitting method--------------------
            handles=Fitting_method_Callback(hObject, eventdata, handles);
            handles=Dmin_Callback(hObject, eventdata, handles);
            handles=Dmax_Callback(hObject, eventdata, handles);
            handles=beta_Callback(hObject, eventdata, handles);
            
            %-----------------SFD method--------------------
            handles=SFD_method_Callback(hObject, eventdata, handles);
            handles=PF_Callback(hObject, eventdata, handles);
            handles=ChronologyFunction_Callback(hObject, eventdata, handles);
            %-------------------target_type----------------------------------
            handles=target_type_Callback(hObject, eventdata, handles);
            %-----------------Plotting--------------------
            handles=line_Callback(hObject, eventdata, handles);
            handles=LineWidth_Callback(hObject, eventdata, handles);
            handles=LineColor_Callback(hObject, eventdata, handles);
            
            handles=Marker_Callback(hObject, eventdata, handles);
            handles=MarkerSize_Callback(hObject, eventdata, handles);
            handles=MarkerColor_Callback(hObject, eventdata, handles);
            
            temp=strrep(filename_added{i,1},'_',' ');
            temp=temp(1:end-4);
            if temp(end)=='.'
                temp(end)=[];
            end
            handles.data.legendname{i}=temp;
            
            handles.data.data_table{i,1}=filename_added{i,1};
            handles.data.data_table{i,2}=handles.data.legendname{i};
        end
        set(handles.table,'data',handles.data.data_table);
        handles.data.getsetting=0;
    case 'Delete'
        name=fieldnames(handles.data);
        if length(index_delete)==length(handles.data.pathname)%reset
            for i=length(name):-1:1
                cmd=['handles.data.' name{i} '=[];'];
                eval(cmd);
            end
        else
            index_delete=sort(index_delete,'descend');%delete 3 and then 2,1...
            for i=length(name):-1:1
                if ~strcmp(name{i},'legend')&&~strcmp(name{i},'data_table')&&~strcmp(name{i},'index_selected')
                    for j=1:length(index_delete)
                        cmd=['handles.data.' name{i} sprintf('(%d)',index_delete(j)) '=[];'];
                        eval(cmd);
                    end
                end
            end
            handles.data.data_table(index_delete,:)=[];
        end
        set(handles.table,'data',handles.data.data_table);
end


% Choose default command line output for craterdating
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in pushbutton_fit.
function pushbutton_fit_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(get(handles.listbox1,'string'))
    handles=main(handles);
else
    msgbox('No data is selected','modal');
    return;
end
% Choose default command line output for craterdating
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in Calibration.
function Calibration_Callback(hObject, eventdata, handles)
% hObject    handle to Calibration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
choice = questdlg('Would you have already conduct crater diameter measurements at Apollo 17 landing region using the data in folder ''calibration''?', ...
	'Calibration data', ...
	'Yes','No','');

switch choice
    case 'Yes'
        
    case 'No'
        d=dialog('Position',[300 300 500 150],'Name','Calibration data');
        uicontrol('Parent',d,...
               'Style','text',...
               'Position',[20 80 410 40],...
               'String','Please find the data in the folder ''calibration'' and use them to conduct crater diameter measurements');
        return;   
    otherwise
         return;   
end
if ~isfield(handles,'data')
    handles.data=[];
end

[filenames, pathname] = uigetfile( ...
    {'*.scc;*.shp;*.dbf;*.csfd;*.txt','All Image Files'},...
    'Please select the data (multiple options)', '.\data', ...
    'MultiSelect', 'on');
if ~(ischar(filenames)||iscell(filenames))
    return;
end
reference_datafile='calibration\CRATER_Apollo17.dbf';
if ischar(filenames)
    filenames={reference_datafile,filenames};
elseif length(filenames)>1
    filenames=filenames';
    filenames={reference_datafile,filenames{:}};
end
for i=2:length(filenames)
    filenames{i}=[pathname filenames{i}];
end
[systematic_error_diameter,precision_diameter]=calibration(filenames);
handles.data.calibration.systematic_error_diameter=systematic_error_diameter;
handles.data.calibration.precision_diameter=precision_diameter;

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function AMA1_Callback(hObject, eventdata, handles)
% hObject    handle to AMA1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AMA1 as text
%        str2double(get(hObject,'String')) returns contents of AMA1 as a double


% --- Executes during object creation, after setting all properties.
function AMA1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AMA1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function AMA1_uncertainty_Callback(hObject, eventdata, handles)
% hObject    handle to AMA1_uncertainty (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AMA1_uncertainty as text
%        str2double(get(hObject,'String')) returns contents of AMA1_uncertainty as a double


% --- Executes during object creation, after setting all properties.
function AMA1_uncertainty_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AMA1_uncertainty (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function AMA2_Callback(hObject, eventdata, handles)
% hObject    handle to AMA2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AMA2 as text
%        str2double(get(hObject,'String')) returns contents of AMA2 as a double


% --- Executes during object creation, after setting all properties.
function AMA2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AMA2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function AMA2_uncertainty_Callback(hObject, eventdata, handles)
% hObject    handle to AMA2_uncertainty (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AMA2_uncertainty as text
%        str2double(get(hObject,'String')) returns contents of AMA2_uncertainty as a double


% --- Executes during object creation, after setting all properties.
function AMA2_uncertainty_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AMA2_uncertainty (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function craterNumber_Callback(hObject, eventdata, handles)
% hObject    handle to craterNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of craterNumber as text
%        str2double(get(hObject,'String')) returns contents of craterNumber as a double

handles.data.comparison.craterNumber=str2double(get(handles.craterNumber,'string'));
% Choose default command line output for craterdating
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function craterNumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to craterNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function density_uncertainty_pos_Callback(hObject, eventdata, handles)
% hObject    handle to density_uncertainty_pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of density_uncertainty_pos as text
%        str2double(get(hObject,'String')) returns contents of density_uncertainty_pos as a double
index= get(handles.listbox1,'value'); %index的值代表我们选的是第几个选项
for i=1:length(index)
    handles.data.comparison.density_uncertainty_pos{index(i)}=str2double(get(hObject,'string'));
end
% Choose default command line output for craterdating
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function density_uncertainty_pos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to density_uncertainty_pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function possibility_of_being_the_same_Callback(hObject, eventdata, handles)
% hObject    handle to possibility_of_being_the_same (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of possibility_of_being_the_same as text
%        str2double(get(hObject,'String')) returns contents of possibility_of_being_the_same as a double

index= get(handles.listbox1,'value'); %index的值代表我们选的是第几个选项
if length(index)==1
    
    density_comparison=handles.data.comparison.density(index);
    density_uncertainty_pos_comparison=handles.data.comparison.density_uncertainty_pos(index);
    density_uncertainty_neg_comparison=handles.data.comparison.density_uncertainty_neg(index);
    possibility_of_being_the_same=0;
else
    set(hObject,'string','N/A');
end

% --- Executes during object creation, after setting all properties.
function possibility_of_being_the_same_CreateFcn(hObject, eventdata, handles)
% hObject    handle to possibility_of_being_the_same (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function density_uncertainty_neg_Callback(hObject, eventdata, handles)
% hObject    handle to density_uncertainty_neg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of density_uncertainty_neg as text
%        str2double(get(hObject,'String')) returns contents of density_uncertainty_neg as a double
index= get(handles.listbox1,'value'); %index的值代表我们选的是第几个选项
for i=1:length(index)
    handles.data.comparison.density_uncertainty_neg{index(i)}=str2double(get(hObject,'string'));
end
% Choose default command line output for craterdating
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function density_uncertainty_neg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to density_uncertainty_neg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function handles=beta_Callback(hObject, eventdata, handles)
% hObject    handle to beta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of beta as text
%        str2double(get(hObject,'String')) returns contents of beta as a double
index= get(handles.listbox1,'value'); %index的值代表我们选的是第几个选项
for i=1:length(index)
    handles.data.systematic_error.beta{index(i)}=str2double(get(handles.beta,'string'));
end
% Choose default command line output for craterdating
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function beta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to beta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit22 as text
%        str2double(get(hObject,'String')) returns contents of edit22 as a double
handles.data.comparison.craterNumber2=str2double(get(handles.craterNumber2,'string'));
% Choose default command line output for craterdating
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
