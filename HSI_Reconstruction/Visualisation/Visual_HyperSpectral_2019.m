function varargout = Visual_HyperSpectral_2019(varargin)
% VISUAL_INTERFACE MATLAB code for Visual_HyperSpectral_2019.fig
%      VISUAL_INTERFACE, by itself, creates a new VISUAL_INTERFACE or raises the existing
%      singleton*.
%
%      H = VISUAL_INTERFACE returns the handle to a new VISUAL_INTERFACE or the handle to
%      the existing singleton*.
%
%      VISUAL_INTERFACE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VISUAL_INTERFACE.M with the given input arguments.
%
%      VISUAL_INTERFACE('Property','Value',...) creates a new VISUAL_INTERFACE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Visual_HyperSpectral_2019_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Visual_HyperSpectral_2019_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Visual_HyperSpectral_2019

% Last Modified by GUIDE v2.5 20-Mar-2018 16:42:50

global Gamma; 
Gamma = .2; % Pour augmentation contraste RGB
% Add functions PATH
addpath Visualisation/VHS_functions

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Visual_HyperSpectral_2019_OpeningFcn, ...
                   'gui_OutputFcn',  @Visual_HyperSpectral_2019_OutputFcn, ...
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


% --- Executes just before Visual_HyperSpectral_2019 is made visible.
function Visual_HyperSpectral_2019_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Visual_HyperSpectral_2019 (see VARARGIN)


% Choose default command line output for Visual_HyperSpectral_2019
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Visual_HyperSpectral_2019 wait for user response (see UIRESUME)
% uiwait(handles.figure1);
% set(handles.axes_RGB,'xtick',[],'ytick',[])
% set(handles.axes_panchro,'xtick',[],'ytick',[])
% set(handles.axes_pl,'xtick',[],'ytick',[])
% set(handles.axes_xy,'xtick',[],'ytick',[])
% set(handles.axes_xl,'xtick',[],'ytick',[])
% set(handles.axes_yl,'xtick',[],'ytick',[])
set(handles.axes_visible,'xtick',[],'ytick',[]);
handles.pushbutton_pixelchoice.Visible='off';
handles.pushbutton_wavelengthchoice.Visible='off';
handles.pushbutton_zoom.Visible='off';
handles.axes_pl.Visible='off';
handles.axes_xl.Visible='off';
handles.axes_yl.Visible='off';
handles.axes_xy.Visible='off';
handles.axes_panchro.Visible='off';
handles.axes_RGB.Visible='off';
handles.slider_W.Visible='off';
handles.text2.Visible='off';
handles.edit_W.Visible='off';
handles.axes_visible.Visible='off';
%keyboard

% --- Outputs from this function are returned to the command line.
function varargout = Visual_HyperSpectral_2019_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_Load.
function pushbutton_Load_Callback(hObject, eventdata, handles)
global Gamma; 
% keyboard
% hObject    handle to pushbutton_Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[Object_filename, pathname] = uigetfile('*.mat', 'Select a file');
if isequal(Object_filename,0)
    disp('User selected Cancel')
    handles.pushbutton_pixelchoice.Visible='off';
    handles.pushbutton_wavelengthchoice.Visible='off';
    handles.pushbutton_zoom.Visible='off';
    handles.axes_pl.Visible='off';
    handles.axes_xl.Visible='off';
    handles.axes_yl.Visible='off';
    handles.axes_xy.Visible='off';
    handles.axes_panchro.Visible='off';
    handles.axes_RGB.Visible='off';
    handles.slider_W.Visible='off';
    handles.text2.Visible='off';
    handles.edit_W.Visible='off';
    handles.axes_visible.Visible='off';
else
    cla(handles.axes_RGB)
    cla(handles.axes_panchro)
    cla(handles.axes_xy)
    cla(handles.axes_xl)
    cla(handles.axes_yl)
    cla(handles.axes_pl)
    pathname = pathname(1:end-1);
    str = fullfile(pathname);
    addpath(str)
    load(Object_filename); 
    Object_init = Object;
    panchro_init = sum(Object_init,3);
    RGB = hypercube_to_RGB(Object_init);
    axes(handles.axes_panchro)
    hold off
    imagesc(panchro_init.^Gamma)
    axis on
    xlabel('x'); ylabel('y')
    title('panchromatic ')
    axes(handles.axes_RGB)
    hold off
    imagesc(RGB.^Gamma)
    axis on
    xlabel('x'); ylabel('y')
    title('RGB')
    setappdata(0,'Object_init',Object_init);
    setappdata(0,'Object_cube',Object_init);
    setappdata(0,'RGB_init',RGB);
    setappdata(0,'RGB',RGB);
    handles.pushbutton_pixelchoice.Visible='on';
    handles.pushbutton_zoom.Visible='on';
    handles.pushbutton_wavelengthchoice.Visible='off';
    handles.axes_pl.Visible='off';
    handles.axes_xl.Visible='off';
    handles.axes_yl.Visible='off';
    handles.axes_xy.Visible='off';
    handles.axes_panchro.Visible='off';
    handles.slider_W.Visible='off';
    handles.text2.Visible='off';
    handles.edit_W.Visible='off';
    handles.axes_visible.Visible='off';    
end

% --- Executes on button press in pushbutton_zoom.
function pushbutton_zoom_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Gamma; 

Object_cube = getappdata(0,'Object_init');
RGB = getappdata(0,'RGB_init');
panchro = sum(Object_cube,3);
axes(handles.axes_panchro)
hold off
imagesc(panchro.^Gamma);
axis on
xlabel('x'); ylabel('y')
title('panchromatic ')
axes(handles.axes_RGB)
hold off
%imshow(RGB,[]);
imagesc(RGB.^Gamma)
axis on
xlabel('x'); ylabel('y')
title('RGB ')

rect = getrec_RV;
rect = round(rect);
[nr,nc] = size(panchro);
if rect(3)>5 && rect(4)>5
    cla(handles.axes_xy)
    cla(handles.axes_xl)
    cla(handles.axes_yl)
    cla(handles.axes_pl)
    X = [rect(1) rect(1) + rect(3) ];
    Y = [rect(2) rect(2) + rect(4) ];
 
    for i=1:2
        if X(i)<1
            X(i) = 1;
        elseif X(i)>nc
            X(i) = nc;
        end
        if Y(i)<1
            Y(i) = 1;
        elseif Y(i)>nr
            Y(i) = nr;
        end
    end
    Object_cube = Object_cube(Y(1): Y(2),X(1):X(2),:);
    axes(handles.axes_panchro)
    hold off
    imagesc(sum(Object_cube,3).^Gamma);
    axis on
    xlabel('x'); ylabel('y')
    title('panchromatic')
    RGB = RGB(Y(1): Y(2),X(1):X(2),:);
    axes(handles.axes_RGB)
    hold off
    imagesc(RGB.^Gamma)
    axis on
    xlabel('x'); ylabel('y')
    title('RGB')
    setappdata(0,'Object_cube',Object_cube)
    setappdata(0,'RGB',RGB)
else
    setappdata(0,'Object_cube',Object_cube)
    setappdata(0,'RGB',RGB)
end


% --- Executes on button press in pushbutton_pixelchoice.
function pushbutton_pixelchoice_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_pixelchoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Gamma; 

Object_cube = getappdata(0,'Object_cube');
axes(handles.axes_panchro)
hold off
imagesc(sum(Object_cube,3).^Gamma)
title('panchromatic ')
axis on
xlabel('x'); ylabel('y')
RGB = getappdata(0,'RGB');
axes(handles.axes_RGB)
hold off
imagesc(RGB.^Gamma)
axis on
xlabel('x'); ylabel('y')
title('RGB ')

[c,r] = ginput(1);
c = floor(c); r = floor(r);
cla(handles.axes_xl)
cla(handles.axes_yl)
cla(handles.axes_pl)

axes(handles.axes_pl)
plot(squeeze(Object_cube(r,c,:))); axis tight
title(['x = ' num2str(c) '      y = ' num2str(r)])
%xlabel('W')

axes(handles.axes_visible)
vis=imread('Visible.jpg');
imagesc(vis); axis off
%xlabel('W')

axes(handles.axes_panchro)
hold on
plot([1 size(Object_cube,2)],[r r], 'red')
plot([c c],[1 size(Object_cube,1)], 'red')

axes(handles.axes_RGB)
hold on
plot([1 size(Object_cube,2)],[r r], 'red')
plot([c c],[1 size(Object_cube,1)], 'red')

axes(handles.axes_xl)
A = Object_cube(r,:,:);
A = reshape(A, size(A,2), size(A,3));
imagesc(A'); colorbar
title(['y = ' num2str(r)])
xlabel('x')
ylabel('W')

axes(handles.axes_yl)
A = Object_cube(:,c,:);
A = reshape(A, size(A,1), size(A,3));
imagesc(A); colorbar
title(['x = ' num2str(c)])
xlabel('W')
ylabel('y')

setappdata(0,'r',r)
setappdata(0,'c',c)
handles.pushbutton_wavelengthchoice.Visible='on';
handles.slider_W.Visible='on';
handles.text2.Visible='on';
handles.edit_W.Visible='on';



% --- Executes on button press in pushbutton_wavelengthchoice.
function pushbutton_wavelengthchoice_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_wavelengthchoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Gamma; 
Object_cube = getappdata(0,'Object_cube');
[val_lig,val_col,~,h0] = ginput_HC(1);
h0_title = h0.Title.Text.String;
%keyboard
W=[];
if ~isempty(h0_title)
    if h0_title(1)=='y'  
        W = val_col;
    elseif h0_title(1)=='x'  
        W = val_lig;
    end
end    
if ~isempty(W),
    W = floor(W);
    cla(handles.axes_xy)
    set(handles.slider_W,'Value',W);
    set(handles.edit_W,'String',W);
    handles.edit_W.Visible='on';
    axes(handles.axes_xy)
    imagesc(Object_cube(:,:,W).^Gamma); % colorbar
    axis on
    title(['W = ' num2str(W)])
    xlabel('x')
    ylabel('y')
end






function edit_W_Callback(hObject, eventdata, handles)
% hObject    handle to edit_W (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_W as text
%        str2double(get(hObject,'String')) returns contents of edit_W as a double
global Gamma; 
Object_cube = getappdata(0,'Object_cube');
W = str2double(get(hObject,'String'));
if W==floor(W) && W>=1 && W<=size(Object_cube,3)
  cla(handles.axes_xy)
  set(handles.slider_W,'Value',W);
  axes(handles.axes_xy)
  imagesc(Object_cube(:,:,W).^Gamma); % colorbar
  title(['W = ' num2str(W)])
  xlabel('x')
  ylabel('y')
else
  set(handles.edit_W,'String',[])
end


% --- Executes during object creation, after setting all properties.
function edit_W_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_W (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on slider movement.
function slider_W_Callback(hObject, eventdata, handles)
% hObject    handle to slider_W (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global Gamma
Object_cube = getappdata(0,'Object_cube');
W = floor(get(hObject,'Value'));
cla(handles.axes_xy)
set(handles.edit_W,'String',W);
handles.edit_W.Visible='on';
axes(handles.axes_xy)
imagesc(Object_cube(:,:,W).^Gamma); % colorbar
axis on
title(['W = ' num2str(W)])
xlabel('x')
ylabel('y')

% --- Executes during object creation, after setting all properties.
function slider_W_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_W (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
