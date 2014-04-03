function varargout = SpectrumGUI(varargin)
% SPECTRUMGUI MATLAB code for SpectrumGUI.fig
%      SPECTRUMGUI, by itself, creates a new SPECTRUMGUI or raises the existing
%      singleton*.
%
%      H = SPECTRUMGUI returns the handle to a new SPECTRUMGUI or the handle to
%      the existing singleton*.
%
%      SPECTRUMGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPECTRUMGUI.M with the given input arguments.
%
%      SPECTRUMGUI('Property','Value',...) creates a new SPECTRUMGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SpectrumGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SpectrumGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SpectrumGUI

% Last Modified by GUIDE v2.5 03-Apr-2014 12:26:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @SpectrumGUI_OpeningFcn, ...
    'gui_OutputFcn',  @SpectrumGUI_OutputFcn, ...
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


% --- Executes just before SpectrumGUI is made visible.
function SpectrumGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SpectrumGUI (see VARARGIN)

% Choose default command line output for SpectrumGUI
handles.output = hObject;


handles.N = str2double(get(handles.edit3, 'string'));
handles.M = str2double(get(handles.edit6, 'string'));
handles.L = str2double(get(handles.edit7, 'string'));

handles.meshSize = [handles.N handles.M handles.L];


set(handles.text8, 'string', num2str(2*handles.N*handles.M*handles.L));

val = get(handles.popupmenu1, 'value');
str = get(handles.popupmenu1, 'string');
handles.scheme = str{val};

if handles.scheme ~= 'SLF'
    set(handles.popupmenu2, 'disabled');
end

handles.BC = get(handles.popupmenu2, 'value');

handles.R = str2double(get(handles.edit1, 'string'));
handles.pct_lambda = str2double(get(handles.edit2, 'string'));

if get(handles.popupmenu3, 'value') == 2
    handles.single_precision = true;
else
    handles.single_precision = false;
end

% empties become NaN
handles.hs_pos = [str2double(get(handles.edit8, 'string'))...
    str2double(get(handles.edit9, 'string')) str2double(get(handles.edit10, 'string'))];

% kludge for windows fontsize...
if ispc
    ch = get(hObject, 'children');
    for m = 1:length(ch)
%         disp(m)
        if m ~= 5 % sub-kludge
            set(ch(m), 'fontsize',7, 'fontname', 'Arial');
        end
    end
end

% default z-scale for plotting eigenvectors and geometry
handles.zscale = 4;


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SpectrumGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SpectrumGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

val = get(handles.popupmenu1, 'value');
str = get(handles.popupmenu1, 'string');
handles.scheme = str{val};

if strcmp(handles.scheme, 'SLF')
    set(handles.popupmenu2, 'enable', 'on');
else
    set(handles.popupmenu2, 'enable', 'off', 'value', 1);
end

set(handles.pushbutton1, 'ForegroundColor', [0 0.5 0]);
set(handles.status, 'string', 'Recompute eigenvalues...');

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2

handles.BC = get(handles.popupmenu2, 'value');

set(handles.pushbutton1, 'ForegroundColor', [0 0.5 0]);
set(handles.status, 'string', 'Recompute eigenvalues...');

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

handles.R = str2double(get(handles.edit1, 'string'));

set(handles.pushbutton1, 'ForegroundColor', [0 0.5 0]); 
set(handles.status, 'string', 'Recompute eigenvalues...');

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double

handles.pct_lambda = str2double(get(handles.edit2, 'string'));

set(handles.pushbutton1, 'ForegroundColor', [0 0.5 0]); 
set(handles.status, 'string', 'Recompute eigenvalues...');

guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.status, 'string', 'Setting up matrix...')
drawnow;

switch handles.scheme
    
    case 'IWB'
        sig = handles.pct_lambda;
        b = 1/4;
        c = 1/16;
        
    case 'OCTA'
        sig = handles.pct_lambda;
        b = 1/2;
        c = 1/4;
        
    case 'CCP'
        sig = handles.pct_lambda;
        b = 1/4;
        c = 0;
        
    case 'SLF'
        sig = handles.pct_lambda  / sqrt(3);
        b = 0;
        c = 0;
        
    case 'IDWM'
        sig = handles.pct_lambda / sqrt(3);
        b = 0.2034;
        c = 0.0438;
        
    case 'IISO'
        sig = handles.pct_lambda * sqrt(3/4);
        b = 1/6;
        c = 0;
        
    case 'IISO2'
        sig = handles.pct_lambda * sqrt(3/4);
        b = 1/6;
        c = 1/48;
        
end

s2 = sig*sig;
bs2 = s2; %


handles.sig = sig;
handles.bsig = sqrt(bs2);

use_sparse = [];


handles.A = set_up_matrix(b, c, s2, bs2, handles.hs_pos, handles.BC, handles.R, handles.single_precision, handles);

tic;

set(handles.status, 'string', 'Computing eigenvalues...');
drawnow;

[handles.eigVects, handles.eigVals] = eig(handles.A);
tt = toc;
% disp([num2str(tt) ' sec to compute eigenvalues and eigenvectors'])

% a rounding threshold for plotting
% we color eigenvalues outside this cushion red
handles.roundingThresh = norm(handles.A*handles.eigVects - handles.eigVects*handles.eigVals);

set(handles.pushbutton2, 'enable','on');
set(handles.pushbutton3, 'enable','on');
if exist('eigtoollib') == 7
    set(handles.pushbutton4, 'enable', 'on');
end

set(handles.pushbutton5, 'enable','on');
set(handles.pushbutton6, 'enable','on');
set(handles.pushbutton7, 'enable','on');
set(handles.pushbutton9, 'enable','on');


set(handles.pushbutton1, 'ForegroundColor', 'k');

set(handles.status, 'string', 'Ready for plotting...');
drawnow;

guidata(hObject, handles);



% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.status, 'string', 'Making figure...')
drawnow;

if size(handles.eigVals, 2) == 1
    handles.eigVals = diag(handles.eigVals);
end



% handles.roundingThresh = norm(handles.A*handles.eigVects - handles.eigVects*handles.eigVals);
% roundingThresh = 1e-8;

handles.eigVals = diag(handles.eigVals);

usInds = find(abs(handles.eigVals) > 1 + handles.roundingThresh);

fh = figure;
hold all;


circle_angles = [0:0.001:2*pi 0];

str = get(handles.popupmenu4, 'string');
val = get(handles.popupmenu4, 'value');

if strcmp(get(handles.normal, 'checked'), 'on')
    
    plot(cos(circle_angles), sin(circle_angles), 'color', 0.8*[1 1 1], 'linewidth', 1);
    
    line([-1 1],[0 0], 'color', [0.8 0.8 0.8], 'linewidth', 1);
    line([0 0], [-1 1], 'color', [0.8 0.8 0.8], 'linewidth', 1);
    
    if strcmp(str(val), 'Box')
        r1 = 0.75;
        plot(r1*cos(circle_angles), r1*sin(circle_angles), 'color', 0.8*[1 1 1], 'linewidth', 1);
        PlotFrequencies(fh, handles, r1);
    end
    
else
    
    fill(cos(circle_angles), sin(circle_angles),[0.8 0.8 0.8], 'edgecolor', 'none');
    
    line([-1 1],[0 0], 'color', [1 1 1], 'linewidth', 1);
    line([0 0], [-1 1], 'color', [1 1 1], 'linewidth', 1);
    
    if strcmp(str(val), 'Box')
        r1 = 0.75;
        plot(r1*cos(circle_angles), r1*sin(circle_angles), 'color', 1*[1 1 1], 'linewidth', 1);
        PlotFrequenciesWhite(fh, handles, r1);
    end
    
end

plot(real(handles.eigVals), imag(handles.eigVals), 'ko', 'markersize', 3, 'markerfacecolor', 'k', 'markeredgecolor', 'none');
plot(real(handles.eigVals(usInds)), imag(handles.eigVals(usInds)), 'ko', 'markersize', 3, 'markerfacecolor', 'k', 'markeredgecolor', 'r');


hold off;

set(gca, 'fontsize', 10, 'fontname', 'Times New Roman');
axis 'equal';
ylim([-1.1 1.1]);
xlim([-1.1 1.1]);

box on;

set(handles.status, 'string', 'Ready for plotting...')
drawnow;



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double

handles.N = str2double(get(handles.edit3, 'string'));
handles.meshSize = [handles.N handles.M handles.L];

if get(handles.popupmenu4, 'value') == 1
    set(handles.text8, 'string', num2str(2*handles.N*handles.M*handles.L));
else
    set(handles.text8, 'string', num2str( 2*(handles.N+2)*(handles.M+2)*(handles.L+2)) );
end

set(handles.pushbutton1, 'ForegroundColor', [0 0.5 0]); 
set(handles.status, 'string', 'Recompute eigenvalues...');

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double

handles.M = str2double(get(handles.edit6, 'string'));
handles.meshSize = [handles.N handles.M handles.L];

if get(handles.popupmenu4, 'value') == 1
    set(handles.text8, 'string', num2str(2*handles.N*handles.M*handles.L));
else
    set(handles.text8, 'string', num2str( 2*(handles.N+2)*(handles.M+2)*(handles.L+2)) );
end

set(handles.pushbutton1, 'ForegroundColor', [0 0.5 0]); 
set(handles.status, 'string', 'Recompute eigenvalues...');

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double

handles.L = str2double(get(handles.edit7, 'string'));
handles.meshSize = [handles.N handles.M handles.L];

if get(handles.popupmenu4, 'value') == 1
    set(handles.text8, 'string', num2str(2*handles.N*handles.M*handles.L));
else
    set(handles.text8, 'string', num2str( 2*(handles.N+2)*(handles.M+2)*(handles.L+2)) );
end

set(handles.pushbutton1, 'ForegroundColor', [0 0.5 0]); 
set(handles.status, 'string', 'Recompute eigenvalues...');

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


set(handles.status, 'string', 'Making figure...')
drawnow;

if size(handles.eigVals, 2) > 1
    handles.eigVals = diag(handles.eigVals);
end

usInds = find(abs(handles.eigVals) > 1 + handles.roundingThresh);
sInds = find(abs(handles.eigVals) <= 1 + handles.roundingThresh);

circle_angles = [0:0.0001:2*pi 0];

figure;
plot(cos(circle_angles), sin(circle_angles), 'color', 0.8*[1 1 1], 'linewidth', 1);
hold all;
plot(real(handles.eigVals(sInds)), imag(handles.eigVals(sInds)), 'ko', 'markersize', 6, 'markerfacecolor', 'none', 'markeredgecolor', 'k');
plot(real(handles.eigVals(usInds)), imag(handles.eigVals(usInds)), 'ko', 'markersize', 6, 'markerfacecolor', 'none', 'markeredgecolor', 'r');
title({'Click on eigenvalues to plot corresponding eigenvectors',...
    'It will find the closest eigenvalue to your mouse click'})
xlim([-1.1 1.1])
ylim([-1.1 1.1])
axis equal;


if strcmp(get(handles.slices, 'checked'), 'off')
    
    % now attach the function to the axes
    set(gca,'ButtonDownFcn', {@mouseclick_callback, handles})
    
    % and we also have to attach the function to the children, in this
    % case that is the line in the axes.
    set(get(gca,'Children'),'ButtonDownFcn', {@mouseclick_callback, handles})
    
else
    
    % now attach the function to the axes
    
    set(gca,'ButtonDownFcn', {@alt_mouseclick_callback, handles})
    
    % and we also have to attach the function to the children, in this
    % case that is the line in the axes.
    set(get(gca,'Children'),'ButtonDownFcn', {@alt_mouseclick_callback, handles})
    
    
end

set(handles.status, 'string', 'Ready for plotting...')
drawnow;


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


cd eigtoollib
set_eigtool_prefs('default')
opts.ax = [-1.5 1.5 -1.5 1.5];
opts.unit_circle = 1;
eigtool(handles.A,opts);
set(gcf,'DeleteFcn',@cd_on_delete);


function cd_on_delete(src,evnt)
cd ../

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.status, 'string', 'Making figure...')
drawnow;

if size(handles.eigVals, 2) > 1
    handles.eigVals = diag(handles.eigVals);
end

kappa = condeig(handles.A);
ang = angle(handles.eigVals);
kappa = kappa(ang >= 0); % only the non-negative frequencies
ang = ang(ang >= 0);


figure;
if handles.single_precision
    plot(ang, log10(0.5*eps(single(1))*(kappa)), 'ok');
else
    plot(ang, log10(0.5*eps*(kappa)), 'xk', 'markersize', 8);
end
set(gca, 'xtick', [0 pi/4 pi/2 3*pi/4 pi], 'xticklabel', {'0' , '\pi/4', '\pi/2', '3\pi/4', '\pi'})
xlim([-0.05 pi+0.05]);
xlabel('\omega k');
ylabel('log_{10} [\epsilon_M/2 \kappa(\lambda, A)]');
title('Estimated Perturbations', 'fontweight', 'bold');

set(handles.status, 'string', 'Ready for plotting...')
drawnow;

guidata(hObject, handles);

function mask = GetBoundaryMask(handles)

str = get(handles.popupmenu4, 'string');
val = get(handles.popupmenu4, 'value');

if strcmp(str(val), 'Box')
    
    mask = zeros(handles.N+2, handles.M+2, handles.L+2);
    mask(2:handles.N+1, 2:handles.M+1, 2:handles.L+1) = 1;
    
    handles.meshSize = [handles.N handles.M handles.L];
    
else
    
    handles.scheme = 'SLF';
    
    if strcmp(str(val), 'L1 (unstable)')
        
        N = handles.meshSize;
        
        mask = zeros(N(1), N(2), N(3));
        mask(2:N(1)-1, 2:N(2)-1, 2:N(3)-1 ) = 1;
        
        mask(1:5, 1:5, :) = 0;
        
        
    elseif strcmp(str(val), 'L2 (stable)' )
        
        N = handles.meshSize;
        
        mask = zeros(N(1), N(2), N(3));
        mask(2:N(1)-1, 2:N(2)-1, 2:N(3)-1 ) = 1;
        
        mask(1:6, 1:5, :) = 0;
        
    elseif strcmp(str(val), 'Corner' )
        
        N = handles.meshSize;
        
        mask = zeros(N(1), N(2), N(3));
        mask(2:N(1)-1, 2:N(2)-1, 2:N(3)-1 ) = 1;
        
        mask(1:4, 1:4, 1:4) = 0;
        
    elseif strcmp(str(val), 'Cylinder' )
        
        N = handles.meshSize;
        
        mask = zeros(N(1), N(2), N(3));
        mask(2:N(1)-1, 2:N(2)-1, 2:N(3)-1 ) = 1;
        
        [X, Y, Z] = meshgrid(linspace(-1, 1, N(1)),...
            linspace(-1, 1, N(2)), linspace(-1, 1, N(3)));
        mask( X.^2 + Y.^2 > 1^2 ) = 0;
        mask(:, :, 1) = 0;
        mask(:, :, N(3)) = 0;
        
    elseif strcmp(str(val), 'Sphere' )
        
        N = handles.meshSize;
        
        mask = zeros(N(1), N(2), N(3));
        mask(2:N(1)-1, 2:N(2)-1, 2:N(3)-1 ) = 1;
        
        [X, Y, Z] = meshgrid(linspace(-1, 1, N(1)),...
            linspace(-1, 1, N(2)), linspace(-1, 1, N(3)));
        mask( X.^2 + Y.^2 + Z.^2 > 1^2 ) = 0;
        mask(:, :, 1) = 0;
        mask(:, :, N(3)) = 0;
        
    end
end

function A = set_up_matrix(b, c, s2, bs2, HSPos, boundaries, reflection_coefficient, single_precision, handles)

sig = sqrt(s2);
bsig = sqrt(bs2);

beta = (reflection_coefficient+1)/(1-reflection_coefficient);

d1 = s2*(1 - 4*b + 4*c);
d2 = s2*(b - 2*c);
d3 = s2*c;
d4 = 2*( 1 - 3*s2 + 6*s2*b - 4*c*s2 );

% have to check that d4 isn't artificially non-zero
if abs(d4) > 0 && abs(d4) < 1e-6
    d4 = 0;
end

BM = GetBoundaryMask(handles);

s = size(BM);

%         if use_sparse
%             A1 = sparse(prod(s));
%             A2 = A1;
%         else
str = get(handles.popupmenu4, 'string');
val = get(handles.popupmenu4, 'value');


if strcmp(str(val), 'Box')
    s = s-2;
    A1 = zeros(prod(s));
    A2 = A1;
else
    A1 = zeros(prod(s));
    A2 = A1;
end

% tic
for n = 1:prod(s)
    
    [i, j, k] = ind2sub(s, n);
    
    if strcmp(str(val), 'Box')
        [nodeType, nodeOrientation] = GetNodeTypeBox(i, j, k, s);
    else
        [nodeType, nodeOrientation] = GetNodeType(i, j, k, BM);
    end
    
    
    switch nodeType
        
        case 'interior'
            
            %%%% 26 nearest neighbors %%%%
            n_ip1 = sub2ind(s, i+1, j, k);
            n_im1 = sub2ind(s, i-1, j, k);
            n_jp1 = sub2ind(s, i, j+1, k);
            n_jm1 = sub2ind(s, i, j-1, k);
            n_kp1 = sub2ind(s, i, j, k+1);
            n_km1 = sub2ind(s, i, j, k-1);
            
            n_ip1_jp1 = sub2ind(s, i+1, j+1, k);
            n_im1_jp1 = sub2ind(s, i-1, j+1, k);
            n_ip1_jm1 = sub2ind(s, i+1, j-1, k);
            n_im1_jm1 = sub2ind(s, i-1, j-1, k);
            
            n_jp1_kp1 = sub2ind(s, i, j+1, k+1);
            n_jm1_kp1 = sub2ind(s, i, j-1, k+1);
            n_jp1_km1 = sub2ind(s, i, j+1, k-1);
            n_jm1_km1 = sub2ind(s, i, j-1, k-1);
            
            n_ip1_kp1 = sub2ind(s, i+1, j, k+1);
            n_im1_kp1 = sub2ind(s, i-1, j, k+1);
            n_ip1_km1 = sub2ind(s, i+1, j, k-1);
            n_im1_km1 = sub2ind(s, i-1, j, k-1);
            
            n_ip1_jp1_kp1 = sub2ind(s, i+1, j+1, k+1);
            n_ip1_jm1_kp1 = sub2ind(s, i+1, j-1, k+1);
            n_im1_jp1_kp1 = sub2ind(s, i-1, j+1, k+1);
            n_im1_jm1_kp1 = sub2ind(s, i-1, j-1, k+1);
            
            n_ip1_jp1_km1 = sub2ind(s, i+1, j+1, k-1);
            n_ip1_jm1_km1 = sub2ind(s, i+1, j-1, k-1);
            n_im1_jp1_km1 = sub2ind(s, i-1, j+1, k-1);
            n_im1_jm1_km1 = sub2ind(s, i-1, j-1, k-1);
            
            %%%%%%% UPDATE %%%%%%%
            A1(n, n) = d4;
            A2(n, n) = -1;
            
            A1(n, n_ip1) = d1;
            A1(n, n_im1) = d1;
            A1(n, n_jp1) = d1;
            A1(n, n_jm1) = d1;
            A1(n, n_kp1) = d1;
            A1(n, n_km1) = d1;
            
            A1(n, n_ip1_jp1) = d2;
            A1(n, n_ip1_jm1) = d2;
            A1(n, n_ip1_kp1) = d2;
            A1(n, n_ip1_km1) = d2;
            
            A1(n, n_jp1_kp1) = d2;
            A1(n, n_jm1_kp1) = d2;
            A1(n, n_jp1_km1) = d2;
            A1(n, n_jm1_km1) = d2;
            
            A1(n, n_im1_jp1) = d2;
            A1(n, n_im1_jm1) = d2;
            A1(n, n_im1_kp1) = d2;
            A1(n, n_im1_km1) = d2;
            
            A1(n, n_ip1_jp1_kp1) = d3;
            A1(n, n_ip1_jm1_kp1) = d3;
            A1(n, n_ip1_jp1_km1) = d3;
            A1(n, n_ip1_jm1_km1) = d3;
            
            A1(n, n_im1_jp1_kp1) = d3;
            A1(n, n_im1_jm1_kp1) = d3;
            A1(n, n_im1_jp1_km1) = d3;
            A1(n, n_im1_jm1_km1) = d3;
        
    
    case 'face'
        if beta == 0
            A1(n, :) = 0;
            A2(n, :) = 0;
        elseif boundaries == 1 % proper boundaries
            clear n_in_axial n_in_corners n_out n_out_axial n_out_corners
            
            % central node outside the plane
            iii = i+nodeOrientation(1);
            jjj = j+nodeOrientation(2);
            kkk = k+nodeOrientation(3);
            
            n_out = sub2ind(s, iii, jjj, kkk);
            
            % unit vectors along the surface
            ii = find( nodeOrientation == 0 );
            sV1 = [0 0 0]; sV2 = sV1; % surface vectors
            sV1(ii(1)) = 1; % sign doesn't matter here since we +/-
            sV2(ii(2)) = 1;
            
            n_in_axial(1) = sub2ind(s, i+sV1(1), j+sV1(2), k+sV1(3));
            n_in_axial(2) = sub2ind(s, i-sV1(1), j-sV1(2), k-sV1(3));
            n_in_axial(3) = sub2ind(s, i+sV2(1), j+sV2(2), k+sV2(3));
            n_in_axial(4) = sub2ind(s, i-sV2(1), j-sV2(2), k-sV2(3));
            
            n_in_corners(1) = sub2ind(s, i+sV1(1)+sV2(1), j+sV1(2)+sV2(2), k+sV1(3)+sV2(3));
            n_in_corners(2) = sub2ind(s, i-sV1(1)+sV2(1), j-sV1(2)+sV2(2), k-sV1(3)+sV2(3));
            n_in_corners(3) = sub2ind(s, i+sV1(1)-sV2(1), j+sV1(2)-sV2(2), k+sV1(3)-sV2(3));
            n_in_corners(4) = sub2ind(s, i-sV1(1)-sV2(1), j-sV1(2)-sV2(2), k-sV1(3)-sV2(3));
            
            n_out_axial(1) = sub2ind(s, iii+sV1(1), jjj+sV1(2), kkk+sV1(3));
            n_out_axial(2) = sub2ind(s, iii-sV1(1), jjj-sV1(2), kkk-sV1(3));
            n_out_axial(3) = sub2ind(s, iii+sV2(1), jjj+sV2(2), kkk+sV2(3));
            n_out_axial(4) = sub2ind(s, iii-sV2(1), jjj-sV2(2), kkk-sV2(3));
            
            n_out_corners(1) = sub2ind(s, iii+sV1(1)+sV2(1), jjj+sV1(2)+sV2(2), kkk+sV1(3)+sV2(3));
            n_out_corners(2) = sub2ind(s, iii-sV1(1)+sV2(1), jjj-sV1(2)+sV2(2), kkk-sV1(3)+sV2(3));
            n_out_corners(3) = sub2ind(s, iii+sV1(1)-sV2(1), jjj+sV1(2)-sV2(2), kkk+sV1(3)-sV2(3));
            n_out_corners(4) = sub2ind(s, iii-sV1(1)-sV2(1), jjj-sV1(2)-sV2(2), kkk-sV1(3)-sV2(3));
            
            
            A1(n, n) = d4 / (1 + sig / beta);
            A2(n, n) = (sig/beta - 1) / (1 + sig / beta);
            
            A1(n, n_out) = 2*d1 / (1 + sig/beta);
            
            for k = 1:4
                A1(n, n_in_axial(k)) = d1 / (1 + sig/beta);
                A1(n, n_out_axial(k)) = 2*d2 / (1 + sig/beta);
                A1(n, n_in_corners(k)) = d2 / (1 + sig/beta);
                A1(n, n_out_corners(k)) = 2*d3 / (1 + sig/beta);
            end
            
        else
            
            clear n_in n_out
            n_out = sub2ind(s, i+nodeOrientation(1), j+nodeOrientation(2), k+nodeOrientation(3));
            
            ind = 1;
            for ii = find( nodeOrientation == 0 )
                temp = [0 0 0];
                temp(ii) = 1;
                n_in(ind) = sub2ind(s, i+temp(1), j+temp(2), k+temp(3));
                n_in(ind+1) = sub2ind(s, i-temp(1), j-temp(2), k-temp(3));
                
                ind = ind + 2;
            end
            
            if boundaries == 2 % SLF boundaries
                
                A1(n, n) = (2-5*bs2) / (1 + bsig/beta/2);
                A2(n, n) = (bsig/beta/2 - 1) / (1 + bsig/beta/2);
                A1(n, n_out) = bs2 / (1 + bsig/beta/2);
                A1(n, n_in(1)) = bs2 / (1 + bsig/beta/2);
                A1(n, n_in(2)) = bs2 / (1 + bsig/beta/2);
                A1(n, n_in(3)) = bs2 / (1 + bsig/beta/2);
                A1(n, n_in(4)) = bs2 / (1 + bsig/beta/2);
                
            end
        end
        
        case 'edge'
            if beta == 0
                A1(n, :) = 0;
                A2(n, :) = 0;
                
            elseif boundaries == 1
                
                clear n_axial n_edge n_wall_corners n_interior
                
                ii = find( nodeOrientation ~= 0 );
                eV = ~nodeOrientation; % unit vector along edge
                uV1 = [0 0 0]; uV1(ii(1)) = nodeOrientation(ii(1));
                uV2 = [0 0 0]; uV2(ii(2)) = nodeOrientation(ii(2));
                
                n_axial(1) = sub2ind(s, i+uV1(1), j+uV1(2), k+uV1(3));
                n_axial(2) = sub2ind(s, i+uV2(1), j+uV2(2), k+uV2(3));
                n_edge(1) = sub2ind(s, i+eV(1), j+eV(2), k+eV(3));
                n_edge(2) = sub2ind(s, i-eV(1), j-eV(2), k-eV(3));
                
                n_wall_corners(1) = sub2ind(s, i+uV1(1)+eV(1), j+uV1(2)+eV(2), k+uV1(3)+eV(3));
                n_wall_corners(2) = sub2ind(s, i+uV1(1)-eV(1), j+uV1(2)-eV(2), k+uV1(3)-eV(3));
                
                n_wall_corners(3) = sub2ind(s, i+uV2(1)+eV(1), j+uV2(2)+eV(2), k+uV2(3)+eV(3));
                n_wall_corners(4) = sub2ind(s, i+uV2(1)-eV(1), j+uV2(2)-eV(2), k+uV2(3)-eV(3));
                
                n_interior(1) = sub2ind(s, i+uV1(1)+uV2(1), j+uV1(2)+uV2(2), k+uV1(3)+uV2(3));
                n_interior(2) = sub2ind(s, i+uV1(1)+uV2(1)+eV(1), j+uV1(2)+uV2(2)+eV(2), k+uV1(3)+uV2(3)+eV(3));
                n_interior(3) = sub2ind(s, i+uV1(1)+uV2(1)-eV(1), j+uV1(2)+uV2(2)-eV(2), k+uV1(3)+uV2(3)-eV(3));
                
                %%%%% update %%%%%
                A1(n, n) = d4 / (1 + 2*sig / beta);
                A2(n, n) = (2*sig/beta - 1) / (1 + 2*sig / beta);
                
                A1(n, n_axial(1)) = 2*d1 / (1 + 2*sig/beta);
                A1(n, n_axial(2)) = 2*d1 / (1 + 2*sig/beta);
                A1(n, n_edge(1)) = d1 / (1 + 2*sig/beta);
                A1(n, n_edge(2)) = d1 / (1 + 2*sig/beta);
                
                A1(n, n_interior(1)) = 4*d2 / (1 + 2*sig/beta);
                A1(n, n_interior(2)) = 4*d3 / (1 + 2*sig/beta);
                A1(n, n_interior(3)) = 4*d3 / (1 + 2*sig/beta);
                
                for k = 1:4
                    A1(n, n_wall_corners(k)) = 2*d2 / (1 + 2*sig/beta);
                end
                
            else
                
                clear n_in n_out
                ind = 1;
                for ii = find( nodeOrientation ~= 0 )
                    temp = [0 0 0];
                    temp(ii) = nodeOrientation(ii);
                    
                    n_out(ind) = sub2ind(s, i+temp(1), j+temp(2), k+temp(3));
                    ind = ind + 1;
                end
                
                
                ind = 1;
                for ii = find( nodeOrientation == 0 )
                    temp = [0 0 0];
                    temp(ii) = 1;
                    n_in(ind) = sub2ind(s, i+temp(1), j+temp(2), k+temp(3));
                    n_in(ind+1) = sub2ind(s, i-temp(1), j-temp(2), k-temp(3));
                    
                    ind = ind + 2;
                end
                
                if boundaries == 2
                    
                    A1(n, n) = (2-4*bs2) / (1 + 2*bsig/beta/2);
                    A2(n, n) = (2*bsig/beta/2 - 1) / (1 + 2*bsig/beta/2);
                    A1(n, n_out(1)) = bs2 / (1 + 2*bsig/beta/2);
                    A1(n, n_out(2)) = bs2 / (1 + 2*bsig/beta/2);
                    A1(n, n_in(1)) = bs2 / (1 + 2*bsig/beta/2);
                    A1(n, n_in(2)) = bs2 / (1 + 2*bsig/beta/2);
                    
                end
            end
            
            
            case 'corner'
                if beta == 0
                    A1(n, :) = 0;
                    A2(n, :) = 0;
                    
                elseif boundaries == 1
                    
                    uV1 = [nodeOrientation(1) 0 0];
                    uV2 = [0 nodeOrientation(2) 0];
                    uV3 = [0 0 nodeOrientation(3)];
                    
                    n_axial(1) = sub2ind(s, i+nodeOrientation(1), j, k);
                    n_axial(2) = sub2ind(s, i, j+nodeOrientation(2), k);
                    n_axial(3) = sub2ind(s, i, j, k+nodeOrientation(3));
                    
                    n_corners(1) = sub2ind(s, i+uV1(1)+uV2(1), j+uV1(2)+uV2(2), k+uV1(3)+uV2(3));
                    n_corners(2) = sub2ind(s, i+uV1(1)+uV3(1), j+uV1(2)+uV3(2), k+uV1(3)+uV3(3));
                    n_corners(3) = sub2ind(s, i+uV3(1)+uV2(1), j+uV3(2)+uV2(2), k+uV3(3)+uV2(3));
                    
                    n_opposite = sub2ind(s, i+uV1(1)+uV3(1)+uV2(1), j+uV1(2)+uV3(2)+uV2(2), k+uV1(3)+uV3(3)+uV2(3));
                    
                    A1(n, n) = d4 / (1 + 3*sig / beta);
                    A2(n, n) = (3*sig/beta - 1) / (1 + 3*sig / beta);
                    
                    for k = 1:3
                        A1(n, n_axial(k)) = 2*d1 / (1 + 3*sig/beta);
                        A1(n, n_corners(k)) = 4*d2 / (1 + 3*sig/beta);
                    end
                    
                    A1(n, n_opposite) = 8*d3 / (1 + 3*sig/beta);
                    
                elseif boundaries == 2
                    
                    n_ip1 = sub2ind(s, i+nodeOrientation(1), j, k);
                    n_jp1 = sub2ind(s, i, j+nodeOrientation(2), k);
                    n_kp1 = sub2ind(s, i, j, k+nodeOrientation(3));
                    
                    
                    
                        
                        
                        A1(n, n) = (2 - 3*bs2) / (1 + 3*bsig/beta/2);
                        A2(n, n) = (3*bsig/beta/2 - 1) / (1 + 3*bsig/beta/2);
                        A1(n, n_ip1) = bs2 / (1 + 3*bsig/beta/2);
                        A1(n, n_jp1) = bs2 / (1 + 3*bsig/beta/2);
                        A1(n, n_kp1) = bs2 / (1 + 3*bsig/beta/2);
                   
                end
    end
end



% now fix hard source node
if ~any(isnan(handles.hs_pos))
    hsRow = sub2ind(s, handles.hs_pos(1), handles.hs_pos(2), handles.hs_pos(3));
    A1(hsRow, :) = 0;
    A2(hsRow, :) = 0;
end



% t1 = toc;
% disp([num2str(t1) ' sec to set up matrix'])

%         if use_sparse
%             A = [[A1 A2]; [speye(size(A1)) sparse(size(A1))]];
%         else
A = [[A1 A2]; [eye(size(A1)) zeros(size(A1))]];
%         end

clear A1 A2

if single_precision
    A = single(A);
end

function mouseclick_callback(gcbo,eventdata, handles)
% the arguments are not important here, they are simply required for
% a callback function. we don't even use them in the function,
% but Matlab will provide them to our function, we we have to
% include them.
%
% first we get the point that was clicked on

eVals = handles.eigVals;
eVects = handles.eigVects;
dims = handles.meshSize;

% define colormap
r = [50 50 50]/255;       %# start
w = [1 1 1];    %# middle
b = [0 206 209]/255;       %# end

% colormap of size 64-by-3, ranging from red -> white -> blue
c1 = zeros(32,3); c2 = zeros(32,3);
for i=1:3
    c1(:,i) = linspace(r(i), w(i), 32);
    c2(:,i) = linspace(w(i), b(i), 32);
end
colors = [c1(1:end-1,:);c2];

cP = get(gca,'Currentpoint');
x = cP(1,1);
y = cP(1,2);

if size(eVals, 2) > 1
    eVals = diag(eVals);
end


[~, index] = min(abs(eVals - x+1i*y));
if ~isempty(index)
    V = reshape(eVects(1:end/2, index), dims);
    V = real(0.999 * V / max(max(max(abs(V)))));
    
    figure;
    s1 = slice(V, dims(2), dims(1) , 1);
    hold all;
    s2 = slice(V, round(dims(2)/2), dims(1) , round(dims(3)/2) );
    set(s2, 'facealpha', 0.4, 'EdgeColor', 'none');
    colormap(colors);
    caxis([-1 1])
    set(gca, 'fontname', 'Times New Roman')
    
    l(1) = line([dims(2) dims(2)], [1 dims(1)], [1 1]);
    l(2) = line([1 dims(2)], [dims(1) dims(1)], [1 1]);
    l(3) = line([dims(2) dims(2)], [dims(1) dims(1)], [1 dims(3)]);
    
    set(l, 'color', 'k', 'linewidth', 1);
    axis equal;
    
end

function alt_mouseclick_callback(gcbo,eventdata, handles)
% the arguments are not important here, they are simply required for
% a callback function. we don't even use them in the function,
% but Matlab will provide them to our function, we we have to
% include them.
%
% first we get the point that was clicked on

eVals = handles.eigVals;
eVects = handles.eigVects;
dims = handles.meshSize;

% define colormap
r = [50 50 50]/255;       %# start
w = [1 1 1];    %# middle
b = [0 206 209]/255;       %# end

% colormap of size 64-by-3, ranging from red -> white -> blue
c1 = zeros(32,3); c2 = zeros(32,3);
for i=1:3
    c1(:,i) = linspace(r(i), w(i), 32);
    c2(:,i) = linspace(w(i), b(i), 32);
end
colors = [c1(1:end-1,:);c2];

cP = get(gca,'Currentpoint');
x = cP(1,1);
y = cP(1,2);

if size(eVals, 2) > 1
    eVals = diag(eVals);
end

[~, index] = min(abs(eVals(:) - x+1i*y));
if ~isempty(index)
    V = reshape(eVects(1:end/2, index), dims);
    VR = real(0.999 * V / max(abs(real(V(:)))));
    
    
    figure;
    sh = slice(1:dims(2), 1:dims(1), handles.zscale*(1:dims(3)), VR, [], [], (2:dims(3)-1)*handles.zscale );
    set(sh, 'linewidth', 0.5);
    colormap(colors);
    caxis([-1 1]);
    zlim([1 dims(3)]);
    set(gca, 'fontname', 'Times New Roman');
    axis off;
    axis equal;
    view([-38 22])
    
end

function [node_type, node_orientation] = GetNodeType(i, j, k, mesh)

% convention is True(1) is in the domain and False (0) is out

if mesh(i, j, k) % if interior
    
    nn = [mesh(i+1, j, k), -mesh(i-1, j, k),...
        mesh(i, j+1, k), -mesh(i, j-1, k), ...
        mesh(i, j, k+1), -mesh(i, j, k-1) ];
    
    
    if all(nn) % all neighbors are inside
        
        node_type = 'interior';
        node_orientation = [sum(nn(1:2)) sum(nn(3:4)) sum(nn(5:6))];
        
    elseif sum(abs(nn)) == 5 % one is out
        node_type = 'face';
        node_orientation = [sum(nn(1:2)) sum(nn(3:4)) sum(nn(5:6))];
        
    elseif sum(abs(nn)) == 4 % 2 are out
        node_type = 'edge';
        node_orientation = [sum(nn(1:2)) sum(nn(3:4)) sum(nn(5:6))];
        
    elseif sum(abs(nn)) == 3 % 3 are out
        node_type = 'corner';
        node_orientation = [sum(nn(1:2)) sum(nn(3:4)) sum(nn(5:6))];
    else
        disp('Node on interior but not a type we expect...');
        disp([i, j, k])
        node_type = NaN;
        node_orientation = NaN;
    end
    
    
else
    node_type = NaN;
    node_orientation = NaN;
    
end

function [node_type, node_orientation] = GetNodeTypeBox(i, j, k, s)

N = s(1);
M = s(2);
L = s(3);

if i > 1 && i < N && j > 1 && j < M && k > 1 && k < L
    node_type = 'interior';
    node_orientation = [];
    
elseif i == 1
    if j > 1 && j < M
        if k > 1 && k < L
            node_type = 'face';
            node_orientation = [1 0 0];
        elseif k == 1
            node_type = 'edge';
            node_orientation = [1 0 1];
        elseif k == L
            node_type = 'edge';
            node_orientation = [1 0 -1];
        end
    elseif j == 1
        if k > 1 && k < L
            node_type = 'edge';
            node_orientation = [1 1 0];
        elseif k == 1
            node_type = 'corner';
            node_orientation = [1 1 1];
        elseif k == L
            node_type = 'corner';
            node_orientation = [1 1 -1];
        end
        
    elseif j == M
        if k > 1 && k < L
            node_type = 'edge';
            node_orientation = [1 -1 0];
        elseif k == 1
            node_type = 'corner';
            node_orientation = [1 -1 1];
        elseif k == L
            node_type = 'corner';
            node_orientation = [1 -1 -1];
        end
    end
    
elseif i == N
    if j > 1 && j < M
        if k > 1 && k < L
            node_type = 'face';
            node_orientation = [-1 0 0];
        elseif k == 1
            node_type = 'edge';
            node_orientation = [-1 0 1];
        elseif k == L
            node_type = 'edge';
            node_orientation = [-1 0 -1];
        end
    elseif j == 1
        if k > 1 && k < L
            node_type = 'edge';
            node_orientation = [-1 1 0];
        elseif k == 1
            node_type = 'corner';
            node_orientation = [-1 1 1];
        elseif k == L
            node_type = 'corner';
            node_orientation = [-1 1 -1];
        end
        
    elseif j == M
        if k > 1 && k < L
            node_type = 'edge';
            node_orientation = [-1 -1 0];
        elseif k == 1
            node_type = 'corner';
            node_orientation = [-1 -1 1];
        elseif k == L
            node_type = 'corner';
            node_orientation = [-1 -1 -1];
        end
    end
    
elseif j == 1
    if k == 1
        node_type = 'edge';
        node_orientation = [0 1 1];
    elseif k == L
        node_type = 'edge';
        node_orientation = [0 1 -1];
    else
        node_type = 'face';
        node_orientation = [0 1 0];
    end
    
elseif j == M
    if k == 1
        node_type = 'edge';
        node_orientation = [0 -1 1];
    elseif k == L
        node_type = 'edge';
        node_orientation = [0 -1 -1];
    else
        node_type = 'face';
        node_orientation = [0 -1 0];
    end
    
elseif k == 1
    node_type = 'face';
    node_orientation = [0 0 1];

elseif k == L
    node_type = 'face';
    node_orientation = [0 0 -1];
    
else
    fprintf('Something went wrong...\n (%d, %d, %d)', i, j, k);
    return;
    
end


function PlotFrequencies(fh, s, r1)

i = 0:s.N-1; j = 0:s.M-1; k = 0:s.L-1;
if s.R < 0
    i = i+1; j = j+1; k = k+1;
end
[ii, jj, kk] = meshgrid(i, j, k);

if s.BC == 2
    norm_modal_frequencies = pi * s.sig * sqrt( (ii/(s.N)).^2 + (jj/(s.M)).^2 + (kk/(s.L)).^2 );
elseif s.BC == 1
    norm_modal_frequencies = pi * s.sig *sqrt( (ii/(s.N-1)).^2 + (jj/(s.M-1)).^2 + (kk/(s.L-1)).^2 );
end

angles = sort(norm_modal_frequencies(:));
angles = angles(angles < pi );

clear ii jj kk norm_modal_frequencies

angles = [angles; -flipud(angles)];

figure(fh);
hold all;
for k = 1:length(angles)
    line([r1*cos(angles(k)) cos(angles(k))],[r1*sin(angles(k)) sin(angles(k))], 'color', [0.8 0.8 0.8]);
end

function PlotFrequenciesWhite(fh, s, r1)

i = 0:s.N-1; j = 0:s.M-1; k = 0:s.L-1;
if s.R < 0
    i = i+1; j = j + 1; k = k+1;
end
[ii, jj, kk] = meshgrid(i, j, k);

if s.BC == 2
    norm_modal_frequencies = pi * sig *sqrt( (ii/(N-3 + 2*sig/bsig)).^2 + (jj/(M-3 + 2*sig/bsig)).^2 + (kk/(L-3 + 2*sig/bsig)).^2 );
elseif s.BC == 1
    norm_modal_frequencies = pi * s.sig *sqrt( (ii/(s.N-1)).^2 + (jj/(s.M-1)).^2 + (kk/(s.L-1)).^2 );
end

angles = sort(norm_modal_frequencies(:));
angles = angles(angles < pi );

clear ii jj kk norm_modal_frequencies

angles = [angles; -flipud(angles)];

figure(fh);
hold all;
for k = 1:length(angles)
    line([r1*cos(angles(k)) cos(angles(k))],[r1*sin(angles(k)) sin(angles(k))], 'color', [1 1 1]);
end


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if get(handles.popupmenu3, 'value') == 2
    handles.single_precision = true;
else
    handles.single_precision = false;
end


set(handles.pushbutton1, 'ForegroundColor', [0 0.5 0]); 
set(handles.status, 'string', 'Recompute eigenvalues...');

guidata(hObject, handles)



% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles.hs_pos = [str2double(get(handles.edit8, 'string'))...
    str2double(get(handles.edit9, 'string')) str2double(get(handles.edit10, 'string'))];

set(handles.pushbutton1, 'ForegroundColor', [0 0.5 0]); 
set(handles.status, 'string', 'Recompute eigenvalues...');

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.hs_pos = [str2double(get(handles.edit8, 'string'))...
    str2double(get(handles.edit9, 'string')) str2double(get(handles.edit10, 'string'))];

set(handles.pushbutton1, 'ForegroundColor', [0 0.5 0]); 
set(handles.status, 'string', 'Recompute eigenvalues...');

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles.hs_pos = [str2double(get(handles.edit8, 'string'))...
    str2double(get(handles.edit9, 'string')) str2double(get(handles.edit10, 'string'))];

set(handles.pushbutton1, 'ForegroundColor', [0 0.5 0]); 
set(handles.status, 'string', 'Recompute eigenvalues...');

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.status, 'string', 'Making figure...')
drawnow;

if size(handles.eigVals, 2) > 1
    handles.eigVals = diag(handles.eigVals);
end

figure;

circle_angles = linspace(-0.05, 0.05, 10000);
plot(cos(circle_angles), sin(circle_angles), 'color', [0.8 0.8 0.8], 'linewidth', 1);


hold all;

% zero axes
line([-1 1],[0 0], 'color', [0.8 0.8 0.8], 'linewidth', 1);

tempEV = handles.eigVals(abs(imag(handles.eigVals)) < 0.01);
us = tempEV(abs(tempEV) > 1);
st = tempEV(abs(tempEV) <= 1);

plot(real(st), imag(st), 'ko', 'markersize', 6, 'markerfacecolor', 0.75*[1 1 1], 'markeredgecolor', 'none');
plot(real(us), imag(us), 'ko', 'markersize', 6, 'markerfacecolor', 'k', 'markeredgecolor', 'none');

if handles.single_precision
    extra = 7.5e-4;
else
    extra = 7.5e-8;
end
xlim([1-extra 1+extra]);
ylim([-extra extra]);

title(handles.scheme, 'fontweight', 'bold')

axis square;

set(handles.status, 'string', 'Ready for plotting...')
drawnow;


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.status, 'string', 'Making figure...')
drawnow;

if size(handles.eigVals, 2) > 1
    handles.eigVals = diag(handles.eigVals);
end

figure;

circle_angles = linspace(pi-0.05, pi+0.05, 10000);
plot(cos(circle_angles), sin(circle_angles), 'color', [0.8 0.8 0.8], 'linewidth', 1);

hold all;

% zero axes
line([-1 1],[0 0], 'color', [0.8 0.8 0.8], 'linewidth', 1);

tempEV = handles.eigVals( abs(imag(handles.eigVals)) < 0.01 );
us = tempEV(abs(tempEV) > 1);
st = tempEV(abs(tempEV) <= 1);

plot(real(st), imag(st), 'ko', 'markersize', 6, 'markerfacecolor', 0.75*[1 1 1], 'markeredgecolor', 'none');
plot(real(us), imag(us), 'ko', 'markersize', 6, 'markerfacecolor', 'k', 'markeredgecolor', 'none');

if handles.single_precision
    extra = 1.6e-3;
else
    extra = 7.5e-8;
end
xlim([-1-extra -1+extra]);
ylim([-extra extra]);

title(handles.scheme, 'fontweight', 'bold')

axis square;

set(handles.status, 'string', 'Ready for plotting...')
drawnow;

% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.status, 'string', 'Making figure...')
drawnow;

if size(handles.eigVals, 2) > 1
    handles.eigVals = diag(handles.eigVals);
end

ang = angle(handles.eigVals);

figure; plot(ang(abs(handles.eigVals) > 1e-5), abs(handles.eigVals(abs(handles.eigVals) > 1e-5)),...
    'ok', 'markerfacecolor', 'k', 'markeredgecolor', 'none');
xlim([-0.1 pi+0.1]);
set(gca, 'xtick', [0 pi/4 pi/2 3*pi/4 pi], 'xticklabel', {'0' , '\pi/4', '\pi/2', '3\pi/4', '\pi'});
pos = get(gcf, 'position');
set(gcf, 'position', [pos(1:3) pos(4)/2]);

SR = max(abs(handles.eigVals));
[~, mreInd] = max(real(handles.eigVals));
[~, minReInd] = min(real(handles.eigVals));


fprintf('\nSpectral Radius is 1 + %e \n\n', SR-1)

fprintf('\nDC Radius is 1 + %e \n\n', abs(handles.eigVals(mreInd))-1);

fprintf('\nNyquist Radius is 1 + %e \n\n', abs(handles.eigVals(minReInd))-1);

set(handles.status, 'string', 'Ready for plotting...')
drawnow;

% --------------------------------------------------------------------
function options_Callback(hObject, eventdata, handles)
% hObject    handle to options (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function EVPlotType_Callback(hObject, eventdata, handles)
% hObject    handle to EVPlotType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --------------------------------------------------------------------
function slices_Callback(hObject, eventdata, handles)
% hObject    handle to slices (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.slices,'Checked','on');
set(handles.alpha_slices,'Checked','off');

guidata(hObject, handles);


% --------------------------------------------------------------------
function alpha_slices_Callback(hObject, eventdata, handles)
% hObject    handle to alpha_slices (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


set(handles.slices,'Checked','off');
set(handles.alpha_slices,'Checked','on');

guidata(hObject, handles);


% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

val = get(handles.popupmenu4, 'value');
str = get(handles.popupmenu4, 'string');

if strcmp(str(val), 'Box')
    set(handles.text8, 'string', num2str(2*handles.N*handles.M*handles.L));
    
    set(handles.edit3, 'enable', 'on');
    set(handles.edit6, 'enable', 'on');
    set(handles.edit7, 'enable', 'on');
    
    set(handles.popupmenu1, 'enable', 'on');
    
    handles.meshSize = [handles.N handles.M handles.L];
    
else
    
    handles.scheme = 'SLF';
    
    if strcmp(str(val),'L1 (unstable)')
        handles.meshSize = [13 13 5];
    elseif strcmp(str(val), 'L2 (stable)')
        handles.meshSize = [13 13 5];
    elseif strcmp(str(val), 'Corner')
        handles.meshSize = [9 9 9];
    elseif strcmp(str(val), 'Cylinder')
        handles.meshSize = [11 11 6];
    elseif strcmp(str(val), 'Sphere')
        handles.meshSize = [9 9 9];
    end
    
    
    set(handles.text8, 'string', num2str( 2*prod(handles.meshSize)));
    
    set(handles.edit3, 'enable', 'off');
    set(handles.edit6, 'enable', 'off');
    set(handles.edit7, 'enable', 'off');
    
    set(handles.popupmenu1, 'enable', 'off', 'value', 1);
    set(handles.popupmenu2, 'enable', 'on');
    
end

set(handles.pushbutton1, 'ForegroundColor', [0 0.5 0]); 
set(handles.status, 'string', 'Recompute eigenvalues...');

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function EVPlotSubMenu_Callback(hObject, eventdata, handles)
% hObject    handle to EVPlotSubMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function normal_Callback(hObject, eventdata, handles)
% hObject    handle to normal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.normal,'Checked','on');
set(handles.inverted,'Checked','off');

guidata(hObject, handles);


% --------------------------------------------------------------------
function inverted_Callback(hObject, eventdata, handles)
% hObject    handle to inverted (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.normal,'Checked','off');
set(handles.inverted,'Checked','on');

guidata(hObject, handles);


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


m = GetBoundaryMask(handles);

dims = size(m);

xi = 2:dims(1)-1;
yi = 2:dims(2)-1;
zi = 2:dims(3)-1;

mask = zeros(size(m));

mask(xi, yi, zi) = mask(xi, yi, zi) + m(xi+1, yi, zi) + m(xi-1, yi, zi)...
    + m(xi, yi+1, zi) + m(xi, yi-1, zi)...
    + m(xi, yi, zi+1) + m(xi, yi, zi-1);

mask(m == 0) = 0;


figure;
slice(1:dims(2), 1:dims(1), handles.zscale*(1:dims(3)), mask, [], [], (2:dims(3)-1)*handles.zscale );
zlim([1 dims(3)]);
axis equal;
hold off;
view([-38 22])
axis off;
title('Number of interior neighbors', 'fontweight', 'bold' )
% colormap(gray)
caxis([1 6]);
colorbar;

guidata(hObject, handles)


% --------------------------------------------------------------------
function zscale_Callback(hObject, eventdata, handles)
% hObject    handle to zscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.zscale = str2double(inputdlg('Set the z-scale for spacing slices of eigenvectors:',...
    'Set Z-Scale' ,1, {num2str(handles.zscale)}));


guidata(hObject, handles);
