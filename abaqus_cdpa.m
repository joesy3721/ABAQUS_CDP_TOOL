function varargout = abaqus_cdpa(varargin)
% ABAQUS_CDPA MATLAB code for abaqus_cdpa.fig
%      ABAQUS_CDPA, by itself, creates a new ABAQUS_CDPA or raises the existing
%      singleton*.
%
%      H = ABAQUS_CDPA returns the handle to a new ABAQUS_CDPA or the handle to
%      the existing singleton*.
%
%      ABAQUS_CDPA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ABAQUS_CDPA.M with the given input arguments.
%
%      ABAQUS_CDPA('Property','Value',...) creates a new ABAQUS_CDPA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before abaqus_cdpa_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to abaqus_cdpa_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help abaqus_cdpa

% Last Modified by GUIDE v2.5 25-Jun-2021 07:41:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @abaqus_cdpa_OpeningFcn, ...
                   'gui_OutputFcn',  @abaqus_cdpa_OutputFcn, ...
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


% --- Executes just before abaqus_cdpa is made visible.
function abaqus_cdpa_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to abaqus_cdpa (see VARARGIN)

% Choose default command line output for abaqus_cdpa
handles.output = hObject;

%program name
handles.title = 'ABAQUS CDP Calculator';

% Update handles structure
guidata(hObject, handles);

execute(handles)





% UIWAIT makes abaqus_cdpa wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = abaqus_cdpa_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_fcm_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fcm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_fcm as text
%        str2double(get(hObject,'String')) returns contents of edit_fcm as a double
execute(handles)

% --- Executes during object creation, after setting all properties.
function edit_fcm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fcm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ecm_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ecm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ecm as text
%        str2double(get(hObject,'String')) returns contents of edit_ecm as a double
execute(handles)

% --- Executes during object creation, after setting all properties.
function edit_ecm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ecm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% compressive behavior
function execute(handles)

% read input
    % Compressive strength
    fcm = str2num(handles.edit_fcm.String);
    
    % Young's Modulus
    Ecm = str2num(handles.edit_ecm.String);
    
    % maximal compression damage
    o_c_max = str2num(handles.edit_d_c.String);
    
    % compression damage at ec
    o_c_ec = str2num(handles.edit_decu.String);
    
    % maximal tensile damage
    o_t_max = str2num(handles.edit_d_t.String);
    
    % data points compression
    poi_c_left = str2num(handles.edit_points_c.String);
    
    % data points tension
    poi_t = str2num(handles.edit_points_t.String);
    
    %softening
    soft_check = handles.checkbox_soft.Value;
    
    % gfi ruler
    ratio_gfi = str2num(handles.edit_gfi_ratio.String);
    
    %characterisitc length
    lchar = str2num(handles.edit_lchar.String);
    
    % gfi tension
    gfi_t = str2num(handles.edit_gfi.String);
    
    % tensile strength
    fctm = str2num(handles.edit_fctm.String);
    
    % beta
    beta = str2num(handles.edit_beta.String);
    % beta Warning
    if beta<0.9
        set(handles.text_beta_warning,'string','Warning: rapidly increasing damage values impair convergence. Increase beta.')
        set(handles.text_beta_warning,'ForegroundColor',[1 1 0.07])
        handles.export.Enable = 'off';
    else
        set(handles.text_beta_warning,'string','Suitable beta value.')
        set(handles.text_beta_warning,'ForegroundColor',[0.47 0.67 0.19])
        handles.export.Enable = 'on';
    end
 
%% tension

%data with with initial gfi_t

% inelastic displacement at 1% of tensile strength
iu001 = -log(0.01)*gfi_t/fctm;

% inelastic displacement-values
iu = linspace(0,iu001,100);  

%normalized stress curve
st_n = exp(-iu*fctm/gfi_t);

% output data indices, logarithmic spacing for better capture of damage at
% the beginning
idx_t = [round(logspace(0,2,poi_t))];

% stress, inelastic displacement, damage output data
st_out = [st_n(idx_t)]*fctm; %in MPa
iu_out = [iu(idx_t)]; %in mm


%% gfi_t iteration
    
    %initial gfi_t due to multilinear approximation
    gfi_t_adjust=trapz(iu_out,st_out); %in Nmm/mm^2 ; always bigger than desired one
    
    %initial gfi_t factor
    gfi_factor=0;
    
    %formula for gfi adjustment to reach desired GFIt with multilinear approximation
    gfi_t_adjust = gfi_t_adjust - gfi_t_adjust.*gfi_factor;
        
    %calculate ratio prior to gfi adjustment
    ratio=gfi_t_adjust/gfi_t; % >1
    n=0;
    %initial factor split if ratio got exceeded
    m=2;
    %iteration: adjust gfi_factor until desired ratio is achieved (1/100
    %tolerance, max 500 iterations)
    while ratio > 1.01 && n<200 % 1% tolerance
        %increase of factor in this step
        plus=ratio/20;
        %save old parameters, if increase was too much
        factor_old = gfi_factor;
        gfi_t_adjust_old = gfi_t_adjust;
        
        %new weight factor for this step
        gfi_factor=gfi_factor+plus;
           
        %gfi adjustment to reach desired GFIt with multilinear approximation
        gfi_t_adjust = gfi_t_adjust - gfi_t_adjust.*gfi_factor;
        

        % new stress-displacement data

            % inelastic displacement at 1% of tensile strength
            iu001 = -log(0.01)*gfi_t_adjust/fctm;

            % inelastic displacement-values
            iu = linspace(0,iu001,100);  

            %normalized stress curve
            st_n = exp(-iu*fctm/gfi_t_adjust);

            % stress, inelastic displacement, damage output data
            st_out = [st_n(idx_t)]*fctm; %in MPa
            iu_out = [iu(idx_t)]; %in mmt      
         
        %new gfi_t
        gfi_t_adjust=trapz(iu_out,st_out);
        
        %calculate new ratio
        ratio=gfi_t_adjust/gfi_t;
        %counter
        n=n+1;
        
        while ratio < .99 && n<200
            %last value of upper loop that produces a ratio>1
            gfi_t_adjust = gfi_t_adjust_old;
            %reduce gfi_factor
            gfi_factor=factor_old+plus/m;
            
            %gfi adjustment to reach desired GFIt with multilinear approximation
            gfi_t_adjust = gfi_t_adjust - gfi_t_adjust.*gfi_factor;

            % new stress-displacement data

                % inelastic displacement at 1% of tensile strength
                iu001 = -log(0.01)*gfi_t_adjust/fctm;

                % inelastic displacement-values
                iu = linspace(0,iu001,100);  

                %normalized stress curve
                st_n = exp(-iu*fctm/gfi_t_adjust);

                % stress, inelastic displacement, damage output data
                st_out = [st_n(idx_t)]*fctm; %in MPa
                iu_out = [iu(idx_t)]; %in mmt 
                
            %new gfi_t
            gfi_t_adjust=trapz(iu_out,st_out);

            %calculate new ratio
            ratio=gfi_t_adjust/gfi_t;
            %counter
            n=n+1;

            %increase of factor split if ratio is still exceeded
            m=m+2;
        end
        %reset factor split to initial value
        m=2;
    end

%%

%damage
o_t=(1-st_n./(st_n+Ecm*iu*(1-beta)))*o_t_max;

%plastic displacement

pe_t = iu - o_t./(1-o_t).*st_n/Ecm*fctm;

%check if pe data is correct
for i=1:length(pe_t)-1 
    if pe_t(i)>pe_t(i+1) | pe_t(i+1)<0
        flag=1;
        set(handles.text_warning_plastic_t,'string','Error: Invalid TENSION data. Negative and/or decreasing plastic displacement. Increase beta or decrease damage.')
        set(handles.text_warning_plastic_t,'ForegroundColor',[1 0 0])
        handles.export.Enable = 'off';
        break        
    else
        set(handles.text_warning_plastic_t,'string','Valid plastic strain data.')
        set(handles.text_warning_plastic_t,'ForegroundColor',[0.47 0.67 0.19])
        handles.export.Enable = 'on';
        flag=0;
    end    
    if(flag==1)
        break
    end
end

% damage output data
o_t_out = round([o_t(idx_t)],2);
    
%make output data available for other functions   
% stress
    global st_iu_out
    st_iu_out =[st_out;iu_out];
% damage
    global omt_iu_out
    omt_iu_out =[o_t_out;iu_out];
        
%% compression
    
%total strain at compressive strength gem tab.3.1, EN1992-1-1
ec1 = min([0.7*fcm^0.31,2.8]);

%total strain at failure
if fcm-8 < 50
    ecu1 = 3.5;
else
    ecu1 = 2.8 + 27*((98-fcm)/100)^4
end
    
%plasticity number
k = 1.05*Ecm*(ec1/1000/fcm);

%strain at 0.4 * compressive strengt, end of elastic domain
e04f = 1/10 * (-sqrt(3)*sqrt(3*k^2+8*k-8)+3*k+4)*ec1;

%strain at 0 bzw. 0.01* compressive strength
e001 = 1/200 * (3*sqrt(11)*sqrt(99*k^2+4*k-4)+99*k+2)*ec1;

%conversion of ec1 to inelastic strain
inec1 = ec1-e04f-(1-0.4)*fcm/Ecm*1000; %bei ec1

%conversion of e001 to inelastic strain
ine001 = e001-e04f-(0.01-0.4)*fcm/Ecm*1000; %bei 0.01*fcm

%initializing number of data points, equal for softening and hardening
poi_c_right = poi_c_left;

%if compression softening is neglected
if soft_check ~= 1
    
    %disable lchar and gfi ratio field, color adjustment
    handles.edit_lchar.Enable = 'off';
    set(handles.text_lchar,'BackgroundColor',[0.64 0.83 0.51])
    handles.slider_gfi_ratio.Enable = 'off';
    handles.edit_gfi_ratio.Enable = 'off';
    set(handles.text_gfi_ratio,'BackgroundColor',[0.64 0.83 0.51])
    
    % strain-values
    ec = linspace(e04f,ecu1,100);  

    eta = ec/ec1;   

    %normalized stress curve
    sc_n = (k*eta-eta.^2)./(1+(k-2)*eta);

    %Plateau after fcm
    n=1;
    for i = sc_n
        n = min([n+1,100]);
        if i>sc_n(n)
            sc_n(n)=1;
        end
    end  
    %overwrite data point number for plateau
    poi_c_right = 2;
    %no regularization warning
    set(handles.text_softening,'string','')
    %no compression gfi
    set(handles.text_gfi_c,'string','-')
    
    %inelastic strain 
    ine = ec-e04f-(sc_n-sc_n(1))*fcm/Ecm*1000;
    
    %conversion of ecu1 to inelastic strain
    inecu1 = ecu1-e04f-(1-0.4)*fcm/Ecm*1000; %bei ecu1
    
    %Damage equation system system; curve goes through 0,0; inec1,o_c_ec;
    %inecu1,o_c_max and has a horizontal tangent at inecu1,o_c_max
    % o_c_max is the damage limit, o_c_ec is the damage when reaching fcm
    Mat = [(inecu1/1000)^3 (inecu1/1000)^2 (inecu1/1000); (inec1/1000)^3 (inec1/1000)^2 (inec1/1000); 3*(inecu1/1000)^2 2*(inecu1/1000) 1];
    Vec = [o_c_max; o_c_ec; 0];
    %solve equation
    Sol = Mat\Vec;
    
    %damage strain curve; polynomial function third order
    o_c = Sol(1)*(ine/1000).^3 + Sol(2)*(ine/1000).^2 + Sol(3)*(ine/1000);
    
    %check if damage data is correct
    for i=1:length(o_c)-1 
        if o_c(i)>o_c(i+1) | o_c(i+1)<0
            flag=1;
            set(handles.text_warning_dam,'string','Error: Invalid COMPRESSIVE DAMAGE data. Negative and/or decreasing damage values occur. Adjust damage.')
            set(handles.text_warning_dam,'ForegroundColor',[1 0 0])
            handles.export.Enable = 'off';
            break
        else
            set(handles.text_warning_dam,'string','Valid damage data.')
            set(handles.text_warning_dam,'ForegroundColor',[0.47 0.67 0.19])
            handles.export.Enable = 'on';
            flag=0;
        end    
        if(flag==1)
            break
        end
    end
        
    elim=max(ine);
    
    % outputdata indices, equal spacing
    %find fcm index I
    [M,I] = max(sc_n);   
    %from start to fcm, generate #poi_c_left indices and round them
    left=round(linspace(1,I,poi_c_left));
    %only start and end index for plateau
    right=round(linspace(I,length(sc_n),poi_c_right));
    %combine hardening and softening branch
    idx = [left,right(2:end)];
    
    % stress, inelastic displacement, damage output data
    sc_out = [sc_n(idx)]*fcm;
    ine_out = [ine(idx)];
    o_c_out = [o_c(idx)];
          
else % if compression softening is accounted for
    
    %enable lchar and gfi ratio field, color adjustment
    handles.edit_lchar.Enable = 'on';
    set(handles.text_lchar,'BackgroundColor',[0.39 0.83 0.07])
    handles.slider_gfi_ratio.Enable = 'on';
    handles.edit_gfi_ratio.Enable = 'on';
    set(handles.text_gfi_ratio,'BackgroundColor',[0.39 0.83 0.07])
    
    % strain-values
    ec = linspace(e04f,ecu1,100);  
   
    eta = ec/ec1;   

    %normalized stress curve
    sc_n = (k*eta-eta.^2)./(1+(k-2)*eta);
    
    %inelastic strain 
    ine = ec-e04f-(sc_n-sc_n(1))*fcm/Ecm*1000;
    
    %determination of inelastic ecu1
    inecu1 = ine(end);
    
    %Damage equation system system;
    Mat = [(inecu1/1000)^3 (inecu1/1000)^2 (inecu1/1000); (inec1/1000)^3 (inec1/1000)^2 (inec1/1000); 3*(inecu1/1000)^2 2*(inecu1/1000) 1];
    Vec = [o_c_max; o_c_ec; 0];
    Sol = Mat\Vec;
    
    %damage strain curve
    o_c = Sol(1)*(ine/1000).^3 + Sol(2)*(ine/1000).^2 + Sol(3)*(ine/1000);
    
    %check if damage data is correct
    for i=1:length(o_c)-1 
        if o_c(i)>o_c(i+1) | o_c(i+1)<0
            flag=1;
            set(handles.text_warning_dam,'string','Error: Invalid COMPRESSIVE DAMAGE data. Negative and/or decreasing damage values occur. Adjust damage.')
            set(handles.text_warning_dam,'ForegroundColor',[1 0 0])
            %export only possible when no warning occurs
            handles.export.Enable = 'off';
            break
        else
            set(handles.text_warning_dam,'string','Valid damage data.')
            set(handles.text_warning_dam,'ForegroundColor',[0.47 0.67 0.19])
            %export only possible when no warning occurs
            handles.export.Enable = 'on';
            flag=0;
        end    
        if(flag==1)
            break
        end
    end
    
    % outputdata indices, equal spacing
    %find fcm index I
    [M,I] = max(sc_n);   
    %from start to fcm, generate #poi_c_left indices and round them
    left=round(linspace(1,I,poi_c_left));
    %from fcm to end, generate #poi_c_right indices and round them
    right=round(linspace(I,length(sc_n),poi_c_right));
    %combine
    idx = [left,right(2:end),101];
    
    %weighting/shifting of softening branch for accurate gfi ratio
    %no shifting in hardening branch
    hardening=zeros(1,I-1);
    %linear shifting in softening branch (zero at fcm, 1 at last value)
    softening=linspace(0,1,100-length(hardening));
    %due to simply appending the 0.01 stress value to the stress array,
    %values are not evenly spaced, last value needs adaptions
    weight=[hardening,softening, (ine001-ine(end))/(ine(end)-ine(end-1))*(softening(2)-softening(1))+1];
    
    %curve with 101 data points
    %append 0.01 stress value
    sc_n = [sc_n,0.01]; %normalized
    ine = [ine,ine001];
    o_c = [o_c,o_c_max];
    
    %initialization of weight/shift factor
    gfi_factor=0;
    
    %formula for strain shift to reach certain GFIc
    ine0 = ine + ine.*weight*gfi_factor;
   
    %curve with customized data points for ABAQUS input
    sc_out = [sc_n(idx)]*fcm; %MPa
    ine_out = [ine(idx)];
    
    %gfi calculation with output data
    gfi_c=trapz(ine_out*lchar/1000,sc_out); 
    
    %gfic output
%     set(handles.text_gfi_c,'string',gfi_c)   
    
    %calculate ratio prior to strain shift
    ratio=gfi_c/gfi_t_adjust;
    %set minimum to gfi ratio before strain shift
    set(handles.slider_gfi_ratio, 'min', ratio);
    set(handles.slider_gfi_ratio, 'max', ratio+300);
    n=0;
    %initial factor split if ratio got exceeded
    m=2;
    %iteration: adjust gfi_factor until desired ratio is achieved (1/100
    %tolerance, max 500 iterations)
    while ratio < ratio_gfi && n<500 && (ratio-ratio_gfi)<ratio_gfi/100
        %increase of weight factor in this step
        plus=(abs(ratio - ratio_gfi))/10;
        %save old factor, if increase was too much
        factor_old = gfi_factor;
        %new weight factor for this step
        gfi_factor=gfi_factor+plus;
        
        %strain shift
        ine = ine0 + ine0.*weight*gfi_factor;

        %curve with customized data points for ABAQUS input
        sc_out = [sc_n(idx)]*fcm; %MPa
        ine_out = [ine(idx)];
        %calculate gfi
        gfi_c=trapz(ine_out*lchar/1000,sc_out); 
        
        %show result gfic
        set(handles.text_gfi_c,'string',gfi_c)   
        
        %calculate new ratio
        ratio=gfi_c/gfi_t_adjust;
        %counter
        n=n+1;
        
        while ratio > ratio_gfi
            %reduce step
            gfi_factor=factor_old+plus/m;
            
            %strain shift
            ine = ine0 + ine0.*weight*gfi_factor;

            %curve with customized data points for ABAQUS input
            sc_out = [sc_n(idx)]*fcm; %MPa
            ine_out = [ine(idx)];
            %calculate gfi
            gfi_c=trapz(ine_out*lchar/1000,sc_out); 

            %show result gfic
            set(handles.text_gfi_c,'string',gfi_c)   
            
            %new ratio
            ratio=gfi_c/gfi_t_adjust;
            
            %counter
            n=n+1;
            %increase of factor split if ratio is still exceeded
            m=m+1;
        end
        %reset factor split to initial value
        m=2;
    end
    
    %compression damage output
    o_c_out = [o_c(idx)];   
    
    %plot xlim
    elim=max(ine);   
    
    %regularization warning
    set(handles.text_softening,'string','Warning: compressive softening causes mesh dependent results')
    
    %reset iteration counter
    n=0;
end

%plastic strain
pe = ine/1000 - o_c./(1-o_c).*sc_n/Ecm*fcm;

%check plastic strain
for i=1:length(pe)-1 
    if pe(i)>pe(i+1) | pe(i+1)<0
        flag=1;
        set(handles.text_warning_plastic,'string','Error: Invalid COMPRESSION data. Negative and/or decreasing plastic displacement occur. Reduce damage.')
        set(handles.text_warning_plastic,'ForegroundColor',[1 0 0])
        handles.export.Enable = 'off';
        break
    else
        set(handles.text_warning_plastic,'string','Valid plastic strain data.')
        set(handles.text_warning_plastic,'ForegroundColor',[0.47 0.67 0.19])
        handles.export.Enable = 'on';
        flag=0;
    end    
    if(flag==1)
        break
    end
end
%global availability
global sc_iu_out
sc_iu_out =[sc_out;ine_out/1000];
%global availability
global omc_iu_out
omc_iu_out =[o_c_out;ine_out/1000];

%% plot compression

% plastic strain
    %plot(handles.axes2, ine, pe, [0,max(ine)],[0,max(ine)],'k:')
    plot(handles.axes2, pe)
    % xlim(handles.axes2,[0,elim]); 
    % ylim(handles.axes2,[0,elim]);
    %xlabel(handles.axes2,'inelastic strain (1/1000)');
    xlabel(handles.axes2,'index');
    ylabel(handles.axes2,'plastic strain (1/1000)');
    
%stress strain
        yyaxis(handles.axes1,'left')
        plot(handles.axes1, ine, sc_n*fcm)
        ylabel(handles.axes1,'compressive stress (MPa)');
        ylim(handles.axes1,[0,round(fcm+5,-1)]);
        hold(handles.axes1, 'on' )
        plot(handles.axes1,ine_out, sc_out, 'k--+')
        hold(handles.axes1, 'off' )
        
%damage strain
        yyaxis(handles.axes1,'right')
        plot(handles.axes1, ine, o_c)
        ylabel(handles.axes1,'compressive damage (-)');
        ylim(handles.axes1,[0,1]);
        hold(handles.axes1, 'on' )
        plot(handles.axes1,ine_out, o_c_out, 'k--+')
        hold(handles.axes1, 'off' )

        xlim(handles.axes1,[0,elim]); 
        xlabel(handles.axes1,'inelastic strain (1/1000)');

%show ec1 & ecu1 parameters
    set(handles.text_ec1,'string',ec1)
    set(handles.text_ecu1,'string',ecu1)

%% plot tension

% plastic displacement
    %plot(handles.axes3, iu, pe_t, [0,max(iu)],[0,max(iu)],'k:')
    plot(handles.axes3, pe_t)
    %xlim(handles.axes3,[0,iu001]); 
    %ylim(handles.axes3,[0,iu001]);
    %xlabel(handles.axes3,'inelastic displacement (mm)');
    xlabel(handles.axes3,'index');
    ylabel(handles.axes3,'plastic displacement (mm)');
    
%stress displacement
        yyaxis(handles.axes4,'left')
        handles.axes4.YColor = [0, 0.4470, 0.7410];
        plot(handles.axes4, iu, st_n*fctm,'Color',[0, 0.4470, 0.7410])
        ylabel(handles.axes4,'tensile stress (MPa)');
        ylim(handles.axes4,[0,round(fctm+.5,0)]);
        hold(handles.axes4, 'on' )
        plot(handles.axes4,iu_out, st_out,'k--+')
        hold(handles.axes4, 'off' )

%damage displacement
        yyaxis(handles.axes4,'right')
        handles.axes4.YColor = [0.8500, 0.3250, 0.0980];
        plot(handles.axes4, iu, o_t,'Color',[0.8500, 0.3250, 0.0980])
        ylabel(handles.axes4,'tensile damage (-)');
        ylim(handles.axes4,[0,1]);
        hold(handles.axes4, 'on' )
        plot(handles.axes4,iu_out, o_t_out,'k--+')
        hold(handles.axes4, 'off' )

        xlim(handles.axes4,[0,iu001]); 
        xlabel(handles.axes4,'inelastic displacement (mm)');



% --- Executes on button press in checkbox_soft.
function checkbox_soft_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_soft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_soft
execute(handles)



% --------------------------------------------------------------------
function file_Callback(hObject, eventdata, handles)
% hObject    handle to file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function new_Callback(hObject, eventdata, handles)
% hObject    handle to new (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close(gcbf)
abaqus_cdpa


% --------------------------------------------------------------------
function open_Callback(hObject, eventdata, handles)
% hObject    handle to open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Funktionsdatei auswählen
[file,path] = uigetfile({'*.cdp','CDPA Datei (*.cdp)'},handles.title);

% falls nicht abgebrochen
if (file)
    % vollständiger Dateiname
    handles.file = fullfile(path,file);
    
    %file öffnen
    fid = fopen(handles.file,'r');

    %falls Öffnen nicht funktioniert
    if (fid==-1)
        warndig(['Die Datei ', handle.file, 'kann nicht geöffnet werden!'],...
            handles.title);
        return
    end
    
    % und lesen
        %data
        handles.edit_psi.String = fscanf(fid,'%s',1);
        handles.edit_ecc.String = fscanf(fid,'%s',1);
        handles.edit_sbo.String = fscanf(fid,'%s',1);
        handles.edit_Kc.String = fscanf(fid,'%s',1);
        handles.edit_mu.String = fscanf(fid,'%s',1);
        handles.edit_fcm.String = fscanf(fid,'%s',1);
        handles.edit_ecm.String = fscanf(fid,'%s',1);
        handles.edit_d_c.String = fscanf(fid,'%s',1);
        handles.edit_decu.String = fscanf(fid,'%s',1);
        handles.edit_points_c.String = fscanf(fid,'%s',1);
        handles.edit_lchar.String = fscanf(fid,'%s',1);
        handles.edit_gfi_ratio.String = fscanf(fid,'%s',1);
        handles.checkbox_soft.Value = fscanf(fid,'%i',1);
        handles.edit_fctm.String = fscanf(fid,'%s',1);
        handles.edit_gfi.String = fscanf(fid,'%s',1);
        handles.edit_d_t.String = fscanf(fid,'%s',1);
        handles.edit_points_t.String = fscanf(fid,'%s',1);
        handles.edit_beta.String = fscanf(fid,'%s',1);
        
        fclose(fid)
    % Menüpunkt save einschalten
    handles.save.Enable = 'on';
         
    % handles save
    guidata(hObject,handles);
    
    execute(handles)  
end

% --------------------------------------------------------------------
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%file öffnen
fid = fopen(handles.file,'wt');

%falls Öffnen nicht funktioniert
if (fid==-1)
    warndig(['Die Datei ', handle.file, 'kann nicht geöffnet werden!'],...
        handles.title);
    return
end

%data
fprintf(fid,'%s\n',handles.edit_psi.String);
fprintf(fid,'%s\n',handles.edit_ecc.String);
fprintf(fid,'%s\n',handles.edit_sbo.String);
fprintf(fid,'%s\n',handles.edit_Kc.String);
fprintf(fid,'%s\n',handles.edit_mu.String);
fprintf(fid,'%s\n',handles.edit_fcm.String);
fprintf(fid,'%s\n',handles.edit_ecm.String);
fprintf(fid,'%s\n',handles.edit_d_c.String);
fprintf(fid,'%s\n',handles.edit_decu.String);
fprintf(fid,'%s\n',handles.edit_points_c.String);
fprintf(fid,'%s\n',handles.edit_lchar.String);
fprintf(fid,'%s\n',handles.edit_gfi_ratio.String);
fprintf(fid,'%i\n',handles.checkbox_soft.Value);
fprintf(fid,'%s\n',handles.edit_fctm.String);
fprintf(fid,'%s\n',handles.edit_gfi.String);
fprintf(fid,'%s\n',handles.edit_d_t.String);
fprintf(fid,'%s\n',handles.edit_points_t.String);
fprintf(fid,'%s\n',handles.edit_beta.String);

fclose(fid)

% Menüpunkt save einschalten
handles.save.Enable = 'on';

% --------------------------------------------------------------------
function save_as_Callback(hObject, eventdata, handles)
% hObject    handle to save_as (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Funktionsdatei auswählen
[file,path] = uiputfile({'*.cdp','CDPA Datei (*.cdp)'},handles.title);

% falls nicht abgebrochen
if (file)
    % vollständiger Dateiname
    handles.file = fullfile(path,file);
    % und save
    save_Callback(handles.save,eventdata,handles);
    % handles save
    guidata(hObject,handles);
end

% --------------------------------------------------------------------
function export_Callback(hObject, eventdata, handles)
% hObject    handle to export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global sc_iu_out

global omc_iu_out

global st_iu_out

global omt_iu_out

% Funktionsdatei auswählen
[file,path] = uiputfile({'*.txt','Textdatei (*.txt)'},handles.title);

% falls nicht abgebrochen
if (file)
    % vollständiger Dateiname
    handles.file = fullfile(path,file);
          
    %file öffnen
    fid = fopen(handles.file,'wt');

    %falls Öffnen nicht funktioniert
    if (fid==-1)
        warndig(['The file ', handle.file, 'can not be opended!'],...
            handles.title);
        return
    end
    
    %data
    fprintf(fid,'*Material, name=CDPA\n');
    fprintf(fid,'*Elastic\n');
    fprintf(fid,'%s, ',handles.edit_ecm.String);
    fprintf(fid,'0.2\n');
    fprintf(fid,'*Concrete Damaged Plasticity\n');
    fprintf(fid,'%s, ',handles.edit_psi.String);
    fprintf(fid,'%s, ',handles.edit_ecc.String);
    fprintf(fid,'%s, ',handles.edit_sbo.String);
    fprintf(fid,'%s, ',handles.edit_Kc.String);
    fprintf(fid,'%s\n',handles.edit_mu.String);      
    fprintf(fid,'*Concrete Compression Hardening\n');
    fprintf(fid,'%4f, %4f\n',sc_iu_out);
    fprintf(fid,'*Concrete Tension Stiffening, type=DISPLACEMENT\n');
    fprintf(fid,'%4f, %4f\n',st_iu_out);
    fprintf(fid,'*Concrete Compression Damage\n');
    fprintf(fid,'%4f, %4f\n',omc_iu_out);
    fprintf(fid,'*Concrete Tension Damage, type=DISPLACEMENT\n');
    fprintf(fid,'%4f, %4f\n',omt_iu_out);

    fclose(fid);
end

% --------------------------------------------------------------------
function close_Callback(hObject, eventdata, handles)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

button = questdlg('Do you really want to close the program?',...
    handles.title,'Yes','No','No');

switch button
    case 'Yes'
        delete(handles.figure1)
    case 'No'
end
    

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
close_Callback(handles.close, eventdata, handles)



function edit_d_c_Callback(hObject, eventdata, handles)
% hObject    handle to edit_d_c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_d_c as text
%        str2double(get(hObject,'String')) returns contents of edit_d_c as a double

val = str2num(hObject.String);
mini = handles.slider_d_c.Min;
maxi = handles.slider_d_c.Max;
if (val>maxi)
    val = maxi;
end
if (val<mini)
    val = mini;
end
hObject.String = num2str(val,2);
% handles.line.LineWidth = val;
handles.slider_d_c.Value = val;

execute(handles)



% --- Executes during object creation, after setting all properties.
function edit_d_c_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_d_c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider_d_c_Callback(hObject, eventdata, handles)
% hObject    handle to slider_d_c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

val = hObject.Value;
% handles.line.LineWidth = val;
handles.edit_d_c.String = num2str(val,2);

execute(handles)


% --- Executes during object creation, after setting all properties.
function slider_d_c_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_d_c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_points_c_Callback(hObject, eventdata, handles)
% hObject    handle to edit_points_c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_points_c as text
%        str2double(get(hObject,'String')) returns contents of edit_points_c as a double

val = round(str2num(hObject.String));
mini = handles.slider_points_c.Min;
maxi = handles.slider_points_c.Max; 
if (val>maxi)
    val = maxi;
end
if (val<mini)
    val = mini;
end
hObject.String = num2str(val,0);
% handles.line.LineWidth = val;
handles.slider_points_c.Value = val;

execute(handles)

% --- Executes during object creation, after setting all properties.
function edit_points_c_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_points_c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider_points_c_Callback(hObject, eventdata, handles)
% hObject    handle to slider_points_c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

val = hObject.Value;
% handles.line.LineWidth = val;
handles.edit_points_c.String = num2str(val,0);

execute(handles)


% --- Executes during object creation, after setting all properties.
function slider_points_c_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_points_c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_gfi_ratio_Callback(hObject, eventdata, handles)
% hObject    handle to edit_gfi_ratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_gfi_ratio as text
%        str2double(get(hObject,'String')) returns contents of edit_gfi_ratio as a double

val = round(str2num(hObject.String));

mini = handles.slider_gfi_ratio.Min;
maxi = handles.slider_gfi_ratio.Max; 
if (val>maxi)
    val = maxi;
end
if (val<mini)
    val = mini;
end
hObject.String = num2str(val,3);
% handles.line.LineWidth = val;
handles.slider_gfi_ratio.Value = val;

execute(handles)

% --- Executes during object creation, after setting all properties.
function edit_gfi_ratio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_gfi_ratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider_gfi_ratio_Callback(hObject, eventdata, handles)
% hObject    handle to slider_gfi_ratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

val = hObject.Value;
% handles.line.LineWidth = val;
handles.edit_gfi_ratio.String = num2str(val,3);

execute(handles)

% --- Executes during object creation, after setting all properties.
function slider_gfi_ratio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_gfi_ratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_lchar_Callback(hObject, eventdata, handles)
% hObject    handle to edit_lchar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_lchar as text
%        str2double(get(hObject,'String')) returns contents of edit_lchar as a double
execute(handles)

% --- Executes during object creation, after setting all properties.
function edit_lchar_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_lchar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_gfi_Callback(hObject, eventdata, handles)
% hObject    handle to edit_gfi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_gfi as text
%        str2double(get(hObject,'String')) returns contents of edit_gfi as a double
execute(handles)

% --- Executes during object creation, after setting all properties.
function edit_gfi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_gfi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_fctm_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fctm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_fctm as text
%        str2double(get(hObject,'String')) returns contents of edit_fctm as a double
execute(handles)


% --- Executes during object creation, after setting all properties.
function edit_fctm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fctm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_points_t_Callback(hObject, eventdata, handles)
% hObject    handle to edit_points_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_points_t as text
%        str2double(get(hObject,'String')) returns contents of edit_points_t as a double

val = round(str2num(hObject.String));
mini = handles.slider_points_t.Min;
maxi = handles.slider_points_t.Max; 
if (val>maxi)
    val = maxi;
end
if (val<mini)
    val = mini;
end
hObject.String = num2str(val,0);
handles.slider_points_t.Value = val;


val = hObject.Value;
handles.edit_points_t.String = num2str(val,0);

execute(handles)

% --- Executes during object creation, after setting all properties.
function edit_points_t_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_points_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider_points_t_Callback(hObject, eventdata, handles)
% hObject    handle to slider_points_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

val = hObject.Value;
handles.edit_points_t.String = num2str(val,0);

execute(handles)


% --- Executes during object creation, after setting all properties.
function slider_points_t_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_points_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_d_t_Callback(hObject, eventdata, handles)
% hObject    handle to edit_d_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_d_t as text
%        str2double(get(hObject,'String')) returns contents of edit_d_t as a double

val = str2num(hObject.String);
mini = handles.slider_d_t.Min;
maxi = handles.slider_d_t.Max;
if (val>maxi)
    val = maxi;
end
if (val<mini)
    val = mini;
end
hObject.String = num2str(val,2);
handles.slider_d_t.Value = val;

execute(handles)

% --- Executes during object creation, after setting all properties.
function edit_d_t_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_d_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider_d_t_Callback(hObject, eventdata, handles)
% hObject    handle to slider_d_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

val = hObject.Value;
handles.edit_d_t.String = num2str(val,2);

execute(handles)

% --- Executes during object creation, after setting all properties.
function slider_d_t_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_d_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_beta_Callback(hObject, eventdata, handles)
% hObject    handle to edit_beta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_beta as text
%        str2double(get(hObject,'String')) returns contents of edit_beta as a double

val = str2num(hObject.String);
mini = handles.slider_beta.Min;
maxi = handles.slider_beta.Max;
if (val>maxi)
    val = maxi;
end
if (val<mini)
    val = mini;
end
hObject.String = num2str(val,4);
handles.slider_beta.Value = val;

execute(handles)

% --- Executes during object creation, after setting all properties.
function edit_beta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_beta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider_beta_Callback(hObject, eventdata, handles)
% hObject    handle to slider_beta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

val = hObject.Value;
handles.edit_beta.String = num2str(val,4);

execute(handles)


% --- Executes during object creation, after setting all properties.
function slider_beta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_beta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_psi_Callback(hObject, eventdata, handles)
% hObject    handle to edit_psi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_psi as text
%        str2double(get(hObject,'String')) returns contents of edit_psi as a double


% --- Executes during object creation, after setting all properties.
function edit_psi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_psi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ecc_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ecc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ecc as text
%        str2double(get(hObject,'String')) returns contents of edit_ecc as a double


% --- Executes during object creation, after setting all properties.
function edit_ecc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ecc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_sbo_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sbo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_sbo as text
%        str2double(get(hObject,'String')) returns contents of edit_sbo as a double


% --- Executes during object creation, after setting all properties.
function edit_sbo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sbo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Kc_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Kc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Kc as text
%        str2double(get(hObject,'String')) returns contents of edit_Kc as a double


% --- Executes during object creation, after setting all properties.
function edit_Kc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Kc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_mu_Callback(hObject, eventdata, handles)
% hObject    handle to edit_mu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_mu as text
%        str2double(get(hObject,'String')) returns contents of edit_mu as a double


% --- Executes during object creation, after setting all properties.
function edit_mu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_mu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_decu_Callback(hObject, eventdata, handles)
% hObject    handle to edit_decu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_decu as text
%        str2double(get(hObject,'String')) returns contents of edit_decu as a double

val = str2num(hObject.String);
mini = handles.slider_decu.Min;
maxi = handles.slider_decu.Max;
if (val>maxi)
    val = maxi;
end
if (val<mini)
    val = mini;
end
hObject.String = num2str(val,2);
handles.slider_decu.Value = val;

execute(handles)

% --- Executes during object creation, after setting all properties.
function edit_decu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_decu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider_decu_Callback(hObject, eventdata, handles)
% hObject    handle to slider_decu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

val = hObject.Value;
handles.edit_decu.String = num2str(val,2);

execute(handles)

% --- Executes during object creation, after setting all properties.
function slider_decu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_decu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
