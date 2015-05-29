function [] = PitchAnalyzer()

% Tool to analyze the pitch from a speaker in realtime.
%
% This program takes usage of the function msound for the real time
% signal processing. You can get it from the jade-hochschule-site:
% http://tgm.jade-hs.de/web/file/Institut_fr_Hrtechnik_und_Audiologie/Software.php
%
% Usage: to start the tool, you have to run the script (F5).
%
%------------------------------------------------------------------------
% Authors: Jan Heimann, Sven Herrmann, Tanja Stüttchen
%          (c) TGM @ Jade Hochschule applied licence see EOF
%
% Version History:
% Ver. 0.01            13-Dez-2014       Heimann, Herrmann, Stüttchen
% Ver. 1.00            12-Jan-2015       Heimann, Herrmann, Stüttchen
%------------------------------------------------------------------------
clc;
close all;
clear all;

% gui
handles.colorOn = [0.8 0.8 0.8];
handles.colorOff = [0.6 0.6 0.6];
handles.colorMain = [224 102 255]/255;

% figure
handles.fig1 = figure(1);
set(handles.fig1,'menubar','none','toolbar','none','color',handles.colorMain,...
    'units','normalized','position',[0.15 0.2 0.7 0.6],'name','Realtime PitchAnalyzer',...
    'numbertitle','off');

% panels and axes for plots
handles.panel1 = uipanel(handles.fig1,'position',[0.4 0 0.6 0.5]);
handles.panel2 = uipanel(handles.fig1,'position',[0.4 0.5 0.6 0.5]);

handles.axes_unten = axes('parent',handles.panel1,'position',[0 0.1 1 0.9]);
axis off;
handles.axes_oben = axes('parent',handles.panel2,'position',[0 0.1 1 0.9]);
axis off;

% start/ stop button
handles.button1 = uicontrol(handles.fig1,'units','normalized','position',...
    [0.025 0.475 0.075 0.075],'backgroundcolor',handles.colorOff,'callback',...
    {@buttonCallbackFcn1},'string','Start','value',1);

% readme button
handles.button2 = uicontrol(handles.fig1,'units','normalized','position',...
    [0.945 0.94 0.05 0.05],'backgroundcolor',handles.colorOn,'callback',...
    {@readme},'string','Readme');

% slider for smoothing the spectrum
handles.slider_speed = uicontrol(handles.fig1,'style','slider','units',...
    'normalized','position',[0.045 0.1 0.3 0.03],'backgroundcolor',...
    handles.colorMain,'value',0.3,'min',0,'max',1,'sliderstep',[0.1 0.3]/1,...
    'callback',{@sliderCallbackFcn,3});

%slider for f0
handles.slider_f0 = uicontrol(handles.fig1,'style','slider','units',...
    'normalized','position',[0.045 0.9 0.3 0.03],'backgroundcolor',...
    handles.colorMain,'value',2500,'min',500,'max',5000,'sliderstep',[500 1000]/4500,...
    'callback',{@sliderCallbackFcn,1});

%slider for bw
handles.slider_bw = uicontrol(handles.fig1,'style','slider','units',...
    'normalized','position',[0.045 0.8 0.3 0.03],'backgroundcolor',...
    handles.colorMain,'value',5000,'min',100,'max',10000,'sliderstep',[100 1000]/9900,...
    'callback',{@sliderCallbackFcn,2});

%checkboxes
handles.check1 = uicontrol(handles.fig1,'style','checkbox','units',...
    'normalized','position',[0.125 0.6 0.05 0.05],'backgroundcolor',...
    handles.colorMain);
handles.check2 = uicontrol(handles.fig1,'style','checkbox','units',...
    'normalized','position',[0.125 0.5 0.05 0.05],'backgroundcolor',...
    handles.colorMain);
handles.check3 = uicontrol(handles.fig1,'style','checkbox','units',...
    'normalized','position',[0.125 0.4 0.05 0.05],'backgroundcolor',...
    handles.colorMain);

% pulldownmenue for fs
handles.pulldownFs = uicontrol(handles.fig1,'style','popupmenu','units',...
    'normalized','position',[0.045 0.2 0.1 0.1],'backgroundcolor',...
    handles.colorMain,'string',{'8000','16000','22050','33075','44100'},...
    'value',3);

% pulldownmenue for bz
handles.pulldownBz = uicontrol(handles.fig1,'style','popupmenu','units',...
    'normalized','position',[0.243 0.2 0.1 0.1],'backgroundcolor',...
    handles.colorMain,'string',{'512','1024','2048','4096','8192'},...
    'value',3);

% diverse textfields
handles.text_f0 = uicontrol(handles.fig1,'style','text','units','normalized',...
    'position',[0.35 0.9 0.05 0.03],'backgroundcolor',handles.colorMain,'foregroundcolor',...
    [1 1 1],'string',get(handles.slider_f0,'value'));
handles.text_f0_2 = uicontrol(handles.fig1,'style','text','units','normalized',...
    'position',[0.002 0.9 0.03 0.03],'backgroundcolor',handles.colorMain,'foregroundcolor',...
    [1 1 1],'string','F0');
handles.text_bw = uicontrol(handles.fig1,'style','text','units','normalized',...
    'position',[0.35 0.8 0.05 0.03],'backgroundcolor',handles.colorMain,'foregroundcolor',...
    [1 1 1],'string',get(handles.slider_bw,'value'));
handles.text_bw_2 = uicontrol(handles.fig1,'style','text','units','normalized',...
    'position',[0.002 0.8 0.03 0.03],'backgroundcolor',handles.colorMain,'foregroundcolor',...
    [1 1 1],'string','BW');
handles.text_fast = uicontrol(handles.fig1,'style','text','units','normalized',...
    'position',[0.004 0.1 0.03 0.03],'backgroundcolor',handles.colorMain,'foregroundcolor',...
    [1 1 1],'string','Fast');
handles.text_slow = uicontrol(handles.fig1,'style','text','units','normalized',...
    'position',[0.348 0.1 0.05 0.03],'backgroundcolor',handles.colorMain,'foregroundcolor',...
    [1 1 1],'string','Slow');
handles.text_result = uicontrol(handles.fig1,'style','text','units','normalized',...
    'position',[0.28 0.475 0.1 0.075],'backgroundcolor',handles.colorMain,'foregroundcolor',...
    [0 0 0],'fontsize',18);
handles.cepstrum_based = uicontrol(handles.fig1,'style','text','units','normalized',...
    'position',[0.15 0.61 0.05 0.03],'backgroundcolor',handles.colorMain,'foregroundcolor',...
    [1 1 1],'string','Cepstrum');
handles.akf_based = uicontrol(handles.fig1,'style','text','units','normalized',...
    'position',[0.14 0.51 0.05 0.03],'backgroundcolor',handles.colorMain,'foregroundcolor',...
    [1 1 1],'string','AKF+CC');
handles.amdf_based = uicontrol(handles.fig1,'style','text','units','normalized',...
    'position',[0.14 0.41 0.05 0.03],'backgroundcolor',handles.colorMain,'foregroundcolor',...
    [1 1 1],'string','AMDF');
handles.fpitch = uicontrol(handles.panel1,'style','text','units','normalized',...
    'position',[0 0.9 0.2 0.1],'string','','fontsize',14,'foregroundcolor',[0.125 0.125 0.125]);
handles.time = uicontrol(handles.panel2,'style','text','units','normalized',...
    'position',[0 0 0.2 0.1],'string','','fontsize',14,'foregroundcolor',[0.125 0.125 0.125]);
handles.text_fs = uicontrol(handles.fig1,'style','text','units','normalized',...
    'position',[0.005 0.26 0.03 0.03],'backgroundcolor',handles.colorMain,'foregroundcolor',...
    [1 1 1],'string','FS');
handles.text_bz = uicontrol(handles.fig1,'style','text','units','normalized',...
    'position',[0.2 0.26 0.03 0.03],'backgroundcolor',handles.colorMain,'foregroundcolor',...
    [1 1 1],'string','BZ');
handles.xachse1 = uicontrol(handles.panel1,'style','text','units','normalized',...
    'position',[0.08 0 0.1 0.08],'string','100','fontsize',12,'foregroundcolor',[0.125 0.125 0.125]);
handles.xachse2 = uicontrol(handles.panel1,'style','text','units','normalized',...
    'position',[0.385 0 0.1 0.08],'string','500','fontsize',12,'foregroundcolor',[0.125 0.125 0.125]);
handles.xachse3 = uicontrol(handles.panel1,'style','text','units','normalized',...
    'position',[0.52 0 0.1 0.08],'string','1k','fontsize',12,'foregroundcolor',[0.125 0.125 0.125]);
handles.xachse4 = uicontrol(handles.panel1,'style','text','units','normalized',...
    'position',[0.64 0 0.1 0.08],'string','2k','fontsize',12,'foregroundcolor',[0.125 0.125 0.125]);
handles.xachse1 = uicontrol(handles.panel1,'style','text','units','normalized',...
    'position',[0.77 0 0.1 0.08],'string','4k','fontsize',12,'foregroundcolor',[0.125 0.125 0.125]);
handles.xachse5 = uicontrol(handles.panel1,'style','text','units','normalized',...
    'position',[0.93 0 0.05 0.08],'string','10K','fontsize',12,'foregroundcolor',[0.125 0.125 0.125]);

guidata(handles.button1,handles);
end

% start button
function [] = buttonCallbackFcn1(handle,event)

handles = guidata(handle);

set(handles.button1,'string','Stop','callback',{@buttonCallbackFcn2},...
    'backgroundcolor',handles.colorOn);

% get fs and bz
fs_value = get(handles.pulldownFs,'value');
fs_string = get(handles.pulldownFs,'string');
set(handles.pulldownFs,'enable','off');

bz_value = get(handles.pulldownBz,'value');
bz_string = get(handles.pulldownBz,'string');
set(handles.pulldownBz,'enable','off');

% initial msound
input_devices = 0;
fs = str2double(fs_string(fs_value));
bz = str2double(bz_string(bz_value));
chan = 1;

% preallokation
spectrum = zeros((bz/2),1);
rec2 = zeros(bz,1);
value2 = 0;

% starting msound
msound('close');
msound('openread',input_devices,fs, bz, chan);

% start time measurement
tic;

% play-loop
while get(handles.button1,'value')
    
    handles = guidata(handle);
    
    % calculation of a and b coefficients
    f0 = get(handles.text_f0,'string');
    bw = get(handles.text_bw,'string');
    [b,a] = bpf_filter_coeff(f0,fs,bw);
    
    % calculation of smoothing
    tau = get(handles.slider_speed,'value');
    alpha = exp(-1/(tau*fs/bz));
    
    % recording and filtering
    rec = msound('getsamples');
    rec = filter(b,a,rec);
    
    % smoothing the spectrum
    rec2 = rec2*0.1 + rec*0.9;
    spec = abs(fft(rec));
    spec = spec(1:end/2);
    spectrum = spectrum*alpha + spec*(1-alpha);
    
    % plot spectral signal
    plotline1 = plot(handles.axes_unten,linspace(1,fs,length(spectrum)),spectrum,'linewidth',2);
    xlim(handles.axes_unten,[50 10000]);
    set(handles.axes_unten,'xscale','log');
    axis (handles.axes_unten,'off');
    
    % plot time signal
    plotline2 = plot(handles.axes_oben,linspace(0,1,length(rec)),rec2,'linewidth',2);
    ylim([-1 1]);
    axis (handles.axes_oben,'off');
    
    drawnow;
    
    % use pitch detectors and smoothing the results
    value = pitch_detektor(handles,rec,fs);
    value2 = value2*0.85 + value*0.15;
    
    % only show values up to 750Hz
    if value2 < fs/750
        value2 = 0;
    end
    
    % calculate frequency from the result and show this in textfield
    if value2 ~= 0
        f_pitch = round(fs/value2);
        set(handles.fpitch,'string',[num2str(f_pitch)  ' Hz']);
    else
        set(handles.fpitch,'string','');
        f_pitch = 0;
    end
    
    % decision between man-, women- and child-voice and show this in tf
    if f_pitch == 0
        set(handles.text_result,'string','');
    elseif f_pitch >= 300 && f_pitch <=450
        set(handles.text_result,'string','Kind');
    elseif f_pitch >= 150 && f_pitch < 300
        set(handles.text_result,'string','Frau');
    elseif f_pitch >= 75 && f_pitch < 150
        set(handles.text_result,'string','Mann');
    else
        set(handles.text_result,'string','-');
    end
    
    % get duration and show this in tf
    time = round(toc*10)/10;
    set(handles.time,'string',[num2str(time) ' sec']);
end

msound('close');
% delete(plotline1,plotline2);
end

% stop button
function [] = buttonCallbackFcn2(handle,event)

handles = guidata(handle);

% set value from button1 to 0 to stop the while(play)loop
set(handles.button1,'value',0,'string','Start','callback',{@buttonCallbackFcn1},...
    'backgroundcolor',handles.colorOff);

set(handles.pulldownFs,'enable','on');
set(handles.pulldownBz,'enable','on');
end

% slider callbacks
function [] = sliderCallbackFcn(handle,event,idx)

handles = guidata(handle);

% round slider values and show the values in tf
switch idx
    case 1 %f0 slider
        value = get(handles.slider_f0,'value');
        value = round(value/500)*500;
        set(handles.text_f0,'string',value);
    case 2 % bw slider
        value = get(handles.slider_bw,'value');
        value = round(value/100)*100;
        set(handles.text_bw,'string',value);
    case 3 % speed slider
        % not necessary ;)
end
end

% bandpass (rbj-cookbook)
function [b,a] = bpf_filter_coeff(f0,fs,bw)

% calculate the coefficents
bw = str2double(bw);
f0 = str2double(f0);
w0 = 2*pi*f0/fs;
q = f0/bw;
alpha = sin(w0)/(2*q);

b0 = alpha;
b1 = 0;
b2 = -alpha;
a0 = 1 + alpha;
a1 = -2*cos(w0);
a2 = 1 - alpha;

a = [a0 a1 a2]./a0;
b = [b0 b1 b2]./a0;
end

% pitch detectors
function [value] = pitch_detektor(handles,data,fs)

% lp-filter @400hz (good for all)
data = filter([0.0008 0.0016 0.0008],[1.0000 -1.9194 0.9226],data);

idx = 0;
value = 0;

% cepstrum based
if 	get(handles.check1,'value') == 1
    idx = idx + 1;
    han_win = hamming(length(data));
    y = fft(data.*han_win);
    cepstrum=fft(log(abs(y)+eps));
    cepstrum_part = cepstrum(round(fs/400):round(fs/55));
    pk = max(abs(cepstrum_part));
    loc = find(pk==abs(cepstrum));
    value = value + loc(1);
end

% akf + center clipping
if 	get(handles.check2,'value') == 1
    idx = idx + 1;
    clip = 0.45;
    thresh = max(abs(data)) * clip;
    x_clip = abs(data) < thresh;
    data(x_clip) = 0;
    data(~x_clip) = data(~x_clip) - (thresh*sign(data(~x_clip)));
    akf = xcorr(data);
    akf = akf(floor(length(akf)/2):end);
    han_win = hann(90);
    han_win_long = [han_win(1:45);ones(length(akf)-90,1);han_win(46:90)];
    akf = akf.*han_win_long;
    n_high = floor(fs/500);
    [~, max_pos] = max(akf(n_high:end));
    max_pos = max_pos + n_high - 1;
    value = value + max_pos(1);
end

% amdf
if 	get(handles.check3,'value') == 1
    idx = idx + 1;
    delay_n = round(fs/80);
    diff_fcn = zeros(delay_n,1);
    data_long = [data;diff_fcn];
    for kk = 1:delay_n
        data_delay = [zeros(kk,1); data; zeros(delay_n-kk,1)];
        diff_fcn(kk) = mean(abs(data_long - data_delay));
    end
    han_win = hann(10);
    han_win_long = [han_win(1:5);ones(delay_n-10,1);han_win(6:10)];
    diff_fcn = diff_fcn.*han_win_long;
    [pks,loc] = findpeaks(diff_fcn);
    loc1 = loc(1);
    loc2 = loc(end);
    result = find(min(diff_fcn(loc1:loc2))==diff_fcn(loc1:loc2)) + loc1;
    value = value + result(1);
end

if value ~= 0
    value = value/idx;
end

end

% readme
function [] = readme(handle,event)

% figure
handles.fig2 = figure(2);
set(handles.fig2,'menubar','none','toolbar','none','color',[0.8 0.8 0.8],...
    'units','normalized','position',[0.47 0.2 0.33 0.6],'name','Readme',...
    'numbertitle','off');

% content
s1 = 'Mit Hilfe dieses Programmes lässt sich die Grundfrequenz eines Sprechers ermitteln und verschieden Pitcherkenner können verglichen werden.';
s2 = 'Zusätzlich können spektrale Eigenschaften der Sprache aufgezeigt werden.';
s3 = 'Mit den beiden Slidern für F0 und BW lässt sich ein Bandpassfilter über das aufgenommene Signal legen.';
s4 = 'Mit den jeweiligen Checkboxen kann ein oder mehrere Erkenner gewählt werden. Falls mehrere gewählt werden, wird das Ergebnis zwischen diesen gemittelt.';
s5 = 'Mit dem Start/Stop Button startet/stopt das Program.';
s6 = 'Mit den beiden Pulldownmenüs kann die Abtastrate und die Blocklänge eingestellt werden. Dies muss vor dem Start geschehen.';
s7 = 'Der Fast/Slow-Slider dient einer Glättung des Spektrums.';
s8 = 'In den rechten beiden Panels wird oben das Zeitsignal, unten das Frequenzspektrum dargestellt.';
s9 = 'In der unteren linken Ecke des Zeitsignals wird die Dauer in Sekunden angegeben.';
s10 = 'In der oberen linken Ecke des Frequenzspektrums wird die ermittelte Grundfrequenz des Sprechers angezeigt.';
s11 = 'Anmerkung: Bei einer sehr großen Bz können die pitcher-Erkenner mehrere Sekunden zum "Einschwingen" benötigen.';

% print
handles.text_readme = uicontrol(handles.fig2,'style','text','units','normalized',...
    'position',[0 0 1 1],'backgroundcolor',[0.8 0.8 0.8],'foregroundcolor',...
    [0 0 0],'fontsize',10,'horizontalalignment','left','string',...
    sprintf('%s\n%s\n\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n\n%s',s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11));

% close button
handles.buttonClose = uicontrol(handles.fig2,'units','normalized','position',...
    [0.4 0.1 0.17 0.07],'backgroundcolor',[0.8 0.8 0.8],'callback',...
    'close(2)','string','close');
end

%--------------------Licence ---------------------------------------------
% Copyright (c) <2015> Jan Heimann, Sven Herrmann, Tanja Stüttchen
% Jade University of Applied Sciences
% Permission is hereby granted, free of charge, to any person obtaining
% a copy of this software and associated documentation files
% (the "Software"), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish,
% distribute, sublicense, and/or sell copies of the Software, and to
% permit persons to whom the Software is furnished to do so, subject
% to the following conditions:
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
% IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
% CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
% TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
% SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.