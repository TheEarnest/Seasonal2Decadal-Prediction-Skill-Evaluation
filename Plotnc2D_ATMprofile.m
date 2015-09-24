% read ECHAM data//
% by Mao-Lin Shen
format long;close all;
clearvars -except Cvname Svname Cfilename Sfilename CCLevelStepVar1 SCLevelStepVar1 CMinVal CMaxVal SMinVal SMaxVal


if ~isempty(strfind(computer,'MACI64')) || ~isempty(strfind(computer,'Darwin'))
	load /Users/earnestshen/Cloudsky/SkyDrive/Program/Matlab/work/Z_MatFiles/MacBook_Workaround
	addpath /Users/earnestshen/Cloudsky/SkyDrive/Program/Matlab/SelfDefine/export_fig
	addpath /Users/earnestshen/Cloudsky/SkyDrive/Program/Matlab/SelfDefine/m_map
	eval(['addpath ' SDpath])
elseif ~isempty(strfind(computer,'PCWIN64'))
	load D:\MyCloud\SkyDrive\Program\Matlab\work\Z_MatFiles\VAIO_Workaround
	addpath D:\MyCloud\SkyDrive\Program\Matlab\SelfDefine\export_fig
	addpath D:\MyCloud\SkyDrive\Program\Matlab\SelfDefine\m_map
	addpath D:\MyCloud\SkyDrive\Program\Matlab\SelfDefine\export_fig
	eval(['addpath ' SDpath])
elseif ~isempty(strfind(computer,'GLNXA64'))
  addpath /home/uib/earnest/tools/Matlab_SelfDefine
  addpath /home/uib/earnest/tools/Matlab_SelfDefine/m_map
  addpath /home/uib/earnest/tools/Matlab_SelfDefine/mlsDefined
  addpath /home/uib/earnest/tools/Matlab_SelfDefine/strings
  addpath /home/uib/earnest/Analysis/matlab/SelfDefine/seawater_ver3_3_1
else
	error(' No workaround for THIS Machine!!!!!!!!!!!!');
end
%% parameters
% =======================================================
DSep = '/';

imgpref = 'TSProfile';
% -----------------------------------------------------------------------
if ~exist('Tvname','var'), Tvname='templvl'; end
if ~exist('Svname','var'), Svname='salnlvl'; end
if ~exist('Tfilename','var'), Cfilename='/work/earnest/temp/SPG_templvl_NorCPM_F19_tn21_ensmean_1985_1994.nc'; end
if ~exist('Sfilename','var'), Sfilename='/work/earnest/temp/SPG_salnlvl_NorCPM_F19_tn21_ensmean_1985_1994.nc'; end

Tmax = 30; Tmin = 0; meshN = 256;
Smax = 36; Smin = 30; 
% -----------------------------------------------------------------------
ScreenPos = [68 56 1111 443];

temp = strsplit(Tfilename, DSep);
CaseName = temp{end}

Data = readnetcdf_one(Tfilename, Cvname); Tvar = Data.data; Tvar(Tvar < -100) = NaN;
Data = readnetcdf_one(Sfilename, Svname); Svar = Data.data; Svar(Svar < -100) = NaN;

temperatureAxis = linspace(Tmin, Tmax, meshN);
salinityAxis = linspace(Smin, Smax, meshN);
[mS, mT] = meshgrid(salinityAxis, temperatureAxis);

mTvar = squeeze(Tvar)';
mSvar = squeeze(Svar)';
return
% -----------------------------------------------------------------------
%% plot figure 
%
FigID = figure; hold on; set(FigID,'position', ScreenPos, 'color', [1 1 1], 'Visible', 'off');hold on;


[CC, MHC] = contour(mx, mdepth, mCvar, 'k-', 'LevelList', CMinVal:CCLevelStepVar1:CMaxVal, 'ShowText','on'); 
% -----------------------------------------------------------------------
nLines = numel(get(MHC,'Children'));  % get the number of contour lines
pH = get(MHC,'Children'); % control of each contour line
clear CDV; CDV = zeros(1).*NaN;
for iL = 1:nLines
  %CDV(iL) = {[(min(get(pH(iL),'Cdata'))+max(get(pH(iL),'Cdata')))/2]}; % take the value of contour
  if strcmp(get(pH(iL),'type'), 'patch')
    CDV(iL) = (min(get(pH(iL),'Cdata'))+max(get(pH(iL),'Cdata')))/2;  % take the value of contour
  else
    set(pH(iL),'fontsize', 10)
  end
end
% find CP == 0  ---------------------------------------------------
ZeroHD = pH(CDV>=0);
for izh = 1:length(ZeroHD)
  %zhLo = get(ZeroHD(izh),'XData');zhLa = get(ZeroHD(izh),'YData');
  %HP = plot(zhLo,zhLa,'k');
  %set(HP, 'LineStyle','--', 'LineWidth', 3)
  set(ZeroHD(izh), 'LineWidth', 1.5, 'linestyle', '-')
end
% find CP < 0  ---------------------------------------------------
ZeroHD = pH(CDV<0);
for izh = 1:length(ZeroHD)
  %zhLo = get(ZeroHD(izh),'XData');zhLa = get(ZeroHD(izh),'YData');
  %HP = plot(zhLo,zhLa,'k');
  %set(HP, 'LineStyle','--', 'LineWidth', 3)
  set(ZeroHD(izh), 'LineWidth', 1.5, 'linestyle', '--')
end
% -----------------------------------------------------------------------

pcolor(mx, mdepth, mSvar); shading interp
HB = colorbar;

caxis([SMinVal SMaxVal])
shadingIntDiv = (SMaxVal-SMinVal)/SCLevelStepVar1;
lim = get(gca,'clim');
%mycmap = redblue((lim(2)-lim(1))*4*shadingIntDiv);
mycmap = redblue(4*shadingIntDiv);
acmap = mycmap(1:4:end,:); %acmap(end/2:(end/2+1),:) = 1;
colormap(acmap);

yHL = ylabel('Depth (m)', 'fontsize', 16);
if (length(long) > 1)
  set(gca, 'XTick', -170:10:170)
  set(gca, 'XTickLabel', ['170W';'160W';'150W';'140W';'130W';'120W';'110W'; ...
                          '110W';'100W';' 90W';' 80W';' 70W';' 60W';' 50W'; ...
                          ' 40W';' 30W';' 20W';' 10W';'   0';' 10E';' 20E'; ...
                          ' 30E';' 40E';' 50E';' 60E';' 70E';' 80E';' 90E'; ...
                          '100E';'110E';'120E';'130E';'140E';'150E';'160E'; ...
                          '170E'])
elseif (length(lati) > 1)
  set(gca, 'XTick', -70:10:70)
  set(gca, 'XTickLabel', ['70S';'60S';'50S';'40S';'30S';'20S';'10S';'  0'; ...
                          '10N';'20N';'30N';'40N';'50N';'60N';'70N'; ...
                          ])
end

axis([min(mx(:)) max(mx(:)) 0 MaxDepth])
box on;
set(gca, 'fontsize', 14, 'YDir', 'reverse', 'linewidth', 2)
%return

title(Cfilename);


figfilename = [imgpref '_' Cvname '_' CaseName '.png']
%mlstightmargin(gca, 0.5);
%IMG_now = getframe(gcf); imwrite(IMG_now.cdata, filename);
set(gcf,'PaperPositionMode','auto')
eval(['print -dpng ' figfilename ])

%eval(['export_fig(''' num2str(FigID) ''', ''' filename(1:end-3)  'pdf'' )']);
%eval(['print -f' num2str(1) '  -dpdf -r300 -opengl ''' filename(1:end-3)  'pdf'' ']); 


