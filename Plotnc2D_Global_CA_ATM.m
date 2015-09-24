% read ECHAM data//
% by Mao-Lin Shen
format long;close all;
clearvars -except vname filename
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
  addpath /home/uib/earnest/Analysis/matlab/SelfDefine/export_fig
else
	error(' No workaround for THIS Machine!!!!!!!!!!!!');
end
%% parameters
% =======================================================
DSep = '/';

imgpref = 'Glabal2D_NA';
if ~exist('vname','var'), vname='sst'; end
if ~exist('filename','var'), filename='/work/earnest/UPData/FF01_Seasonal_pM01_noDA_pn3_pL6/ACC_SST/ACC_SST_r005_05.nc'; end

temp = strsplit(filename, DSep);
CaseName = temp{end};

%MinVal = -5; MaxVal = 5; CLevelStepVar1 = 2;  ICMAP = 2;
MinVal = -10; MaxVal = 10; CLevelStepVar1 = 2;  ICMAP = 2;
%MinVal = -100; MaxVal = 100; CLevelStepVar1 = 20;  ICMAP = 2/10;
%MinVal = -8; MaxVal = 8; CLevelStepVar1 = 2;  ICMAP = 2;
%MinVal = -50; MaxVal = 50; CLevelStepVar1 = 10;  ICMAP = 0.4;

%RLW = -180; RLE = 180; RLS = -75; RLN = 75; ScreenPos = [68    56   765   452];
RLW = -100; RLE = 40; RLS = 10; RLN = 75; ScreenPos = [68    56   765   452];

%filename = ['/work/earnest/UPData' DSep expDiffStr DSep FileStr ];

Data = readnetcdf_one(filename, vname); var = Data.data; %var(var==Data.att__FillValue) = NaN;
long = Data.lon; lati = Data.lat;
[~, iLW] = min(abs(long - RLW)); [~, iLE] = min(abs(long - RLE));
[~, iLS] = min(abs(lati - RLS)); [~, iLN] = min(abs(lati - RLN));
%iLong = iLW:iLE; iLati = iLS:iLN;
iLong = 1:length(long); iLati = 1:length(lati);
temp = long; 
long(1:end/2) = temp((end/2+1):end); 
long((end/2+1):end) = temp(1:end/2);
long(long>=180) = long(long>=180) -360;
clear temp
temp = var;
var(1:end/2,:) = temp((end/2+1):end,:);
var((end/2+1):end,:) = temp(1:end/2,:);
var(var > 110) = 110;
var(var < -110) = -110;
clear temp


[mlon, mlat] = meshgrid(long, lati);
mmvar = var(iLong,iLati)';
%% plot figure 
FigID = figure; hold on; set(FigID,'position', ScreenPos, 'color', [1 1 1], 'Visible', 'off');hold on;
m_proj('mercator','lat',[RLS RLN],'long',[RLW RLE]);

%[CC, MHC] = m_contour(mlon, mlat, mmvar, 'k-', 'LevelList', MinVal:CLevelStepVar1:MaxVal, 'ShowText','on'); 
[CC, MHC] = m_contour(mlon, mlat, mmvar, 'LevelList', MinVal:CLevelStepVar1:MaxVal, 'ShowText','on');
%[CC, MHC] = m_contour(mlon, mlat, mmvar, 'k-', 'LevelList', MinVal:CLevelStepVar1:MaxVal);

nLines = numel(get(MHC,'Children'));	% get the number of contour lines
pH = get(MHC,'Children');	% control of each contour line
clear CDV; CDV = zeros(1).*NaN;
for iL = 1:nLines
	%CDV(iL) = {[(min(get(pH(iL),'Cdata'))+max(get(pH(iL),'Cdata')))/2]};	% take the value of contour
	if strcmp(get(pH(iL),'type'), 'patch') 
		CDV(iL) = (min(get(pH(iL),'Cdata'))+max(get(pH(iL),'Cdata')))/2;	% take the value of contour
	else
		set(pH(iL),'fontsize', 13)
	end
end
% find CP == 0  ---------------------------------------------------
ZeroHD = pH(CDV>=0);
for izh = 1:length(ZeroHD)
	%zhLo = get(ZeroHD(izh),'XData');zhLa = get(ZeroHD(izh),'YData');
	%HP = plot(zhLo,zhLa,'k');
	%set(HP, 'LineStyle','--', 'LineWidth', 3)
	set(ZeroHD(izh), 'LineWidth', 0.8, 'linestyle', '-')
end
% find CP < 0  ---------------------------------------------------
ZeroHD = pH(CDV<0);
for izh = 1:length(ZeroHD)
	%zhLo = get(ZeroHD(izh),'XData');zhLa = get(ZeroHD(izh),'YData');
	%HP = plot(zhLo,zhLa,'k');
	%set(HP, 'LineStyle','--', 'LineWidth', 3)
	set(ZeroHD(izh), 'LineWidth', 0.8, 'linestyle', '--')
end

m_pcolor(mlon, mlat, mmvar); shading interp
HB = colorbar;

caxis([MinVal MaxVal])
lim = get(gca,'clim');
mycmap = redblue((lim(2)-lim(1))*2*ICMAP);
acmap = mycmap(1:4:end,:); %acmap(end/2:(end/2+1),:) = 1;
colormap(acmap);

set(gca, 'fontsize', 14)
costfile = ['Global' '_domain_gshhs_c'];	
if ~exist([costfile '.mat'], 'file')
	m_gshhs_c('save',costfile);
end
%return

%m_usercoast(costfile,'patch',[0.8 0.8 0.8]);
m_usercoast(costfile,'color','k');
m_grid('box','on','LineStyle','none','tickdir','in','FontSize',20,'fontweight','bold');
% axis normal
title(filename);


figfilename = [imgpref '_' vname '_' CaseName '.tif'];
%mlstightmargin(gca, 0.5);
%IMG_now = getframe(gcf); imwrite(IMG_now.cdata, filename);
%set(gcf,'PaperPositionMode','auto')
eval(['print -dpng ' figfilename ])

eval(['export_fig(''' num2str(FigID) ''', ''' figfilename(1:end-3)  'pdf'' )']);
%eval(['print -f' num2str(1) '  -dpdf -r300 -opengl ''' filename(1:end-3)  'pdf'' ']); 

