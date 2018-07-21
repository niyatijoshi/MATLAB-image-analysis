clc
clearvars -except axis_intercept device_dimension theta
close all
%% Get avi video file to be processed for image analysis.
filename = uigetfile;
%% Inserting parameters to be read during image analysis. 
prompt = {'Start Frame: ','Stop Frame: ','Number: ','Start Analysis at: '};
dlg_title = 'Input';
num_lines=1;
answer = inputdlg(prompt, dlg_title,num_lines);
%% A dialog box is created with above mentioned commands and input to this box is accepted using the commands below.
start_frame= str2num(answer{1}); %  This frame is used as a reference frame for ellipse detection using image subtraction. 
stop_frame = str2num(answer{2});% This is mainly to limit the number of frames processed.
drop_num = str2num(answer{3});% Ellipse number is used to make sure that center of device coordinates are consistent throughout the video 
%file. The coordinates are calculated only for first.
start_at = str2num(answer{4});
%% The user input of the frame numbers are required from the above commands to tell the code the number of frames that are to be read. 
%% Reading avi file
v = VideoReader(filename); %% Reading avi video file
start = start_at; % starting analysis from this frame
fstart =start_frame; 
fstop = stop_frame;
Num = fstop-fstart+1; 
%% Creating 3D image array
for i= fstart:fstop  
ind = i-fstart+1;
img(:,:,ind) = read(v,i);
i
end
%%
if mod(Num,2)==0
    A=Num;
else
    A=Num+1;
end
%%
i=A/2
I  = img(:,:,i);
figure(100);uiwait(msgbox('Adjust to get low_in and high_in values')); 
imshow(I)
imcontrast
%%
lowin = 18;
highin = 64;
%% Image reference and changing contrast for ref image
Imref = img(:,:,1); % reference image
imref = imadjust(Imref,[lowin highin]/255);
figure();imshow(imref)
%% User input on first frame ellipse centroid
uiwait(msgbox('Click on the approximate centroid of the ellipse in the background subtracted image'));
i = start-fstart+1;
I = img(:,:,i)-imref; 
figure(100);imshow(I);
pts = ginput(1);
centroidx_select_prev = pts(1);
centroidy_select_prev = pts(2);
%% Changing contrast for all images in the image array
for i=start-fstart+1:Num 
 K = img(:,:,i);
 Kadj(:,:,i) = imadjust(K,[lowin highin]/255);
end
%% Image processing and ellipse
center detection
stat_num = 1;
for i=start-fstart+1:Num
    i
   I = Kadj(:,:,i)-imref;
%  I = img(:,:,i)-imref; %% Image subtraction
%  Iadj = imadjust(I,[lowin highin]/255);
 level = graythresh(I);
 Ibw2 = im2bw(I,level);
 stats= regionprops(Ibw2,'Centroid','MajorAxisLength','MinorAxisLength','Orientation','PixelList','Area');
%  stats= regionprops(Iadj,'Centroid','MajorAxisLength','MinorAxisLength','Orientation','PixelList','Area');
 figure(200);
 hold off
 imshow(img(:,:,i))
%    imshow(I)  
     cell_data = struct2cell(stats);
     centroid_cell = cell2mat(cell_data(2,:));
     centroidx_select = centroid_cell(1:2:length(centroid_cell)-1);
     centroidy_select = centroid_cell(2:2:length(centroid_cell));
     area_select = cell2mat(cell_data(1,:));
          
     
     if i>start-fstart+1 % Assuming that regionprops managed to only detect the ellipse
      distance_diff = sqrt((centroidx_select_prev-centroidx_select).^2+ (centroidy_select_prev-centroidy_select).^2);
      area_diff  = abs(area_select_prev-area_select);
      Error = sqrt(distance_diff.^2+area_diff.^2);
      stat_num = find(Error==min(Error));
     else
          distance_diff = sqrt((centroidx_select_prev-centroidx_select).^2+ (centroidy_select_prev-centroidy_select).^2);
           Error = distance_diff;
           stat_num = find(Error==min(Error));
     end
     centroidx_select_prev = centroidx_select(stat_num);
     centroidy_select_prev = centroidy_select(stat_num);
     area_select_prev = area_select(stat_num);

     
      
 stats_final(i) = stats(stat_num);
 pixel_list = stats_final(i).PixelList
 hold on;plot(pixel_list(:,1),pixel_list(:,2),'r*')
 pause(0.0001)
end
%% Extracting centroid trajectory from structure
count = 1
for i=start-fstart+1:Num
    i
    centroid(count,1:2) = stats_final(i).Centroid;
    
    major_axis(count) = stats_final(i).MajorAxisLength;
    minor_axis(count) = stats_final(i).MinorAxisLength;
    area(count) = stats_final(i).Area;
    Orientation(count) = stats_final(i).Orientation;
    count = count+1;
end
centroidx = centroid(:,1);
centroidy = centroid(:,2);
%% Plotting the trajectory
i=1;
I = img(:,:,i);
figure(100);imshow(I)
hold on;
plot(centroidx,centroidy,'r*')
num = length(centroidx);
resoln = 1.75; % pixel resolution
strain_radius = 200/resoln; %% in microns

checks = sqrt((centroidx-axis_intercept(1)).^2+(centroidy-axis_intercept(2)).^2);
ind = (checks<strain_radius);
new_centroidx = centroidx(ind);
new_centroidy = centroidy(ind);
new_major_axis = major_axis(ind);
new_minor_axis = minor_axis(ind);
new_area = area(ind);
new_orientation = Orientation(ind);
figure(100);hold on;plot(centroidx,centroidy,'r*');
plot(new_centroidx,new_centroidy,'bo','MarkerFaceColor','b')
%% Straight line fit to the horizontal portion of the ellipse movement
new_num = length(new_centroidx);
xtry = new_centroidx(1:floor(new_num/10)); % only first 10% of data points
ytry = new_centroidy(1:floor(new_num/10));
[xData, yData] = prepareCurveData(xtry,ytry);
ft = fittype('poly1');
opts = fitoptions('Method','LinearLeastSquares');
[fitresult,gof] = fit(xData,yData,ft,opts);
p1_entry = fitresult.p1; % Slope of line
p2_entry = fitresult.p2; % X intercept
%% Straight line fit to the vertical portion of the trajectory 
new_centroidx_rev = new_centroidx(end:-1:1); % flip the data set
new_centroidy_rev = new_centroidy(end:-1:1);
xtry_rev = new_centroidx_rev(1:floor(new_num/10)); % only first 10% of data points
ytry_rev = new_centroidy_rev(1:floor(new_num/10));
[xData, yData] = prepareCurveData(xtry_rev,ytry_rev);
ft = fittype('poly1');
opts = fitoptions('Method','LinearLeastSquares');
[fitresult,gof] = fit(xData,yData,ft,opts);
p1_exit= fitresult.p1;  % Slope of line
p2_exit = fitresult.p2; % Y intercept
%% Plotting the fitted straight lines 
[r c d] = size(img);
x_check = 1:c;
y_checkent = p1_entry.*x_check+p2_entry;
y_checkext = p1_exit.*x_check+p2_exit;
figure(101);hold off;imshow(imref);
hold on;plot(new_centroidx,new_centroidy,'bo','MarkerFaceColor','b')
plot(x_check,y_checkent,x_check,y_checkext)

distances = sqrt((xcenter(1)-new_centroidx).^2+(xcenter(2)-new_centroidy).^2);
ind_short = find(distances==min(distances));
xshort = new_centroidx(ind_short);
yshort = new_centroidy(ind_short);
figure(101); hold on;plot(xshort,yshort,'s','MarkerFaceColor','y')
%%  Rotation and exponential fit
theta1 = atan(p1_entry); 
theta1 = theta1*-1; 
centroidx_rotate = new_centroidx.*cos(theta1)-new_centroidy.*sin(theta1); % rotating the trajectory
centroidy_rotate = new_centroidx.*sin(theta1)+new_centroidy.*cos(theta1);
xshort_rot = xshort.*cos(theta1)-yshort.*sin(theta1); % Rotating the bifurcating point
yshort_rot = xshort.*sin(theta1)+yshort.*cos(theta1);
xcenter_rot_approach = xcenter(1).*cos(theta1)-xcenter(2).*sin(theta1); % Rotating the intersectio pointg
ycenter_rot_approach = xcenter(1).*sin(theta1)+xcenter(2).*cos(theta1);
figure(101);plot(centroidx_rotate,centroidy_rotate,xshort_rot,yshort_rot,'r*',xcenter_rot_approach,ycenter_rot_approach,'d','MarkerFaceColor','k')
%%  line1
ind_approach = 1:(ind_short-1);  % Changed
x_approach = centroidx_rotate(ind_approach);
y_approach = centroidy_rotate(ind_approach);
major_approach = new_major_axis(ind_approach);
minor_approach = new_minor_axis(ind_approach);
area_approach = new_area(ind_approach);
orientation_approach = new_orientation(ind_approach);
figure(101);hold on;plot(x_approach,y_approach,'co')
%% Exponetial data prep for line1
time_approach = (0:length(x_approach)-1)./v.FrameRate; % time series in seconds
x_exp =abs(x_approach-xcenter_rot_approach);  % position series wrt the intersection point
figure();plot(time_approach,x_exp)
[xData, yData] = prepareCurveData( time_approach, x_exp );
% Set up fittype and options.
ft = fittype( 'exp1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
% Fit model to data.
[fitresult_approach, gof_approach] = fit( xData, yData, ft, opts );
% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult_approach, xData, yData );
legend( h, 'x vs. time', 'untitled fit 1', 'Location', 'NorthEast' );
% Label axes
xlabel time
ylabel x_exp
title(['Rsquare = ' num2str(gof_approach.rsquare) ',    Entering Strain rate,G: ' num2str(fitresult_approach.b)])
grid on
%%  Rotation for line2
theta2 = atan(p1_exit);
theta2 = theta2*-1;  % Changed
centroidx_rotate_exit = new_centroidx.*cos(theta2)-new_centroidy.*sin(theta2);
centroidy_rotate_exit = new_centroidx.*sin(theta2)+new_centroidy.*cos(theta2);
xshort_rot_exit = xshort.*cos(theta2)-yshort.*sin(theta2);
yshort_rot_exit = xshort.*sin(theta2)+yshort.*cos(theta2);
xcenter_rot_exit = xcenter(1).*cos(theta2)-xcenter(2).*sin(theta2);
ycenter_rot_exit = xcenter(1).*sin(theta2)+xcenter(2).*cos(theta2);
figure(102); imshow(imref); hold on;plot(new_centroidx,new_centroidy,'r');
plot(centroidx_rotate_exit,centroidy_rotate_exit,'b',xshort_rot_exit,yshort_rot_exit,'r*')
%% Selecting the line2
% ind_exit = (centroidy_rotate>yshort_rot);
ind_exit = (ind_short):length(centroidx_rotate_exit); % Changed
x_exit = centroidx_rotate_exit(ind_exit);
y_exit= centroidy_rotate_exit(ind_exit);
major_exit = new_major_axis(ind_exit);
minor_exit = new_minor_axis(ind_exit);
area_exit = new_area(ind_exit);
orientation_exit = new_orientation(ind_exit);
figure(102);hold on;plot(x_exit,y_exit,'co')
%%  Exponetial data prep for vertical/departure trajectory
time_exit = (0:length(x_exit)-1)./v.FrameRate;
y_exp =abs(x_exit-xcenter_rot_exit); % Changed % Since we rotated the departure trajectory parallel to x axis
figure();plot(time_exit,y_exp)
[xData, yData] = prepareCurveData( time_exit, y_exp );
% Set up fittype and options.
ft = fittype( 'exp1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
% Fit model to data.
[fitresult_exit, gof_exit] = fit( xData, yData, ft, opts );
% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult_exit, xData, yData );
legend( h, 'y vs. time', 'untitled fit 1', 'Location', 'NorthEast' );
% Label axes
xlabel time
ylabel y_exp
title(['Rsquare = ' num2str(gof_exit.rsquare) ',    Exit Strain rate,G: ' num2str(fitresult_exit.b)])
grid on
%%  Excel export
 parameterz = {'Centroidx','Centroidy','Approach: t','Approach: x', 'Major Axis Length', ' Minor Axis Length','Area','Orientation','R-Approach', 'Exit: t','Exit: y', 'Major Axis Length', ' Minor Axis Length','Area','Orientation','R-Exit'};
 drop_num;
 xlswrite('data',parameterz,drop_num,'A1');
 % Trajectory
 xlswrite('data',centroidx,drop_num,'A2');
 xlswrite('data',centroidy,drop_num,'B2');
 % Approach
 xlswrite('data',time_approach',drop_num,'C2');
 xlswrite('data',x_exp,drop_num,'D2');
 xlswrite('data',((1.75*major_approach)/2)',drop_num,'E2');
 xlswrite('data',((1.75*minor_approach)/2)',drop_num,'F2');
 xlswrite('data',area_approach',drop_num,'G2');
 xlswrite('data',orientation_approach',drop_num,'H2');
 xlswrite('data',1.75*(sqrt((major_approach.*minor_approach)/4))',drop_num,'I2');
 % Exit
 xlswrite('data',time_exit',drop_num,'J2');
 xlswrite('data',y_exp,drop_num,'K2');
 xlswrite('data',((1.75*major_exit)/2)',drop_num,'L2');
 xlswrite('data',((1.75*minor_exit)/2)',drop_num,'M2');
 xlswrite('data',area_exit',drop_num,'N2');
 xlswrite('data',orientation_exit',drop_num,'O2');
 xlswrite('data',1.75*(sqrt((major_exit.*minor_exit)/4))',drop_num,'P2');

clear I Ibw2 img imref 
save('data')
clear I Ibw2 img imref 
save(['data' num2str(drop_num)])











