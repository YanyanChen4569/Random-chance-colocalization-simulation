function [Num_Colo_multiple_simulation]=ColocolizationAnalysis(MaskName,Num_Ch1,Num_Ch2,Num_iteration,i)
%%% Batch simulation for random chance colocalization 
%%% Yanyan Chen (08/27/2024; Columbia University; yc4569@cumc.columbia.edu)

% Read the Mask image 

% Folderpath='/Users/chenyanyan/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Confocal Facility YC/Image analysis/Image analysis for Alondra Schweizer Burguete/masks/';
Folderpath='/Users/chenyanyan/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Confocal Facility YC/Image analysis/Image analysis for Alondra Schweizer Burguete/323-2 images/082224_323-2 Images 1,2,3,4,6 masks/';

[A,map]=imread([Folderpath MaskName]); 

% Convert the 0-1 binary image
threshold=0.5;  
binary_image=A>threshold;
A1=single(binary_image);
%figure()
%imshow(A1)

% Simulate the radom points in the mask area

maskPixelIndex=find(binary_image==1);
Num_Pixel=length(maskPixelIndex);

% Assume Channel1: 35 random dots, Channel2: 7 random dots 
% (need to acquire the density of the foci, or directly acquire
% the number of foci in each mask)

for k=1:1:Num_iteration

% Acquire the random position in the Mask 
Ch1_dots_Num=Num_Ch1;
Ch2_dots_Num=Num_Ch2;
Ch1_Id=randperm(Num_Pixel,Ch1_dots_Num); % Positions are unique
Ch2_Id=randperm(Num_Pixel,Ch2_dots_Num);

% Convert to the Pixel Index 
Ch1_PixelId=maskPixelIndex(Ch1_Id);
Ch2_PixelId=maskPixelIndex(Ch2_Id);

% Convert to the Pixel coordinates on the image
[RowNum,ColNum]=size(A1);
Ch1_Xcoord=mod(Ch1_PixelId,RowNum); % get the reminder after division
Ch1_Xcoord(Ch1_Xcoord==0)=RowNum;   % If the reminder is 0, then the rowNum is the the Max row number
Ch1_Ycoord=fix(Ch1_PixelId./RowNum)+1;
% When the number is dividible
Reminder_1=mod(Ch1_PixelId,RowNum);
Id_Ch1_Re0=find(Reminder_1==0);
Ch1_Ycoord(Id_Ch1_Re0)=Ch1_PixelId(Id_Ch1_Re0)/RowNum;

Ch2_Xcoord=mod(Ch2_PixelId,RowNum); % get the reminder after division
Ch2_Xcoord(Ch2_Xcoord==0)=RowNum;   % If the reminder is 0, then the rowNum is the the Max row number
Ch2_Ycoord=fix(Ch2_PixelId./RowNum)+1;
% When the number is dividible
Reminder_2=mod(Ch2_PixelId,RowNum);
Id_Ch2_Re0=find(Reminder_2==0);
Ch2_Ycoord(Id_Ch2_Re0)=Ch2_Ycoord(Id_Ch2_Re0)/RowNum;

% % Plot the coordinates in the mask area 
% if i==17
% figure()
% [B,L] = bwboundaries(A1,'noholes');
% imshow(label2rgb(L, @jet, [.5 .5 .5]))
% hold on
% 
% for kk = 1:length(B)
%    boundary = B{kk};
%    plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
%    hold on 
%    plot(Ch1_Ycoord,Ch1_Xcoord,'r+','MarkerSize', 4,'LineWidth', 1);
%    hold on 
%    plot(Ch2_Ycoord,Ch2_Xcoord,'k+','MarkerSize', 4,'LineWidth', 1);
% end
% saveas(gcf,'Random simulation output_Complete dataset/Random simulation images/simulation_'+string(k)+'_'+MaskName);
% close
% end

% Precent of colocalization analysis

Ch1_Coord=[Ch1_Xcoord Ch1_Ycoord];
Ch2_Coord=[Ch2_Xcoord Ch2_Ycoord];
MinDist_Colo=2;    % Distance for colocalization

Num_Colo=0;
for j=1:1:Num_Ch1
    if Num_Ch1~=0
       Coordinates_Ch1=Ch1_Coord(j);  % get one X,Y coord from Ch1 
       Diff_Ch1_Ch2=abs(Ch2_Coord-Coordinates_Ch1);
       Delta_X=Diff_Ch1_Ch2(:,1);
       Delta_Y=Diff_Ch1_Ch2(:,2);
       Distance=sqrt(Delta_X.^2+Delta_Y.^2); % sqrt((Δx)^2+abs(Δy)^2)<=2
       YorN_withinMin_dist=Distance<=MinDist_Colo; % Return the logical array
       YorN_withinMin_dist_temp=double(YorN_withinMin_dist); % Convert the logical array to numbers
       NumofCh2_ColoWith_Ch1=sum(YorN_withinMin_dist_temp); % 0:no spot from Ch2 colocalize with one of the Coord from Ch1; 1/2/3: 1/2/3 spot(s) from Ch2 colocalize with one of the Coord from Ch1
       
       if NumofCh2_ColoWith_Ch1>=1
          Num_Ch1_Colocalized=1;
       else
          Num_Ch1_Colocalized=0; 
       end

        else
          Num_Ch1_Colocalized=0;
        end
          Num_Colo=Num_Colo+Num_Ch1_Colocalized;
end 

Num_Colo_multiple_simulation(k)=Num_Colo;

end

end
