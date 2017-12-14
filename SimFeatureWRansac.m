clear;clc;close all
imagefiles = dir('*.jpg');
nfiles = length(imagefiles);
M = cell(nfiles,1);

for i =1:nfiles
    currentfilename = imagefiles(i).name;
    currentimage = imread(currentfilename);
    images{i} = currentimage;
    M2{i} = rgb2gray(images{i});
    M{i} = double(M2{i});   
end

A = size(M{1});  
B = A(1); C = A(2);  %stores Pixel Sizes
EL = 3;             %Prewitt mask size
PrewX(1:EL,1) = -1;
PrewX(1:EL,2) = 0;
PrewX(1:EL,3) = 1;
PrewY = PrewX';
SS = floor(EL/2);
SST = ceil(EL/2);
%%
tic
for k = 1:nfiles
    Ix{k} = zeros(B,C);
    Iy{k} = zeros(B,C);
    Mag{k} = zeros(B,C);                % Edge Detection Storage
    Mag2{k} = zeros(B,C);
end
for k = 1:nfiles                         % conv2 is dead to me
    for i = SST:B-SST
        for j = SST:C-SST
            Ix{k}(i,j) = sum(sum(PrewX.*(M{k}(i-SS:i+SS,j-SS:j+SS))));
            Iy{k}(i,j) = sum(sum(PrewY.*(M{k}(i-SS:i+SS,j-SS:j+SS))));
        end
    end
end

%%
Thresh = 100;
for k = 1:nfiles
    for i = 1:B
        for j = 1:C
            Mag{k}(i,j) = sqrt(Ix{k}(i,j)^2+Iy{k}(i,j)^2); %X,Y Magnitude
            if Mag{k}(i,j) > Thresh
                Mag2{k}(i,j) = 1;       %0,1 Map Magnitudes
            end           
        end
    end
end

Gauss = fspecial('gaussian',7,1.4);
for k = 1:nfiles
    Ix2{k} = conv2(Ix{k}.^2,Gauss,'same');   
    Iy2{k} = conv2(Iy{k}.^2,Gauss,'same');
    IxIy{k} = conv2(Ix{k}.*Iy{k},Gauss,'same');
    R{k} = zeros(B,C);
    R2{k} = zeros(B,C);
    
end
%% R Values
Thresh2 = 1e8; % Corner Threshold, Higher = less corners
count = 1;
count2 = 1;
Matrix2 = zeros(3);
Mat2 = zeros(2);
for k = 1:nfiles
    for i = 2:B-1
        
        for j = 2:C-1
            for l = -1:1
                for m = -1:1
                    Matrix1{count} = ([Ix2{k}(i+l,j+m) IxIy{k}(i+l,j+m); IxIy{k}(i+l,j+m) Iy2{k}(i+l,j+m)]);
                    count = count + 1; %Stores values centered at i,j and 8 neighbors
                    
                end
            end
            for ii = 1:9            
                Res = (Matrix1{ii} - Matrix1{5}).^2;  % Gradient difference from center and neighboards
                Mat2 = Mat2 + Res;              
            end                         
             Matrix3 = sqrt(Mat2);
             R{k}(i,j) = det(Matrix3) - .05*trace(Matrix3)^2;  
            if R{k}(i,j) > Thresh2
                R2{k}(i,j) = 1;
               
            end
            count = 1;
            count2 = 1;
            Mat2 = zeros(2);
        end
    end
end
                
%% Suppress Corners
for k = 1:nfiles
    count = 1;
    for i = 2:B-1
        for j = 2:C-1
            if R2{k}(i,j) == 1
             Window = R{k}(i-1:i+1,j-1:j+1);
             Max = Window(2,2);
             Window(2,2) = 0;
             
             if Max<max(max(Window))
                 R2{k}(i,j) = 0;
             else
                 R3{k}(count,1) = i;
                 R3{k}(count,2) = j;
                 count = count + 1;
             end
            end
        end
    end 
end
%% NCC
Xset = R3{1}(:,1);
Yset = R3{1}(:,2);
count3 = 0;
TT = 160; % Sets Pixel search limit for second image. % keep under 175
%Corner at i,j image1, image2 i+/- TT, j +/- TT
% 
NCCW = 9; %NCC Window, %%%% Significant impact on accuracy %%%% 5 min, 3 = bad results
NCCMin = floor(NCCW/2);
for ii = 1:length(Xset)
    if Xset(ii) < ceil(NCCW/2)
        Xset(ii) = ceil(NCCW/2);
    end
    if Xset(ii) > B - floor(NCCW/2)
        Xset(ii) = B-floor(NCCW/2);
    end
    if Yset(ii) < ceil(NCCW/2)
        Yset(ii) = ceil(NCCW/2);
    end
    if Yset(ii) > C - floor(NCCW/2)
        Yset(ii) = C - floor(NCCW/2);
    end
      
    count3 = count3 + 1;
    Result = 0;
    FFF = M{1}(Xset(ii)-NCCMin:Xset(ii)+NCCMin,Yset(ii)-NCCMin:Yset(ii)+NCCMin);
    FMat{ii} = FFF./sqrt(sum(sum(FFF.^2)));
    
   
    XX = Xset(ii);
    YY = Yset(ii);
    %Ensures the search window won't go outside actual image 
    
    if XX <= TT+NCCMin
        XXMin = ceil(NCCW/2);
        XXMax = XX + TT;
    elseif XX >= B - TT-NCCMin
        XXMax = B-NCCMin;
        XXMin = XX;
    else
        XXMax = XX + TT;
        XXMin = XX - TT;
    end
    if YY <= TT+NCCMin
        YMin = ceil(NCCW/2);
        YMax = YY;
    elseif YY >=C - TT-NCCMin
        YMax = C-NCCMin;
        YMin = YY-TT;
    else
        YMin = YY-TT;
        YMax = YY+TT;
    end
    % Only calculates between image1 and image2
    % i j search windows for image2
    for i = XXMin:XXMax
        for j = YMin:YMax
            GGG = M{2}(i-NCCMin:i+NCCMin,j-NCCMin:j+NCCMin);
            GGG = GGG./sqrt(sum(sum(GGG.^2)));
           
            NCC = sum(sum(FMat{ii}.*GGG));
           
            if NCC > Result
         
                Result = NCC;
                Image2{1}(count3,1) = i;
                Image2{1}(count3,2) = j;
                %i,j values of image2 wrt corner in image1
            end
        end
    end
end
%% Homography Matrix

Trials = 2500;  %Significant impact on run time
RANLen = length(R3{1}(:,1));

Set = 1:RANLen-rem(RANLen,4);
count = 1;
for i = 1:Trials  
   I = randperm(RANLen,4);
   count = 1;
   
       
   AA = [R3{1}(I(1),1) R3{1}(I(1),2) 1 0 0 0 -R3{1}(I(1),1)...
   *Image2{1}(I(1),1) -R3{1}(I(1),2)*Image2{1}(I(1),1);...
   0 0 0 R3{1}(I(1),1) R3{1}(I(1),2) 1 -R3{1}(I(1),1)*Image2{1}(I(1),2) ...
   -R3{1}(I(1),2)*Image2{1}(I(1),2);...
   R3{1}(I(2),1) R3{1}(I(2),2) 1 0 0 0 -R3{1}(I(2),1)*Image2{1}(I(2),1) ...
   -R3{1}(I(2),2)*Image2{1}(I(2),1);...
   0 0 0 R3{1}(I(2),1) R3{1}(I(2),2) 1 -R3{1}(I(2),1)*Image2{1}(I(2),2) ...
   -R3{1}(I(2),2)*Image2{1}(I(2),2);...
   R3{1}(I(3),1) R3{1}(I(3),2) 1 0 0 0 -R3{1}(I(3),1)*Image2{1}(I(3),1) ...
   -R3{1}(I(3),2)*Image2{1}(I(3),1);...
   0 0 0 R3{1}(I(3),1) R3{1}(I(3),2) 1 -R3{1}(I(3),1)*Image2{1}(I(3),2) ...
   -R3{1}(I(3),2)*Image2{1}(I(3),2);...
   R3{1}(I(4),1) R3{1}(I(4),2) 1 0 0 0 -R3{1}(I(4),1)*Image2{1}(I(4),1) ...
   -R3{1}(I(4),2)*Image2{1}(I(4),1);...
   0 0 0 R3{1}(I(4),1) R3{1}(I(4),2) 1 -R3{1}(I(4),1)*Image2{1}(I(4),2) ...
   -R3{1}(I(4),2)*Image2{1}(I(4),2)];

 BB = [Image2{1}(I(1),1) Image2{1}(I(1),2) Image2{1}(I(2),1) Image2{1}(I(2),2) ...
     Image2{1}(I(3),1) Image2{1}(I(3),2) Image2{1}(I(4),1) Image2{1}(I(4),2)]';
 Res2 = AA\BB;
 XX2{i} = Res2;
 
end
   % More succint but theres a bug
% for i = 1:Trials  
%    I = randperm(RANLen,4);
%    count = 1;
%    for jj = 1:2:7
%        AA(jj,1:8) = [R3{1}(I(count),1) R3{1}(I(count),2) 1 0 0 0 -R3{1}(I(count),1)...
%            *Image2{1}(I(count),1) -R3{1}(I(count),2)*Image2{1}(I(count),1)];
%        AA(jj+1,1:8) = [0 0 0 R3{1}(I(count),1) R3{1}(I(count),2) 1 -R3{1}(I(count),1)...
%            *Image2{1}(I(count),2) -R3{1}(I(count),2)*Image2{1}(I(count),2)];
% 
%        BB(jj,1) = Image2{1}(I(count),1);
%        BB(jj+1,1) = Image2{1}(I(count),2);
%    count = count + 1;
%    Res2 = AA\BB;
%    XX2{i} = Res2;
%    end
% end
                                                                                                     
%% Best Homography Matrix
for i = 1:Trials
    RS1{i} = zeros(1,RANLen*2);
end
%%
%Trials*RANLen operations

for i = 1:Trials
    XMat = [XX2{i}(1:3)';XX2{i}(4:6)';XX2{i}(7:8)',1];
    store1 = 1;
    store2 = 1;
    for j = 1:RANLen
    L1 = [R3{1}(j,1);R3{1}(j,2);1];
    ANS = XMat*L1;
    RS1{i}(store1,1) = ANS(1);
    RS1{i}(store1,2) = ANS(2);   
    store1 = store1 + 1;
    
    end
end

%%
for i = 1:Trials
    for j = 1:store1-1
         DifferencesX{i}(j,1) = (RS1{i}(j,1)-Image2{1}(j,1));
    end
end
for i = 1:Trials
    for j = 1:store1-1
         DifferencesY{i}(j,1) = (RS1{i}(j,2)-Image2{1}(j,2));
   
    end
end

%%
Winner = 0;
for i = 1:Trials
    MINS = 0;
    for j = 1:store1-1
        if abs(DifferencesX{i}(j,1))<3  && abs(DifferencesY{i}(j,1))<3
            MINS = MINS + 1;
        end
           
        if MINS>Result
            Result = MINS;
            WinnerX = i;
        end
    end
end

%% Winning X!
XWin = [XX2{WinnerX}(1:3)';XX2{WinnerX}(4:6)';XX2{WinnerX}(7:8)',1];
for ii = 1:RANLen
    XPoint1(ii) = R3{1}(ii,1);
    YPoint1(ii) = R3{1}(ii,2);
    MAT = [XPoint1(ii);YPoint1(ii);1];
    AXB = XWin*MAT;
    XPoint2(ii) = AXB(1);
    YPoint2(ii) = AXB(2);
end
toc
%%

% figure(1)
% subplot(1,2,1); imshow(images{1}); hold on; plot(R3{1}(:,2),R3{1}(:,1),'*k')
% subplot(1,2,2); imshow(images{2}); hold on; plot(R3{2}(:,2),R3{2}(:,1),'*r')

%% Close enough plots
% Mostly good, can be filtered further
% figure()
% imshowpair(images{1},images{2},'montage');
% for i = 1:RANLen
%     hold on;plot(YPoint1(i),XPoint1(i),'*r'); hold on; 
%     plot(round(YPoint2(i))+C,round(XPoint2(i)),'*r')
%     hold on;
%     plot([YPoint1(i),(YPoint2(i))+C],[XPoint1(i),(XPoint2(i))],'LineWidth',3)
% end
%% Extreme Criteria Plot
% Very few results, only exact correspondences will show
% If poor results, lower Thresh2, not enough corners for proper homography
% matrix
%Second option increase the size of NCCW
% 3rd lower TH3
%More correspondences, lower Thresh2, longer run times
TH3 = .25;
figure(2)
imshowpair(images{1},images{2},'montage');
for i = 1:RANLen
    if abs(YPoint2(i) - Image2{1}(i,2))<TH3 && abs(XPoint2(i)-Image2{1}(i,1))<TH3
    hold on;plot(YPoint1(i),XPoint1(i),'*r'); hold on;
    plot(round(YPoint2(i))+C,round(XPoint2(i)),'*r')
    hold on;
    plot([YPoint1(i),round(YPoint2(i))+C],[XPoint1(i),round(XPoint2(i))],'LineWidth',3)
    end
end

 

