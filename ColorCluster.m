% 统计图像在绣线颜色库中的直方图，取直方图中最多的几个为初始聚类中心
% 对原图进行聚类，聚类的每次迭代时，新的聚类中心在绣线颜色库中选择
% License
% This software is made publicly for research use only. It may be modified 
% and redistributed under the terms of the GNU General Public License. 
%
% Author: hujiagao@gmail.com
% 

addpath(genpath('./GCO'));
%matlabpool('open','local',3); 
clear;
tic;

img_dir = './data/';
img_name = 'hua2.jpg';
K = 3;
lamdaWeight = 8;   % 数据项与平滑项的权重因子，越大越平滑

%% Read the color list
filename = 'threadcolor.txt';
delimiter = '\t';
startRow = 2;
fileID = fopen(filename,'r');
formatSpec = '%f%f%f%[^\n\r]';
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);
threadcolor = [dataArray{1:end-1}];
threadcolor = unique(threadcolor, 'rows');  % 删除重复的颜色
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

%% read the image and calculating the histogram
disp('calculating the histogram of the image...');
img = imread([img_dir img_name]);
figure('Name','The Source Image');
imshow(img);
[imgRow, imgCol, ~] = size(img);
h = fspecial('gaussian', [7 7], 0.5);  % 对原图进行平滑，认为邻域像素颜色会影响该像素观感颜色
imgSmooth = imfilter(img, h, 'replicate', 'conv');
imgDataSmooth = reshape(imgSmooth, imgRow*imgCol, 3);
% 将RGB数据转换为可以用来计算颜色距离的数据（Lab）
imgDataCalDis = myColorDisMeasures(imgDataSmooth);%imgData;
threadColorDis = myColorDisMeasures(threadcolor);%threadcolor;
% 统计颜色直方图，计算每个像素与绣线颜色最近的那个
[~, minIndex] = pdist2(threadColorDis, imgDataCalDis, 'euclidean', 'smallest', 1);
figure('Name','The quantized Image');
imshow(uint8(reshape(threadcolor(minIndex,:), imgRow, imgCol, 3)));
histogram = tabulate(minIndex); % histogram中不一定包含所有的颜色，有些没出现在原图中的颜色可能不在这里面
[sorted, sInx] = sort(histogram(:,3), 'descend'); 
iniColorIndex = histogram(sInx(1:4*K), 1);   % 选取前4K种颜色
CenterColorCalDis = threadColorDis(iniColorIndex, :);
Z = linkage(CenterColorCalDis, 'complete', 'euclidean');
C = cluster(Z, 'maxclust', K);  % 换成'cutoff',sigma，sigma为颜色是否相似的阈值，则可自动确定聚类数
color = zeros(K, size(threadColorDis, 2));  % 计算聚类后每类的中心色。可用该类中数量最多的颜色代替，
for i=1:K                                   % 或按各颜色像素比例计算中心色                
	color(i,:) = mean(CenterColorCalDis(C==i, :), 1); 
end
CenterCalDis = color;   % 以直方图聚类中心色作为原图聚类的初始绣线中心色
% clearvars h imgSmooth imgDataSmooth minIndex histogram iniColorIndex CenterColor ...
%     Z C i color dis indexmin;
toc;

%% cluster the image
disp('clustering ...');
MAXITERATOR = 100;
iterator = 0;
% 计算每个像素所属类别
[~, labels] = pdist2(CenterCalDis, imgDataCalDis, 'euclidean', 'smallest', 1);
labels2 = labels;
while iterator < MAXITERATOR
    % 重新计算聚类中心，在绣线库中找聚类中心,保证类内方差最小（通过求均值，然后在绣线库中找与均值最接近的颜色实现，正确性待证明）
    color = zeros(K, size(threadColorDis, 2));  % 计算聚类后每类的中心色
    for i=1:K                                              
        color(i,:) = mean(imgDataCalDis(labels2==i, :), 1); 
    end
    [~, CenterColorIndex] = pdist2(threadColorDis, color, 'euclidean', 'smallest', 1); % 利用色板中与区域平均色最接近的颜色替换原色
    CenterCalDis(:,:)  = threadColorDis(CenterColorIndex,:);    % 原图的初始绣线中心颜色

    % 重新计算每个像素所属类别
    [~, labels2] = pdist2(CenterCalDis, imgDataCalDis, 'euclidean', 'smallest', 1);
    
    if isequal(labels2, labels)  % labels不变，聚类结束
        break;
    end

    labels = labels2;
    iterator = iterator+1;
end

disp(['  total iterators:' num2str(iterator)]);
CenterColor = threadcolor(CenterColorIndex, :);
imgCluster = CenterColor(labels, :);
imgCluster = int32(imgCluster);
binaryImgSmooth = bitshift(imgCluster(:,1),16)+bitshift(imgCluster(:,2),8)+imgCluster(:,3);
imgCluster = reshape(uint8(imgCluster), size(img));
figure('Name','The Cluster Result');
imshow(imgCluster);
% clearvars MAXITERATOR labels2 i indexmin color imgCluster; %iterator CenterColor
toc;

%% smooth
disp('smooth ...');
h = GCO_Create(imgRow*imgCol, K);
disp('setting data cost ...');
dataDis = pdist2(CenterCalDis, imgDataCalDis);
dataDis = 1./(dataDis+0.0000001);   % 距离越远，概率越小
dataSum = sum(dataDis, 1);
dataProb = dataDis ./ dataSum(ones(K,1),:);    % 这里的dataProb为每个像素属于各类的概率
GCO_SetDataCost(h, int32(-log(dataProb)));
toc;

smoothCost = int32(ones(K, K));
smoothCost(1:(K+1):K*K) = 0;
GCO_SetSmoothCost(h, int32(smoothCost));
disp('setting neighbors ...');
neighbors = myGetNeighbors(imgDataCalDis, imgRow, imgCol);
GCO_SetNeighbors(h, lamdaWeight*neighbors);
toc;

GCO_SetLabeling(h, labels);

disp('expansion ...');
GCO_SetVerbosity(h, 1);
Energy = GCO_Expansion(h);
newLabel = GCO_GetLabeling(h);
GCO_Delete(h);
newColorTable = tabulate(newLabel);
toc;

disp('display the resuslt...');
imgSmoothShow = zeros(imgRow*imgCol, 3);
smoothChoosedColor = zeros(K, 3);   % 所选绣线颜色
for i=1:size(newColorTable, 1)
    index = find(newLabel==newColorTable(i,1));
    color = calRegionColor(threadcolor, imgDataCalDis(index, :));   % 利用迭代计算区域颜色
    smoothChoosedColor(i, :) = color;
    imgSmoothShow(index, :) = color(ones(size(index,1), 1), :); 
end
imgSmoothShow = int32(imgSmoothShow);
imgSmoothShow = reshape(uint8(imgSmoothShow), size(img));
figure('Name','The Smooth Result In The ColorList');
imshow(imgSmoothShow);
imwrite(imgSmoothShow, [img_dir img_name '_result.bmp'])
dlmwrite([img_dir img_name '.txt'], [imgRow imgCol], '-append', 'delimiter', ' ')
dlmwrite([img_dir img_name '.txt'], reshape(newLabel, imgRow, imgCol) , '-append', 'delimiter', ' ')
dlmwrite([img_dir img_name '.txt'], smoothChoosedColor , '-append', 'delimiter', ' ')

toc;