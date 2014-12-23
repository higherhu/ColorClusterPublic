function [ imgSmoothShow, newLabel, smoothChoosedColor, runTime ] = ColorClusterFunc( img, threadcolor, K )
%COLORCLUSTERFUNC Summary of this function goes here
%   Detailed explanation goes here
%   输入：
%   img：            imread得到的输入图像
%   threadcolor:     绣线颜色库
%   K：              聚类类别数
%   
%   输出：
%       平滑后利用最接近的绣线颜色替换的结果imgSmoothShow，
%       平滑后个像素的标记newLabel
%       平滑后各类像素所选绣线颜色smoothChoosedColor
%       各阶段运行时间
%
% License
% This software is made publicly for research use only. It may be modified 
% and redistributed under the terms of the GNU General Public License. 
%
% Author: hujiagao@gmail.com
% 
    tic;
    lamdaWeight = 8;   % 数据项与平滑项的权重因子，越大越平滑
    %% calculating the histogram
    disp('calculating the histogram of the image...');
    [imgRow, imgCol, ~] = size(img);
    h = fspecial('gaussian', [7 7], 0.5);  % 对原图进行平滑，认为领域像素颜色会影响该像素观感颜色
    imgSmooth = imfilter(img, h, 'replicate', 'conv');%img;%        % 是否对输入图像进行初始平滑，即是否考虑邻域像素对像素观感颜色的影响
    imgDataSmooth = reshape(imgSmooth, imgRow*imgCol, 3);
    % 将RGB数据转换为可以用来计算颜色距离的数据（Lab）
    imgDataCalDis = myColorDisMeasures(imgDataSmooth);%imgData;
    threadColorDis = myColorDisMeasures(threadcolor);%threadcolor;
    % 统计颜色直方图，计算每个像素与绣线颜色最近的那个
    [~, minIndex] = pdist2(threadColorDis, imgDataCalDis, 'euclidean', 'smallest', 1);
    histogram = tabulate(minIndex); % histogram中不一定包含所有的颜色，有些没出现在原图中的颜色可能不在这里面
    [~, sInx] = sort(histogram(:,3), 'descend');
    iniColorIndex = histogram(sInx(1:4*K), 1);  
    CenterColorCalDis = threadColorDis(iniColorIndex, :);
    Z = linkage(CenterColorCalDis, 'complete', 'euclidean');
    C = cluster(Z, 'maxclust', K);  % 换成'cutoff',sigma，sigma为颜色是否相似的阈值，则可自动确定聚类数
    color = zeros(K, size(threadColorDis, 2));  % 计算聚类后每类的中心色。可用该类中数量最多的颜色代替，
    for i=1:K                                   % 或按各颜色像素比例计算中心色                
        color(i,:) = mean(CenterColorCalDis(C==i, :), 1); 
    end
    CenterCalDis = color;  
    t1=toc;
    
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
    disp(['  cluster: total iterators:' num2str(iterator)]);
    t2=toc;

    %% smooth
    disp('smooth ...');
    h = GCO_Create(imgRow*imgCol, K);
    dataDis = pdist2(CenterCalDis, imgDataCalDis);
    dataDis = 1./(dataDis+0.0000001);   % 距离越远，概率越小
    dataSum = sum(dataDis, 1);
    dataProb = dataDis ./ dataSum(ones(K,1),:);    % 这里的dataProb为每个像素属于各类的概率
    GCO_SetDataCost(h, int32(-log(dataProb)));
    smoothCost = ones(K, K);
    smoothCost(1:(K+1):K*K) = 0;
    GCO_SetSmoothCost(h, int32(smoothCost));
    neighbors = myGetNeighbors(imgDataCalDis, imgRow, imgCol);
    GCO_SetNeighbors(h, lamdaWeight*neighbors);
    t3 = toc;
    
    GCO_Expansion(h);
    newLabel = GCO_GetLabeling(h);
    newColorTable = tabulate(newLabel);
    t4 = toc;

    imgSmoothShow = zeros(imgRow*imgCol, 3);
    smoothChoosedColor = zeros(K, 3);   % 所选绣线颜色
    for i=1:size(newColorTable, 1)
        index = find(newLabel==newColorTable(i,1));
        color = calRegionColor(threadcolor, imgDataCalDis(index, :));   % 利用迭代计算区域颜色
        smoothChoosedColor(i, :) = color;
        imgSmoothShow(index, :) = color(ones(size(index,1), 1), :); 
    end
    imgSmoothShow = reshape(uint8(imgSmoothShow), size(img));
    GCO_Delete(h);
    t5=toc;
    
    % 运行时间，依次为：颜色聚类耗时、设定datacost和smoothcost耗时、设置neighbors耗时
    % 运行图割耗时、总耗时
    runTime = [t1, t2-t1, t3-t2, t4-t3, t5];
    
    disp('Finish!');

end

