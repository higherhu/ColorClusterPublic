function resultColor = calRegionColor( threadcolor, regionColor )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% 计算区域平均色以选取绣线，可以考虑迭代算法：
% 1、	计算区域平均色
% 2、	计算区域内每个像素与平均色之间的差距
% 3、	根据上一步计算的差距，确定各像素对平均色的权重，距离越近权重越大
% 4、	将各像素颜色加权，计算平均色，若此轮平均色变化不大（所选绣线不变），停止；否则，返回步骤2.

    threadcolorlab = myColorDisMeasures(threadcolor);
    resultColor = mean(regionColor, 1);
    [~, indexmin] = pdist2(threadcolorlab, resultColor, 'euclidean', 'smallest', 1);
    iterator = 0;
    while iterator < 0    % 最多迭代100次  
        dis = pdist2(threadcolorlab, resultColor);

        weight = exp(-dis); % 距离越远，权重越小
        weight = weight / sum(weight);
        resultColor = sum(regionColor .* weight(:,ones(1, size(regionColor, 2))), 1);
        
        dis = pdist2(threadcolorlab, resultColor);
        [~,indexmin2] = min(dis);
        if indexmin == indexmin2    % 所选绣线颜色不变
            break;
        end
        indexmin = indexmin2;
        iterator = iterator + 1;
    end
    if iterator > 0
        disp(['  calRegionColor: iterator:' num2str(iterator)]);
    end
    resultColor = threadcolor(indexmin, :);     % 返回绣线颜色
end

