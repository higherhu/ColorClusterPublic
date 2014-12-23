function [ neighbors ] = myGetNeighbors( imgData, rowNum, colNum )
%MYGETNEIGHBORS Summary of this function goes here
%   Detailed explanation goes here
%   imgData:    numPixel*3 
%
    imgSize = size(imgData, 1);
    
    index = -1*ones(24,1);
    nowPos= 1;
    varianceSquared = getvs(imgData, rowNum, colNum );
    pixelDis = [1 1.4142 1 1.4142 1 1.4142 1 1.4142 2 2.2361 2.8284 2.2361 2 ...
        2.2361 2.8284 2.2361 2 2.2361 2.8284 2.2361 2 2.2361 2.8284 2.2361];

    for i=1:imgSize
        row = mod(i-1, rowNum)+1;
        col = ceil(i/rowNum);
        index(1) = (col-2)*rowNum + row; % LL
        index(2) = index(1) + 1;         % LD
        index(3) = i + 1;                % DD
        index(4) = (col)*rowNum + row+1; % DR
        index(5) = index(4) - 1;         % RR
        index(6) = index(5) - 1;         % RT
        index(7) = i - 1;                % TT
        index(8) = index(1) - 1;         % TL
        
        index(9) = (col-3)*rowNum + row; 
        index(10) = index(9) + 1;
        index(11) = index(10) + 1;
        index(12) = index(2) + 1;
        index(13) = index(3) + 1;
        index(14) = index(4) + 1;
        index(15) = (col+1)*rowNum + row+2;
        index(16) = index(15) - 1;
        index(17) = index(16) - 1;
        index(18) = index(17) - 1;
        index(19) = index(18) - 1;
        index(20) = index(6) - 1;
        index(21) = index(7) - 1;
        index(22) = index(8) - 1;
        index(23) = index(9) - 2;
        index(24) = index(23) + 1;
    
        if row == rowNum    % 最下一行下方不存在相邻像素
            index([2 3 4 10 16 11 12 13 14 15]) = -1;
        elseif row == 1    % 最上一行上方不存在相邻像素
            index([6 7 8 18 19 20 21 22 23 24]) = -1;
        elseif row==2
            index([19 20 21 22 23]) = -1;
        elseif row==rowNum-1
            index([11 12 13 14 15]) = -1;
        end
        if col ==  colNum  % 最右一列右方不存在相邻像素
            index([4 5 6 14 15 16 17 18 19 20]) = -1;
        elseif col == 1    % 最左一列左方不存在相邻像素
            index([1 2 8 9 10 11 12 22 23 24]) = -1;
        elseif col==2
            index([9 10 11 23 24]) = -1;
        elseif col==colNum-1
            index([15 16 17 18 19]) = -1;
        end
        
        indexExist = index(index > 0);
        neighDis = pixelDis(index > 0)';
        neighDis = neighDis(indexExist>i);
        indexExist = indexExist(indexExist>i);  % 只计算上三角阵
        existNum = size(indexExist,1);
        
       %% 两组像素编号，以及对应两组像素之间的距离，为了加速稀疏矩阵的建立过程
        pixel1(nowPos:nowPos+existNum-1) = i*ones(existNum,1);      % 
        pixel2(nowPos:nowPos+existNum-1) = indexExist;                        % 
        X = imgData(i, :);
        Y = imgData(indexExist, :);    
    % exp( -(C1-C2)^2/(2*varianceSquared) )/dis
        dis(nowPos:nowPos+existNum-1) = ...
            exp( -sum((X(ones(1,existNum), :)-Y).^2, 2)/(2*varianceSquared))./neighDis;% 
        nowPos = nowPos + existNum;        
    end
    neighbors = sparse(pixel1, pixel2, dis, imgSize, imgSize);

end

function varianceSquared = getvs(imgData, rowNum, colNum )
    counter = 0;
    varianceSquared = 0;
    index = -1*ones(24,1);
    for i=1:size(imgData,1)
        row = mod(i-1, rowNum)+1;
        col = ceil(i/rowNum);
        index(1) = (col-2)*rowNum + row; % LL
        index(2) = index(1) + 1;         % LD
        index(3) = i + 1;                % DD
        index(4) = (col)*rowNum + row+1; % DR
        index(5) = index(4) - 1;         % RR
        index(6) = index(5) - 1;         % RT
        index(7) = i - 1;                % TT
        index(8) = index(1) - 1;         % TL
        
        index(9) = (col-3)*rowNum + row; 
        index(10) = index(9) + 1;
        index(11) = index(10) + 1;
        index(12) = index(2) + 1;
        index(13) = index(3) + 1;
        index(14) = index(4) + 1;
        index(15) = (col+1)*rowNum + row+2;
        index(16) = index(15) - 1;
        index(17) = index(16) - 1;
        index(18) = index(17) - 1;
        index(19) = index(18) - 1;
        index(20) = index(6) - 1;
        index(21) = index(7) - 1;
        index(22) = index(8) - 1;
        index(23) = index(9) - 2;
        index(24) = index(23) + 1;
     
        if row == rowNum    % 最下一行下方不存在相邻像素
            index([2 3 4 10 16 11 12 13 14 15]) = -1;
        elseif row == 1    % 最上一行上方不存在相邻像素
            index([6 7 8 18 19 20 21 22 23 24]) = -1;
        elseif row==2
            index([19 20 21 22 23]) = -1;
        elseif row==rowNum-1
            index([11 12 13 14 15]) = -1;
        end
        if col ==  colNum  % 最右一列右方不存在相邻像素
            index([4 5 6 14 15 16 17 18 19 20]) = -1;
        elseif col == 1    % 最左一列左方不存在相邻像素
            index([1 2 8 9 10 11 12 22 23 24]) = -1;
        elseif col==2
            index([9 10 11 23 24]) = -1;
        elseif col==colNum-1
            index([15 16 17 18 19]) = -1;
        end
        
        indexExist = index(index >= 1);
        indexExist = indexExist(indexExist>i);  % 只计算上三角阵
        existNum = size(indexExist,1);
        
        X = imgData(i, :);
        Y = imgData(indexExist, :);
        varianceSquared = varianceSquared + sum(sum((X(ones(1,existNum), :)-Y).^2, 2));
        counter = counter+size(indexExist, 1);
    end
    varianceSquared = varianceSquared/counter;
end

