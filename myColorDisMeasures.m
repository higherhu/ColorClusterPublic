function [ output_args ] = myColorDisMeasures( input_args )
%MYCOLORDISMEASURES Summary of this function goes here
%   Detailed explanation goes here

    input_args = double(input_args);

    output_args = colorspace('rgb->lab', input_args/255.0);
    
end

