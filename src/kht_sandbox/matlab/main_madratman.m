% LINES = KHT(IMAGE,...) function usage. The complete description of the
% implemented techinique can be found at the paper below. If you use this
% implementation, please reference the paper.
%
%    Leandro A. F. Fernandes, Manuel M. Oliveira
%    Real-time line detection through an improved Hough transform voting scheme
%    Pattern Recognition (PR), Elsevier, 41:1, 2008, 299-314.
%    DOI: http://dx.doi.org/10.1016/j.patcog.2007.04.003
%    Project Page: http://www.inf.ufrgs.br/~laffernandes/kht.html

% Copyright (C) 2008 Leandro A. F. Fernandes and Manuel M. Oliveira

close all
clear
clc

basedir = '/home/ratneshmadaan/data/turkey_512/pcnn_output'
q = dir(basedir)
filename = {q.name}'

for i=3:length(filename)
    i
    I = imread(strcat(basedir, '/',filename{i}));
    I = rgb2gray(I);
    medfilt2(I);
    se = strel('disk',3,8);
    size(I);
%     I = bwmorph(I,'skel',Inf);
%     I = bwareaopen(I,2);

    I_dil = imdilate(I,se);

    figure;
%     I_dil = medfilt2(I_dil);

    lines = kht(I_dil);
    no_of_lines_to_consider = 10;

    if length(lines)<no_of_lines_to_consider
        no_of_lines_to_consider = length(lines);
    end
    % do kmeans on theta       
    %     clusters = kmeans(lines(:,2),2);

    [height width c] = size(I);
    
    axes;
    set(gca,'XLim',[0 height-1]);
    set(gca,'YLim',[0 width-1]);
    hold on;
    imshow(I_dil);


    lines = lines(1:no_of_lines_to_consider, :);
    for j=1:no_of_lines_to_consider
        if 
        
    for j=1:no_of_lines_to_consider
        if sind(lines(j,2)) ~= 0
            x = [-width/2 width/2-1];
            y = (lines(j,1) - x*cosd(lines(j,2)))/sind(lines(j,2));
        else
            x = [lines(j,1) lines(j,1)];
            y = [-height/2 height/2-1];
        end
        patch(x+width/2,y+height/2,[1 1 0],'EdgeColor','r','LineWidth',2);
    end
    
    hold off;
    saveas(gcf, strcat('/home/ratneshmadaan/data/turkey_512/pcnn_kht_lines/', int2str(i)), 'png');
    close all;
end
