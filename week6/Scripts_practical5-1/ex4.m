picture = imread('mop.png ');
picture = imresize(picture, [224, 224]);
[labels, scores ] = classify(GoogLeNet, picture);
figure 
imshow( picture)
title(string(label)+ "," + num2str(100*scores(classNames == labels), 3) + "%");
[~, idx ]= sort(scores, 'descend');
idx= idx(3:-1:1);
figure
barh(scores)
xlim([0 1])
xlabel("probabilities")
title("top 3 predictopns")
yticketlabels(classNamesTop)