function plot_spectra_neighbors(image)
[R,C,W] = size(image);

if W > 1
    image_2D = sum(image,3);
end

imagesc(image_2D);
[x, y] = ginput(1);

x = round(x); y = round(y);

figure
subplot 331
plot(squeeze(image(x,y,:)))
title(sprintf('Subplot 1: X = %d, Y = %d', max(x-1,1),max(y-1,1)))
grid on
ax = gca;
ax.GridColor = [0 .5 .5];
ax.GridLineStyle = '--';
ax.GridAlpha = 0.5;
ax.Layer = 'top';

subplot 332
plot(squeeze(image(max(x-1,1),y,:)))
title(sprintf('Subplot 2: X = %d, Y = %d', x,max(y-1,1)))
grid on
ax = gca;
ax.GridColor = [0 .5 .5];
ax.GridLineStyle = '--';
ax.GridAlpha = 0.5;
ax.Layer = 'top';

subplot 333
plot(squeeze(image(max(x-1,1),max(y-1,1),:)))
title(sprintf('Subplot 3: X = %d, Y = %d', min(x+1,C),max(y-1,1)))
grid on
ax = gca;
ax.GridColor = [0 .5 .5];
ax.GridLineStyle = '--';
ax.GridAlpha = 0.5;
ax.Layer = 'top';

subplot 334
plot(squeeze(image(max(x-1,1),min(y+1,R),:)))
title(sprintf('Subplot 4: X = %d, Y = %d', max(x-1,1), y))
grid on
ax = gca;
ax.GridColor = [0 .5 .5];
ax.GridLineStyle = '--';
ax.GridAlpha = 0.5;
ax.Layer = 'top';

subplot 335
plot(squeeze(image(min(x+1,C),y,:)))
title(sprintf('Subplot 5: X = %d, Y = %d', x,y))
grid on
ax = gca;
ax.GridColor = [0 .5 .5];
ax.GridLineStyle = '--';
ax.GridAlpha = 0.5;
ax.Layer = 'top';

subplot 336
plot(squeeze(image(min(x+1,C),max(y-1,1),:)))
title(sprintf('Subplot 6: X = %d, Y = %d', min(x+1,C),y))
grid on
ax = gca;
ax.GridColor = [0 .5 .5];
ax.GridLineStyle = '--';
ax.GridAlpha = 0.5;
ax.Layer = 'top';

subplot 337
plot(squeeze(image(min(x+1,C),min(y+1,R),:)))
title(sprintf('Subplot 7: X = %d, Y = %d', max(x-1,1),min(y+1,R)))
grid on
ax = gca;
ax.GridColor = [0 .5 .5];
ax.GridLineStyle = '--';
ax.GridAlpha = 0.5;
ax.Layer = 'top';

subplot 338
plot(squeeze(image(x,max(y-1,1),:)))
title(sprintf('Subplot 8: X = %d, Y = %d', x,min(y+1,R)))
grid on
ax = gca;
ax.GridColor = [0 .5 .5];
ax.GridLineStyle = '--';
ax.GridAlpha = 0.5;
ax.Layer = 'top';

subplot 339
plot(squeeze(image(x,min(y+1,R),:)))
title(sprintf('Subplot 9: X = %d, Y = %d', min(x+1,C),min(y+1,R)))
grid on
ax = gca;
ax.GridColor = [0 .5 .5];
ax.GridLineStyle = '--';
ax.GridAlpha = 0.5;
ax.Layer = 'top';