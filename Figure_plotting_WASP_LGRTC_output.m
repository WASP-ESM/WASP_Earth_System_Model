% This is example MATLAB code for plotting output from the WASP/LGRTC model presented in A computationally efficient model for probabilistic local warming projections constrained by history matching and pattern scaling

% Manuscript for Geosceintific Model Development Discussions


% Philip Goodwin1,*, Martin Leduc2, Antti-Ilari Partanen3, H. Damon Matthews4 and Alex Rogers5


% 1 Ocean and Earth Science, University of Southampton, National Oceanography Centre Southampton, Southampton, SO14 3ZH, UK
% 2 Ouranos, Montreal, Canada
% 3 Climate System Research, Finnish Meteorological Institute, Helsinki, Finland
% 4 Department of Geography, Planning and Environment, Concordia University, Montreal, Canada
% 5 Department of Computer Science, University of Oxford, Oxford, UK

% * Corresponding author, email: p.a.goodwin@soton.ac.uk

%This code plots a figure in MATLAB from the WASP/LGRTC output
%The figure shows three panels: (top) a median spatial warming projection from the WASP/LGRTC ensemble, along with the (middle) 83rd percentile and (bottom) the 17th percentile


lat=[-87.8638
     -85.096527
     -82.312912
     -79.525612
     -76.7369
     -73.947517
     -71.157753
     -68.36776
     -65.577606
     -62.787354
     -59.997021
     -57.206635
     -54.416203
     -51.625736
     -48.835243
     -46.044727
     -43.254196
     -40.46365
     -37.673092
     -34.882523
     -32.091946
     -29.301363
     -26.510773
     -23.720177
     -20.929577
     -18.138973
     -15.348368
     -12.557758
     -9.767148
     -6.976536
     -4.185923
     -1.395309
     1.395309
     4.185923
     6.976536
     9.767148
     12.557758
     15.348368
     18.138973
     20.929577
     23.720177
     26.510773
     29.301363
     32.091946
     34.882523
     37.673092
     40.46365
     43.254196
     46.044727
     48.835243
     51.625736
     54.416203
     57.206635
     59.997021
     62.787354
     65.577606
     68.36776
     71.157753
     73.947517
     76.7369
     79.525612
     82.312912
     85.096527
     87.8638];

lon=[-180
     -177.1875
     -174.375
     -171.5625
     -168.75
     -165.9375
     -163.125
     -160.3125
     -157.5
     -154.6875
     -151.875
     -149.0625
     -146.25
     -143.4375
     -140.625
     -137.8125
     -135
     -132.1875
     -129.375
     -126.5625
     -123.75
     -120.9375
     -118.125
     -115.3125
     -112.5
     -109.6875
     -106.875
     -104.0625
     -101.25
     -98.4375
     -95.625
     -92.8125
     -90
     -87.1875
     -84.375
     -81.5625
     -78.75
     -75.9375
     -73.125
     -70.3125
     -67.5
     -64.6875
     -61.875
     -59.0625
     -56.25
     -53.4375
     -50.625
     -47.8125
     -45
     -42.1875
     -39.375
     -36.5625
     -33.75
     -30.9375
     -28.125
     -25.3125
     -22.5
     -19.6875
     -16.875
     -14.0625
     -11.25
     -8.4375
     -5.625
     -2.8125
     0
     2.8125
     5.625
     8.4375
     11.25
     14.0625
     16.875
     19.6875
     22.5
     25.3125
     28.125
     30.9375
     33.75
     36.5625
     39.375
     42.1875
     45
     47.8125
     50.625
     53.4375
     56.25
     59.0625
     61.875
     64.6875
     67.5
     70.3125
     73.125
     75.9375
     78.75
     81.5625
     84.375
     87.1875
     90
     92.8125
     95.625
     98.4375
     101.25
     104.0625
     106.875
     109.6875
     112.5
     115.3125
     118.125
     120.9375
     123.75
     126.5625
     129.375
     132.1875
     135
     137.8125
     140.625
     143.4375
     146.25
     149.0625
     151.875
     154.6875
     157.5
     160.3125
     163.125
     165.9375
     168.75
     171.5625
     174.375
     177.1875];



for n = 1:64
    for m = 1:128

        B = T_2081_2100(n,m,:);
        C = mean(B);
        C2 = prctile(B,82);
        C3 = prctile(B,18);


        D(n,m)= C;
        D2(n,m)= C2;
        D3(n,m)= C3;
    end
end




colormap(flipud(hot));

subplot(3,1,1)
worldmap world

load coastlines
whos

caxis([0,15])
colorbar

surfacem(lat,lon,D);

plotm(coastlat,coastlon,'k');

subplot(3,1,2)
worldmap world

load coastlines
whos

caxis([0,15])
colorbar

surfacem(lat,lon,D2);

plotm(coastlat,coastlon,'k');

subplot(3,1,3)
worldmap world

load coastlines
whos

caxis([0,15])
colorbar

surfacem(lat,lon,D3);

plotm(coastlat,coastlon,'k');









