% clear all;
%db1,4,8,12
%level 1,2,3,



dwtmode('per')
startup;
waves = [];
points = [1];

wavelet = 'db1';

depth =3;
n = 1;
N = 2048;

for leaf = 0:(2^depth) - 1
% leaf = min(leaf,(2^depth) - 1);

% ns=floor(linspace(1,(N/(2^depth))/2,8)); %father
% % ns=floor(linspace(1+(N/(2^depth))/2,(N/(2^depth)),8)); %mother
for i = 1:n

    
    %initialize tree
    y=zeros(N,1);
    T = wpdec(y,depth,wavelet);
    
    %set all to zero
    Nleaves = length(read(T,'data'));
    y=zeros(Nleaves,1);
    T = write(T,'data',y);
    
    %set a singe one
    Nleaf = length(read(T,'data',[depth 0]));
    y = zeros(Nleaf,1);
    
%     sp = floor(linspace(1,Nleaf,n));
    sp = round(Nleaf/2*ones(n,1));
    
    y(sp(i)) = 1;
    T = write(T,'data',[depth leaf], y);
    
    hold on;    
   
    
end



%% select

threshold = 0.01;

 x= wprec(T);
 Nhigh = N;
 Nlow = 1;
 while x(Nlow)  < threshold
    Nlow = Nlow +1; 
 end
 Nlow = Nlow -20;
 while x(Nhigh)  < threshold
    Nhigh = Nhigh -1; 
 end
 Nhigh = Nhigh +20;
 
 
 %% and add to waves
 
 waves = [waves; x(Nlow:Nhigh)];
 points = [points; points(end)+ Nhigh - Nlow+1];

 
end 

%bypass thinking about indexes
 points = [points(1:end-1); length(waves)];

 %% 
 

 %% look at the waves
 
 xmax = 800;
 xmin= 0;
 
 shift = (xmax - points(end))/2;
 
 
 % Create figure
figure1 = figure('Position', [100, 100, 2000, 600]);

% Create axes
axes1 = axes('Parent',figure1);

 xlim([xmin xmax])
hold(axes1,'on');
 
% Create ylabel
% ylabel({'height of components'});

% Create xlabel
xlabel({'position in original signal'});

 for i = 1:(length(points) -1)
     plot((points(i):points(i+1))+shift,waves(points(i):points(i+1)),'DisplayName',strcat('(',num2str(depth),',',num2str(i-1),')'));
     
 end
 
 % Create legend
legend1 = legend(axes1,'show');
set(legend1,'FontSize',9);

%% save it
% waveletDec = strcat(wavelet, ' with depth ');
% saveas(figure1, strcat(waveletDec,num2str( depth )) , 'png')
 
% 