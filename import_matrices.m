clear 
close all
%int2str
%strcat
dims = [20,40,80];
pxdim=10.*dims;
for i=1:length(dims)
    dimstr = int2str(dims(i));
    pxstr = int2str(pxdim(i));
    
    filename = strcat('analytical_',dimstr,'_',dimstr,'.dat');
    orig=readmatrix(filename);
    filename = strcat('schwefel_',dimstr,'_',dimstr,'_',pxstr,'_',pxstr,'.dat');
    bilin=readmatrix(filename);
    filename = strcat('analytical_',pxstr,'_',pxstr,'.dat');
    ana = readmatrix(filename);
    
    [bilin,ana] = reshapematrices(pxdim(i),pxdim(i),size(bilin,2),bilin,ana);
    orig=reshape(orig,dims(i),dims(i),3);
    plotall(bilin,'bilinear')
    plotall(ana,'exact')
    plotall(orig,'input')
    
end

%[X,Y,Z] = peaks(25);
%CO(:,:,1) = zeros(25); % red
%CO(:,:,2) = ones(25).*linspace(0.5,0.6,25); % green
%CO(:,:,3) = ones(25).*linspace(0,1,25); % blue
%surf(X,Y,Z,CO)

function plotall(matrix,title_)
    x=matrix(:,:,1);y=matrix(:,:,2);z=matrix(:,:,3);
    figure()
    surf(x,y,z,z,'edgecolor','none');
    title(title_)
    colorbar
    figure()
    contourf(x,y,z,20);
    title(title_)
end

function varargout = reshapematrices(d1,d2,d3,varargin)
    %varargout{1} = zeros(d1*d2,d3);
    for ii = 1:length(varargin)
        varargout{ii}=varargin{ii};
        varargout{ii}=reshape(varargout{ii},d1,d2,d3);
    end
end