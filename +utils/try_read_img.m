function [img, imgname] = try_read_img(imgpaths, imgnames, readimg)
    if nargin < 3
        readimg = true;
    end
    img = [];
    if ~iscell(imgpaths) && ~iscell(imgnames)
        imgdir = [imgpaths, filesep, imgnames];
        if ~exist(imgdir, 'file')
            imgname = [];
            return
        end
        imgname = struct('folder',imgpaths,'name',imgnames);
        if readimg
            img = imread(imgdir);
        end
        return
    end
    if ~iscell(imgpaths)
        imgpaths = {imgpaths};
    end
    if ~iscell(imgnames)
        imgnames = {imgnames};
    end
    v = combvec(1:numel(imgpaths), 1:numel(imgnames));
    for t = 1:size(v,2)
        p0 = v(1,t);
        n0 = v(2,t);
        imgdir = [imgpaths{p0}, filesep, imgnames{n0}];
        if exist(imgdir, 'file')
            imgname = struct('folder',imgpaths{p0},'name',imgnames{n0});
            if readimg
                img = imread(imgdir);
            end
            return
        end
    end
    imgname = [];
end