clc; clear
addpath('align_config')
options = parse_yaml(['align_config', filesep, 'optimize.yaml']);
parentpath = get_parent_path;
ext = '.png';

anti_alias_scl = 1;

shareM0 = false;
rigidFirst = true;
affineFirst = true;
toRender = false;
elastic = true;
full_cover = true;


imgpath = [parentpath, filesep, 'stitched'];
maskrpath = [parentpath, filesep, 'masks_region'];
maskwpath = [parentpath, filesep, 'masks_wrinkle'];
ptpairpath = [parentpath, filesep, 'ptpairs'];
outmeshpath = [parentpath, filesep, 'meshes'];
outimgpath = [parentpath, filesep, 'output'];
outxypath = [parentpath, filesep, 'displacement'];

imglist = dir([imgpath, filesep, '*', ext]);
align_endpt = {};
reference_names = {'control_section57'};
if isempty(imglist)
    return
end
if ~exist(outmeshpath, 'dir')
    mkdir(outmeshpath);
end
if ~exist(outxypath, 'dir')
    mkdir(outxypath);
end
if ~exist(outimgpath, 'dir')
    mkdir(outimgpath);
end

[~, refidx] = ismember(reference_names, {imglist.name});
if isempty(align_endpt)
    imglist = imglist(min(refidx):max(refidx));
else
    [~, tidx] = ismember(align_endpt, {imglist.name});
    imglist = imglist(min(tidx):max(tidx));
end
Nimg = numel(imglist);
matlist = dir([ptpairpath, filesep, '*.mat']);
matnames = strrep({matlist.name},'.mat', '');
matnames = utils.split_pair_name(matnames(:));
imagenames = strrep({imglist.name}, ext, '');
idxt = contains(matnames(:,1), imagenames,'IgnoreCase',true) & ...
    contains(matnames(:,2), imagenames,'IgnoreCase',true);
matlist = matlist(idxt);
matnames = strrep({matlist.name},'.mat', '');
matnames = utils.split_pair_name(matnames(:));

Nlink = numel(matlist);
%%
if shareM0
    IMG = imread([imgpath, filesep, imglist(1).name]);
    M0 = elastic_mesh.gen_eqtriang_mesh(size(IMG), options.optimizer.mesh_space);
end
Ms = cell(Nimg, 1);
fixed_sec = false(Nimg, 1);
[~, matidxs] = ismember(matnames, imagenames);
imagenames = {imglist.name};
A = repmat(eye(3), 1,1,Nimg);

for k = 1:Nimg
    if any(strcmpi(reference_names, imglist(k).name))
        load([outmeshpath, filesep, strrep(imglist(k).name, ext, '.mat')], 'M')
        Ms{k} = M;
        fixed_sec(k) = true;
    else
        IMG = imread([imgpath, filesep, imglist(k).name]);
        if exist([maskrpath, filesep, imglist(k).name], 'file')
            maskr = imread([maskrpath, filesep, imglist(k).name]);
        else
            maskr = ones(size(IMG), 'uint8');
        end
        if exist([maskwpath, filesep, imglist(k).name], 'file')
            maskw = imread([maskwpath, filesep, imglist(k).name]);
        else
            maskw = zeros(size(IMG), 'uint8');
        end
        if ~shareM0
            M0 = elastic_mesh.gen_eqtriang_mesh_mask(maskw<255, options.optimizer.mesh_space, 1);
        end
        Ms{k} = elastic_mesh.init_mesh_subregion(M0, maskr, maskw == 200, options.optimizer);
    end
    if rigidFirst
        idxt = find(any(matidxs == k, 2) & (sum(matidxs,2) < 2*k));
        if ~isempty(idxt)
            n_ref = numel(idxt(:));
            XY0 = cell(n_ref,1);
            XY1 = cell(n_ref,1);
            for t = 1:n_ref
                load([matlist(idxt(t)).folder, filesep, matlist(idxt(t)).name], 'ptpairs','imgnames');
                k1 = find(strcmpi(imagenames, imgnames{1}));
                k2 = find(strcmpi(imagenames, imgnames{2}));
                kt = min(k1, k2);
                if k1 == k
                    ptpairs = geometries.flip_ptpairs(ptpairs);
                end
                XY1{t} = fliplr(ptpairs.yx2);
                XY0{t} = fliplr(ptpairs.yx1) * A(1:2,1:2,kt) + A(3,1:2,kt);
            end
            xy1 = cat(1, XY1{:});
            xy0 = cat(1, XY0{:});
            [~,R] = geometries.fit_affine(xy0, xy1);
            A(:,:,k) = R;
        end 
    end
end

if rigidFirst
    idxt = find(fixed_sec, 1);
    xy0 = Ms{idxt}.TR0.Points;
    xy1 = Ms{idxt}.TR.Points;
    pt_num0 = Ms{idxt}.pt_num0;
    region_num = Ms{idxt}.region_num;
    xy1 = permute(reshape(xy1',2, pt_num0, region_num),[2,1,3]);
    xy1 = nanmean(xy1, 3);
    R0 = geometries.fit_affine(xy1, xy0);
    R0 = A(:,:,idxt) \ R0;
    for t = 1:numel(Ms)
        if fixed_sec(t)
            continue
        end
        R = A(:,:,t) * R0;
        Ms{t}.TR.Points = Ms{t}.TR.Points * R(1:2,1:2) + R(3,1:2);
    end
end

links =  cell(Nlink,1);
for k = 1:Nlink
    load([matlist(k).folder, filesep, matlist(k).name], 'ptpairs','imgnames');
    k1 = find(strcmpi(imagenames, imgnames{1}));
    k2 = find(strcmpi(imagenames, imgnames{2}));
    lnk = elastic_mesh.ptpairs_to_links(ptpairs, Ms{k1}, Ms{k2});
    lnk_struct = elastic_mesh.link_struct(lnk, [k1, k2]);
    links{k} = lnk_struct;
end
links = vertcat(links{:});

%%
Nrep = 100;
cost_history = cell(Nrep+1, 1);
h0 = options.optimizer.params.learning_rate;
options.optimizer.params.learning_rate = 15 * h0;
[Ms_out, cost_history0] = elastic_mesh.BGD(Ms, links, options.optimizer.params, fixed_sec, 1);
if affineFirst
    for k = 1:30
        tic;
        [Ms_out, cost_history0] = elastic_mesh.BGD(Ms_out, links, options.optimizer.params, fixed_sec, 3);
        toc;
    end
    cost_history{1} = cost_history0(~isnan(cost_history0));
end
%%
options.optimizer.params.learning_rate = h0;
if any(~fixed_sec)
    disp(' ')    
    for k = 1:Nrep
        tic;
        [Ms_out, cost_history0] = elastic_mesh.BGD(Ms_out, links, options.optimizer.params, fixed_sec);
        cost_history{k+1} = cost_history0(~isnan(cost_history0));
        toc;
    end
    plot(cost_history0)
end

imglist0 = dir([outimgpath, filesep, '*', ext]);
refIMG = imread([outimgpath, filesep, imglist0(1).name]);
outsz = size(refIMG);

for k = 1:Nimg
    if fixed_sec(k) && any(~fixed_sec)
        continue
    end
    M = Ms_out{k};
    img = imread([imgpath, filesep, imglist(k).name]);
    if exist([maskrpath, filesep, imglist(k).name], 'file')
        maskr = imread([maskrpath, filesep, imglist(k).name]);
    else
        maskr = ones(size(img), 'uint8');
    end
    if exist([maskwpath, filesep, imglist(k).name], 'file')
        maskw = uint8(imread([maskwpath, filesep, imglist(k).name]) == 200);
    else
        maskw = zeros(size(img), 'uint8');
    end
    save([outmeshpath, filesep, strrep(imglist(k).name, ext, '.mat')], 'M', 'maskr', 'maskw','outsz');
end
%%

if ~toRender
    return
end


for k = 1:Nimg
    tic;
    if fixed_sec(k) && any(~fixed_sec)
        continue
    end
    M = Ms_out{k};
    img = imread([imgpath, filesep, imglist(k).name]);
    if exist([maskrpath, filesep, imglist(k).name], 'file')
        maskr = imread([maskrpath, filesep, imglist(k).name]);
    else
        maskr = ones(size(img), 'uint8');
    end
    if exist([maskwpath, filesep, imglist(k).name], 'file')
        maskw = uint8(imread([maskwpath, filesep, imglist(k).name]) == 200);
    else
        maskw = zeros(size(img), 'uint8');
    end
    [mxy, wt] = elastic_mesh.subregion_movingcoord(Ms_out{k}, cat(3,maskr, maskw), options.render, [], outsz);
    D = imgaussfilt(elastic_mesh.movcoord2displ(mxy),3);
    imgt = imwarp(img, D);
    % imgt = imresize(imwarp(imresize(img,anti_alias_scl), imresize(D*anti_alias_scl,anti_alias_scl), 'cubic'), outsz);
    imagesc(imgt); axis equal; colormap(gray);shg
    imwrite(imgt, [outimgpath, filesep, imglist(k).name]);
    save([outmeshpath, filesep, strrep(imglist(k).name, ext, '.mat')], 'M', 'maskr', 'maskw','outsz');
    % save([outxypath, filesep, strrep(imglist(k).name, ext, '.mat')],
    % 'mxy');
    toc;
end
