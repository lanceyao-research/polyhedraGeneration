function main()
    addpath('astra\astra-1.9.0.dev11\tools')
    addpath('astra\astra-1.9.0.dev11\mex')
    simPath='CScube/tilt/';
    GTPath='CScube/stack/';
    nGrid = 256;
    [X, Y, Z] = meshgrid(linspace(-127.5,127.5, nGrid));
    QP = [X(:) Y(:) Z(:)];
    rng(999)
    betas = -90:3:90;
    %%%%%%%%%%%%%%%
    sampleVolume = 400;
    [x0,y0,z0] = sphere(100);x0 = x0(:);y0 = y0(:);z0 = z0(:);
    for sample = 1:sampleVolume
        disp(['current sample: ',num2str(sample),'/',num2str(sampleVolume)])
        simname=strcat(simPath,num2str(sample),'.tif');
        GTname=strcat(GTPath,num2str(sample),'.tif');
        
        t_rad = 5+rand()*15; %5~20 = 10
        x1 = x0 *t_rad;
        y1 = y0 *t_rad;
        z1 = z0 *t_rad;
        sz = 60+rand()*40; %50~150 =100
        T = [x1+sz/2,     y1+sz/2,     z1+sz/2;
                x1+sz/2,     y1-sz/2,     z1+sz/2;
                x1-sz/2,     y1-sz/2,     z1+sz/2;
                x1-sz/2,     y1+sz/2,     z1+sz/2;
             
                x1+sz/2,     y1+sz/2,    z1-sz/2;
                x1+sz/2,     y1-sz/2,    z1-sz/2;
                x1-sz/2,     y1-sz/2,    z1-sz/2;
                x1-sz/2,     y1+sz/2,    z1-sz/2];
        
        t_rad = 1+rand()*4; %1~5 = 10
        x1 = x0 *t_rad;
        y1 = y0 *t_rad;
        z1 = z0 *t_rad;
        sz = sz*(0.5+rand()*0.3);
        T2 = [x1+sz/2,     y1+sz/2,     z1+sz/2;
                x1+sz/2,     y1-sz/2,     z1+sz/2;
                x1-sz/2,     y1-sz/2,     z1+sz/2;
                x1-sz/2,     y1+sz/2,     z1+sz/2;
             
                x1+sz/2,     y1+sz/2,    z1-sz/2;
                x1+sz/2,     y1-sz/2,    z1-sz/2;
                x1-sz/2,     y1-sz/2,    z1-sz/2;
                x1-sz/2,     y1+sz/2,    z1-sz/2];

        %scatter3(T(:,1),T(:,2),T(:,3))
        T_stack = rand_rot_3D([T;T2]);
        T_stack = rand_trans_3D(T_stack);
        T = T_stack(1:size(T_stack,1)/2,:);
        T2 = T_stack(size(T_stack,1)/2+1:end,:);

        S1 = delaunayTriangulation(T);
        indexIntersect1 = (~isnan(pointLocation(S1, QP)));
        mask1 = double(reshape(indexIntersect1, [nGrid nGrid nGrid]));

        S2 = delaunayTriangulation(T2);
        indexIntersect2 = (~isnan(pointLocation(S2, QP)));
        mask2 = reshape(indexIntersect2, [nGrid nGrid nGrid]);
        mask1(mask2) = 0.5;


        imSaver(mask1,GTname)
        projs = projector(mask1,betas,nGrid);
        imSaver(projs/max(projs(:)),simname)
    end
    
    
end

function vertices = rand_rot_3D(vertices)
    t = rand()*360;
    RZ1 = [cosd(t) -sind(t) 0;
          sind(t)  cosd(t) 0;
          0 0 1];  
    t = rand()*360;
    RX = [1           0            0;
           0 cosd(t) -sind(t);
          0 sind(t)  cosd(t)];
    t = rand()*360;
    RZ2 = [cosd(t) -sind(t) 0;
          sind(t)  cosd(t) 0;
          0 0 1];      
    vertices=vertices*RZ1;
    vertices=vertices*RX;
    vertices=vertices*RZ2;   
end

function vertices = rand_trans_3D(vertices)
    vertices = vertices + (rand([1,3])*10-5);
end

function proj = projector(mask,betas,nGrid)
    vol_geom = astra_create_vol_geom(nGrid,nGrid,nGrid);
    proj_geom = astra_create_proj_geom('parallel3d', 1.0, 1.0, nGrid, nGrid, betas/180*pi+pi/2); %carefully cabrilated, don't change
    [proj_id, proj_data] = astra_create_sino3d_cuda(permute(mask,[3,2,1]), proj_geom, vol_geom);%carefully cabrilated, don't change
    proj_data = permute(proj_data,[3,1,2]);%carefully cabrilated, don't change
    proj = proj_data;
    astra_mex_data3d('delete', proj_id);
end

function imSaver(im,path)
    for i = 1:size(im,3)
        if i==1
            imwrite(im(:, :, i), path, 'WriteMode', 'overwrite',  'Compression','none');
        else
            imwrite(im(:, :, i), path, 'WriteMode', 'append',  'Compression','none');
        end 
    end

end