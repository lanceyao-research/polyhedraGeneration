function main()
    addpath('astra\astra-1.9.0.dev11\tools')
    addpath('astra\astra-1.9.0.dev11\mex')
    simPath='prism/tilt/';
    GTPath='prism/stack/';
    nGrid = 256;
    [X, Y, Z] = meshgrid(linspace(-127.5,127.5, nGrid));
    QP = [X(:) Y(:) Z(:)];
    rng(999)
    betas = -90:3:90;
    %%%%%%%%%%%%%%%
    sampleVolume = 1000;
    for sample = 1:sampleVolume
        disp(['current sample: ',num2str(sample),'/',num2str(sampleVolume)])
        simname=strcat(simPath,num2str(sample),'.tif');
        GTname=strcat(GTPath,num2str(sample),'.tif');
        P0 = rand()*40+80;
        P1 = rand()*10+10;
        [tipX,tipY,tipZ] = sphere(100);
        tipX = tipX(:)*P1;
        tipY = tipY(:)*P1;
        tipZ = tipZ(:)*P1;
        T = [tipX - P0/2, tipY, tipZ;
            tipX + P0/2, tipY, tipZ;
            tipX     , tipY + P0/2*sqrt(3), tipZ];
        T(:,2) = T(:,2) - P0/2*sqrt(3)/3;
        scatter3(T(:,1),T(:,2),T(:,3))
        T = rand_rot_3D(T);
        T = rand_trans_3D(T);
        S = delaunayTriangulation(T);
        indexIntersect = (~isnan(pointLocation(S, QP)));
        mask = double(reshape(indexIntersect, [nGrid nGrid nGrid]));
        imSaver(mask,GTname)
        projs = projector(mask,betas,nGrid);
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