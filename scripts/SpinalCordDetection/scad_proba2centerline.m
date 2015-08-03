function scad_proba2centerline(fname)
% scad_proba2centerline(fname)
img=load_nii(fname);
centerline=false(img.dims(1:3));

%% search for maximum slice by slice
for iz=1:img.dims(3)
    slice=img.img(:,:,iz,1);
    slice=imgaussian(slice,3/mean(img.scales(1:2))); % smooth (3mm) for higher robustness in maxium 
    [~,I]=max(slice(:));
    [X,Y]=ind2sub(img.dims(1:2),I);
    centerline(X,Y,iz)=true;
end

%% todo : z regulation
%% save centerline
save_nii_v2(centerline,'centerline.nii.gz',fname,8)

disp(['sct_unix(''fslview ' fname ' centerline.nii.gz -l Red'')'])