function scad_symetry(fname)
% scad_symetry(fname)
% example: scad_symetry HC_SC_002.nii.gz
u=load_nii(fname);
clear out
p=squeeze(sum(u.img,2));
for iz=1:size(p,2)
    out(:,iz)=xcorr(p(:,iz),p(end:-1:1,iz));
    out(:,iz)=out(:,iz)./max(out(:,iz));
end

% halfing
index1=round(linspace(1,size(out,1)-1,u.dims(1)));
index2=round(linspace(2,size(out,1),u.dims(1)));
out=1/2*(out(index1,:)+out(index2,:));

save_nii_v2(repmat(permute(out,[1 3 2]),[1 u.dims(2) 1]),[sct_tool_remove_extension(fname,1),'_symetry.nii.gz'],fname,64)
disp(['unix(''fslview ' fname ' ' sct_tool_remove_extension(fname,1) '_symetry_matlab -l "Red" -t 0.5 -b 0.8,1'')'])