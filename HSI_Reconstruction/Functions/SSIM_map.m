function [ssim_val,map]= SSIM_map(GroundTruth,Reconstructed_Cube)
        [x,y,z]=size(GroundTruth);
        map=zeros(x,y,z);
        ssim_val=0;
         for i =1:z
            [ssimval,ssimmap]=ssim(Reconstructed_Cube(:,:,i),GroundTruth(:,:,i));
            ssim_val=ssim_val+abs(ssimval);
            map(:,:,i)=ssimmap;
         end
        ssim_val=ssim_val/z;
end