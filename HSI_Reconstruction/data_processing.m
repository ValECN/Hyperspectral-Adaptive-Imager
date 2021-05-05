function [I,FCS] = data_processing()

addpath A_FC 
addpath A_CCD

path_FC  = [pwd, '/A_FC'];
path_CCD = [pwd, '/A_CCD'];
FC       = dir(path_FC);
CCD        = dir(path_CCD);
CCD(1:2) = [];
FC(1:2) = [];
l_FC     = length(FC);
l_CCD    = length(CCD);
I = zeros(26,26,l_CCD);
FCS = zeros(26,26,264,l_FC);

if l_FC == l_CCD
    for k = 1:l_FC
        name_I     = CCD(k).name;
        name_FC    = FC(k).name;
        a          = load(name_I);
        fc         = load(name_FC);
        I(:,:,k)   = a.data;
        FCS(:,:,:,k) = fc.data;
    end
end
        


        
        
       
    
    
