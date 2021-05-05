function [map]= SAM_map(GroundTruth,Reconstructed_Cube)
        [X,Y,Z]=size(GroundTruth);
        map=zeros(X,Y);
         for x =1:X
             for y =1:Y
                s1=squeeze(GroundTruth(x,y,:));
                s2=squeeze(Reconstructed_Cube(x,y,:));
                SAM=acos((s1'*s2)/(norm(s1)*norm(s2)));
                map(x,y)=SAM;
             end
         end
        
end