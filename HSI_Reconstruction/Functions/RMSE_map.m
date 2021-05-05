function [map]= RMSE_map(GroundTruth,Reconstructed_Cube)
        [X,Y,Z]=size(GroundTruth);
        map=zeros(X,Y);
         for x =1:X
             for y =1:Y
                so=squeeze(GroundTruth(x,y,:));
                sr=squeeze(Reconstructed_Cube(x,y,:));
                RMSE=norm(sr-so)/norm(so);
                map(x,y)=RMSE;
             end
         end
        
end