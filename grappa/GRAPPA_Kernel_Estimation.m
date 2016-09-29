

function coef = GRAPPA_Kernel_Estimation(Ref_Img, option)

if option.2D
    coef = GRAPPA_Kernel_2D(Ref_Img, option);
    
elseif option.3D
    coef = GRAPPA_Kernel_3D(Ref_Img, option);

else
    coef = GRAPPA_Kernel_2D(Ref_Img, option);

end

return