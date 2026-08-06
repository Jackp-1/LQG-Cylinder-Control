function eta=Eta_calc(Base_mean_drag,Controlled_mean_drag,min_cd)

eta=(Base_mean_drag-Controlled_mean_drag)/(Base_mean_drag-min_cd)*100;

end