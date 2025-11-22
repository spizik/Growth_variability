## calculates RWI
for(i in names(prepared.data.tree.means)){
  
  # print(i)
  
  write.table(create.rwi(prepared.data.tree.means[[i]], detrendig_method), paste0("Calculated_datasets/rwi_data/",i,".txt"), dec = ".", sep = ";")
}


