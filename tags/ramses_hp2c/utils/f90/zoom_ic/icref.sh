

/home/rteyssie/icref/geticref -gid true -xc 0.5698 -yc 0.5905 -zc 0.5615 -rad 0.006 -inp output_00385
/home/rteyssie/icref/geticref -fid true -inp output_00001 -smt 4 -rsm 12 -gfc ic_files/AqC_boxlen12p5_n256

#mv ic_refmap ic_files/AqC_boxlen12p5_n256/.
#/home/rteyssie/icref/icdegrade ic_files/AqC_boxlen12p5_n256 ic_files/AqC_boxlen12p5_n128
#/home/rteyssie/icref/icdegrade ic_files/AqC_boxlen12p5_n128 ic_files/AqC_boxlen12p5_n64
#/home/rteyssie/icref/icinject  ic_files/AqC_boxlen12p5_n64  ic_files/AqC_boxlen25_n128
#/home/rteyssie/icref/icdegrade ic_files/AqC_boxlen25_n128   ic_files/AqC_boxlen25_n64
#/home/rteyssie/icref/icinject  ic_files/AqC_boxlen25_n64    ic_files/AqC_boxlen50_n128
#/home/rteyssie/icref/icdegrade ic_files/AqC_boxlen50_n128   ic_files/AqC_boxlen50_n64
#/home/rteyssie/icref/icinject  ic_files/AqC_boxlen50_n64    ic_files/AqC_boxlen100_n128





