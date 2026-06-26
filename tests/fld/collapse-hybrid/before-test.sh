# Warning, it's assumed that you have all the files needed (EOS, opacities) in a directory lib/collapse in the test_suite directory

lib_path=../../lib/collapse;

ln -s ${lib_path}/* .
ls ${lib_path} > to_be_removed;
