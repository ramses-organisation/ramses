# This is a little script to build the diffs for all files in a given patch directory. The
# diffs should be put under version control rather than the actual source code files. This 
# makes the maintainance of RAMSES patches easier when the corresponding source code files
# change. However, it keeps alive the RAMSES patching strategy where various patches can be
# combined, etc.

# Andreas

# USAGE : sh build_diffs.sh path/to/patch
# Should be run from ramses directory (the one that contains /bin /amr etc)


PATCHDIR=$1
SOLVER=hydro
tmp_dir="/tmp/"`date +%s`
mkdir $tmp_dir

cp $SOLVER/*.f90 $tmp_dir
cp aton/*.f90 $tmp_dir
cp hydro/*.f90 $tmp_dir
cp pm/*.f90 $tmp_dir
cp poisson/*.f90 $tmp_dir
cp amr/*.f90 $tmp_dir

for file in `find $PATCHDIR -name "*.f90" -printf "%f\n"`
do 
    echo $file
    diff  $tmp_dir/$file $PATCHDIR$file > $PATCHDIR/$file"_diff"
done

rm -r $tmp_dir