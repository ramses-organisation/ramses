# This is a little script to apply the diffs in a patch directory to the corresponding file
# in your local working tree. The idea is to keep only the diffs in the patch directories under 
# version control and then apply this script to the diffs to create the actual source code files 
# which are then compiled the usual way. This makes the maintainance of RAMSES patches easier when 
# the corresponding source code files change. However, it keeps alive the RAMSES patching strategy 
# where various patches can be combined, etc.

# Andreas

# USAGE : sh apply_diffs.sh path/to/patch
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

wdir=`pwd`

cd $tmp_dir

for file in `find $wdir/$PATCHDIR -name "*.f90_diff" -printf "%f\n"`
do 
    echo $wdir/$PATCHDIR/$file
    patch `echo $file | sed -e 's/_diff//g'` < $wdir/$PATCHDIR/$file
    cp `echo $file | sed -e 's/_diff//g'` $wdir/$PATCHDIR
done

cd $wdir

rm -r $tmp_dir