if [ -f to_be_removed ]; then
    rm_list=$(cat to_be_removed);
    rm -f $rm_list;
    rm to_be_removed;
fi
