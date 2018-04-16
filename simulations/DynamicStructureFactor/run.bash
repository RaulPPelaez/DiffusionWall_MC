
ROOT=../..


cd $ROOT

make 

if ! ls mc
then
    echo "ERROR: Compilation failed! I cant find mc! Check Makefile"
    cd -
    exit 1
else

    cd -
    cp $ROOT/mc .
    cp $ROOT/libCall*so .
    
    ./mc > uammd.log 2>&1

    cd tools

    bash analysis.bash

    cd ..
fi
