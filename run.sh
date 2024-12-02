
if [ $# -ne 2 ]
then
    echo "Usage: $0 <CPU/GPU> <normalCT/offsetCT>"
    exit -1
fi

tmp=`dirname $0`
PROJECT_ROOT=`cd $tmp; pwd`
cd ${PROJECT_ROOT}

PLATFORM=$1
MODE=$2
#cd ${PLATFORM}

#build_path="build"
#if [[ -d $build_path ]]; then
#    rm -rf $build_path
#fi
#mkdir -p $build_path

output_path="output"
if [[ -d $output_path ]]; then
    rm -rf $output_path
fi
mkdir -p ${output_path}
mkdir -p ${output_path}/$PLATFORM

#cd $build_path
#cmake .. 
#make -j4

NB=128

if [[ "$PLATFORM" == "CPU" ]]; then
    exe=./PBCT/x64/Release/cpuCT.exe
else
    exe=./PBCT/x64/Release/gpuCT.exe
fi


if [ $MODE == 'normalCT' ]; then
    cmd="$exe ./data/raw ./$output_path/${PLATFORM} 300 740 900 0 391.0 740 3 $NB"
elif [ $MODE == 'offsetCT' ]; then
    cmd="$exe ./data/offsetCT ./$output_path/${PLATFORM} 32 1000 3600 1 909.0 1818 4 $NB"
else
    echo "error MODE <normalCT/offsetCT>"
    exit
fi
echo cmd=$cmd; eval $cmd;

