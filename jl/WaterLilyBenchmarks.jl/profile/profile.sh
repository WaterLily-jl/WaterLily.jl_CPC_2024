#!/bin/bash
# Usage example
#
# sh profile.sh -c "tgv sphere cylinder" -p "8 5 6"


THIS_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Utils
# Grep current julia version
julia_version () {
    julia_v=($(julia -v))
    echo "${julia_v[2]}"
}
# Grep current julia version
waterlily_profile_branch () {
    cd $WATERLILY_DIR
    git checkout kernel_profiling
    julia --project -e "using Pkg; Pkg.update();"
    cd $THIS_DIR
}
update_environment () {
    if check_if_juliaup; then
        echo "Updating environment to Julia $version"
        julia --project=$THIS_DIR -e "using Pkg; Pkg.develop(PackageSpec(path=get(ENV, \"WATERLILY_DIR\", \"\"))); Pkg.update();"
    fi
}

# Run profiling
run_profiling () {
    full_args=(--project=${THIS_DIR} --startup-file=no $args)
    echo "Running: nsys profile --force-overwrite true -o $THIS_DIR/data/$case/$case julia ${full_args[@]}"
    nsys profile --force-overwrite true -o $THIS_DIR/data/$case/$case julia "${full_args[@]}"
}

# Print benchamrks info
display_info () {
    echo "--------------------------------------"
    echo "Running benchmark tests for:
 - Julia:        $VERSION
 - Backends:     $BACKEND
 - Cases:        ${CASES[@]}
 - Size:         ${LOG2P[@]:0:$NCASES}
 - Sim. steps:   $MAXSTEPS
 - Data type:    $FTYPE"
    echo "--------------------------------------"; echo
}

# Default backends
JULIA_USER_VERSION=$(julia_version)
VERSION=$JULIA_USER_VERSION
BACKEND='CuArray'
# Default cases. Arrays below must be same length (specify each case individually)
CASES=() # ('tgv' 'sphere' 'cylinder')
LOG2P=() # ('5' '5' '5')
MAXSTEPS='250'
FTYPE='Float32'

# Parse arguments
while [ $# -gt 0 ]; do
case "$1" in
    --versions|-v)
    VERSION=($2)
    shift
    ;;
    --backends|-b)
    BACKEND=($2)
    shift
    ;;
    --cases|-c)
    CASES=($2)
    shift
    ;;
    --log2p|-p)
    LOG2P=($2)
    shift
    ;;
    --max_steps|-s)
    MAXSTEPS=($2)
    shift
    ;;
    --float_type|-ft)
    FTYPE=($2)
    shift
    ;;
    *)
    printf "ERROR: Invalid argument %s\n" "${1}" 1>&2
    exit 1
esac
shift
done

# Assert all case arguments have equal size
NCASES=${#CASES[@]}
NLOG2P=${#LOG2P[@]}
st=0
for i in $NLOG2P; do
    [ "$NCASES" = "$i" ]
    st=$(( $? + st ))
done
if [ $st != 0 ]; then
    echo "ERROR: cases and log2p arrays of different sizes."
    exit 1
fi

# Display information
display_info

# Checkout to WaterLily profiling branch and update it
# waterlily_profile_branch

# Update this environment
# update_environment

# Profiling
args_case="--backend=$BACKEND --max_steps=$MAXSTEPS --ftype=$FTYPE"
for ((i = 0; i < ${#CASES[@]}; ++i)); do
    case=${CASES[$i]}
    mkdir -p $THIS_DIR/data/$case
    args="${THIS_DIR}/profile.jl --case=$case --log2p=${LOG2P[$i]} $args_cases"
    run_profiling
done

# Postprocessing results



echo "All done!"
exit 0
