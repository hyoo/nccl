set -x

NPKIT_FLAGS_CPU_PREFIX="-DENABLE_NPKIT"
NPKIT_FLAGS_GPU_PREFIX="-DENABLE_NPKIT -DENABLE_NPKIT_EVENT_TIME_SYNC_CPU -DENABLE_NPKIT_EVENT_TIME_SYNC_GPU"
NPKIT_NCCL_TEST_BIN_ALLREDUCE="/nccl-tests-master/build/all_reduce_perf"
NPKIT_NCCL_TEST_BIN_ALLTOALL="/nccl-tests-master/build/alltoall_perf"

export NPKIT_NCCL_TEST_BIN=${NPKIT_NCCL_TEST_BIN_ALLREDUCE}
# export NPKIT_NCCL_TEST_BIN=${NPKIT_NCCL_TEST_BIN_ALLTOALL}

export NPKIT_NCCL_PROTO="Simple"
# export NPKIT_NCCL_PROTO="LL"
# export NPKIT_NCCL_PROTO="LL128"

export NPKIT_NCCL_ALGO="Ring"
# export NPKIT_NCCL_ALGO="Tree"

export NPKIT_FLAGS=${NPKIT_FLAGS_GPU_PREFIX}" -DENABLE_NPKIT_EVENT_ALL_REDUCE_RING_ENTRY -DENABLE_NPKIT_EVENT_ALL_REDUCE_RING_EXIT"
# export NPKIT_FLAGS=${NPKIT_FLAGS_GPU_PREFIX}" -DENABLE_NPKIT_EVENT_ALL_REDUCE_TREE_UPDOWN_ENTRY -DENABLE_NPKIT_EVENT_ALL_REDUCE_TREE_UPDOWN_EXIT"
# export NPKIT_FLAGS=${NPKIT_FLAGS_GPU_PREFIX}" -DENABLE_NPKIT_EVENT_ALL_REDUCE_TREE_SPLIT_ENTRY -DENABLE_NPKIT_EVENT_ALL_REDUCE_TREE_SPLIT_EXIT"

# export NPKIT_FLAGS=${NPKIT_FLAGS_GPU_PREFIX}" -DENABLE_NPKIT_EVENT_COPY_SEND_ENTRY -DENABLE_NPKIT_EVENT_COPY_SEND_EXIT"
# export NPKIT_FLAGS=${NPKIT_FLAGS_GPU_PREFIX}" -DENABLE_NPKIT_EVENT_DIRECT_COPY_SEND_ENTRY -DENABLE_NPKIT_EVENT_DIRECT_COPY_SEND_EXIT"
# export NPKIT_FLAGS=${NPKIT_FLAGS_GPU_PREFIX}" -DENABLE_NPKIT_EVENT_DIRECT_RECV_ENTRY -DENABLE_NPKIT_EVENT_DIRECT_RECV_EXIT"
# export NPKIT_FLAGS=${NPKIT_FLAGS_GPU_PREFIX}" -DENABLE_NPKIT_EVENT_DIRECT_RECV_COPY_SEND_ENTRY -DENABLE_NPKIT_EVENT_DIRECT_RECV_COPY_SEND_EXIT"
# export NPKIT_FLAGS=${NPKIT_FLAGS_GPU_PREFIX}" -DENABLE_NPKIT_EVENT_DIRECT_RECV_REDUCE_COPY_SEND_ENTRY -DENABLE_NPKIT_EVENT_DIRECT_RECV_REDUCE_COPY_SEND_EXIT"
# export NPKIT_FLAGS=${NPKIT_FLAGS_GPU_PREFIX}" -DENABLE_NPKIT_EVENT_DIRECT_SEND_ENTRY -DENABLE_NPKIT_EVENT_DIRECT_SEND_EXIT"
# export NPKIT_FLAGS=${NPKIT_FLAGS_GPU_PREFIX}" -DENABLE_NPKIT_EVENT_DIRECT_SEND_FROM_OUTPUT_ENTRY -DENABLE_NPKIT_EVENT_DIRECT_SEND_FROM_OUTPUT_EXIT"
# export NPKIT_FLAGS=${NPKIT_FLAGS_GPU_PREFIX}" -DENABLE_NPKIT_EVENT_RECV_ENTRY -DENABLE_NPKIT_EVENT_RECV_EXIT"
# export NPKIT_FLAGS=${NPKIT_FLAGS_GPU_PREFIX}" -DENABLE_NPKIT_EVENT_RECV_COPY_SEND_ENTRY -DENABLE_NPKIT_EVENT_RECV_COPY_SEND_EXIT"
# export NPKIT_FLAGS=${NPKIT_FLAGS_GPU_PREFIX}" -DENABLE_NPKIT_EVENT_RECV_REDUCE_COPY_ENTRY -DENABLE_NPKIT_EVENT_RECV_REDUCE_COPY_EXIT"
# export NPKIT_FLAGS=${NPKIT_FLAGS_GPU_PREFIX}" -DENABLE_NPKIT_EVENT_RECV_REDUCE_COPY_SEND_ENTRY -DENABLE_NPKIT_EVENT_RECV_REDUCE_COPY_SEND_EXIT"
# export NPKIT_FLAGS=${NPKIT_FLAGS_GPU_PREFIX}" -DENABLE_NPKIT_EVENT_RECV_REDUCE_SEND_ENTRY -DENABLE_NPKIT_EVENT_RECV_REDUCE_SEND_EXIT"
# export NPKIT_FLAGS=${NPKIT_FLAGS_GPU_PREFIX}" -DENABLE_NPKIT_EVENT_SEND_ENTRY -DENABLE_NPKIT_EVENT_SEND_EXIT"
# export NPKIT_FLAGS=${NPKIT_FLAGS_GPU_PREFIX}" -DENABLE_NPKIT_EVENT_SEND_FROM_OUTPUT_ENTRY -DENABLE_NPKIT_EVENT_SEND_FROM_OUTPUT_EXIT"

# export NPKIT_FLAGS=${NPKIT_FLAGS_GPU_PREFIX}" -DENABLE_NPKIT_EVENT_PRIM_SIMPLE_WAIT_PEER_ENTRY -DENABLE_NPKIT_EVENT_PRIM_SIMPLE_WAIT_PEER_EXIT"
# export NPKIT_FLAGS=${NPKIT_FLAGS_GPU_PREFIX}" -DENABLE_NPKIT_EVENT_PRIM_SIMPLE_REDUCE_OR_COPY_MULTI_ENTRY -DENABLE_NPKIT_EVENT_PRIM_SIMPLE_REDUCE_OR_COPY_MULTI_EXIT"

# export NPKIT_FLAGS=${NPKIT_FLAGS_GPU_PREFIX}" -DENABLE_NPKIT_EVENT_PRIM_LL_WAIT_SEND_ENTRY -DENABLE_NPKIT_EVENT_PRIM_LL_WAIT_SEND_EXIT"
# export NPKIT_FLAGS=${NPKIT_FLAGS_GPU_PREFIX}" -DENABLE_NPKIT_EVENT_PRIM_LL_DATA_PROCESS_ENTRY -DENABLE_NPKIT_EVENT_PRIM_LL_DATA_PROCESS_EXIT"

# export NPKIT_FLAGS=${NPKIT_FLAGS_GPU_PREFIX}" -DENABLE_NPKIT_EVENT_PRIM_LL128_WAIT_SEND_ENTRY -DENABLE_NPKIT_EVENT_PRIM_LL128_WAIT_SEND_EXIT"
# export NPKIT_FLAGS=${NPKIT_FLAGS_GPU_PREFIX}" -DENABLE_NPKIT_EVENT_PRIM_LL128_DATA_PROCESS_ENTRY -DENABLE_NPKIT_EVENT_PRIM_LL128_DATA_PROCESS_EXIT"

# export NPKIT_FLAGS=${NPKIT_FLAGS_CPU_PREFIX}" -DENABLE_NPKIT_EVENT_NET_SEND_ENTRY -DENABLE_NPKIT_EVENT_NET_SEND_EXIT -DENABLE_NPKIT_EVENT_NET_RECV_ENTRY -DENABLE_NPKIT_EVENT_NET_RECV_EXIT"

bash run_nccl_tests_with_npkit.sh