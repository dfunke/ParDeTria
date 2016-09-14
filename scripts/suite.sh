bin/getRun
run=$?
echo "Run number: $run"

REPS=3
THREADS=32

#uniform
bin/Benchmark --minN 10000000 --maxN 50000000 --occupancy 100 --no-parallel-base --reps ${REPS} --maxThreads ${THREADS} --dist u --run-number ${run}

#normal
bin/Benchmark --minN 1000000  --maxN 5000000  --occupancy 100 --no-parallel-base --reps ${REPS} --maxThreads ${THREADS} --dist n --run-number ${run}

#bubble
bin/Benchmark --minN 10000000 --maxN 50000000 --occupancy 100 --no-parallel-base --reps ${REPS} --maxThreads ${THREADS} --dist b --run-number ${run}

#ellipsoid
bin/Benchmark --minN 100000   --maxN 500000   --occupancy 100 --no-parallel-base --reps ${REPS} --maxThreads ${THREADS} --dist e --run-number ${run}

#lines
bin/Benchmark --minN 10000    --maxN 50000    --occupancy 100 --no-parallel-base --reps ${REPS} --maxThreads ${THREADS} --dist l --run-number ${run}
