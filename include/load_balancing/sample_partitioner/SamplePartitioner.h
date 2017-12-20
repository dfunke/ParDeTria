#pragma once
#include "load_balancing/Partitioner.h"
#include "Sampler.h"

namespace LoadBalancing
{
    template <uint D, typename Precision>
    struct SamplePartitioner : Partitioner<D, Precision>
    {
        virtual const Sampling<D, Precision>& sampling() const = 0;
    };
}
