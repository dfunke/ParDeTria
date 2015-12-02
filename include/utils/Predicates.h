#pragma once

template<uint D, typename Precision>
class Predicates {

    public:
        static Precision orient(const Precision *pa, const Precision *pb, const Precision *pc);

        static Precision orient(const Precision *pa, const Precision *pb, const Precision *pc, const Precision *pd);

        static Precision insphere(const Precision *pa, const Precision *pb, const Precision *pc, const Precision *pd);

        static Precision insphere(const Precision *pa, const Precision *pb, const Precision *pc, const Precision *pd, const Precision *pe);

};