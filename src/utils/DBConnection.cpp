//
// Created by dfunke on 4/7/15.
//

#include "DBConnection.h"
#include "ASSERT.h"

DBConnection::DBConnection(const std::string &file, const std::string &collection,
                           const int mode) {
    m_ejdb = ejdbnew();

    if (!m_ejdb) {
        _handleError("EJDB Error: Could not create EJDB instance");
    }

    if (!ejdbopen(m_ejdb, file.c_str(), mode))
        _handleError();

    m_coll = ejdbgetcoll(m_ejdb, collection.c_str());
    if (!m_coll && mode & JBOWRITER)
        m_coll = ejdbcreatecoll(m_ejdb, collection.c_str(), nullptr);
    if (!m_coll)
        _handleError();

}

DBConnection::~DBConnection() {
    if (m_ejdb) {
        if (ejdbisopen(m_ejdb))
            ejdbclose(m_ejdb);

        ejdbdel(m_ejdb);
    }
}

std::string DBConnection::_getMaximum(const std::string &field) const {

    /* Execute the following query:
     * Query: *
     * Hint: { "$orderby" : [ (field, -1) ],
               "$fields" :  { field : 1 } }
     */

    std::string result;

    bson bsquery;
    bson_init_as_query(&bsquery);
    bson_finish(&bsquery);

    if (bsquery.err)
        _handleError(bsquery.errstr);

    bson bshints;
    bson_init_as_query(&bshints);
    bson_append_start_object(&bshints, "$orderby");
    bson_append_int(&bshints, field.c_str(), -1); //DESC order on field
    bson_append_finish_object(&bshints);
    bson_append_start_object(&bshints, "$fields");
    bson_append_int(&bshints, field.c_str(), 1); //include field
    bson_append_finish_object(&bshints);
    bson_finish(&bshints);

    if (bshints.err)
        _handleError(bshints.errstr);

    EJQ *q1 = ejdbcreatequery(m_ejdb, &bsquery, nullptr, 0, &bshints);
    if (!q1)
        _handleError();

    uint32_t count;
    TCLIST *res = ejdbqryexecute(m_coll, q1, &count, JBQRYFINDONE, NULL);
    if (!res)
        _handleError();

    if (count < 1)
        result = "";
    else {

        const char *bsdata = (const char *) TCLISTVALPTR(res, 0);
        //bson_print_raw(bsdata, 0);

        bson_iterator it;
        bson_iterator_from_buffer(&it, bsdata);

#ifdef NDEBUG
        __attribute__((unused))
#endif
        bson_type bt = bson_find_from_buffer(&it, bsdata, field.c_str());
        ASSERT(bt == BSON_STRING);

        const char *value = bson_iterator_string(&it);
        result = std::string(value);

    }

    //Dispose result set
    ejdbqresultdispose(res);

    //Dispose query
    ejdbquerydel(q1);
    bson_destroy(&bshints);
    bson_destroy(&bsquery);

    return result;
}

#include "Timings.h"

template<>
void DBConnection::bsonify(bson *bobj, const ExperimentRun &o) {

    for (const auto &t : o.traits()) {
        bson_append_string(bobj, t.first.c_str(), t.second.c_str());
    }

    bson_append_start_array(bobj, "times");
    for (uint i = 0; i < o.times().size(); ++i) {
        bson_append_long(bobj, std::to_string(i).c_str(), o.times()[i].count());
    }
    bson_append_finish_array(bobj);

    bson_append_start_array(bobj, "memory");
    for (uint i = 0; i < o.mem().size(); ++i) {
        bson_append_long(bobj, std::to_string(i).c_str(), o.mem()[i]);
    }
    bson_append_finish_array(bobj);

}
