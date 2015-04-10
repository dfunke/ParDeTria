//
// Created by dfunke on 4/7/15.
//

#include "DBConnection.h"

DBConnection::DBConnection(const std::string & file, const std::string & collection,
             const int mode) {
    m_ejdb = ejdbnew();

    if(!m_ejdb){
        _handleError("EJDB Error: Could not create EJDB instance");
    }

    if(!ejdbopen(m_ejdb, file.c_str(), mode))
        _handleError();

    m_coll = ejdbgetcoll(m_ejdb, collection.c_str());
    if(!m_coll && mode & JBOWRITER)
        m_coll = ejdbcreatecoll(m_ejdb, collection.c_str(), nullptr);
    if(!m_coll)
        _handleError();

}

DBConnection::~DBConnection() {
    if(m_ejdb) {
        if (ejdbisopen(m_ejdb))
            ejdbclose(m_ejdb);

        ejdbdel(m_ejdb);
    }
}

#include "Timings.h"

template <>
void DBConnection::bsonify(bson *bobj, const ExperimentRun & o) {

    for(const auto & t : o.traits()){
        bson_append_string(bobj, t.first.c_str(), t.second.c_str());
    }

    bson_append_start_array(bobj, "times");
    for(uint i = 0; i < o.times().size(); ++i){
        bson_append_long(bobj, std::to_string(i).c_str() ,o.times()[i].count());
    }
    bson_append_finish_array(bobj);

    bson_append_start_array(bobj, "memory");
    for(uint i = 0; i < o.mem().size(); ++i){
        bson_append_long(bobj, std::to_string(i).c_str() ,o.mem()[i]);
    }
    bson_append_finish_array(bobj);

}
