//
// Created by dfunke on 4/7/15.
//

#pragma once

#include <exception>
#include <string>

#include <boost/noncopyable.hpp>

#include <ejdb/ejdb.h>
#include <ejdb/bson.h>

#include "Logger.h"

class EJDBException : public std::exception {

public:
    EJDBException() { };

    EJDBException(const std::string &what)
            : _what(what) { }

    const char *what() const noexcept { return _what.c_str(); }

private:
    std::string _what;
};

class DBConnection : private boost::noncopyable {

public:
    DBConnection(const std::string &file, const std::string &collection,
                 const int mode = JBOWRITER | JBOCREAT);

    ~DBConnection();

    template<class T>
    void save(const T &o) {
        bson bobj;
        bson_init(&bobj);

        //BSONify object
        bsonify(&bobj, o);

        bson_finish(&bobj);

        if (IS_VERBOSE)
            bson_print_raw(bobj.data, 2);

        bson_oid_t oid;
        if (!ejdbsavebson(m_coll, &bobj, &oid))
            _handleError();

        bson_destroy(&bobj);
    }

    template<typename T>
    T getMaximum(const std::string &field) const {
        std::string max = _getMaximum(field);
        T o = T();

        if (max != "") {
            std::istringstream ss(max);
            ss >> o;
        }

        return o;
    }

private:
    inline void _handleError() const {
        int ec = ejdbecode(m_ejdb);
        if (ec)
            throw EJDBException(ejdberrmsg(ec));
    }

    inline void _handleError(const std::string &msg) const {
        throw EJDBException(msg);
    }

    std::string _getMaximum(const std::string &field) const;

    template<class T>
    void bsonify(bson *bobj, const T &o);

private:
    EJDB *m_ejdb;
    EJCOLL *m_coll;

};

