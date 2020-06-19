#ifndef Production_EntryUtils_hpp_
#define Production_EntryUtils_hpp_

#include <Hadrons/SqlEntry.hpp>

using namespace Grid;
using namespace Hadrons;

struct BaseEntry: public SqlEntry
{
    HADRONS_SQL_FIELDS(SqlNotNull<std::vector<double>>, k,
                       SqlNotNull<std::vector<double>>, p,
                       SqlNotNull<std::vector<double>>, q,
                       SqlNotNull<int>,                 tK,
                       SqlNotNull<int>,                 tJ,
                       SqlNotNull<int>,                 tPi);
};

BaseEntry makeBaseEntry(std::vector<double> k,
                        std::vector<double> p,
                        std::vector<double> q,
                        int tK, int tJ, int tPi)
{
    BaseEntry baseEntry;
    baseEntry.k   = k;
    baseEntry.p   = p;
    baseEntry.q   = q;
    baseEntry.tK  = tK;
    baseEntry.tJ  = tJ;
    baseEntry.tPi = tPi;
    return baseEntry;
}

struct Entry2pt: public SqlEntry
{
    HADRONS_SQL_FIELDS(SqlNotNull<std::string>, diag,
                       SqlNotNull<std::string>, sink,
                       SqlNotNull<std::string>, mom);
};

MergedSqlEntry<BaseEntry, Entry2pt> makeEntry2pt(BaseEntry   baseEntry,
                                                 std::string diag,
                                                 std::string sink,
                                                 std::string mom)
{
    Entry2pt entry2pt;
    entry2pt.diag = diag;
    entry2pt.sink = sink;
    entry2pt.mom  = mom;

    auto mergedEntry(mergeSqlEntries(baseEntry, entry2pt));

    return mergedEntry; 
}

struct Entry3ptMeson: public SqlEntry
{
    HADRONS_SQL_FIELDS(SqlNotNull<std::string>, diag,
                       SqlNotNull<std::string>, sink);
};

MergedSqlEntry<BaseEntry, Entry3ptMeson> makeEntry3ptMeson(BaseEntry   baseEntry,
                                                           std::string diag,
                                                           std::string sink)
{
    Entry3ptMeson entry3pt;
    entry3pt.diag  = diag;
    entry3pt.sink  = sink;

    auto mergedEntry(mergeSqlEntries(baseEntry, entry3pt));

    return mergedEntry; 
}

struct Entry3pt: public SqlEntry
{
    HADRONS_SQL_FIELDS(SqlNotNull<std::string>, diag,
                       SqlNotNull<std::string>, mom);
};

MergedSqlEntry<BaseEntry, Entry3pt> makeEntry3pt(BaseEntry   baseEntry,
                                                 std::string diag,
                                                 std::string mom)
{
    Entry3pt entry3pt;
    entry3pt.diag  = diag;
    entry3pt.mom   = mom;

    auto mergedEntry(mergeSqlEntries(baseEntry, entry3pt));

    return mergedEntry; 
}

struct Entry3ptHwE: public SqlEntry
{
    HADRONS_SQL_FIELDS(SqlNotNull<std::string>, diag,
                       SqlNotNull<std::string>, sink,
                       SqlNotNull<std::string>, mom,
                       SqlNotNull<std::string>, gim,
                       SqlNotNull<int>,         noise);
};

MergedSqlEntry<BaseEntry, Entry3ptHwE> makeEntry3ptHwE(BaseEntry   baseEntry,
                                                       std::string diag,
                                                       std::string sink,
                                                       std::string mom,
                                                       std::string gim,
                                                       int         noise)
{
    Entry3ptHwE entry3pt;
    entry3pt.diag  = diag;
    entry3pt.sink  = sink;
    entry3pt.mom   = mom;
    entry3pt.gim   = gim;
    entry3pt.noise = noise;

    auto mergedEntry(mergeSqlEntries(baseEntry, entry3pt));

    return mergedEntry; 
}

struct Entry4pt: public SqlEntry
{
    HADRONS_SQL_FIELDS(SqlNotNull<std::string>, diag);
};

MergedSqlEntry<BaseEntry, Entry4pt> makeEntry4pt(BaseEntry   baseEntry,
                                                 std::string diag)
{
    Entry4pt entry4pt;
    entry4pt.diag = diag;

    auto mergedEntry(mergeSqlEntries(baseEntry, entry4pt));

    return mergedEntry; 
}

struct Entry4ptGim: public SqlEntry
{
    HADRONS_SQL_FIELDS(SqlNotNull<std::string>, diag,
                       SqlNotNull<std::string>, gim,
                       SqlNotNull<int>,         noise);
};

MergedSqlEntry<BaseEntry, Entry4ptGim> makeEntry4ptGim(BaseEntry   baseEntry,
                                                       std::string diag,
                                                       std::string gim,
                                                       int         noise)
{
    Entry4ptGim entry4pt;
    entry4pt.diag  = diag;
    entry4pt.gim   = gim;
    entry4pt.noise = noise;

    auto mergedEntry(mergeSqlEntries(baseEntry, entry4pt));

    return mergedEntry; 
}

struct DiscLoopEntry: public SqlEntry
{
    HADRONS_SQL_FIELDS(SqlNotNull<std::string>, diag,
                       SqlNotNull<std::string>, mom,
                       SqlNotNull<std::string>, flavor,
                       SqlNotNull<int>,         noise);
};

MergedSqlEntry<BaseEntry, DiscLoopEntry> makeDiscLoopEntry(BaseEntry   baseEntry,
                                                           std::string diag,
                                                           std::string mom,
                                                           std::string flavor,
                                                           int         noise)
{
    DiscLoopEntry entryDL;
    entryDL.diag   = diag;
    entryDL.mom    = mom;
    entryDL.flavor = flavor;
    entryDL.noise  = noise;

    auto mergedEntry(mergeSqlEntries(baseEntry, entryDL));

    return mergedEntry; 
}

#endif // Production_EntryUtils_hpp_