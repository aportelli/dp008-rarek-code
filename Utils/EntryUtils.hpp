/*
 * EntryUtils.hpp, part of RareK (https://github.com/aportelli/dp008-rarek-code)
 *
 * Copyright (C) 2015 - 2022
 *
 * Authors: Felix Erben, Fionn Ó hÓgáin, and Antonin Portelli
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution 
 * directory.
 */

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
    HADRONS_SQL_FIELDS(SqlNotNull<int>,                 tSrc,
                       SqlNotNull<std::string>,         diag,
                       SqlNotNull<std::string>,         sink,
                       SqlNotNull<std::vector<double>>, mom);
};

Entry2pt makeEntry2pt(int         tSrc,
                      std::string diag,
                      std::string sink,
                      std::vector<double> mom)
{
    Entry2pt entry2pt;
    entry2pt.tSrc = tSrc;
    entry2pt.diag = diag;
    entry2pt.sink = sink;
    entry2pt.mom  = mom;

    return entry2pt; 
}

struct Entry3ptVC: public SqlEntry
{
    HADRONS_SQL_FIELDS(SqlNotNull<std::vector<double>>, k,
                       SqlNotNull<std::vector<double>>, p,
                       SqlNotNull<std::vector<double>>, q,
                       SqlNotNull<int>,                 tSrc,
                       SqlNotNull<int>,                 tJ,
                       SqlNotNull<std::string>,         diag,
                       SqlNotNull<std::string>,         sink);
};

Entry3ptVC makeEntry3ptVC(std::vector<double> k,
                          std::vector<double> p,
                          std::vector<double> q,
                          int                 tSrc,
                          int                 tJ,
                          std::string         diag,
                          std::string         sink)
{
    Entry3ptVC entry3ptVC;
    entry3ptVC.k    = k;
    entry3ptVC.p    = p;
    entry3ptVC.q    = q;
    entry3ptVC.tSrc = tSrc;
    entry3ptVC.tJ   = tJ;
    entry3ptVC.diag = diag;
    entry3ptVC.sink = sink;

    return entry3ptVC; 
}

struct Entry3ptHw: public SqlEntry
{
    HADRONS_SQL_FIELDS(SqlNotNull<int>,                 tK,
                       SqlNotNull<int>,                 tPi,
                       SqlNotNull<std::string>,         diag,
                       SqlNotNull<std::vector<double>>, mom,
                       std::string,                     gim,
                       int,                             hit,
                       int,                             noise);
};

Entry3ptHw makeEntry3ptHw(int                 tK,
                          int                 tPi,
                          std::string         diag,
                          std::vector<double> mom,
                          std::string         gim,
                          int                 hit,
                          int                 noise,
                          bool                nullHit=false,
                          bool                nullNoise=false)
{
    Entry3ptHw entry3pt;
    entry3pt.tK    = tK;
    entry3pt.tPi   = tPi;
    entry3pt.diag  = diag;
    entry3pt.mom   = mom;
    entry3pt.gim   = gim;
    entry3pt.hit   = hit;
    entry3pt.noise = noise;
    entry3pt.nullify.hit = nullHit;
    entry3pt.nullify.noise = nullNoise;

    return entry3pt; 
}

struct Entry4pt: public SqlEntry
{
    HADRONS_SQL_FIELDS(SqlNotNull<std::string>, diag,
                       std::string,             gim,
                       int,                     hit,
                       int,                     noise);
};

MergedSqlEntry<BaseEntry, Entry4pt> makeEntry4pt(BaseEntry   baseEntry,
                                                 std::string diag,
                                                 std::string gim,
                                                 int         hit,
                                                 int         noise,
                                                 bool        nullHit=false,
                                                 bool        nullNoise=false)
{
    Entry4pt entry4pt;
    entry4pt.diag  = diag;
    entry4pt.gim   = gim;
    entry4pt.hit   = hit;
    entry4pt.noise = noise;
    entry4pt.nullify.hit = nullHit;
    entry4pt.nullify.noise = nullNoise;

    auto mergedEntry(mergeSqlEntries(baseEntry, entry4pt));

    return mergedEntry; 
}

struct DiscLoopEntry: public SqlEntry
{
    HADRONS_SQL_FIELDS(SqlNotNull<std::vector<double>>, mom,
                       SqlNotNull<std::string>,         flavor,
                       int,                             hit,
                       int,                             noise);
};

DiscLoopEntry makeDiscLoopEntry(std::vector<double> mom,
                                std::string         flavor,
                                int                 hit,
                                int                 noise,
                                bool                nullHit=false)
{
    DiscLoopEntry entryDL;
    entryDL.mom    = mom;
    entryDL.flavor = flavor;
    entryDL.hit    = hit;
    entryDL.noise  = noise;
    entryDL.nullify.hit = nullHit;

    return entryDL; 
}

#endif // Production_EntryUtils_hpp_
