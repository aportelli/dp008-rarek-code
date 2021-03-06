/*
 * ApplicationUtils.hpp, part of RareK (https://github.com/aportelli/dp008-rarek-code)
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

#ifndef Production_ApplicationUtils_hpp_
#define Production_ApplicationUtils_hpp_

#include "EntryUtils.hpp"
#include "NamingUtils.hpp"
#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>
#include <Hadrons/VirtualMachine.hpp>

using namespace Grid;
using namespace Hadrons;

void makeScalarPointSink(Application &application, const std::string mom,
                         const std::string name)
{
    MSink::ScalarPoint::Par scalarPointSinkPar;
    scalarPointSinkPar.mom = mom;
    application.createModule<MSink::ScalarPoint>(name, scalarPointSinkPar);
}

void makePointSink(Application &application, const std::string mom,
                         const std::string name)
{
    MSink::Point::Par pointSinkPar;
    pointSinkPar.mom = mom;
    application.createModule<MSink::Point>(name, pointSinkPar);
}

void makeZGaugeProp(Application &application, const std::string solver,
                    const std::string source, const std::string name)
{
    if (!(VirtualMachine::getInstance().hasModule(name)))
    {
      MFermion::ZGaugeProp::Par zPropPar;
      zPropPar.solver = solver;
      zPropPar.source = source;
      application.createModule<MFermion::ZGaugeProp>(name, zPropPar);
    }
}

void makeGaugeProp(Application &application, const std::string solver,
                   const std::string source, const std::string name)
{
    if (!(VirtualMachine::getInstance().hasModule(name)))
    {
      MFermion::GaugeProp::Par propPar;
      propPar.solver = solver;
      propPar.source = source;
      application.createModule<MFermion::GaugeProp>(name, propPar);
    }
}

void makeGaugeProps(Application &application, const std::string solver,
                    const std::vector<std::string> sources, const std::string name,
                    std::vector<std::string> &props, bool zProp)
{
    std::string propName, srcName;
    for (unsigned int i = 0; i < props.size(); ++i)
    {
      propName = name + "_" + std::to_string(i);
      srcName  = sources[i];
      if(zProp)
      {	      
        makeZGaugeProp(application, solver, srcName, propName);
      }
      else
      {
        makeGaugeProp(application, solver, srcName, propName);
      }
      props[i] = propName;
    }
}

void makeWallSource(Application &application, const std::string mom,
                    const unsigned int tW, const std::string name)
{
    if (!(VirtualMachine::getInstance().hasModule(name)))
    {
        MSource::Wall::Par wallPar;
        wallPar.tW = tW;
        wallPar.mom = mom;
        application.createModule<MSource::Wall>(name, wallPar);
    }
}

void makeWallProp(Application &application, const std::string solver,
                  const std::string mom, const unsigned int tW,
                  const std::string propName, bool zProp)
{
    std::string srcName = makeWallSourceName(tW, mom);
    makeWallSource(application, mom, tW, srcName);
    if(zProp)
    {
      makeZGaugeProp(application, solver, srcName, propName);
    }
    else
    {
      makeGaugeProp(application, solver, srcName, propName);
    }
}

void makeSmearedProp(Application &application, const std::string prop,
                     const std::string sink, const std::string name)
{
    if (!(VirtualMachine::getInstance().hasModule(name)))
    {
      MSink::Smear::Par smearPropPar;
      smearPropPar.q = prop;
      smearPropPar.sink = sink;
      application.createModule<MSink::Smear>(name, smearPropPar);
    }
}

void makeZSequentialSource(Application &application, const std::string q,
                           const std::string source,
                           const std::string action, const unsigned int tJ,
                           const std::string mom, const std::string name,
                           const unsigned int mu = 3, const Current curr_type = Current::Vector)
{
    if (!(VirtualMachine::getInstance().hasModule(name)))
    {
      std::string q5d = q + "_5d";
      MSource::ZSeqConserved::Par zSeqSrcPar;
      zSeqSrcPar.q = q5d;
      zSeqSrcPar.source = source;
      zSeqSrcPar.action = action;
      zSeqSrcPar.tA = tJ;
      zSeqSrcPar.tB = tJ;
      zSeqSrcPar.curr_type = curr_type;
      zSeqSrcPar.mu_min = mu;
      zSeqSrcPar.mu_max = mu;
      zSeqSrcPar.mom = mom;
      zSeqSrcPar.photon = "";
      application.createModule<MSource::ZSeqConserved>(name, zSeqSrcPar);
    }
}

void makeSequentialSource(Application &application, const std::string q,
                           const std::string source,
                           const std::string action, const unsigned int tJ,
                           const std::string mom, const std::string name,
                           const unsigned int mu = 3, const Current curr_type = Current::Vector)
{
    if (!(VirtualMachine::getInstance().hasModule(name)))
    {
      std::string q5d = q + "_5d";
      MSource::SeqConserved::Par seqSrcPar;
      seqSrcPar.q = q5d;
      seqSrcPar.source = source;
      seqSrcPar.action = action;
      seqSrcPar.tA = tJ;
      seqSrcPar.tB = tJ;
      seqSrcPar.curr_type = curr_type;
      seqSrcPar.mu_min = mu;
      seqSrcPar.mu_max = mu;
      seqSrcPar.mom = mom;
      seqSrcPar.photon = "";
      application.createModule<MSource::SeqConserved>(name, seqSrcPar);
    }
}

void makeSeqProp(Application &application, const std::string solver,
                 const std::string action, const unsigned int tJ,
                 const std::string mom, const std::string propName,
                 const std::string source, const std::string seqPropName,
                 bool zProp)
{
    std::string srcName = makeSeqSourceName(tJ, mom, seqPropName);
    if(zProp)
    {
      makeZSequentialSource(application, propName, source, action, tJ, mom, srcName);
      makeZGaugeProp(application, solver, srcName, seqPropName);
    }
    else
    {
      makeSequentialSource(application, propName, source, action, tJ, mom, srcName);
      makeGaugeProp(application, solver, srcName, seqPropName);
    }
}

void makeSeqProps(Application &application, const std::string solver,
                   const std::string action, const unsigned int tJ,
                   const std::string mom, const std::vector<std::string> props,
                   const std::vector<std::string> sources, const std::string name,
                   std::vector<std::string> &seqProps, bool zProp)
{
    std::string propName, srcName, seqPropName;
    for (unsigned int i = 0; i < props.size(); ++i)
    {
      seqPropName = name + "_" + std::to_string(i);
      propName = props[i];
      srcName  = sources[i];
      makeSeqProp(application, solver, action, tJ, mom, propName,
                   srcName, seqPropName, zProp);
      seqProps[i] = seqPropName;
    }
}

void makeZSparseSpinColorDiagonal(Application &application, const unsigned int nsrc,
                                  const unsigned int nsparse, const std::string name)
{
    if (!(VirtualMachine::getInstance().hasModule(name)))
    {
        MNoise::SparseSpinColorDiagonal::Par sparsePar;
        sparsePar.nsrc = nsrc;
        sparsePar.nsparse = nsparse;
        application.createModule<MNoise::SparseSpinColorDiagonal>(name, sparsePar);
    }
}

void makeZ2Diluted(Application &application, const std::string noise,
                   const std::string name)
{
    if (!(VirtualMachine::getInstance().hasModule(name)))
    {
        MSource::Z2Diluted::Par dilutedPar;
        dilutedPar.noise = noise;
        application.createModule<MSource::Z2Diluted>(name, dilutedPar);
    }
}

void makeZ2SparseSources(Application &application, const unsigned int nsrc,
                         const unsigned int nsparse, const std::string name)
{
    std::string scdName = name + "_spinColordiagonal";
    makeZSparseSpinColorDiagonal(application, nsrc, nsparse, scdName);
    makeZ2Diluted(application, scdName, name);
}

void makePropagatorVectorUnpack(Application &application, const std::string input,
                                const unsigned int size, const std::string name)
{
    if (!(VirtualMachine::getInstance().hasModule(name)))
    {
        MUtilities::PropagatorVectorUnpack::Par unpackPar;
        unpackPar.input = input;
        unpackPar.size = size;
        application.createModule<MUtilities::PropagatorVectorUnpack>(name, unpackPar);
    }
}

void unpackProps(Application &application, const std::string input,
                 const unsigned int size, const std::string name,
                 std::vector<std::string> &unpackedProps)
{
    makePropagatorVectorUnpack(application, input, size, name);
    for (unsigned int i = 0; i < unpackedProps.size(); ++i)
    {
      unpackedProps[i] = name + "_"  + std::to_string(i);
    }
}

void makeLoop(Application &application, const std::string q,
              const std::string eta, const std::string name)
{
    MContraction::Loop::Par loopPar;
    loopPar.q = q;
    loopPar.eta = eta;
    application.createModule<MContraction::Loop>(name, loopPar);
}

void makeLoops(Application &application, const std::vector<std::string> props,
               const std::vector<std::string> srcs, std::vector<std::string> &loops)
{
    std::string loopName;
    for (unsigned int i = 0; i < loops.size(); ++i)
    {
      loopName = props[i] + "_loop";
      makeLoop(application, props[i], srcs[i], loopName);
      loops[i] = loopName;
    }
}

template <typename E>
void makeMeson(Application &application, const std::string prop1,
               const std::string prop2, const std::string sink,
               const std::string gammas, const std::string output,
               const std::string name, const E entry,
               const std::string entryName)
{
    MContraction::Meson::Par mesonPar;
    mesonPar.q1 = prop1;
    mesonPar.q2 = prop2;
    mesonPar.sink = sink;
    mesonPar.gammas = gammas;
    mesonPar.output = output;
    application.createModule<MContraction::Meson>(name, mesonPar);
    application.setResultMetadata(name, entryName, entry);
}

template <typename E>
void makeGamma3pt(Application &application, const std::string prop1,
                  const std::string prop2, const std::string prop3,
                  const unsigned int tSnk, const std::string gamma,
                  const std::string output, const std::string name,
                  const E entry, const std::string entryName)
{
    MContraction::Gamma3pt::Par gamma3ptPar;
    gamma3ptPar.q1 = prop1;
    gamma3ptPar.q2 = prop2;
    gamma3ptPar.q3 = prop3;
    gamma3ptPar.tSnk = tSnk;
    gamma3ptPar.gamma = gamma;
    gamma3ptPar.output = output;
    application.createModule<MContraction::Gamma3pt>(name, gamma3ptPar);
    application.setResultMetadata(name, entryName, entry);
}

template <typename E>
void makeWeakNonEye(Application &application, const std::string qLeft,
                    const std::string qBarLeft, const std::string qRight,
                    const std::string qBarRight, const Gamma::Algebra gammaIn,
                    const Gamma::Algebra gammaOut, const std::string output,
                    const std::string name, const E entry,
                    const std::string entryName)
{
    MContraction::WeakNonEye3pt::Par weakNonEye3ptPar;
    weakNonEye3ptPar.qLeft = qLeft;
    weakNonEye3ptPar.qBarLeft = qBarLeft;
    weakNonEye3ptPar.qRight = qRight;
    weakNonEye3ptPar.qBarRight = qBarRight;
    weakNonEye3ptPar.gammaIn  = gammaIn;
    weakNonEye3ptPar.gammaOut = gammaOut;
    weakNonEye3ptPar.output = output;
    application.createModule<MContraction::WeakNonEye3pt>(name, weakNonEye3ptPar);
    application.setResultMetadata(name, entryName, entry);
}

template <typename E>
void makeWeakEye(Application &application, const std::string qBarLeft,
                 const std::string qBarRight, const std::string qSpectator,
                 const std::string loop, const Gamma::Algebra gammaIn,
                 const Gamma::Algebra gammaOut, const unsigned int tOut,
                 const std::string output, const std::string name,
                 const E entry, const std::string entryName)
{
    MContraction::WeakEye3pt::Par weakEye3ptPar;
    weakEye3ptPar.qBarLeft = qBarLeft;
    weakEye3ptPar.qBarRight = qBarRight;
    weakEye3ptPar.qSpectator = qSpectator;
    weakEye3ptPar.loop = loop;
    weakEye3ptPar.gammaIn  = gammaIn;
    weakEye3ptPar.gammaOut = gammaOut;
    weakEye3ptPar.tOut = tOut;
    weakEye3ptPar.output = output;
    application.createModule<MContraction::WeakEye3pt>(name, weakEye3ptPar);
    application.setResultMetadata(name, entryName, entry);
}

template <typename E>
void makeDiscLoop(Application &application, const std::string q_loop,
                  const std::vector<std::string> mom, const std::string output,
                  const std::string name, const E entry,
                  const std::string entryName, std::string gammas = "all")
{
    MContraction::DiscLoop::Par discLoopPar;
    discLoopPar.q_loop = q_loop;
    discLoopPar.gammas = gammas;
    discLoopPar.mom    = mom;
    discLoopPar.output = output;
    application.createModule<MContraction::DiscLoop>(name, discLoopPar);
    application.setResultMetadata(name, entryName, entry);
}

template <typename E>
void makeRareKaonNeutralDisc(Application &application, const std::string q1,
                             const std::string q2, const std::string q3,
                             const std::string q4,
                             const std::string output, const std::string name,
                             const E entry, const std::string entryName)
{
    MContraction::RareKaonNeutralDisc::Par RareKaonNeutralDiscPar;
    RareKaonNeutralDiscPar.q1 = q1;
    RareKaonNeutralDiscPar.q2 = q2;
    RareKaonNeutralDiscPar.q3 = q3;
    RareKaonNeutralDiscPar.q4 = q4;
    RareKaonNeutralDiscPar.output = output;
    application.createModule<MContraction::RareKaonNeutralDisc>(name, RareKaonNeutralDiscPar);
    application.setResultMetadata(name, entryName, entry);
}

#endif // Production_ApplicationUtils_hpp_
