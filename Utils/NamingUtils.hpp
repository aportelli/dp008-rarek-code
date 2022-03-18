/*
 * NamingUtils.hpp, part of RareK (https://github.com/aportelli/dp008-rarek-code)
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

#ifndef Production_NamingUtils_hpp_
#define Production_NamingUtils_hpp_

#include <Grid/Grid.h>
#include <Hadrons/Modules.hpp>
using namespace Grid;
using namespace Hadrons;

std::string sanitizeMom(std::string mom)
{
    std::string sMom;
    std::vector<double> dMom = strToVec<double>(mom);
    int len = dMom.size();
    for (int i = 0; i < len; ++i)
    {
      int iMom = static_cast<int>(dMom[i]);
      sMom += std::to_string(iMom);
    }
    return sMom;
}

std::vector<std::vector<std::string>> makeDiscMom(std::string kMom, std::string pMom)
{
    std::vector<std::vector<std::string>> discMom;
    std::vector<std::string> discKMom, discPMom;
    std::string sKMom, sPMom;
    std::vector<double> dKMom = strToVec<double>(kMom);
    std::vector<double> dPMom = strToVec<double>(pMom);
    int len = dKMom.size();
    for (int i = 0; i < len; ++i)
    {
      int iKMom = static_cast<int>(dKMom[i]);
      sKMom += std::to_string(iKMom) + " ";
      int iPMom = static_cast<int>(dPMom[i]);
      sPMom += std::to_string(iPMom) + " ";
    }
    discKMom.push_back(sKMom);    discPMom.push_back(sPMom);
    discMom.push_back(discKMom);  discMom.push_back(discPMom);
    return discMom;
}

std::string makeWallSourceName(const unsigned int tW, const std::string mom)
{
    std::string wallSourceName = "wall_" + std::to_string(tW) + "_" + mom;
    return  wallSourceName;
}

std::string makeZ2SparseSourceName(const std::string src, const unsigned int i)
{
    std::string z2SparseSourceName = src + "_" + std::to_string(i);
    return  z2SparseSourceName;
}

std::string makeSeqSourceName(const unsigned int tJ, const std::string mom,
                              const std::string propName, const Gamma::Algebra gamma = Gamma::Algebra::GammaT)
{
    std::string Gmu = std::to_string(gamma);
    std::string seqSourceName = "G" + Gmu +
                                 "_tJ_" + std::to_string(tJ) + "_" + mom +
                                 "_" + propName;
    return  seqSourceName;
}

#endif // Production_NamingUtils_hpp_
