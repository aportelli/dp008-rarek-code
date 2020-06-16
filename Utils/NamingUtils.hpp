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