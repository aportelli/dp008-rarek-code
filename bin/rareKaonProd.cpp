#include "Utils.hpp"
#include <Grid/Grid.h>
#include <Hadrons/A2AVectors.hpp>
#include <Hadrons/Modules/MContraction/A2AMesonField.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MContraction;

namespace RareKaonInputs
{
    class TrajRange : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(TrajRange,
                                        unsigned int, start,
                                        unsigned int, end,
                                        unsigned int, step);
    };

    class TimePar : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(TimePar,
                                        unsigned int, dt,
                                        unsigned int, dtK,
                                        unsigned int, dtJ,
                                        unsigned int, dtP);
    };

    class MomPar : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(MomPar,
                                        std::string, kmom,
                                        std::string, pmom,
                                        std::string, qmom,
                                        std::string, mqmom,
                                        std::string, sinkkmom,
                                        std::string, sinkpmom);
    };

    class GammaPar : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(GammaPar,
                                        std::string,    twoPtGammas,
                                        Gamma::Algebra, localGamma,
                                        Gamma::Algebra, threePtGammaIn,
                                        Gamma::Algebra, threePtGammaOut,
                                        std::string   , threePtGammaInsertions);
    };

    class IOPar : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(IOPar,
                                        std::string,  gaugeFile,
                                        std::string,  epackFile,
                                        unsigned int, epackSize,
                                        unsigned int, epackLs,
                                        bool,         epackMultiFile,
                                        std::string,  scheduleFile,
                                        std::string,  resultPStem,
                                        std::string,  xmlFileName);
    };

    class ZMobiusPar : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(ZMobiusPar,
                                        unsigned int, Ls,
                                        double,       M5,
                                        double,       b,
                                        double,       c,
                                        std::vector<ComplexD>, omega);
    };

    class LightActionPar : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(LightActionPar,
                                        double, mass,
                                        double, residual);
    };

    class CharmActionPar : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(CharmActionPar,
                                        double, mass,
                                        double, residual);
    };

    class StrangeActionPar : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(StrangeActionPar,
                                        double,       mass,
                                        unsigned int, Ls,
                                        double,       M5,
                                        double,       scale,
                                        double,       residual);
    };
} // namespace RareKaonInputs

struct RareKaonPar
{
    RareKaonInputs::StrangeActionPar strangeActionPar;
    RareKaonInputs::CharmActionPar   charm1ActionPar;
    RareKaonInputs::CharmActionPar   charm2ActionPar;
    RareKaonInputs::CharmActionPar   charm3ActionPar;
    RareKaonInputs::LightActionPar   lightActionPar;
    RareKaonInputs::ZMobiusPar       zMobiusPar;
    RareKaonInputs::TrajRange        trajRange;
    RareKaonInputs::GammaPar         gammaPar;
    RareKaonInputs::TimePar          timePar;
    RareKaonInputs::MomPar           momPar;
    RareKaonInputs::IOPar            ioPar;
};
int main(int argc, char *argv[])
{
    // parse command line
    std::string parFilename;

    if (argc < 2)
    {
        std::cerr << "usage: " << argv[0] << " <parameter file>";
        std::cerr << std::endl;

        return EXIT_FAILURE;
    }
    parFilename = argv[1];

    // parse parameter file
    RareKaonPar par;
    XmlReader reader(parFilename);

    read(reader, "strangeActionPar", par.strangeActionPar);
    read(reader,  "charm1ActionPar",  par.charm1ActionPar);
    read(reader,  "charm2ActionPar",  par.charm2ActionPar);
    read(reader,  "charm3ActionPar",  par.charm3ActionPar);
    read(reader,   "lightActionPar",   par.lightActionPar);
    read(reader,       "zMobiusPar",       par.zMobiusPar);
    read(reader,        "trajRange",        par.trajRange);
    read(reader,         "gammaPar",         par.gammaPar);
    read(reader,          "timePar",          par.timePar);
    read(reader,           "momPar",           par.momPar);
    read(reader,            "ioPar",            par.ioPar);

    unsigned int trajStart = par.trajRange.start;
    unsigned int trajEnd   = par.trajRange.end;
    unsigned int trajStep  = par.trajRange.step;

    unsigned int dt  = par.timePar.dt;
    unsigned int dtK = par.timePar.dtK;
    unsigned int dtJ = par.timePar.dtJ;
    unsigned int dtP = par.timePar.dtP;

    std::string kmom  = par.momPar.kmom;
    std::string pmom  = par.momPar.pmom;
    std::string qmom  = par.momPar.qmom;
    std::string mqmom = par.momPar.mqmom;
    std::string sinkkmom = par.momPar.sinkkmom;
    std::string sinkpmom = par.momPar.sinkpmom;

    // initialization //////////////////////////////////////////////////////////
    Grid_init(&argc, &argv);
    HadronsLogError.Active(GridLogError.isActive());
    HadronsLogWarning.Active(GridLogWarning.isActive());
    HadronsLogMessage.Active(GridLogMessage.isActive());
    HadronsLogIterative.Active(GridLogIterative.isActive());
    HadronsLogDebug.Active(GridLogDebug.isActive());
    LOG(Message) << "Grid initialized" << std::endl;

    // global parameters
    Application application;
    Application::GlobalPar globalPar;
    globalPar.runId = "Rare_K";
    globalPar.trajCounter.start = trajStart;
    globalPar.trajCounter.end = trajEnd;
    globalPar.trajCounter.step = trajStep;
    application.setPar(globalPar);

    // action parameters ///////////////////////////////////////////////////////////////
    std::vector<std::string> flavour = {"s", "l", "c1", "c2", "c3"};
    std::vector<std::string> epack = {"", "epack_l", "", "", ""};
    std::vector<double> mass = {par.strangeActionPar.mass,
                                par.lightActionPar.mass,
                                par.charm1ActionPar.mass,
                                par.charm2ActionPar.mass,
                                par.charm3ActionPar.mass};
    std::vector<double> residual = {par.strangeActionPar.residual,
                                    par.lightActionPar.residual,
                                    par.charm1ActionPar.residual,
                                    par.charm2ActionPar.residual,
                                    par.charm3ActionPar.residual};
    std::vector<unsigned int> Ls = {par.strangeActionPar.Ls,
                                    par.zMobiusPar.Ls};
    std::vector<double> M5 = {par.strangeActionPar.M5,
                              par.zMobiusPar.M5};
    std::string boundary = "1 1 1 -1";
    std::string twist = "0. 0. 0. 0.";
    unsigned int nt = GridDefaultLatt()[Tp];

    // gauge field
    if (!(par.ioPar.gaugeFile.empty()))
    {
        MIO::LoadNersc::Par loadPar;
        loadPar.file = par.ioPar.gaugeFile;
        application.createModule<MIO::LoadNersc>("gauge", loadPar);
    }
    else
    {
        application.createModule<MGauge::Unit>("gauge");
    }
    

    // gauge field cast
    MUtilities::GaugeSinglePrecisionCast::Par gaugefPar;
    gaugefPar.field = "gauge";
    application.createModule<MUtilities::GaugeSinglePrecisionCast>("gaugef", gaugefPar);

    // gauge fixing
    MGauge::GaugeFix::Par gaugeFixPar;
    gaugeFixPar.gauge = "gauge";
    gaugeFixPar.alpha = 0.05;
    gaugeFixPar.maxiter = 1000000;
    gaugeFixPar.Omega_tol = 1e-8;
    gaugeFixPar.Phi_tol = 1e-8;
    gaugeFixPar.gaugeFix = Grid::Hadrons::MGauge::Fix::coulomb;
    gaugeFixPar.Fourier = true;
    application.createModule<MGauge::GaugeFix>("gaugeFix", gaugeFixPar);

    // gauge field cast
    MUtilities::GaugeSinglePrecisionCast::Par gaugefFixPar;
    gaugefFixPar.field = "gaugeFix";
    application.createModule<MUtilities::GaugeSinglePrecisionCast>("gaugefFix", gaugefFixPar);
    
    // epack

    if (!(par.ioPar.epackFile.empty()))
    {
        MIO::LoadFermionEigenPackIo32::Par epackPar;
        epackPar.filestem = par.ioPar.epackFile;
        epackPar.multiFile = par.ioPar.epackMultiFile;
        epackPar.size = par.ioPar.epackSize;
        epackPar.Ls = par.ioPar.epackLs;
        epackPar.gaugeXform = "gaugeFix_xform";
        application.createModule<MIO::LoadFermionEigenPackIo32>("epack_l", epackPar);
    }
    else
    {
        epack[1] = "";
    }

    // DWF action: strange
    MAction::ScaledDWF::Par actionStrangePar;
    actionStrangePar.gauge = "gaugeFix";
    actionStrangePar.Ls = Ls[0];
    actionStrangePar.M5 = M5[0];
    actionStrangePar.mass = mass[0];
    actionStrangePar.boundary = boundary;
    actionStrangePar.scale = par.strangeActionPar.scale;
    actionStrangePar.twist = twist;
    application.createModule<MAction::ScaledDWF>("dwf_" + flavour[0], actionStrangePar);

    // actionF strange
    MAction::ScaledDWFF::Par actionFStrangePar;
    actionFStrangePar.gauge = "gaugefFix";
    actionFStrangePar.Ls = Ls[0];
    actionFStrangePar.M5 = M5[0];
    actionFStrangePar.mass = mass[0];
    actionFStrangePar.boundary = boundary;
    actionFStrangePar.scale = par.strangeActionPar.scale;
    actionFStrangePar.twist = twist;
    application.createModule<MAction::ScaledDWFF>("dwff_" + flavour[0], actionFStrangePar);

    // solver strange
    MSolver::MixedPrecisionRBPrecCG::Par solverParStrange;
    solverParStrange.innerAction = "dwff_" + flavour[0];
    solverParStrange.outerAction = "dwf_" + flavour[0];
    solverParStrange.residual = residual[0];
    solverParStrange.maxInnerIteration = 30000;
    solverParStrange.maxOuterIteration = 100;
    solverParStrange.eigenPack = epack[0];
    application.createModule<MSolver::MixedPrecisionRBPrecCG>("mcg_" + flavour[0], solverParStrange);


    // Zmobius action: light, charm1, charm2, charm3
    for (unsigned int i = 1; i < flavour.size(); ++i)
    {
        MAction::ZMobiusDWF::Par ZMobActionPar;
        ZMobActionPar.gauge = "gaugeFix";
        ZMobActionPar.Ls = Ls[1];
        ZMobActionPar.M5 = M5[1];
        ZMobActionPar.mass = mass[i];
        ZMobActionPar.boundary = boundary;
        ZMobActionPar.b = par.zMobiusPar.b;
        ZMobActionPar.c = par.zMobiusPar.c;
        ZMobActionPar.omega = par.zMobiusPar.omega;
        ZMobActionPar.twist = twist;
        application.createModule<MAction::ZMobiusDWF>("dwf_" + flavour[i], ZMobActionPar);

        // actionF light
        MAction::ZMobiusDWFF::Par ZMobFAction;
        ZMobFAction.gauge = "gaugefFix";
        ZMobFAction.Ls = Ls[1];
        ZMobFAction.M5 = M5[1];
        ZMobFAction.mass = mass[i];
        ZMobFAction.boundary = boundary;
        ZMobFAction.b = par.zMobiusPar.b;
        ZMobFAction.c = par.zMobiusPar.c;
        ZMobFAction.omega = par.zMobiusPar.omega;
        ZMobFAction.twist = twist;
        application.createModule<MAction::ZMobiusDWFF>("dwff_" + flavour[i], ZMobFAction);

        // solver light
        MSolver::ZMixedPrecisionRBPrecCG::Par ZMobSolverPar;
        ZMobSolverPar.innerAction = "dwff_" + flavour[i];
        ZMobSolverPar.outerAction = "dwf_" + flavour[i];
        ZMobSolverPar.residual = residual[i];
        ZMobSolverPar.maxInnerIteration = 30000;
        ZMobSolverPar.maxOuterIteration = 100;
        ZMobSolverPar.eigenPack = epack[i];
        application.createModule<MSolver::ZMixedPrecisionRBPrecCG>("mcg_" + flavour[i], ZMobSolverPar);

    }

    // Sparse Noise Sources
    std::string sparseNoise = "sparseNoise";
    std::string sparseNoises = "sparseNoises";
    std::string solver, sparseProp, sparseProps;
    unsigned int nsrc = 1, nsparse = 2, nds = pow(nsparse, 4);
    makeZ2SparseSources(application, nsrc, nsparse, sparseNoises);
    unpackProps(application, sparseNoises, sparseNoise);

    std::vector<std::vector<std::string>> sparseLoops, seqSparseLoops;//, seqSparseProps;//, sparseProps, 
    // sparseProps.resize(flavour.size());
    sparseLoops.resize(flavour.size(), std::vector<std::string>(nds));
    seqSparseLoops.resize(flavour.size(), std::vector<std::string>(nds));
    // seqSparseLoops.resize(flavour.size());
    for (unsigned int i = 1; i < flavour.size(); ++i)
    {
        solver = "mcg_" + flavour[i];
        sparseProps = "sparseProps_" + flavour[i];
        sparseProp = "sparseProp_" + flavour[i];

        makeZGaugeProp(application, solver, sparseNoises, sparseProps);
        unpackProps(application, sparseProps, sparseProp);
        makeLoops(application, sparseProp, sparseNoise, sparseLoops[i]);
    }

    // sinks
    std::string sink2ptZeromom = "2ptZeromomsink";
    makeScalarPointSink(application, sinkkmom, sink2ptZeromom);
    std::string sink2ptPmom = "2ptPmomsink";
    makeScalarPointSink(application, sinkpmom, sink2ptPmom);

    std::string sinkZeromom = "sinkZeromom";
    makePointSink(application, sinkkmom, sinkZeromom);
    std::string sinkPmom = "sinkPmom";
    makePointSink(application, sinkpmom, sinkPmom);


    std::string  twoPtGammas     = par.gammaPar.twoPtGammas;
    std::string  gammaInsertions = par.gammaPar.threePtGammaInsertions;
    Gamma::Algebra localGamma    = par.gammaPar.localGamma;
    Gamma::Algebra gammaIn       = par.gammaPar.threePtGammaIn;
    Gamma::Algebra gammaOut      = par.gammaPar.threePtGammaOut;

    unsigned int tk, tj, tp;
    std::string stk, stj, stp, timeStamp;
    std::string resultPStem = par.ioPar.resultPStem;
    std::string smu = std::to_string(localGamma);
    std::string strangeSolver = "mcg_" + flavour[0];
    std::string lightSolver   = "mcg_" + flavour[1];
    std::string strangeAction = "dwf_" + flavour[0];
    std::string lightAction   = "dwf_" + flavour[1];
    for(unsigned int t = 0; t < nt; t+=dt)
    {
        //////////////////////////////////////////////////
        // Time translations
        //////////////////////////////////////////////////
        tk = (t + dtK) % nt;
        tj = (t + dtJ) % nt;
        tp = (t + dtP) % nt;
        stk = std::to_string(t);
        stj = std::to_string(tj);
        stp = std::to_string(tp);
        timeStamp = stk + "_tJ_" + stj + "_tP_" + stp;

        //////////////////////////////////////////////////
        // Propagators
        //////////////////////////////////////////////////
        // Zero momentum
        std::string qWallksZeromom = "QWall_s_0mom_" + stk;
        makeWallProp(application, strangeSolver, kmom, tk, qWallksZeromom);
        std::string qWallklZeromom = "QWall_l_0mom_" + stk;
        makeWallZProp(application, lightSolver, kmom, tk, qWallklZeromom);
        std::string qWallplZeromom = "QWall_l_0mom_" + stp;

        // Smeared propagators: sink zero mom
        std::string smearedqWallksZeromom = "smearedQWall_s_0mom_" + stk;
        makeSmearedProp(application, qWallksZeromom, sinkZeromom, smearedqWallksZeromom);
        std::string smearedqWallklZeromom = "smearedQWall_l_0mom_" + stk;
        makeSmearedProp(application, qWallklZeromom, sinkZeromom, smearedqWallklZeromom);
        std::string smearedqWallplZeromom = "smearedQWall_l_0mom_" + stp;

        // Nonzero momentum
        std::string qWallksPmom = "QWall_s_Pmom_" + stk;
        makeWallProp(application, strangeSolver, pmom, tk, qWallksPmom);
        std::string qWallklPmom = "QWall_l_Pmom_" + stk;
        makeWallZProp(application, lightSolver, pmom, tk, qWallklPmom);
        std::string qWallplbarPmom = "QWall_l_Pmom_" + stp;
        std::string qWallpsPmom = "QWall_s_Pmom_" + stp;

        // Smeared propagators
        std::string smearedqWallksPmom = "smearedQWall_s_Pmom" + stk;
        makeSmearedProp(application, qWallksPmom, sinkZeromom, smearedqWallksPmom);
        std::string smearedqWallklPmom = "smearedQWall_l_Pmom" + stk;
        makeSmearedProp(application, qWallklPmom, sinkZeromom, smearedqWallklPmom);
        std::string smearedqWallplPmom = "smearedQWall_l_Pmom" + stp;

        //////////////////////////////////////////////////
        // Sequential propagators
        //////////////////////////////////////////////////
        std::string seqGmuKLQmom = "G" + smu + "_KL_" + qWallklZeromom;
        makeSeqZProp(application, lightSolver, lightAction, tj, qmom, qWallklZeromom, seqGmuKLQmom);
        std::string seqGmuKSMqmom = "G" + smu + "_KS_" + qWallksZeromom;
        makeSeqProp(application, strangeSolver, strangeAction, tj, mqmom, qWallksZeromom, seqGmuKSMqmom);
        std::string seqGmuPLMqmom = "G" + smu + "_PL_" + qWallplZeromom;
        makeSeqZProp(application, lightSolver, lightAction, tj, mqmom, qWallplZeromom, seqGmuPLMqmom);
        std::string seqGmuPLbarQmom = "G" + smu + "_PLbar_" + qWallplbarPmom;
        makeSeqZProp(application, lightSolver, lightAction, tj, qmom, qWallplbarPmom, seqGmuPLbarQmom);

        // Smeared sequential propagators
        std::string smearedqWallGmuKLQmom = "smearedQWall_Kl_VC" + smu + "_Qmom_" + stk;
        makeSmearedProp(application, seqGmuKLQmom, sinkPmom, smearedqWallGmuKLQmom);
        std::string smearedqWallGmuKspecQmom = "smearedQWall_kspec_VC" + smu + "_Qmom_" + stk;
        makeSmearedProp(application, seqGmuKLQmom, sinkZeromom, smearedqWallGmuKspecQmom);
        std::string smearedqWallGmuKSmQmom = "smearedQWall_Ks_VC" + smu + "_mQmom_" + stk;
        makeSmearedProp(application, seqGmuKSMqmom, sinkZeromom, smearedqWallGmuKSmQmom);
        std::string smearedqWallGmuPLmQmom = "smearedQWall_PL_VC" + smu + "_mQmom_" + stk;
        makeSmearedProp(application, seqGmuPLMqmom, sinkZeromom, smearedqWallGmuPLmQmom);
        std::string smearedqWallGmuPLbarQmom = "smearedQWall_PL_VC" + smu + "_Qmom_" + stk;
        makeSmearedProp(application, seqGmuPLbarQmom, sinkZeromom, smearedqWallGmuPLbarQmom);

        //////////////////////////////////////////////////
        // Contractions
        //////////////////////////////////////////////////

        //////////////////////////////////////////////////
        // Kaon Meson Contractions: Wall-point
        //////////////////////////////////////////////////
        // Kaon momentum 
        std::string wallPointMesonKZeromomRes = resultPStem + "/2pt/WallPoint/2ptK_mom_K_tK_" + timeStamp;
        std::string wallPointMesonKZeromom = "wallPointMesonKZeromom_" + stk;
        makeMeson(application, qWallklZeromom, qWallksZeromom, sink2ptZeromom,
                  twoPtGammas, wallPointMesonKZeromomRes, wallPointMesonKZeromom);
        // Pion momentum
        std::string wallPointMesonKPmomRes = resultPStem + "/2pt/WallPoint/2ptK_mom_P_tK_" + timeStamp;
        std::string wallPointMesonKPmom = "wallPointMesonKPmom_" + stk;
        makeMeson(application, qWallklZeromom, qWallksPmom, sink2ptPmom,
                  twoPtGammas, wallPointMesonKPmomRes, wallPointMesonKPmom);
        // Light insetion: K->pi momentum
        std::string wallPointMesonKGmuLRes = resultPStem + "/3pt/WallPoint/3ptK_VC" + smu + "L_tK_" + timeStamp;
        std::string wallPointMesonKGmuL = "wallPointMesonKGmuL_" + stk;
        makeMeson(application, seqGmuKLQmom, qWallksZeromom, sink2ptPmom,
                  twoPtGammas, wallPointMesonKGmuLRes, wallPointMesonKGmuL);
        // Strange insetion: K->pi momentum
        std::string wallPointMesonKGmuSRes = resultPStem + "/3pt/WallPoint/3ptK_VC" + smu + "S_tK_" + timeStamp;
        std::string wallPointMesonKGmuS = "wallPointMesonKGmuS_" + stk;
        makeMeson(application, qWallklZeromom, seqGmuKSMqmom, sink2ptPmom,
                  twoPtGammas, wallPointMesonKGmuSRes, wallPointMesonKGmuS);

        //////////////////////////////////////////////////
        // Kaon Meson Contractions: Wall-wall
        //////////////////////////////////////////////////
        // Kaon momentum 
        std::string wallWallMesonKZeromomRes = resultPStem + "/2pt/WallWall/2ptK_mom_K_tK_" + timeStamp;
        std::string wallWallMesonKZeromom = "wallWallMesonKZeromom_" + stk;
        makeMeson(application, smearedqWallklZeromom, smearedqWallksZeromom, sink2ptZeromom,
                  twoPtGammas, wallWallMesonKZeromomRes, wallWallMesonKZeromom);
        // Pion momentum
        std::string wallWallMesonKPmomRes = resultPStem + "/2pt/WallWall/2ptK_mom_P_tK_" + timeStamp;
        std::string wallWallMesonKPmom = "wallWallMesonKPmom_" + stk;
        makeMeson(application, smearedqWallklPmom, smearedqWallksPmom, sink2ptZeromom,
                  twoPtGammas, wallWallMesonKPmomRes, wallWallMesonKPmom);
        // Light insertion: K->pi momentum
        std::string wallWallMesonKLRes = resultPStem + "/3pt/WallWall/3ptK_VC" + smu + "L_tK_" + timeStamp;
        std::string wallWallMesonKL = "wallWallMesonKL_" + stk;
        makeMeson(application, smearedqWallGmuKLQmom, smearedqWallksZeromom, sink2ptPmom,
                  twoPtGammas, wallWallMesonKLRes, wallWallMesonKL);
        // Strange instertion: K->pi momentum
        std::string wallWallMesonKSRes = resultPStem + "/3pt/WallWall/3ptK_VC" + smu + "S_tK_" + timeStamp;
        std::string wallWallMesonKS = "wallWallMesonKS_" + stk;
        makeMeson(application, smearedqWallklPmom, smearedqWallGmuKSmQmom, sink2ptPmom,
                  twoPtGammas, wallWallMesonKSRes, wallWallMesonKS);

        //////////////////////////////////////////////////
        // Pion Meson Contractions: Wall-point
        //////////////////////////////////////////////////
        // Kaon momentum
        std::string wallPointPiZeromomRes = resultPStem + "/2pt/WallPoint/2ptPi_mom_K_tK_" + timeStamp;
        std::string wallPointPiZeromom = "wallPointPiZeromom_" + stk;
        makeMeson(application, qWallplZeromom, qWallplZeromom, sink2ptZeromom,
                  twoPtGammas, wallPointPiZeromomRes, wallPointPiZeromom);
        // Pion momentum
        std::string wallPointPiPmomRes = resultPStem + "/2pt/WallPoint/2ptPi_mom_P_tK_" + timeStamp;
        std::string wallPointPiPmom = "wallPointPiPmom_" + stk;
        makeMeson(application, qWallplZeromom, qWallplbarPmom, sink2ptPmom,
                  twoPtGammas, wallPointPiPmomRes, wallPointPiPmom);
        // Light insetion: K->pi momentum
        std::string wallPointPiGmuLRes = resultPStem + "/3pt/WallPoint/3ptPi_VC" + smu + "L_tK_" + timeStamp;
        std::string wallPointPiGmuL = "wallPointPiGmuL_" + stk;
        makeMeson(application, seqGmuPLMqmom, qWallplbarPmom, sink2ptZeromom,
                  twoPtGammas, wallPointPiGmuLRes, wallPointPiGmuL);
        // Lightbar insetion: K->pi momentum
        std::string wallPointPiGmuLbarRes = resultPStem + "/3pt/WallPoint/3ptPi_VC" + smu + "Lbar_tK_" + timeStamp;
        std::string wallPointPiGmuLbar = "wallPointPiGmuLbar_" + stk;
        makeMeson(application, qWallplZeromom, seqGmuPLbarQmom, sink2ptZeromom,
                  twoPtGammas, wallPointPiGmuLbarRes, wallPointPiGmuLbar);

        //////////////////////////////////////////////////
        // Pion Meson Contractions: Wall-wall
        //////////////////////////////////////////////////
        // Kaon momentum
        std::string wallWallMesonPiZeromomRes = resultPStem + "/2pt/WallWall/2ptPi_mom_K_tK_" + timeStamp;
        std::string wallWallMesonPiZeromom = "wallWallMesonPiZeromom_" + stk;
        makeMeson(application, smearedqWallplZeromom, smearedqWallplZeromom, sink2ptZeromom,
                  twoPtGammas, wallWallMesonPiZeromomRes, wallWallMesonPiZeromom);
        // Pion momentum
        std::string wallWallMesonPiPmomRes = resultPStem + "/2pt/WallWall/2ptPi_mom_P_tK_" + timeStamp;
        std::string wallWallMesonPiPmom = "wallWallMesonPiPmom_" + stk;
        makeMeson(application, smearedqWallplPmom, smearedqWallplPmom, sink2ptPmom,
                  twoPtGammas, wallWallMesonPiPmomRes, wallWallMesonPiPmom);
        // Light insertion: K->pi momentum
        std::string wallWallmesonPiLRes = resultPStem + "/3pt/WallWall/3ptPi_VC" + smu + "l_tK_" + timeStamp;
        std::string wallWallmesonPiL = "wallWallmesonPiL_" + stk;
        makeMeson(application, smearedqWallGmuPLmQmom, smearedqWallplPmom, sink2ptZeromom,
                  twoPtGammas, wallWallmesonPiLRes, wallWallmesonPiL);
        // Lightbar insertion: K->pi momentum
        std::string wallWallmesonPiLbarRes = resultPStem + "/3pt/WallWall/3ptPi_VC" + smu + "lbar_tK_" + timeStamp;
        std::string wallWallmesonPiLbar = "wallWallmesonPiLbar_" + stk;
        makeMeson(application, smearedqWallplZeromom, smearedqWallGmuPLbarQmom, sink2ptZeromom,
                  twoPtGammas, wallWallmesonPiLbarRes, wallWallmesonPiLbar);

        //////////////////////////////////////////////////
        // Gamma 3pt Contractions
        //////////////////////////////////////////////////
        // s -> d
        // Kaon momentum
        std::string kaonSDGamma3ptRes = resultPStem + "/3pt/3pt_sd_mom_K_tK_" + timeStamp;
        std::string kaonSDGamma3pt = "kaonSDGamma3pt_" + stp;
        makeGamma3pt(application, smearedqWallklZeromom, qWallksZeromom, qWallplZeromom, tp,
                  gammaInsertions, kaonSDGamma3ptRes, kaonSDGamma3pt);
        // Pion momentum       
        std::string pionSDGamma3ptRes = resultPStem + "/3pt/3pt_sd_mom_P_tK_" + timeStamp;
        std::string pionSDGamma3pt = "pionSDGamma3pt_" + stp;
        makeGamma3pt(application, smearedqWallklZeromom, qWallksPmom, qWallplbarPmom, tp,
                  gammaInsertions, pionSDGamma3ptRes, pionSDGamma3pt);
        // Spec insertion
        std::string specSDGamma3ptRes = resultPStem + "/4pt/4pt_sd_VC" + smu + "_spec_tK_" + timeStamp;
        std::string specSDGamma3pt = "specSDGamma3pt_" + stp;
        makeGamma3pt(application, smearedqWallGmuKspecQmom, qWallksZeromom, qWallplbarPmom, tp,
                  gammaInsertions, specSDGamma3ptRes, specSDGamma3pt);
        // Strange insertion
        std::string strangeSDGamma3ptRes = resultPStem + "/4pt/4pt_sd_VC" + smu + "_s_tK_" + timeStamp;
        std::string strangeSDGamma3pt = "strangeSDGamma3pt_" + stp;
        makeGamma3pt(application, smearedqWallklZeromom, seqGmuKSMqmom, qWallplbarPmom, tp,
                  gammaInsertions, strangeSDGamma3ptRes, strangeSDGamma3pt);
        // Light insertion
        std::string lightSDGamma3ptRes = resultPStem + "/4pt/4pt_sd_VC" + smu + "_lbar_tK_" + timeStamp;
        std::string lightSDGamma3pt = "lightSDGamma3pt_" + stp;
        makeGamma3pt(application, smearedqWallklZeromom, qWallksZeromom, seqGmuPLbarQmom, tp,
                  gammaInsertions, lightSDGamma3ptRes, lightSDGamma3pt);

        //////////////////////////////////////////////////
        // Weak Hamiltonian Non-Eye Contractions
        //////////////////////////////////////////////////
        // 3pt Contractions
        // Kaon momentum
        std::string resWHNE3ptKmom = resultPStem + "/3pt/3pt_Hw_Non_Eye_mom_K_tK_" + timeStamp;
        std::string WHNE3ptKmom = "WHNE_3pt_Kmom" + stk;
        makeWeakNonEye(application, qWallklZeromom, qWallksZeromom, qWallplZeromom, qWallplZeromom,
                       gammaIn, gammaOut, resWHNE3ptKmom, WHNE3ptKmom);
        // Pion momentum
        std::string resWHNE3ptPmom = resultPStem + "/3pt/3pt_Hw_Non_Eye_mom_Pi_tK_" + timeStamp;
        std::string WHNE3ptPmom = "WHNE_3pt_Pmom_" + stk;
        makeWeakNonEye(application, qWallklZeromom, qWallksPmom, qWallplZeromom, qWallplbarPmom,
                       gammaIn, gammaOut, resWHNE3ptPmom, WHNE3ptPmom);

        // 4pt Contractions
        // Kaon to pion momentum
        std::string resWHNE4ptGmuKL = resultPStem + "/4pt/4pt_Non_Eye_VC" + smu + "KL_tK_" + timeStamp;
        std::string WHNE4ptGmuKL = "4pt_VC" + smu + "_KL_" + stk;
        makeWeakNonEye(application, seqGmuKLQmom, qWallksZeromom, qWallplZeromom, qWallplbarPmom,
                       gammaIn, gammaOut, resWHNE4ptGmuKL, WHNE4ptGmuKL);

        std::string resWHNE4ptGmuKS = resultPStem + "/4pt/4pt_Non_Eye_VC" + smu + "KS_tK_" + timeStamp;
        std::string WHNE4ptGmuKS = "4pt_VC" + smu + "_KS_" + stk;
        makeWeakNonEye(application, qWallklZeromom, seqGmuKSMqmom, qWallplZeromom, qWallplbarPmom,
                       gammaIn, gammaOut, resWHNE4ptGmuKS, WHNE4ptGmuKS);

        std::string resWHNE4ptGmuPL = resultPStem + "/4pt/4pt_Non_Eye_VC" + smu + "PL_tK_" + timeStamp;
        std::string WHNE4ptGmuPL = "4pt_VC" + smu + "_PL_" + stk;
        makeWeakNonEye(application, qWallklZeromom, qWallksZeromom, seqGmuPLMqmom, qWallplbarPmom,
                       gammaIn, gammaOut, resWHNE4ptGmuPL, WHNE4ptGmuPL);

        std::string resWHNE4ptGmuPLbar = resultPStem + "/4pt/4pt_Non_Eye_VC" + smu + "PLbar_tK_" + timeStamp;
        std::string WHNE4ptGmuPLbar = "4pt_VC" + smu + "_PLbar_" + stk;
        makeWeakNonEye(application, qWallklZeromom, qWallksZeromom, qWallplZeromom, seqGmuPLbarQmom,
                       gammaIn, gammaOut, resWHNE4ptGmuPLbar, WHNE4ptGmuPLbar);

        //////////////////////////////////////////////////
        // Weak Hamiltonian Eye Contractions
        //////////////////////////////////////////////////
        for (unsigned int i = 1; i < flavour.size(); ++i)
        {
            for (unsigned int j = 0; j < sparseLoops[i].size(); ++j)
            {
                std::string loop3pt = sparseLoops[i][j];
                std::string sparseName = flavour[i] + "_" + std::to_string(j);
                // 3pt Contractions
                // Kaon momentum
                std::string resWHE3ptKmom = resultPStem + "/3pt/3pt_Hw_Eye_" + sparseName
                                            + "_mom_K_tK_" + timeStamp;
                std::string WHE3ptKmom = "WHE_Kmom_" + sparseName + "_tk_" +  stk;
                makeWeakEye(application, qWallksZeromom, qWallplZeromom, smearedqWallklZeromom, loop3pt,
                            gammaIn, gammaOut, tp, resWHE3ptKmom, WHE3ptKmom);
                // Pion momentum
                std::string resWHE3ptPmom = resultPStem + "/3pt/3pt_Hw_Eye_" + sparseName
                                               + "_mom_Pi_tK_" + timeStamp;
                std::string WHE3ptPmom = "WHE_Pmom_" + sparseName + "_tk_" +  stk;
                makeWeakEye(application, qWallksPmom, qWallplbarPmom, smearedqWallklZeromom, loop3pt,
                            gammaIn, gammaOut, tp, resWHE3ptPmom, WHE3ptPmom);
            
                // 4pt Contractions
                // Kaon to pion momentum
                std::string resWHE4ptGmuspec = resultPStem + "/4pt/4pt_Eye_" + sparseName
                                               + "_VC" + smu + "_spec_tK_" + timeStamp;
                std::string WHE4ptGmuspec = "4pt_VC" + smu + "_Eye_" + sparseName + "_spec_" + stk;
                makeWeakEye(application, qWallksZeromom, qWallplbarPmom, smearedqWallGmuKspecQmom, loop3pt,
                            gammaIn, gammaOut, tp, resWHE4ptGmuspec, WHE4ptGmuspec);

                std::string resWHE4ptGmuKS = resultPStem + "/4pt/4pt_Eye_" + sparseName
                                               + "_VC" + smu + "_KS_tK_" + timeStamp;
                std::string WHE4ptGmuKS = "4pt_VC" + smu + "_Eye_" + sparseName + "_KS_" + stk;
                makeWeakEye(application, seqGmuKSMqmom, qWallplbarPmom, smearedqWallklZeromom, loop3pt,
                            gammaIn, gammaOut, tp, resWHE4ptGmuKS, WHE4ptGmuKS);

                std::string resWHE4ptGmuKLbar = resultPStem + "/4pt/4pt_Eye_" + sparseName
                                               + "_VC" + smu + "_PLbar_tK_" + timeStamp;
                std::string WHE4ptGmuKLbar = "4pt_VC" + smu + "_Eye_" + sparseName + "_KLbar_" + stk;
                makeWeakEye(application, qWallksZeromom, seqGmuPLbarQmom, smearedqWallklZeromom, loop3pt,
                            gammaIn, gammaOut, tp, resWHE4ptGmuKLbar, WHE4ptGmuKLbar);
            }
        }
        //////////////////////////////////////////////////
        // Weak Hamiltonian Eye Loop Current Insertion
        //////////////////////////////////////////////////
        for (unsigned int i = 1; i < flavour.size(); ++i)
        {
            solver = "mcg_" + flavour[i];
            std::string action = "dwf_" + flavour[i];
            std::string seqSparseProps = "seqSparseProps_" + flavour[i] + "_" + timeStamp;
            std::string seqSparseProp = "seqSparseProp_" + flavour[i] + "_" + timeStamp;

            makeSeqZProp(application, solver, action, tj, qmom, sparseProps, seqSparseProps);
            unpackProps(application, seqSparseProps, seqSparseProp);
            makeLoops(application, seqSparseProp, sparseNoise, seqSparseLoops[i]);
        }

        for (unsigned int i = 1; i < flavour.size(); ++i)
        {
            for (unsigned int j = 0; j < sparseLoops[i].size(); ++j)
            {
                std::string seqLoop = seqSparseLoops[i][j];
                std::string sparseName = flavour[i] + std::to_string(tj) + "_" + std::to_string(j);

                std::string resWHE4ptGmuLoop = resultPStem + "/4pt/4pt_Eye" + sparseName
                                          + "_VC" + smu + "_Loop_tK_" + timeStamp;
                std::string WHE4ptGmuLoop = "4pt_VC" + smu + "_Eye_Loop_" + sparseName + "_tK_" + stk;
                makeWeakEye(application, qWallksZeromom, qWallplbarPmom, smearedqWallklZeromom, seqLoop,
                            gammaIn, gammaOut, tp, resWHE4ptGmuLoop, WHE4ptGmuLoop);
            }
        }
        //////////////////////////////////////////////////
        // Disconnected Loop Current Insertion
        //////////////////////////////////////////////////
        for (unsigned int i = 1; i < flavour.size(); ++i)
        {
            for (unsigned int j = 0; j < sparseLoops[i].size(); ++j)
            {
                std::string seqLoop = seqSparseLoops[i][j];
                std::string sparseName = flavour[i] + "_" + std::to_string(j);
                std::string resDiscLoop = resultPStem + "/disc/disc_VC" + smu
                                               + "_" + sparseName + "_tK_" + timeStamp;
                std::string discLoop = "disc_" + flavour[i] + "_" + std::to_string(tj) + "_" + std::to_string(j);
                makeDiscLoop(application, seqLoop, resDiscLoop, discLoop);
            }
        }
    }

    std::string xmlFileName = par.ioPar.xmlFileName;
    // execution
    unsigned int prec = 16;
    application.saveParameterFile(xmlFileName, prec);

    if (par.ioPar.scheduleFile.empty())
    {
        application.schedule();
        application.saveSchedule("rarek.sched");
    }
    else
    {
        application.loadSchedule(par.ioPar.scheduleFile);
    }
    application.run();

    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;

    Grid_finalize();

    return EXIT_SUCCESS;
}