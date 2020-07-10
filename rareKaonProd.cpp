#include "Utils/ApplicationUtils.hpp"
#include <Grid/Grid.h>
#include <Hadrons/A2AVectors.hpp>
#include <Hadrons/Modules/MContraction/A2AMesonField.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MContraction;

namespace RareKaonInputs
{
    class RunPar : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(RunPar,
                                        std::string, runId);
    };

    class TrajRange : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(TrajRange,
                                        unsigned int, start,
                                        unsigned int, end,
                                        unsigned int, step);
    };

    class DbPar : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(DbPar,
                                        std::string, applicationDb,
                                        std::string, resultDb,
                                        bool,        restoreSchedule,
                                        bool,        restoreModules,
                                        bool,        restoreMemoryProfile,
                                        bool,        makeStatDb);
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
                                        std::string,  resultStem,
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
    RareKaonInputs::RunPar           runPar;
    RareKaonInputs::MomPar           momPar;
    RareKaonInputs::DbPar            dbPar;
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
    read(reader,           "runPar",           par.runPar);
    read(reader,           "momPar",           par.momPar);
    read(reader,            "dbPar",            par.dbPar);
    read(reader,            "ioPar",            par.ioPar);

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
    std::string sKMom = sanitizeMom(sinkkmom);
    std::string sPMom = sanitizeMom(sinkpmom);
    std::vector<double> vKMom = strToVec<double>(sinkkmom);
    std::vector<double> vPMom = strToVec<double>(sinkpmom);
    std::vector<double> vQMom = strToVec<double>("-" + sinkpmom);
    
    std::string resultStem = par.ioPar.resultStem;

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
    globalPar.runId                         = par.runPar.runId;
    globalPar.trajCounter.start             = par.trajRange.start;
    globalPar.trajCounter.end               = par.trajRange.end;
    globalPar.trajCounter.step              = par.trajRange.step;
    globalPar.database.applicationDb        = par.dbPar.applicationDb;
    globalPar.database.resultDb             = par.dbPar.resultDb;
    globalPar.database.restoreSchedule      = par.dbPar.restoreSchedule;
    globalPar.database.restoreModules       = par.dbPar.restoreModules;
    globalPar.database.restoreMemoryProfile = par.dbPar.restoreMemoryProfile;
    globalPar.database.makeStatDb           = par.dbPar.makeStatDb;
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

    std::vector<std::vector<std::string>> sparseLoops, seqSparseLoops;
    sparseLoops.resize(flavour.size(), std::vector<std::string>(nds));
    seqSparseLoops.resize(flavour.size(), std::vector<std::string>(nds));
    for (unsigned int i = 1; i < flavour.size(); ++i)
    {
        solver = "mcg_" + flavour[i];
        sparseProps = "sparseProps_" + flavour[i];
        sparseProp = "sparseProp_" + flavour[i];

        makeZGaugeProp(application, solver, sparseNoises, sparseProps);
        unpackProps(application, sparseProps, sparseProp);
        makeLoops(application, sparseProp, sparseNoise, sparseLoops[i]);
    }

    //////////////////////////////////////////////////
    // Disconnected Loop
    //////////////////////////////////////////////////
    std::vector<std::vector<std::string>> discMom = makeDiscMom(sinkkmom, sinkpmom);
    std::vector<std::vector<double>> vDiscMom;
    std::vector<std::string> sDiscMom;
    vDiscMom.push_back(vKMom); vDiscMom.push_back(vPMom);
    sDiscMom.push_back(sKMom); sDiscMom.push_back(sPMom);
    for (unsigned int i = 1; i < flavour.size(); ++i)
    for (unsigned int j = 0; j < sparseLoops[i].size(); ++j)
    for (unsigned int m = 0; m < discMom.size(); ++m)
    {
        std::string loop = sparseLoops[i][j];
        std::string sparseName = flavour[i] + "_" + std::to_string(j);
        std::string resDiscLoop = resultStem + "/disc/disc_"
                                  + sDiscMom[m] + "_" + sparseName;
        std::string discLoop = "disc_" + sDiscMom[m] + "_" + sparseName;
        auto dlVcEntry = makeDiscLoopEntry(vDiscMom[m], flavour[i], j);
        makeDiscLoop(application, loop, discMom[m], resDiscLoop, discLoop,
                     dlVcEntry, "Disconnected");
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
        BaseEntry baseEntry = makeBaseEntry(vKMom, vPMom, vQMom,
                                            tk, tj, tp);

        //////////////////////////////////////////////////
        // Propagators
        //////////////////////////////////////////////////
        // Zero momentum
        std::string qWallksZeromom = "QWall_s_0mom_" + stk;
        makeWallProp(application, strangeSolver, kmom, tk, qWallksZeromom);
        std::string qWallklZeromom = "QWall_l_0mom_" + stk;
        makeWallZProp(application, lightSolver, kmom, tk, qWallklZeromom);
        std::string qWallplZeromom = "QWall_l_0mom_" + stp;
        makeWallZProp(application, lightSolver, kmom, tp, qWallplZeromom);

        // Smeared propagators: sink zero mom
        std::string smearedqWallksZeromom = "smearedQWall_s_0mom_" + stk;
        makeSmearedProp(application, qWallksZeromom, sinkZeromom, smearedqWallksZeromom);
        std::string smearedqWallklZeromom = "smearedQWall_l_0mom_" + stk;
        makeSmearedProp(application, qWallklZeromom, sinkZeromom, smearedqWallklZeromom);
        std::string smearedqWallplZeromom = "smearedQWall_l_0mom_" + stp;
        makeSmearedProp(application, qWallklZeromom, sinkZeromom, smearedqWallplZeromom);

        // Nonzero momentum
        std::string qWallksPmom = "QWall_s_Pmom_" + stk;
        makeWallProp(application, strangeSolver, pmom, tk, qWallksPmom);
        std::string qWallklPmom = "QWall_l_Pmom_" + stk;
        makeWallZProp(application, lightSolver, pmom, tk, qWallklPmom);
        std::string qWallplbarPmom = "QWall_l_Pmom_" + stp;
        makeWallZProp(application, lightSolver, pmom, tp, qWallplbarPmom);

        // Smeared propagators
        std::string smearedqWallksPmom = "smearedQWall_s_Pmom" + stk;
        makeSmearedProp(application, qWallksPmom, sinkZeromom, smearedqWallksPmom);
        std::string smearedqWallklPmom = "smearedQWall_l_Pmom" + stk;
        makeSmearedProp(application, qWallklPmom, sinkZeromom, smearedqWallklPmom);
        std::string smearedqWallplPmom = "smearedQWall_l_Pmom" + stp;
        makeSmearedProp(application, qWallklPmom, sinkZeromom, smearedqWallplPmom);

        //////////////////////////////////////////////////
        // Sequential propagators
        //////////////////////////////////////////////////
        std::string wallSourceKmomK = makeWallSourceName(tk, kmom);
        std::string wallSourcePmomK = makeWallSourceName(tp, kmom);
        std::string wallSourcePmomP = makeWallSourceName(tp, pmom);
        std::string seqVcKLQmom = "VC" + smu + "_KL_" + qWallklZeromom;
        makeSeqZProp(application, lightSolver, lightAction, 
                     tj, qmom, qWallklZeromom, wallSourceKmomK, seqVcKLQmom);
        std::string seqVcKSMqmom = "VC" + smu + "_KS_" + qWallksZeromom;
        makeSeqProp(application, strangeSolver, strangeAction,
                    tj, mqmom, qWallksZeromom, wallSourceKmomK, seqVcKSMqmom);
        std::string seqVcPLMqmom = "VC" + smu + "_PL_" + qWallplZeromom;
        makeSeqZProp(application, lightSolver, lightAction,
                     tj, mqmom, qWallplZeromom, wallSourcePmomK, seqVcPLMqmom);
        std::string seqVcPLbarQmom = "VC" + smu + "_PLbar_" + qWallplbarPmom;
        makeSeqZProp(application, lightSolver, lightAction,
                     tj, qmom, qWallplbarPmom, wallSourcePmomP, seqVcPLbarQmom);

        // Smeared sequential propagators
        std::string smearedqWallVcKLQmom = "smearedQWall_Kl_VC" + smu + "_Qmom_" + stk;
        makeSmearedProp(application, seqVcKLQmom, sinkPmom, smearedqWallVcKLQmom);
        std::string smearedqWallVcKspecQmom = "smearedQWall_kspec_VC" + smu + "_Qmom_" + stk;
        makeSmearedProp(application, seqVcKLQmom, sinkZeromom, smearedqWallVcKspecQmom);
        std::string smearedqWallVcKSmQmom = "smearedQWall_Ks_VC" + smu + "_mQmom_" + stk;
        makeSmearedProp(application, seqVcKSMqmom, sinkZeromom, smearedqWallVcKSmQmom);
        std::string smearedqWallVcPLmQmom = "smearedQWall_PL_VC" + smu + "_mQmom_" + stk;
        makeSmearedProp(application, seqVcPLMqmom, sinkZeromom, smearedqWallVcPLmQmom);
        std::string smearedqWallVcPLbarQmom = "smearedQWall_PL_VC" + smu + "_Qmom_" + stk;
        makeSmearedProp(application, seqVcPLbarQmom, sinkZeromom, smearedqWallVcPLbarQmom);

        //////////////////////////////////////////////////
        // Contractions
        //////////////////////////////////////////////////

        //////////////////////////////////////////////////
        // Kaon Meson Contractions: Point-wall
        //////////////////////////////////////////////////
        // Kaon momentum 
        std::string PointWallMesonKZeromomRes = resultStem + "/2pt/PointWall/2ptKaon_PW_mom" + sKMom +"_tK_" + timeStamp;
        std::string PointWallMesonKZeromom = "PointWallMesonKZeromom_" + stk;
        auto pwKK2ptEntry = makeEntry2pt(tk, "kaon", "point", vKMom);
        makeMeson(application, qWallklZeromom, qWallksZeromom, sink2ptZeromom,
                  twoPtGammas, PointWallMesonKZeromomRes, PointWallMesonKZeromom,
                  pwKK2ptEntry, "Two_point");
        // Pion momentum
        std::string PointWallMesonKPmomRes = resultStem + "/2pt/PointWall/2ptKaon_PW_mom" + sPMom + "_tK_" + timeStamp;
        std::string PointWallMesonKPmom = "PointWallMesonKPmom_" + stk;
        auto pwKP2ptEntry = makeEntry2pt(tk, "kaon", "point", vPMom);
        makeMeson(application, qWallklZeromom, qWallksPmom, sink2ptPmom,
                  twoPtGammas, PointWallMesonKPmomRes, PointWallMesonKPmom,
                  pwKP2ptEntry, "Two_point");
        // Light insetion: K->pi momentum
        std::string PointWallMesonKVcLRes = resultStem + "/3pt/VC" + smu + "/3ptKaon_PW_VC" + smu + "L_tK_" + timeStamp;
        std::string PointWallMesonKVcL = "PointWallMesonKVcL_" + stk;
        auto pwKL3ptEntry = makeEntry3ptVC(vKMom, vQMom, vPMom, tk, tj, "kaon_l", "point");
        makeMeson(application, seqVcKLQmom, qWallksZeromom, sink2ptPmom,
                  twoPtGammas, PointWallMesonKVcLRes, PointWallMesonKVcL,
                  pwKL3ptEntry, "Three_point_Vc");
        // Strange insetion: K->pi momentum
        std::string PointWallMesonKVcSRes = resultStem + "/3pt/VC" + smu + "/3ptKaon_PW_VC" + smu + "S_tK_" + timeStamp;
        std::string PointWallMesonKVcS = "PointWallMesonKVcS_" + stk;
        auto pwKS3ptEntry = makeEntry3ptVC(vKMom, vQMom, vPMom, tk, tj, "kaon_sbar", "point");
        makeMeson(application, qWallklZeromom, seqVcKSMqmom, sink2ptPmom,
                  twoPtGammas, PointWallMesonKVcSRes, PointWallMesonKVcS,
                  pwKS3ptEntry, "Three_point_Vc");

        //////////////////////////////////////////////////
        // Kaon Meson Contractions: Wall-wall
        //////////////////////////////////////////////////
        // Kaon momentum 
        std::string wallWallMesonKZeromomRes = resultStem + "/2pt/WallWall/2ptKaon_WW_mom" + sKMom +"_tK_" + timeStamp;
        std::string wallWallMesonKZeromom = "wallWallMesonKZeromom_" + stk;
        auto wwKK2ptEntry = makeEntry2pt(tk, "kaon", "wall", vKMom);
        makeMeson(application, smearedqWallklZeromom, smearedqWallksZeromom, sink2ptZeromom,
                  twoPtGammas, wallWallMesonKZeromomRes, wallWallMesonKZeromom,
                  wwKK2ptEntry, "Two_point");
        // Pion momentum
        std::string wallWallMesonKPmomRes = resultStem + "/2pt/WallWall/2ptKaon_WW_mom" + sPMom + "_tK_" + timeStamp;
        std::string wallWallMesonKPmom = "wallWallMesonKPmom_" + stk;
        auto wwKP2ptEntry = makeEntry2pt(tk, "kaon", "wall", vPMom);
        makeMeson(application, smearedqWallklPmom, smearedqWallksPmom, sink2ptZeromom,
                  twoPtGammas, wallWallMesonKPmomRes, wallWallMesonKPmom,
                  wwKP2ptEntry, "Two_point");
        // Light insertion: K->pi momentum
        std::string wallWallMesonKLRes = resultStem + "/3pt/VC" + smu + "/3ptKaon_WW_VC" + smu + "L_tK_" + timeStamp;
        std::string wallWallMesonKL = "wallWallMesonKL_" + stk;
        auto wwKL3ptEntry = makeEntry3ptVC(vKMom, vQMom, vPMom, tk, tj, "kaon_l", "wall");
        makeMeson(application, smearedqWallVcKLQmom, smearedqWallksZeromom, sink2ptPmom,
                  twoPtGammas, wallWallMesonKLRes, wallWallMesonKL,
                  wwKL3ptEntry, "Three_point_Vc");
        // Strange instertion: K->pi momentum
        std::string wallWallMesonKSRes = resultStem + "/3pt/VC" + smu + "/3ptKaon_WW_VC" + smu + "S_tK_" + timeStamp;
        std::string wallWallMesonKS = "wallWallMesonKS_" + stk;
        auto wwKS3ptEntry = makeEntry3ptVC(vKMom, vQMom, vPMom, tk, tj, "kaon_sbar", "wall");
        makeMeson(application, smearedqWallklPmom, smearedqWallVcKSmQmom, sink2ptPmom,
                  twoPtGammas, wallWallMesonKSRes, wallWallMesonKS,
                  wwKS3ptEntry, "Three_point_Vc");

        //////////////////////////////////////////////////
        // Pion Meson Contractions: Point-wall
        //////////////////////////////////////////////////
        // Kaon momentum
        std::string PointWallPiZeromomRes = resultStem + "/2pt/PointWall/2ptPion_PW_mom" + sKMom +"_tK_" + timeStamp;
        std::string PointWallPiZeromom = "PointWallPiZeromom_" + stk;
        auto pwPK2ptEntry = makeEntry2pt(tp, "pion", "point", vKMom);
        makeMeson(application, qWallplZeromom, qWallplZeromom, sink2ptZeromom,
                  twoPtGammas, PointWallPiZeromomRes, PointWallPiZeromom,
                  pwPK2ptEntry, "Two_point");
        // Pion momentum
        std::string PointWallPiPmomRes = resultStem + "/2pt/PointWall/2ptPion_PW_mom" + sPMom + "_tK_" + timeStamp;
        std::string PointWallPiPmom = "PointWallPiPmom_" + stk;
        auto pwPP2ptEntry = makeEntry2pt(tp, "pion", "point", vPMom);
        makeMeson(application, qWallplZeromom, qWallplbarPmom, sink2ptPmom,
                  twoPtGammas, PointWallPiPmomRes, PointWallPiPmom,
                  pwPP2ptEntry, "Two_point");
        // Light insetion: K->pi momentum
        std::string PointWallPiVcLRes = resultStem + "/3pt/VC" + smu + "/3ptPion_PW_VC" + smu + "L_tK_" + timeStamp;
        std::string PointWallPiVcL = "PointWallPiVcL_" + stk;
        auto pwPL3ptEntry = makeEntry3ptVC(vKMom, vQMom, vPMom, tp, tj, "pion_l", "point");
        makeMeson(application, seqVcPLMqmom, qWallplbarPmom, sink2ptZeromom,
                  twoPtGammas, PointWallPiVcLRes, PointWallPiVcL,
                  pwPL3ptEntry, "Three_point_Vc");
        // Lightbar insetion: K->pi momentum
        std::string PointWallPiVcLbarRes = resultStem + "/3pt/VC" + smu + "/3ptPion_PW_VC" + smu + "Lbar_tK_" + timeStamp;
        std::string PointWallPiVcLbar = "PointWallPiVcLbar_" + stk;
        auto pwPLbar3ptEntry = makeEntry3ptVC(vKMom, vQMom, vPMom, tp, tj, "pion_lbar", "point");
        makeMeson(application, qWallplZeromom, seqVcPLbarQmom, sink2ptZeromom,
                  twoPtGammas, PointWallPiVcLbarRes, PointWallPiVcLbar,
                  pwPLbar3ptEntry, "Three_point_Vc");

        //////////////////////////////////////////////////
        // Pion Meson Contractions: Wall-wall
        //////////////////////////////////////////////////
        // Kaon momentum
        std::string wallWallMesonPiZeromomRes = resultStem + "/2pt/WallWall/2ptPion_WW_mom" + sKMom +"_tK_" + timeStamp;
        std::string wallWallMesonPiZeromom = "wallWallMesonPiZeromom_" + stk;
        auto wwPK2ptEntry = makeEntry2pt(tp, "pion", "wall", vKMom);
        makeMeson(application, smearedqWallplZeromom, smearedqWallplZeromom, sink2ptZeromom,
                  twoPtGammas, wallWallMesonPiZeromomRes, wallWallMesonPiZeromom,
                  wwPK2ptEntry, "Two_point");
        // Pion momentum
        std::string wallWallMesonPiPmomRes = resultStem + "/2pt/WallWall/2ptPion_WW_mom" + sPMom + "_tK_" + timeStamp;
        std::string wallWallMesonPiPmom = "wallWallMesonPiPmom_" + stk;
        auto wwPP2ptEntry = makeEntry2pt(tp, "pion", "wall", vPMom);
        makeMeson(application, smearedqWallplPmom, smearedqWallplPmom, sink2ptPmom,
                  twoPtGammas, wallWallMesonPiPmomRes, wallWallMesonPiPmom,
                  wwPP2ptEntry, "Two_point");
        // Light insertion: K->pi momentum
        std::string wallWallmesonPiLRes = resultStem + "/3pt/VC" + smu + "/3ptPion_WW_VC" + smu + "l_tK_" + timeStamp;
        std::string wallWallmesonPiL = "wallWallmesonPiL_" + stk;
        auto wwPL3ptEntry = makeEntry3ptVC(vKMom, vQMom, vPMom, tp, tj, "pion_l", "wall");
        makeMeson(application, smearedqWallVcPLmQmom, smearedqWallplPmom, sink2ptZeromom,
                  twoPtGammas, wallWallmesonPiLRes, wallWallmesonPiL,
                  wwPL3ptEntry, "Three_point_Vc");
        // Lightbar insertion: K->pi momentum
        std::string wallWallmesonPiLbarRes = resultStem + "/3pt/VC" + smu + "/3ptPion_WW_VC" + smu + "lbar_tK_" + timeStamp;
        std::string wallWallmesonPiLbar = "wallWallmesonPiLbar_" + stk;
        auto wwPLbar3ptEntry = makeEntry3ptVC(vKMom, vQMom, vPMom, tp, tj, "pion_lbar", "wall");
        makeMeson(application, smearedqWallplZeromom, smearedqWallVcPLbarQmom, sink2ptZeromom,
                  twoPtGammas, wallWallmesonPiLbarRes, wallWallmesonPiLbar,
                  wwPLbar3ptEntry, "Three_point_Vc");

        //////////////////////////////////////////////////
        // Gamma 3pt Contractions
        //////////////////////////////////////////////////
        // s -> d
        // Kaon momentum
        std::string kaonSDGamma3ptRes = resultStem + "/3pt/sd/3pt_sd_mom" + sKMom +"_tK_" + timeStamp;
        std::string kaonSDGamma3pt = "kaonSDGamma3pt_" + stp;
        auto sdK3ptEntry = makeEntry3ptHw(tk, tp, "sd", vKMom, "", -1, true);
        makeGamma3pt(application, smearedqWallklZeromom, qWallksZeromom, qWallplZeromom, tp,
                  gammaInsertions, kaonSDGamma3ptRes, kaonSDGamma3pt,
                  sdK3ptEntry, "Three_point_Hw");
        // Pion momentum       
        std::string pionSDGamma3ptRes = resultStem + "/3pt/sd/3pt_sd_mom" + sPMom + "_tK_" + timeStamp;
        std::string pionSDGamma3pt = "pionSDGamma3pt_" + stp;
        auto sdP3ptEntry = makeEntry3ptHw(tk, tp, "sd", vPMom, "", -1, true);
        makeGamma3pt(application, smearedqWallklZeromom, qWallksPmom, qWallplbarPmom, tp,
                  gammaInsertions, pionSDGamma3ptRes, pionSDGamma3pt,
                  sdP3ptEntry, "Three_point_Hw");
        // Spec insertion
        std::string specSDGamma3ptRes = resultStem + "/4pt/sd/4pt_sd_VC" + smu + "_spec_tK_" + timeStamp;
        std::string specSDGamma3pt = "specSDGamma3pt_" + stp;
        auto sdSpec4ptEntry = makeEntry4pt(baseEntry, "sd_l", "", -1, true);
        makeGamma3pt(application, smearedqWallVcKspecQmom, qWallksZeromom, qWallplbarPmom, tp,
                  gammaInsertions, specSDGamma3ptRes, specSDGamma3pt,
                  sdSpec4ptEntry, "Four_point");
        // Strange insertion
        std::string strangeSDGamma3ptRes = resultStem + "/4pt/sd/4pt_sd_VC" + smu + "_s_tK_" + timeStamp;
        std::string strangeSDGamma3pt = "strangeSDGamma3pt_" + stp;
        auto sdS4ptEntry = makeEntry4pt(baseEntry, "sd_sbar", "", -1, true);
        makeGamma3pt(application, smearedqWallklZeromom, seqVcKSMqmom, qWallplbarPmom, tp,
                  gammaInsertions, strangeSDGamma3ptRes, strangeSDGamma3pt,
                  sdS4ptEntry, "Four_point");
        // Light insertion
        std::string lightSDGamma3ptRes = resultStem + "/4pt/sd/4pt_sd_VC" + smu + "_lbar_tK_" + timeStamp;
        std::string lightSDGamma3pt = "lightSDGamma3pt_" + stp;
        auto sd4ptEntry = makeEntry4pt(baseEntry, "sd_lbar", "", -1, true);
        makeGamma3pt(application, smearedqWallklZeromom, qWallksZeromom, seqVcPLbarQmom, tp,
                  gammaInsertions, lightSDGamma3ptRes, lightSDGamma3pt,
                  sd4ptEntry, "Four_point");

        //////////////////////////////////////////////////
        // Weak Hamiltonian Non-Eye Contractions
        //////////////////////////////////////////////////
        // 3pt Contractions
        // Kaon momentum
        std::string resWHNE3ptKmom = resultStem + "/3pt/HW/3pt_Hw_Non_Eye_mom" + sKMom +"_tK_" + timeStamp;
        std::string WHNE3ptKmom = "WHNE_3pt_Kmom" + stk;
        auto wneK3ptEntry = makeEntry3ptHw(tk, tp, "NE", vKMom, "", -1, true);
        makeWeakNonEye(application, qWallklZeromom, qWallksZeromom, qWallplZeromom, qWallplZeromom,
                       gammaIn, gammaOut, resWHNE3ptKmom, WHNE3ptKmom,
                       wneK3ptEntry, "Three_point_Hw");
        // Pion momentum
        std::string resWHNE3ptPmom = resultStem + "/3pt/HW/3pt_Hw_Non_Eye_mom" + sPMom + "_tK_" + timeStamp;
        std::string WHNE3ptPmom = "WHNE_3pt_Pmom_" + stk;
        auto wneP3ptEntry = makeEntry3ptHw(tk, tp, "NE", vPMom, "", -1, true);
        makeWeakNonEye(application, qWallklZeromom, qWallksPmom, qWallplZeromom, qWallplbarPmom,
                       gammaIn, gammaOut, resWHNE3ptPmom, WHNE3ptPmom,
                       wneP3ptEntry, "Three_point_Hw");

        // 4pt Contractions
        // Kl
        std::string resWHNE4ptVcKL = resultStem + "/4pt/RK/4pt_Non_Eye_VC" + smu + "KL_tK_" + timeStamp;
        std::string WHNE4ptVcKL = "4pt_VC" + smu + "_KL_" + stk;
        auto wneKl4ptEntry = makeEntry4pt(baseEntry, "NE_Kl", "", -1, true);
        makeWeakNonEye(application, seqVcKLQmom, qWallksZeromom, qWallplZeromom, qWallplbarPmom,
                       gammaIn, gammaOut, resWHNE4ptVcKL, WHNE4ptVcKL,
                       wneKl4ptEntry, "Four_point");
        // Ksbar
        std::string resWHNE4ptVcKS = resultStem + "/4pt/RK/4pt_Non_Eye_VC" + smu + "KS_tK_" + timeStamp;
        std::string WHNE4ptVcKS = "4pt_VC" + smu + "_KS_" + stk;
        auto wneKsbar4ptEntry = makeEntry4pt(baseEntry, "NE_Ksbar", "", -1, true);
        makeWeakNonEye(application, qWallklZeromom, seqVcKSMqmom, qWallplZeromom, qWallplbarPmom,
                       gammaIn, gammaOut, resWHNE4ptVcKS, WHNE4ptVcKS,
                       wneKsbar4ptEntry, "Four_point");
        // Pil
        std::string resWHNE4ptVcPL = resultStem + "/4pt/RK/4pt_Non_Eye_VC" + smu + "PL_tK_" + timeStamp;
        std::string WHNE4ptVcPL = "4pt_VC" + smu + "_PL_" + stk;
        auto wnePil4ptEntry = makeEntry4pt(baseEntry, "NE_Pil", "", -1, true);
        makeWeakNonEye(application, qWallklZeromom, qWallksZeromom, seqVcPLMqmom, qWallplbarPmom,
                       gammaIn, gammaOut, resWHNE4ptVcPL, WHNE4ptVcPL,
                       wnePil4ptEntry, "Four_point");
        // Pilbar
        std::string resWHNE4ptVcPLbar = resultStem + "/4pt/RK/4pt_Non_Eye_VC" + smu + "PLbar_tK_" + timeStamp;
        std::string WHNE4ptVcPLbar = "4pt_VC" + smu + "_PLbar_" + stk;
        auto wnePilbar4ptEntry = makeEntry4pt(baseEntry, "NE_Pilbar", "", -1, true);
        makeWeakNonEye(application, qWallklZeromom, qWallksZeromom, qWallplZeromom, seqVcPLbarQmom,
                       gammaIn, gammaOut, resWHNE4ptVcPLbar, WHNE4ptVcPLbar,
                       wnePilbar4ptEntry, "Four_point");

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
                std::string resWHE3ptKmom = resultStem + "/3pt/HW/3pt_Hw_Eye_" + sparseName
                                            + "_mom" + sKMom +"_tK_" + timeStamp;
                std::string WHE3ptKmom = "WHE_Kmom_" + sparseName + "_tk_" +  stk;
                auto weK3ptEntry = makeEntry3ptHw(tk, tp, "E", vKMom, flavour[i], j);
                makeWeakEye(application, qWallksZeromom, qWallplZeromom, smearedqWallklZeromom, loop3pt,
                            gammaIn, gammaOut, tp, resWHE3ptKmom, WHE3ptKmom,
                            weK3ptEntry, "Three_point_Hw");
                // Pion momentum
                std::string resWHE3ptPmom = resultStem + "/3pt/HW/3pt_Hw_Eye_" + sparseName
                                               + "_mom" + sPMom + "_tK_" + timeStamp;
                std::string WHE3ptPmom = "WHE_Pmom_" + sparseName + "_tk_" +  stk;
                auto weP3ptEntry = makeEntry3ptHw(tk, tp, "E", vPMom, flavour[i], j);
                makeWeakEye(application, qWallksPmom, qWallplbarPmom, smearedqWallklZeromom, loop3pt,
                            gammaIn, gammaOut, tp, resWHE3ptPmom, WHE3ptPmom,
                            weP3ptEntry, "Three_point_Hw");
            
                // 4pt Contractions
                // l
                std::string resWHE4ptVcspec = resultStem + "/4pt/RK/4pt_Eye_" + sparseName
                                               + "_VC" + smu + "_spec_tK_" + timeStamp;
                std::string WHE4ptVcspec = "4pt_VC" + smu + "_Eye_" + sparseName + "_spec_" + stk;
                auto weL4ptEntry = makeEntry4pt(baseEntry, "E_l", flavour[i], j);
                makeWeakEye(application, qWallksZeromom, qWallplbarPmom, smearedqWallVcKspecQmom, loop3pt,
                            gammaIn, gammaOut, tp, resWHE4ptVcspec, WHE4ptVcspec,
                            weL4ptEntry, "Four_point");

                // Ksbar
                std::string resWHE4ptVcKS = resultStem + "/4pt/RK/4pt_Eye_" + sparseName
                                               + "_VC" + smu + "_KS_tK_" + timeStamp;
                std::string WHE4ptVcKS = "4pt_VC" + smu + "_Eye_" + sparseName + "_KS_" + stk;
                auto weKS4ptEntry = makeEntry4pt(baseEntry, "E_Ksbar", flavour[i], j);
                makeWeakEye(application, seqVcKSMqmom, qWallplbarPmom, smearedqWallklZeromom, loop3pt,
                            gammaIn, gammaOut, tp, resWHE4ptVcKS, WHE4ptVcKS,
                            weKS4ptEntry, "Four_point");

                // Pilbar
                std::string resWHE4ptVcKLbar = resultStem + "/4pt/RK/4pt_Eye_" + sparseName
                                               + "_VC" + smu + "_PLbar_tK_" + timeStamp;
                std::string WHE4ptVcKLbar = "4pt_VC" + smu + "_Eye_" + sparseName + "_KLbar_" + stk;
                auto wePiLbar4ptEntry = makeEntry4pt(baseEntry, "E_PiLbar", flavour[i], j);
                makeWeakEye(application, qWallksZeromom, seqVcPLbarQmom, smearedqWallklZeromom, loop3pt,
                            gammaIn, gammaOut, tp, resWHE4ptVcKLbar, WHE4ptVcKLbar,
                            wePiLbar4ptEntry, "Four_point");
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

            makeSeqZProp(application, solver, action, tj, qmom, sparseProps, sparseNoises, seqSparseProps);
            unpackProps(application, seqSparseProps, seqSparseProp);
            makeLoops(application, seqSparseProp, sparseNoise, seqSparseLoops[i]);
        }

        for (unsigned int i = 1; i < flavour.size(); ++i)
        {
            for (unsigned int j = 0; j < sparseLoops[i].size(); ++j)
            {
                std::string seqLoop = seqSparseLoops[i][j];
                std::string sparseName = flavour[i] + "_" + std::to_string(j);

                std::string resWHE4ptVcLoop = resultStem + "/4pt/RK/4pt_Eye_" + sparseName
                                          + "_VC" + smu + "_Loop_tK_" + timeStamp;
                std::string WHE4ptVcLoop = "4pt_VC" + smu + "_Eye_Loop_" + sparseName + "_tK_" + stk;
                auto weLoop4ptEntry = makeEntry4pt(baseEntry, "E_loop", flavour[i], j);
                makeWeakEye(application, qWallksZeromom, qWallplbarPmom, smearedqWallklZeromom, seqLoop,
                            gammaIn, gammaOut, tp, resWHE4ptVcLoop, WHE4ptVcLoop,
                            weLoop4ptEntry, "Four_point");
            }
        }

        //////////////////////////////////////////////////
        // Neutral Disconnected Diagram
        //////////////////////////////////////////////////
        for (unsigned int i = 1; i < flavour.size(); ++i)
        {
            for (unsigned int j = 0; j < sparseLoops[i].size(); ++j)
            {
                std::string loop = sparseLoops[i][j];
                std::string sparseName = flavour[i] + "_" + std::to_string(j);
                std::string resDisc0Loop = resultStem + "/4pt/disc0/disc0_VC" + smu
                                           + "_" + sparseName + "mom" + sPMom + "_tK_" + timeStamp;
                std::string disc0Loop = "disc0_VC" + smu + "_" + flavour[i]
                                        + "_" + std::to_string(tj) + "_" + std::to_string(j);
                auto pi04ptEntry = makeEntry4pt(baseEntry, "pi0", flavour[i], j);
                makeRareKaonNeutralDisc(application, qWallklZeromom, qWallksZeromom,
                                        seqVcPLbarQmom, loop, resDisc0Loop, disc0Loop,
                                        pi04ptEntry, "Four_point");
            }
        }
    }

    std::string xmlFileName = par.ioPar.xmlFileName;
    // execution
    unsigned int prec = 16;
    application.saveParameterFile(xmlFileName, prec);
    application.run();

    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;

    Grid_finalize();

    return EXIT_SUCCESS;
}