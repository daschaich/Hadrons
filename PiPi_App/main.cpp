#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

using namespace Grid;
using namespace Hadrons;

struct MesonEntry: public SqlEntry
{
    HADRONS_SQL_FIELDS(SqlNotNull<std::string>, q1,
                       SqlNotNull<std::string>, q2,
                       SqlNotNull<std::string>, source);
};

struct PiPiEntry: public SqlEntry
{
    HADRONS_SQL_FIELDS(SqlNotNull<std::string>, q1,
                        SqlNotNull<std::string>, q2,
                        SqlNotNull<std::string>, q3,
                        SqlNotNull<std::string>, q4
                    );
};

int main(int argc, char *argv[])
{
    // initialise Grid /////////////////////////////////////////////////////////
    Grid_init(&argc, &argv);
    HadronsLogError.Active(GridLogError.isActive());
    HadronsLogWarning.Active(GridLogWarning.isActive());
    HadronsLogMessage.Active(GridLogMessage.isActive());
    HadronsLogIterative.Active(GridLogIterative.isActive());
    HadronsLogDebug.Active(GridLogDebug.isActive());
    LOG(Message) << "Grid initialized" << std::endl;

    // initialise application //////////////////////////////////////////////////
    Application            application;

    //Global parameters
    Application::GlobalPar globalPar;
    globalPar.runId = "pion";
    globalPar.database.resultDb = "pionResult.db";
    globalPar.trajCounter.start = 1500;
    globalPar.trajCounter.end = 1520;
    globalPar.trajCounter.step = 20;


    // global initialisation
    application.setPar(globalPar);

    // create modules //////////////////////////////////////////////////////////

    // gauge field
    application.createModule<MGauge::Random>("gauge");
    std::string flavor = "l";

    // stout smearing
    MGauge::StoutSmearing::Par stoutPar;
    stoutPar.gauge = "gauge";
    stoutPar.steps = 3;
    stoutPar.rho = 0.1;
    application.createModule<MGauge::StoutSmearing>("stout", stoutPar);

    MUtilities::GaugeSinglePrecisionCast::Par stoutfPar;
    stoutfPar.field = "stout";
    application.createModule<MUtilities::GaugeSinglePrecisionCast>("stoutf", stoutfPar);


    //action
    MAction::ScaledDWF::Par actionPar;
    actionPar.gauge = "stout";
    actionPar.Ls = 12;
    actionPar.M5 = 1.8;
    actionPar.scale=2.0;
    actionPar.mass = 0.01;
    actionPar.boundary = "1 1 1 -1";
    application.createModule<MAction::ScaledDWF>("DWF_outer_" + flavor, actionPar);

    actionPar.gauge = "stoutf";
    application.createModule<MAction::ScaledDWFF>("DWF_inner_" + flavor, actionPar);


    //solver
    MSolver::MixedPrecisionRBPrecCG::Par solverPar;
    solverPar.innerAction = "DWF_inner_" + flavor;
    solverPar.outerAction = "DWF_outer_" + flavor;
    solverPar.residual = 1.0e-8;
    solverPar.maxInnerIteration = 10000;
    solverPar.maxOuterIteration = 10;
    application.createModule<MSolver::MixedPrecisionRBPrecCG>("CG_" + flavor, solverPar);

    // source
    MSource::Point::Par ptPar;
    ptPar.position = "0 0 0 0";
    application.createModule<MSource::Point>("pt", ptPar);

    //prop
    MFermion::GaugeProp::Par quarkPar;
    quarkPar.solver = "CG_" + flavor;
    quarkPar.source = "pt";
    application.createModule<MFermion::GaugeProp>("Qpt_" + flavor, quarkPar);


    // sink
    MSink::Point::Par sinkPar;
    sinkPar.mom = "0 0 0";
    application.createModule<MSink::ScalarPoint>("sink", sinkPar);


    // pion contraction
    MContraction::Meson::Par mesPar;
    mesPar.output   = "mesons/pt_ll";
    mesPar.q1       = "Qpt_l";
    mesPar.q2       = "Qpt_l";
    mesPar.gammas   = "all";
    mesPar.sink     = "sink";
    application.createModule<MContraction::Meson>("meson_pt_ll", mesPar);


    // sink to spin-color matrix
    MSink::Point::Par scSinkPar;
    scSinkPar.mom = "0 0 0";
    application.createModule<MSink::SCMatPoint>("scSink1", scSinkPar);
    application.createModule<MSink::SCMatPoint>("scSink2", scSinkPar);

    // my pion contraction
    MContraction::Pion::Par pionPar;
    pionPar.output = "mesons/pion";
    pionPar.q1 = "Qpt_l";
    pionPar.q2 = "Qpt_l";
    pionPar.sink1 = "scSink1";
    application.createModule<MContraction::Pion>("pion_my", pionPar);


    // pi-pi contraction
    MContraction::PiPi::Par pipiPar;
    pipiPar.output = "pipi/pt_llll";
    pipiPar.q1 = "Qpt_l";
    pipiPar.q2 = "Qpt_l";
    pipiPar.q3 = "Qpt_l";
    pipiPar.q4 = "Qpt_l";
    pipiPar.sink1 = "scSink1";
    pipiPar.sink2 = "scSink2";
    application.createModule<MContraction::PiPi>("pipi_pt_llll", pipiPar);








    // execution ///////////////////////////////////////////////////////////////
    try
    {
        //application.saveParameterFile("spectrum.xml");
        application.run();
    }
    catch (const std::exception& e)
    {
        Exceptions::abort(e);
    }

    // epilogue ////////////////////////////////////////////////////////////////
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();

    return EXIT_SUCCESS;
}
