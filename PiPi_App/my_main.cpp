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


    // global initialisation
    application.setPar(globalPar);

    // create modules //////////////////////////////////////////////////////////

    // gauge field
    application.createModule<MGauge::Unit>("gauge");



    // source
    MSource::Point::Par ptPar;
    ptPar.position = "0 0 0 0";
    application.createModule<MSource::Point>("pt", ptPar);
    // sink
    MSink::Point::Par sinkPar;
    sinkPar.mom = "0 0 0";
    application.createModule<MSink::Point>("sink", sinkPar);

    std::string flavor = "l";

    //action
    MAction::DWF::Par actionPar;
    actionPar.gauge = "gauge";
    actionPar.Ls = 12;
    actionPar.M5 = 1.8;
    actionPar.mass = 0.01;
    actionPar.boundary = "1 1 1 -1";
    actionPar.twist = "0. 0. 0. 0.";
    application.createModule<MAction::DWF>("DWF_" + flavor, actionPar);

    //solver
    MSolver::RBPrecCG::Par solverPar;
    solverPar.action = "DWF_" + flavor;
    solverPar.residual = 1.0e-8;
    solverPar.maxIteration = 10000;
    application.createModule<MSolver::RBPrecCG>("CG_" + flavor, solverPar);

    //prop
    MFermion::GaugeProp::Par quarkPar;
    quarkPar.solver = "CG_" + flavor;
    quarkPar.source = "pt";
    application.createModule<MFermion::GaugeProp>("Qpt_" + flavor, quarkPar);

    //mesons
    MContraction::Meson::Par mesPar;
    MesonEntry mesEntry;

    mesPar.output = "mesons/pt_" + flavor + flavor;
    mesPar.q1 = "Qpt_" + flavor;
    mesPar.q2 = "Qpt_" + flavor;
    mesPar.gammas = "Gamma5";
    mesPar.sink = "sink";

    mesEntry.q1 = flavor;
    mesEntry.q2 = flavor;
    mesEntry.source = "pt";
    application.createModule<MContraction::Meson>("meson_pt_" + flavor + flavor, mesPar);
    application.setResultMetadata("meson_pt_" + flavor + flavor, "meson", mesEntry);


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
