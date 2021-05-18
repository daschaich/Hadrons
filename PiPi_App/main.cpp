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
    application.createModule<MGauge::Unit>("gauge");
    std::string flavor = "l";

    //action
    MAction::DWF::Par actionPar;
    actionPar.gauge = "gauge";
    actionPar.Ls = 12;
    actionPar.M5 = 1.8;
    actionPar.mass = 0.01;
    actionPar.boundary = "1 1 1 -1";
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


    // source
    MSource::Point::Par ptPar;
    ptPar.position = "0 0 0 0";
    application.createModule<MSource::Point>("pt", ptPar);
    // scalar sink
    MSink::ScalarPoint::Par scalarSinkPar;
    scalarSinkPar.mom = "0 0 0";
    application.createModule<MSink::ScalarPoint>("scalarSink", scalarSinkPar);


    MSink::Point::Par ptSinkPar;
    ptSinkPar.mom = "0 0 0";
    application.createModule<MSink::Point>("ptSink", ptSinkPar);


    MSink::Smear::Par smearPar;
    smearPar.q = "Qpt_" + flavor;
    smearPar.sink = "ptSink"; //can't smear with scalarSink
    application.createModule<MSink::Smear>("smearSink", smearPar);

    //mesons
    MContraction::Meson::Par mesPar;
    MesonEntry mesEntry;
    //mesons with sliced propagators
    mesPar.output = "mesons/pt_smeared_" + flavor + flavor;
    mesPar.q1 = "smearSink";
    mesPar.q2 = "smearSink";
    mesPar.gammas = "(Gamma5 Gamma5)";  //all works, how to do one type?
    application.createModule<MContraction::Meson>("meson_pt_smear" + flavor + flavor, mesPar);
    application.setResultMetadata("meson_pt_smear" + flavor + flavor, "meson", mesEntry);

/*
    //point source-scalar sink
    mesPar.output = "mesons/pt_" + flavor + flavor;
    mesPar.q1 = "Qpt_" + flavor;
    mesPar.q2 = "Qpt_" + flavor;
    mesPar.gammas = "(Gamma5 Gamma5)";
    mesPar.sink = "scalarSink";

    mesEntry.q1 = flavor;
    mesEntry.q2 = flavor;
    mesEntry.source = "pt";
    application.createModule<MContraction::Meson>("meson_pt_" + flavor + flavor, mesPar);
    application.setResultMetadata("meson_pt_" + flavor + flavor, "meson", mesEntry);
*/




/*
    //pi-pi with sliced propagators, aka already smeared fields
    MContraction::PiPiPar pipiPar;
    PiPiEntry pipiEntry;

    pipiPar.output = "pipi/pt_" + flavor + flavor + flavor + flavor;
    pipiPar.q1 = "smearSink";
    pipiPar.q2 = "smearSink";
    pipiPar.q3 = "smearSink";
    pipiPar.q4 = "smearSink";

    pipiEntry.q1 = flavor;
    pipiEntry.q2 = flavor;
    pipiEntry.q3 = flavor;
    pipiEntry.q4 = flavor;

    application.createModule<MContraction::PiPi>("pipi_ptpt_" + flavor + flavor + flavor + flavor, pipiPar);
//    application.setResultMetadata("pipi_ptpt_" + flavor + flavor + flavor + flavor, "pipi", pipiEntry);
*/


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
