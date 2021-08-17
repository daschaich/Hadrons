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
  std::string parameterFileName;
  if(argc<2)
  {
    std::cerr << "usage: " << argv[0] << " <parameter file> [Grid options]";
    std::cerr << std::endl;
    std::exit(EXIT_FAILURE);
  }
  parameterFileName = argv[1];

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

    {
      XmlReader reader(parameterFileName);
      read(reader, "global", globalPar);
    }

//      read(reader, "inverterMass", inverterMass)


    // global initialisation
    application.setPar(globalPar);
    /// create modules //////////////////////////////////////////////////////////

    //for remaining parameters
    XmlReader reader(parameterFileName);

    // gauge field
    MIO::LoadNersc::Par gaugePar;
    read(reader, "configFilePrefix", gaugePar.file);
    application.createModule<MIO::LoadNersc>("gauge", gaugePar);
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


    //action // TODO: MAKE CORRECT
    MAction::ScaledDWF::Par actionPar;
    actionPar.gauge = "stout";
    actionPar.Ls = 16;
    actionPar.M5 = 1.0;
    actionPar.scale=2.0;
    read(reader, "inverterMass", actionPar.mass);
    actionPar.boundary = "1 1 1 -1";
    application.createModule<MAction::ScaledDWF>("DWF_outer_" + flavor, actionPar);

    actionPar.gauge = "stoutf";
    application.createModule<MAction::ScaledDWFF>("DWF_inner_" + flavor, actionPar);


    //solver
    MSolver::MixedPrecisionRBPrecCG::Par solverPar;
    solverPar.innerAction = "DWF_inner_" + flavor;
    solverPar.outerAction = "DWF_outer_" + flavor;
    read(reader, "solverResidual", solverPar.residual);
    read(reader, "solverMaxInnerIter",solverPar.maxInnerIteration);
    read(reader, "solverMaxOuterIter",solverPar.maxOuterIteration);
    application.createModule<MSolver::MixedPrecisionRBPrecCG>("CG_" + flavor, solverPar);

    /// sink
    MSink::Point::Par sinkPar;
    sinkPar.mom = "0 0 0";
    application.createModule<MSink::ScalarPoint>("sink", sinkPar);

    // sink to spin-color matrix
    MSink::Point::Par scSinkPar;
    scSinkPar.mom = "0 0 0";
    application.createModule<MSink::SCMatPoint>("scSink1", scSinkPar);
    application.createModule<MSink::SCMatPoint>("scSink2", scSinkPar);

    Application::ObjectId id;
    if(!push(reader, "sources"))
    {
        HADRONS_ERROR(Parsing, "Cannot open node 'sources' in parameter file'"
                               + parameterFileName + "'");
    }
    std::vector<std::string> sourceNames;
    if(!push(reader, "source"))
    {
        HADRONS_ERROR(Parsing, "Cannot open node 'sources/source' in parameter file '"
                            + parameterFileName + "'");
    }
    do
    {
        read(reader, "id", id);
        application.createModule(id.name, id.type, reader);
        sourceNames.push_back(id.name);
    } while(reader.nextElement("source"));
    pop(reader);
    pop(reader);

    for(const auto &source: sourceNames)
    {
        MFermion::GaugeProp::Par quarkPar;
        quarkPar.solver = "CG_l";
        quarkPar.source = source;
        auto propName = "Q"+source+"_l";
        application.createModule<MFermion::GaugeProp>(propName, quarkPar);

        MContraction::Meson::Par mesPar;
        mesPar.output   = "mesons/"+source+"_ll";
        mesPar.q1       = propName;
        mesPar.q2       = propName;
        mesPar.gammas   = "all";
        mesPar.sink     = "sink";
        application.createModule<MContraction::Meson>("meson_"+source+"_ll", mesPar);

        MContraction::PiPi::Par pipiWallPar;
        pipiWallPar.output = "pipi/"+source+"_llll";
        pipiWallPar.q1 = propName;
        pipiWallPar.q2 = propName;
        pipiWallPar.q3 = propName;
        pipiWallPar.q4 = propName;
        pipiWallPar.sink1 = "scSink1";
        pipiWallPar.sink2 = "scSink2";
        application.createModule<MContraction::PiPi>("pipi_"+source+"_llll", pipiWallPar);
    }

    // source
    //MSource::Point::Par ptPar;
    //ptPar.position = "0 0 0 0";                                 //TODO read from XML
    //application.createModule<MSource::Point>("pt", ptPar);

    // Wall source
    //MSource::Wall::Par wallPar;
    //wallPar.tW = 0;                                             //TODO read from XML
    //wallPar.mom = "0. 0. 0.";
    //application.createModule<MSource::Wall>("wall", wallPar);

    /// read in option to create pt, wall, or both.
    //prop
    //MFermion::GaugeProp::Par quarkPar;
    //quarkPar.solver = "CG_" + flavor;
    //quarkPar.source = "pt";
    //application.createModule<MFermion::GaugeProp>("Qpt_" + flavor, quarkPar);



    //wall prop
    //MFermion::GaugeProp::Par wallquarkPar;
    //wallquarkPar.solver = "CG_" + flavor;
    //wallquarkPar.source = "wall";
    //application.createModule<MFermion::GaugeProp>("Qwall_" + flavor, wallquarkPar);

/*
    // use flag from above to switch on/off pt wall contractions
    // pion contraction
    MContraction::Meson::Par mesPar;
    mesPar.output   = "mesons/pt_ll";
    mesPar.q1       = "Qpt_l";
    mesPar.q2       = "Qpt_l";
    mesPar.gammas   = "all";
    mesPar.sink     = "sink";
    application.createModule<MContraction::Meson>("meson_pt_ll", mesPar);
*/


/*
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

    // pi-pi wall sources contraction
    MContraction::PiPi::Par pipiWallPar;
    pipiWallPar.output = "pipi/wall_llll";
    pipiWallPar.q1 = "Qwall_l";
    pipiWallPar.q2 = "Qwall_l";
    pipiWallPar.q3 = "Qwall_l";
    pipiWallPar.q4 = "Qwall_l";
    pipiWallPar.sink1 = "scSink1";
    pipiWallPar.sink2 = "scSink2";
    application.createModule<MContraction::PiPi>("pipi_wall_llll", pipiWallPar);
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
