// by Chris Culver
// Copying and modifying as appropriate Meson.cpp
// this will compute the correlation function for pion-pion scattering
// in isospin 2.


#ifndef _pipi_I2_hpp_
#define _pipi_I2_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

BEGIN_MODULE_NAMESPACE(MContraction)

class PiPiPar: Serializable
{
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(PiPiPar,
                                        std::string, q1,
                                        std::string, q2,
                                        std::string, q3,
                                        std::string, q4,
                                        std::string, output);
};


template <typename FImpl1, typename FImpl2, typename FImpl3, typename FImpl4>
class TPiPi: public Module<PiPiPar>
{
    public:
        FERM_TYPE_ALIASES(FImpl1, 1);
        FERM_TYPE_ALIASES(FImpl2, 2);
        FERM_TYPE_ALIASES(FImpl3, 3);
        FERM_TYPE_ALIASES(FImpl4, 4);
        BASIC_TYPE_ALIASES(ScalarImplCR, Scalar);
        SINK_TYPE_ALIASES(Scalar);

        class Result: Serializable
        {
            public:
                GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                                std::vector<Complex>, corr);
        };

    public:
        TPiPi(const std::string name);
        virtual ~TPiPi(void) {};

        virtual std::vector<std::string> getInput(void);
        virtual std::vector<std::string> getOutput(void);
        virtual std::vector<std::string> getOutputFiles(void);
    protected:
        virtual void setup(void);
        virtual void execute(void);
};

MODULE_REGISTER_TMP(PiPi, ARG(TPiPi<FIMPL, FIMPL, FIMPL, FIMPL>), MContraction);

template <typename FImpl1, typename FImpl2, typename FImpl3, typename FImpl4>
TPiPi<FImpl1, FImpl2, FImpl3, FImpl4>::TPiPi(const std::string name)
: Module<PiPiPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3, typename FImpl4>
std::vector<std::string> TPiPi<FImpl1, FImpl2, FImpl3, FImpl4>::getInput(void)
{
    std::vector<std::string> input = {par().q1, par().q2, par().q3, par().q4};

    return input;
}

template <typename FImpl1, typename FImpl2, typename FImpl3, typename FImpl4>
std::vector<std::string> TPiPi<FImpl1, FImpl2, FImpl3, FImpl4>::getOutput(void)
{
    std::vector<std::string> output = {};

    return output;
}

template <typename FImpl1, typename FImpl2, typename FImpl3, typename FImpl4>
std::vector<std::string> TPiPi<FImpl1, FImpl2, FImpl3, FImpl4>::getOutputFiles(void)
{
    std::vector<std::string> output = {resultFilename(par().output)};

    return output;
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3, typename FImpl4>
void TPiPi<FImpl1, FImpl2, FImpl3, FImpl4>::setup(void)
{
    envTmpLat(LatticeComplex, "c");
}



#define d0(q1, q2, q3, q4) \
((g5*q1)*(g5*q2)*(adj(g5)*adj(q3))*(adj(g5)*adj(q4)))

#define d1(q1, q2) \
(g5*q1*adj(g5)*adj(q2))


template <typename FImpl1, typename FImpl2, typename FImpl3, typename FImpl4>
void TPiPi<FImpl1, FImpl2, FImpl3, FImpl4>::execute(void)
{
    LOG(Message) << "Computing pipi contractions '" << getName() << "' using"
                 << " quarks '" << par().q1 << "', '" << par().q2
                 << "', '" << par().q3 << "' and '" << par().q4 << std::endl;


    Result result;
    Gamma g5(Gamma::Algebra::Gamma5);
    int nt=env().getDim(Tp);
    std::vector<TComplex> buf;

    result.corr.resize(nt);

    if( envHasType(SlicedPropagator1, par().q1) and
        envHasType(SlicedPropagator2, par().q2) and
        envHasType(SlicedPropagator3, par().q3) and
        envHasType(SlicedPropagator4, par().q4))
    {
        auto &q1 = envGet(SlicedPropagator1, par().q1);
        auto &q2 = envGet(SlicedPropagator2, par().q2);
        auto &q3 = envGet(SlicedPropagator3, par().q3);
        auto &q4 = envGet(SlicedPropagator4, par().q4);

        LOG(Message) << "[CC](propagator already sinked)" << std::endl;

        for(unsigned int t=0; t<nt; ++t)
        {
            result.corr[t] = TensorRemove(
                               - trace(d0(q1[t], q2[0], q3[t], q4[0]))
                               +trace(d1(q1[t], q2[0    ]))*trace(d1(q3[t], q4[0]))
                                         );
        }

    }
    else
    {
        LOG(Message) << "Error: Only use sliced propagators!";
        /*
        auto &q1 = envGet(PropagatorField1, par().q1);
        auto &q2 = envGet(PropagatorField2, par().q2);
        auto &q3 = envGet(PropagatorField3, par().q3);
        auto &q4 = envGet(PropagatorField4, par().q4);

        envGetTmp(LatticeComplex, c);
        LOG(Message) << "[CC](using sink1 '" << par().sink1 << "')" << std::endl;
        LOG(Message) << "[CC](using sink2 '" << par().sink2 << "')" << std::endl;
        std::string ns1, ns2;

        ns1=vm().getModuleNamespace(env().getObjectModule(par().sink1));
        ns2=vm().getModuleNamespace(env().getObjectModule(par().sink2));

        if( (ns1=="MSource") && (ns2=="MSource") )
        {
            LOG(Message) << "[CC]Error: Haven't implemented what to do if my sinks are sources?";
        }
        else if( (ns1=="MSink") && (ns2=="MSink") )
        {
            SinkFnScalar &sink1 = envGet(SinkFnScalar, par().sink1);
            SinkFnScalar &sink2 = envGet(SinkFnScalar, par().sink2);

            c = -trace(d0(q1, q2, q3, q4))
               +trace(d1(q1, q2))*trace(d1(q3, q4));
            buf = sink1(c);
        }
        LOG(Message) << "[CC]One sink was successfully applied";
        for(unsigned int t=0; t<buf.size(); ++t)
        {
            result.corr[t] = TensorRemove(buf[t]);
        }
        */
    }
    saveResult(par().output, "pipi", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif
