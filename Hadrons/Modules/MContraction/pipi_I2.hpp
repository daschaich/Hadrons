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
                                        std::string, sink1,
                                        std::string, sink2,
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

        typedef Lattice<iScalar<iMatrix<iMatrix<vComplex,Nc>,Ns>>> SCFieldMat;
        typedef std::vector<SCFieldMat::scalar_object> SlicedPropagatorSCMat;
        typedef std::function<SlicedPropagatorSCMat (const SCFieldMat &)> SinkFnSCMat;

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
    std::vector<std::string> input = {par().q1, par().q2, par().q3, par().q4,
                                      par().sink1, par().sink2};

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
    //envTmpLat(SpinColorMatrixField, "tmp1");
    //envTmpLat(SpinColorMatrixField, "tmp1");
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
    std::vector<TComplex> buf1, buf2;
    std::vector<TComplex> p1, p2;

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
                               - trace(d0(q1[t], q2[t], q3[t], q4[t]))
                               +trace(d1(q1[t], q2[t]))*trace(d1(q3[t], q4[t]))
                                         );
        }

    }
    else
    {
        //LOG(Message) << "Error: Only use sliced propagators!";

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
            PropagatorField1 &sink1 = envGet(PropagatorField1, par().sink1);
            PropagatorField2 &sink2 = envGet(PropagatorField1, par().sink2);

            //auto p1 = sliceInnerProductVector(buf1, q1*adj(q2)*sink1, q3*adj(q4)*sink2, Tp);
            //auto p2 = sliceSum(q3*adj(q4)*sink2, buf2, Tp);
            //auto d1 = trace(p1*p2);
            //auto d2 = trace(p1)*trace(p2)

        /*
            //single sum over x - pi(x)pi(x) ?
            auto d1 = trace(q1*adj(q2)*q3*adj(q4)*sink1*sink2);
            auto d2_0 = trace(q1*adj(q2)*sink1);
            auto d2_1 = trace(q3*adj(q4)*sink2);
            c = -d1 + d2_0*d2_1;
            sliceSum(c, buf, Tp);
        */
        }
        else if( (ns1=="MSink") && (ns2=="MSink") )
        {
            SinkFnSCMat &sink1 = envGet(SinkFnSCMat, par().sink1);
            SinkFnSCMat &sink2 = envGet(SinkFnSCMat, par().sink2);

            auto tmp1 = g5*adj(q1)*g5*adj(g5)*q2*g5;
            auto tmp2 = g5*adj(q3)*g5*adj(g5)*q4*g5;
            auto s1 = sink1(tmp1);
            auto s2 = sink2(tmp2);


            for(size_t t=0; t<s1.size(); ++t)
            {
                auto d1 = trace(s1[t]*s2[t]);
                auto d2 = trace(s1[t])*trace(s2[t]);
                //result.d1[t] = TensorRemove(d1);
                //result.d2[t] = TensorRemove(d2);
                result.corr[t] = TensorRemove(-d1+d2);
            }


        }
        //LOG(Message) << "[CC]One sink was successfully applied";
        //for(unsigned int t=0; t<buf.size(); ++t)
        //{
            //auto d1 = trace(p1[t]*p2[t]);
            //auto d2 = trace(p1[t])*trace(p2[t]);
            //result.corr[t] = TensorRemove(buf[t]);
            //result.corr[t] = TensorRemove(-d1+d2);
        //}

    }
    saveResult(par().output, "pipi", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif
