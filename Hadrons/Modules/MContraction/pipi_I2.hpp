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
                                                std::vector<Complex>, corr,
                                                std::vector<Complex>, d0,
                                                std::vector<Complex>, d1,
                                                std::vector<Complex>, s0,
                                                std::vector<Complex>, s1,
                                                std::vector<Complex>, cTst);
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
    result.d0.resize(nt);
    result.d1.resize(nt);
    result.s0.resize(nt);
    result.s1.resize(nt);
    result.cTst.resize(nt);
    auto &q1 = envGet(PropagatorField1, par().q1);
    auto &q2 = envGet(PropagatorField2, par().q2);
    auto &q3 = envGet(PropagatorField3, par().q3);
    auto &q4 = envGet(PropagatorField4, par().q4);
    SinkFnSCMat &sink1 = envGet(SinkFnSCMat, par().sink1);
    SinkFnSCMat &sink2 = envGet(SinkFnSCMat, par().sink2);

    auto tmp0 = g5*adj(q1)*g5*adj(g5)*q2*g5;
    auto tmp1 = g5*adj(q3)*g5*adj(g5)*q4*g5;
    auto s0 = sink1(tmp0);
    auto s1 = sink2(tmp1);


    for(size_t t=0; t<s1.size(); ++t)
    {
        auto d0 = trace(s0[t]*s1[t]);
        auto d1 = trace(s0[t])*trace(s1[t]);
        auto dTst = -trace(s0[t]*s0[t])+trace(s0[t])*trace(s0[t]);
        result.s0[t] = TensorRemove(trace(s0[t]));
        result.s1[t] = TensorRemove(trace(s1[t]));
        result.d0[t] = TensorRemove(d0);
        result.d1[t] = TensorRemove(d1);
        result.corr[t] = TensorRemove(-d0+d1);
        result.cTst[t] = TensorRemove(dTst);
    }

    saveResult(par().output, "pipi", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif
