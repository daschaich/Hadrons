// by Chris Culver
// Copying and modifying as appropriate Meson.cpp
// this will compute the correlation function for pion-pion scattering
// in isospin 2.


#ifndef _pion_hpp_
#define _pion_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

BEGIN_MODULE_NAMESPACE(MContraction)

class PionPar: Serializable
{
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(PionPar,
                                        std::string, q1,
                                        std::string, q2,
                                        std::string, sink1,
                                        std::string, output);
};


template <typename FImpl1, typename FImpl2>
class TPion: public Module<PionPar>
{
    public:
        FERM_TYPE_ALIASES(FImpl1, 1);
        FERM_TYPE_ALIASES(FImpl2, 2);
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
        TPion(const std::string name);
        virtual ~TPion(void) {};

        virtual std::vector<std::string> getInput(void);
        virtual std::vector<std::string> getOutput(void);
        virtual std::vector<std::string> getOutputFiles(void);
    protected:
        virtual void setup(void);
        virtual void execute(void);
};

MODULE_REGISTER_TMP(Pion, ARG(TPion<FIMPL, FIMPL>), MContraction);

template <typename FImpl1, typename FImpl2>
TPion<FImpl1, FImpl2>::TPion(const std::string name)
: Module<PionPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
std::vector<std::string> TPion<FImpl1, FImpl2>::getInput(void)
{
    std::vector<std::string> input = {par().q1, par().q2,
                                      par().sink1};

    return input;
}

template <typename FImpl1, typename FImpl2>
std::vector<std::string> TPion<FImpl1, FImpl2>::getOutput(void)
{
    std::vector<std::string> output = {};

    return output;
}

template <typename FImpl1, typename FImpl2>
std::vector<std::string> TPion<FImpl1, FImpl2>::getOutputFiles(void)
{
    std::vector<std::string> output = {resultFilename(par().output)};

    return output;
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
void TPion<FImpl1, FImpl2>::setup(void)
{
    envTmpLat(LatticeComplex, "c");
    //envTmpLat(SpinColorMatrixField, "tmp1");
    //envTmpLat(SpinColorMatrixField, "tmp1");
}



template <typename FImpl1, typename FImpl2>
void TPion<FImpl1, FImpl2>::execute(void)
{
    LOG(Message) << "Computing pion contractions '" << getName() << "' using"
                 << " quarks '" << par().q1 << "', '" << par().q2
                 << std::endl;


    Result result;
    Gamma g5(Gamma::Algebra::Gamma5);
    int nt=env().getDim(Tp);
    std::vector<TComplex> buf;
    std::vector<TComplex> buf1, buf2;
    std::vector<TComplex> p1, p2;

    result.corr.resize(nt);
    //LOG(Message) << "Error: Only use sliced propagators!";

    auto &q1 = envGet(PropagatorField1, par().q1);
    auto &q2 = envGet(PropagatorField2, par().q2);

    envGetTmp(LatticeComplex, c);
    LOG(Message) << "[CC](using sink1 '" << par().sink1 << "')" << std::endl;
    std::string ns1;

    ns1=vm().getModuleNamespace(env().getObjectModule(par().sink1));

    SinkFnSCMat &sink1 = envGet(SinkFnSCMat, par().sink1);

    auto tmp1 = g5*adj(q1)*g5*adj(g5)*q2*g5;
    auto s1 = sink1(tmp1);

    for(size_t t=0; t<s1.size(); ++t)
    {
        result.corr[t] = TensorRemove(trace(s1[t]));
    }

    saveResult(par().output, "pion", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif
