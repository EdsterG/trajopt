#include "utils/logging.hpp"
#include "qpoases_interface.hpp"
#include "qpOASES.hpp"

using namespace std;

/**
 * TODO: Add description here.
 */
double QPOASES_BIG = 1e+30;

namespace sco {

extern void simplify2(vector<int>& inds, vector<double>& vals);
extern vector<int> vars2inds(const vector<Var>& vars);
extern vector<int> cnts2inds(const vector<Cnt>& cnts);

ModelPtr createQPOASESModel() {
    ModelPtr out(new QPOASESModel());
    return out;
}

QPOASESModel::QPOASESModel(){
}

QPOASESModel::~QPOASESModel() {
}

Var QPOASESModel::addVar(const string& name) {
    m_vars.push_back(new VarRep(m_vars.size(), name, this));
    m_lbs.push_back(-QPOASES_BIG);
    m_ubs.push_back(QPOASES_BIG);
    return m_vars.back();
}
Cnt QPOASESModel::addEqCnt(const AffExpr& expr, const string& name) {
    m_cnts.push_back(new CntRep(m_cnts.size(), this));
    m_cntExprs.push_back(expr);
    m_cntTypes.push_back(EQ);
    return m_cnts.back();
}
Cnt QPOASESModel::addIneqCnt(const AffExpr& expr, const string& name) {
    m_cnts.push_back(new CntRep(m_cnts.size(), this));
    m_cntExprs.push_back(expr);
    m_cntTypes.push_back(INEQ);
    return m_cnts.back();
}
Cnt QPOASESModel::addIneqCnt(const QuadExpr&, const string& name) {
    assert( 0 && "NOT IMPLEMENTED");
    return 0;
}
void QPOASESModel::removeVars(const VarVector& vars) {
    vector<int>inds = vars2inds(vars);
    for (int i=0; i < vars.size(); ++i) vars[i].var_rep->removed = true;
}

void QPOASESModel::removeCnts(const vector<Cnt>& cnts) {
    vector<int>inds = cnts2inds(cnts);
    for (int i=0; i < cnts.size(); ++i) cnts[i].cnt_rep->removed = true;
}

void QPOASESModel::update() {
    {
        unsigned int inew = 0;
        for (int iold=0; iold < m_vars.size(); ++iold) {
            const Var& var = m_vars[iold];
            if (!var.var_rep->removed) {
                m_vars[inew] = var;
                m_lbs[inew] = m_lbs[iold];
                m_ubs[inew] = m_ubs[iold];
                var.var_rep->index = inew;
                ++inew;
            }
            else delete var.var_rep;
        }
        m_vars.resize(inew);
        m_lbs.resize(inew);
        m_ubs.resize(inew);
    }
    {
        unsigned int inew = 0;
        for (int iold = 0; iold < m_cnts.size(); ++iold) {
            const Cnt& cnt = m_cnts[iold];
            if (!cnt.cnt_rep->removed) {
                m_cnts[inew] = cnt;
                m_cntExprs[inew] = m_cntExprs[iold];
                m_cntTypes[inew] = m_cntTypes[iold];
                cnt.cnt_rep->index = inew;
                ++inew;
            }
            else delete cnt.cnt_rep;
        }
        m_cnts.resize(inew);
        m_cntExprs.resize(inew);
        m_cntTypes.resize(inew);
    }
}

void QPOASESModel::setVarBounds(const vector<Var>& vars, const vector<double>& lower, const vector<double>& upper) {
    for (int i=0; i < vars.size(); ++i) {
        int varind = vars[i].var_rep->index;
        m_lbs[varind] = lower[i];
        m_ubs[varind] = upper[i];
    }
}
vector<double> QPOASESModel::getVarValues(const VarVector& vars) const {
    vector<double> out(vars.size());
    for (int i=0; i < vars.size(); ++i) {
        int varind = vars[i].var_rep->index;
        out[i] = m_soln[varind];
    }
    return out;
}

void QPOASESModel::setObjective(const AffExpr& expr) {
    m_objective.affexpr = expr;
}
void QPOASESModel::setObjective(const QuadExpr& expr) {
    m_objective = expr;
}
void QPOASESModel::writeToFile(const string& fname) {
    // assert(0 && "NOT IMPLEMENTED");
}
VarVector QPOASESModel::getVars() const {
    return m_vars;
}

CvxOptStatus QPOASESModel::optimize() {
    update();
    // Load problem data
    int n = m_vars.size();
    int m = m_cnts.size();

    vector<double> H(n*n,0), A(n*m,0), g(n,0), lb(n), ub(n), lbA(m), ubA(m);

    for (int iVar=0; iVar < n; ++iVar) {
        lb[iVar] = fmax(m_lbs[iVar], -QPOASES_BIG);
        ub[iVar] = fmin(m_ubs[iVar], QPOASES_BIG);
    }

    for (int iCnt=0; iCnt < m; ++iCnt) {
        const AffExpr& aff = m_cntExprs[iCnt];
        vector<int> inds = vars2inds(aff.vars);

        for (int i=0; i < aff.vars.size(); ++i) {
            A[iCnt * m + inds[i]] += aff.coeffs[i];
        }

        lbA[iCnt] = (m_cntTypes[iCnt] == INEQ) ? -QPOASES_BIG : -aff.constant;
        ubA[iCnt] = -aff.constant;
    }

    for (int i=0; i < m_objective.size(); ++i) {
        int idx1 = m_objective.vars1[i].var_rep->index, idx2 = m_objective.vars2[i].var_rep->index;
        H[idx1 + idx2 * n] += m_objective.coeffs[i];
        H[idx2 + idx1 * n] += m_objective.coeffs[i];
    }

    for (int i=0; i < m_objective.affexpr.size(); ++i) {
        g[m_objective.affexpr.vars[i].var_rep->index] += m_objective.affexpr.coeffs[i];
    }

    USING_NAMESPACE_QPOASES

    /* Setting up QProblem object. */
    QProblem env(n,m);

    Options options;
    env.setOptions(options);

    /* Solve first QP. */
    int_t nWSR = 10; // Max allowed CPU time
    env.init(H.data(),g.data(),A.data(),lb.data(),ub.data(),lbA.data(),ubA.data(), nWSR);

    /* Get and print solution of first QP. */
    m_soln = vector<double>(n);
    env.getPrimalSolution(m_soln.data());

    if (env.isSolved()) return CVX_SOLVED;
    else if (env.isInfeasible()) return CVX_INFEASIBLE;
    else return CVX_FAILED;
}

}

