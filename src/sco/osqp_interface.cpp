#include "utils/logging.hpp"
#include "osqp_interface.hpp"
#include "osqp.h"

using namespace std;

/**
 * TODO: Add description here.
 */
double OSQP_BIG = 1e+30;

namespace sco {

extern void simplify2(vector<int>& inds, vector<double>& vals);
extern vector<int> vars2inds(const vector<Var>& vars);
extern vector<int> cnts2inds(const vector<Cnt>& cnts);

ModelPtr createOSQPModel() {
    ModelPtr out(new OSQPModel());
    return out;
}

OSQPModel::OSQPModel(){
}

OSQPModel::~OSQPModel() {
}

Var OSQPModel::addVar(const string& name) {
    m_vars.push_back(new VarRep(m_vars.size(), name, this));
    m_lbs.push_back(-OSQP_BIG);
    m_ubs.push_back(OSQP_BIG);
    return m_vars.back();
}
Cnt OSQPModel::addEqCnt(const AffExpr& expr, const string& name) {
    m_cnts.push_back(new CntRep(m_cnts.size(), this));
    m_cntExprs.push_back(expr);
    m_cntTypes.push_back(EQ);
    return m_cnts.back();
}
Cnt OSQPModel::addIneqCnt(const AffExpr& expr, const string& name) {
    m_cnts.push_back(new CntRep(m_cnts.size(), this));
    m_cntExprs.push_back(expr);
    m_cntTypes.push_back(INEQ);
    return m_cnts.back();
}
Cnt OSQPModel::addIneqCnt(const QuadExpr&, const string& name) {
    assert( 0 && "NOT IMPLEMENTED");
    return 0;
}
void OSQPModel::removeVars(const VarVector& vars) {
    vector<int>inds = vars2inds(vars);
    for (int i=0; i < vars.size(); ++i) vars[i].var_rep->removed = true;
}

void OSQPModel::removeCnts(const vector<Cnt>& cnts) {
    vector<int>inds = cnts2inds(cnts);
    for (int i=0; i < cnts.size(); ++i) cnts[i].cnt_rep->removed = true;
}

void OSQPModel::update() {
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

void OSQPModel::setVarBounds(const vector<Var>& vars, const vector<double>& lower, const vector<double>& upper) {
    for (int i=0; i < vars.size(); ++i) {
        int varind = vars[i].var_rep->index;
        m_lbs[varind] = lower[i];
        m_ubs[varind] = upper[i];
    }
}
vector<double> OSQPModel::getVarValues(const VarVector& vars) const {
    vector<double> out(vars.size());
    for (int i=0; i < vars.size(); ++i) {
        int varind = vars[i].var_rep->index;
        out[i] = m_soln[varind];
    }
    return out;
}

void OSQPModel::setObjective(const AffExpr& expr) {
    m_objective.affexpr = expr;
}
void OSQPModel::setObjective(const QuadExpr& expr) {
    m_objective = expr;
}
void OSQPModel::writeToFile(const string& fname) {
    // assert(0 && "NOT IMPLEMENTED");
}
VarVector OSQPModel::getVars() const {
    return m_vars;
}

CvxOptStatus OSQPModel::optimize() {
    update();
    // Load problem data
    int n = m_vars.size();
    int m = m_cnts.size();

    vector<int> P_i, P_p(n+1), A_i, A_p(n+1);
    vector<double> P_x, q(n,0), A_x, lbound(m+n), ubound(m+n);

    for (int iVar=0; iVar < n; ++iVar) {
        lbound[m+iVar] = fmax(m_lbs[iVar], -OSQP_BIG);
        ubound[m+iVar] = fmin(m_ubs[iVar], OSQP_BIG);
    }

    vector< vector<int> > var2cntinds(n);
    vector< vector<double> > var2cntvals(n);
    for (int iCnt=0; iCnt < m; ++iCnt) {
        const AffExpr& aff = m_cntExprs[iCnt];
        // cout << "adding constraint " << aff << endl;
        vector<int> inds = vars2inds(aff.vars);

        for (int i=0; i < aff.vars.size(); ++i) {
            var2cntinds[inds[i]].push_back(iCnt);
            var2cntvals[inds[i]].push_back(aff.coeffs[i]); // xxx maybe repeated/
        }

        lbound[iCnt] = (m_cntTypes[iCnt] == INEQ) ? -OSQP_BIG : -aff.constant;
        ubound[iCnt] = -aff.constant;

    }

    for (int iVar=0; iVar < n; ++iVar) {
        // TODO: check if var2cntvals contains 0's...
        simplify2(var2cntinds[iVar], var2cntvals[iVar]);

        // A_x : A matrix nonzero values
        // A_i : A matrix row index
        // A_p : A matrix col, pointer into vector A_x
        A_p[iVar] = A_x.size();
        A_x.insert(A_x.end(), var2cntvals[iVar].begin(), var2cntvals[iVar].end());
        A_i.insert(A_i.end(), var2cntinds[iVar].begin(), var2cntinds[iVar].end());

        // Add individual variable to A to apply lower and upper bounds
         A_x.push_back(1);
         A_i.push_back(m+iVar);
    }
    A_p[n] = A_x.size();

    vector< vector<double> > var2qcoeffs(n);
    vector< vector<int> > var2qinds(n);
    for (int i=0; i < m_objective.size(); ++i) {
        int idx1 = m_objective.vars1[i].var_rep->index, idx2 = m_objective.vars2[i].var_rep->index;
        var2qinds[idx1].push_back(idx2);
        var2qcoeffs[idx1].push_back(m_objective.coeffs[i]);
        var2qinds[idx2].push_back(idx1);
        var2qcoeffs[idx2].push_back(m_objective.coeffs[i]);
    }

    for (int iVar=0; iVar < n; ++iVar) {
        simplify2(var2qinds[iVar], var2qcoeffs[iVar]);

        // P_x : A matrix nonzero values
        // P_i : A matrix row index
        // P_p : A matrix col, pointer into vector A_x
        P_p[iVar] = P_x.size();
        P_x.insert(P_x.end(), var2qcoeffs[iVar].begin(), var2qcoeffs[iVar].end());
        P_i.insert(P_i.end(), var2qinds[iVar].begin(), var2qinds[iVar].end());
    }
    P_p[n] = P_x.size();

    for (int i=0; i < m_objective.affexpr.size(); ++i) {
        q[m_objective.affexpr.vars[i].var_rep->index] += m_objective.affexpr.coeffs[i];
    }


    // Problem settings
    OSQPSettings * settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));

    // Structures
    OSQPWorkspace * work;  // Workspace
    OSQPData * data;  // OSQPData


    // Populate data
    data = (OSQPData *)c_malloc(sizeof(OSQPData));
    data->n = n;
    data->m = m+n;
    data->P = csc_matrix(data->n, data->n, P_x.size(), P_x.data(), P_i.data(), P_p.data());
    data->q = q.data();
    data->A = csc_matrix(data->m, data->n, A_x.size(), A_x.data(), A_i.data(), A_p.data());
    data->l = lbound.data();
    data->u = ubound.data();


    // Define Solver settings as default
    set_default_settings(settings);
    settings->rho = 1e-3;
    settings->max_iter = 10000;
    if (util::GetLogLevel() < util::LevelDebug) {
        settings->verbose = 0;
    }

    // Setup workspace
    work = osqp_setup(data, settings);

    // Solve Problem
    osqp_solve(work);

    CvxOptStatus retcode;
    if (work->info->status_val == OSQP_SOLVED){
        retcode = CVX_SOLVED;
        m_soln = vector<double>(work->solution->x, work->solution->x+n);
        LOG_DEBUG("solver objective value: %.3e", m_objective.value(m_soln));
    }
    else if (work->info->status_val == OSQP_PRIMAL_INFEASIBLE ||
             work->info->status_val == OSQP_DUAL_INFEASIBLE) retcode = CVX_INFEASIBLE;
    else retcode = CVX_FAILED;

    // Cleanup
    osqp_cleanup(work);
    c_free(data->A);
    c_free(data->P);
    c_free(data);
    c_free(settings);

    return retcode;
}

}

