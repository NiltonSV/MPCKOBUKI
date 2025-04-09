#ifndef PRED_CONTROL_H
#define PRED_CONTROL_H

#include <OsqpEigen/OsqpEigen.h>
#include "math.h"

struct Obstacle {
    double x;
    double y;
    double radius;
};

class PredControl
{
public:

    PredControl();
    ~PredControl() = default;

    virtual void UpdateOptimizationProblem(Eigen::MatrixXd &H,
                                           Eigen::MatrixXd &F,
                                           Eigen::MatrixXd &Ain,
                                           Eigen::VectorXd &lowerBound,
                                           Eigen::VectorXd &upperBound);

    virtual void UpdatePredictionModel();

    virtual void UpdateStates(const Eigen::VectorXd &x);

    virtual void SetInternalVariables();

    virtual void ResizeMatrices();

    virtual void SetWeightMatrices(const Eigen::MatrixXd &Q, const Eigen::MatrixXd &R);

    virtual void SetReference(int i, const Eigen::MatrixXd &ref);

    virtual void DefinePhi();

    virtual void DefineConstraintMtxs();

    virtual void SetConstants(double _ts, int _N, int _M);

    virtual void SetConsBounds(const Eigen::MatrixXd &Lb, const Eigen::MatrixXd &Ub, const std::vector<Obstacle>& obstacles);

    virtual void UpdateObstacleConstraints(const std::vector<Obstacle>& obstacles, int iter);

    void ConfPO();

    void ClearPO();

    void ClearData();

    void ResetPO();

    int SolvePO();

    Eigen::MatrixXd A, B, C, C_cons, C_consV, G, Phi, G_cons, Phi_cons, Q_block, R_block, Uc_block, Lc_block, Ref_block, aux_mdl, aux_cons, Ucv_block, Lcv_block, Lcv_var, Ucv_var, H, F, Ain;
    
    Eigen::MatrixXd C_consStates, G_consStates, Phi_consStates;

    Eigen::VectorXd x, Uc, Lc, Ub, Lb, QPSolution;
    Eigen::SparseMatrix<double> hessian_sparse, linearMatrix;
    OsqpEigen::Solver solver;
    Eigen::Matrix<double, 2, 1> delta_qr;

    int N, M;
    double ts;
    int nx, nxa, ny, nu, nc, ncv, nc_vel, nc_states;
    bool debug = false;
    bool constraints = 0;
    bool error_flag = false;
    int last_op = -1;
    double obj_val;
    bool first_conf = 0;

};
#endif