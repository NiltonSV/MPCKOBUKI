#include "testes/pred_control.h"

PredControl::PredControl() {
}

void PredControl::SetInternalVariables()
{
    this->nx = this->A.rows();                      // # states
    this->nu = this->B.cols();                      // # decision variables
    this->ny = this->C.rows();                      // # outputs for the optimization problem
    this->nxa = this->nx;                           // # aumented states
    this->nc  = this->C_cons.rows();                // # constraints
    this->ncv = this->C_consV.cols();               // # constraints variables
}

void PredControl::ResizeMatrices()
{
    // Model prediction matrices: y̅ = G*u̅ + Phi*x[k]
    this->x.resize(this->nxa, 1);
    this->x.setZero();

    this->G.resize(this->ny * this->N, this->nu * this->M);
    this->G.setZero();

    this->Phi.resize(this->ny * N, this->nxa);
    this->Phi.setZero();

    this->H.resize(this->nu * this->M, this->nu * this->M);
    this->H.setZero();

    this->F.resize(this->nu * this->M, 1);
    this->F.setZero();

    // Auxiliare matrices
    this->aux_mdl.resize(this->ny, this->nu);
    this->aux_mdl.setZero();

    this->aux_cons.resize(this->nc, this->nu);
    this->aux_cons.setZero();

    // Constraints

    this->Ain.resize(this->nc * this->N, this->nu * this->M);
    this->Ain.setZero();

    this->Lb.resize(this->nc * this->N);
    this->Lb.setZero();

    this->Ub.resize(this->nc * this->N);
    this->Ub.setZero();

    this->G_cons.resize(this->nc * this->N, this->nu * this->M);
    this->G_cons.setZero();

    this->Phi_cons.resize(this->nc * N, this->nxa);
    this->Phi_cons.setZero();

    this->Uc_block.resize(this->nc * this->N, 1);
    this->Uc_block.setZero();

    this->Lc_block.resize(this->nc * this->N, 1);
    this->Lc_block.setZero();

    this->Lcv_block.resize(this->nc * this->N, this->ncv);
    this->Ucv_block.resize(this->nc * this->N, this->ncv);
    
    this->Lcv_var.resize(this->ncv, 1);
    this->Ucv_var.resize(this->ncv, 1);

    // Weight matrices: The cost function must be writen as J = (r̅-y̅)'*Q_blcok*(r̅-y̅) + u̅'*R_block*u̅

    // this->Ref_block.resize(this->ny * this->N, 1);
    this->Ref_block.resize(this->ny * this->N, 1);
    this->Ref_block.setZero();

    this->Q_block.resize(this->ny * this->N, this->ny * this->N);
    this->Q_block.setZero();

    this->R_block.resize(this->nu * this->M, this->nu * this->M);
    this->R_block.setZero();
}

void PredControl::SetConstants(double _ts, int _N, int _M)
{
    this->ts = _ts;
    this->N = _N;
    this->M = _M;

    this->ResizeMatrices();
}

void PredControl::SetWeightMatrices(const Eigen::MatrixXd &Q, const Eigen::MatrixXd &R)
{
    if (Q.rows() != this->ny || Q.cols() != this->ny)
    {
        std::cout << "Wrong Q dimension!" << std::endl;
        std::cout << "Expected: " << this->ny << "x" << this->ny << " - Received: " << Q.rows() << "x" << Q.cols() << std::endl;

        return;
    }

    if (R.rows() != this->nu || R.cols() != this->nu)
    {
        std::cout << "Wrong R dimension!" << std::endl;
        std::cout << "Expected: " << this->nu << "x" << this->nu << " - Received: " << R.rows() << "x" << R.cols() << std::endl;
        return;
    }

    for (int i = 0; i < this->N; i++)
    {

        // this->Q_block.block(i * this->ny, i * this->ny, this->ny, this->ny) = Q;
        if (i < N-1){
            this->Q_block.block(i * this->ny, i * this->ny, this->ny, this->ny) = pow(2, i)*Q;
        }
        else{
            this->Q_block.block(i * this->ny, i * this->ny, this->ny, this->ny) = pow(2, i)*30*Q;
        }
        if (i < this->M)
            this->R_block.block(i * this->nu, i * this->nu, this->nu, this->nu) = R;
    }
}

void PredControl::SetConsBounds(const Eigen::MatrixXd &Lb, const Eigen::MatrixXd &Ub, const std::vector<Obstacle>& obstacles)
{
    if (!this->nc)
        return;

    if ((Lb.rows() != this->nc_vel || Lb.cols() != 1))
    {
        std::cout << "Wrong Lower Bound dimension" << std::endl;
        std::cout << "Expected: " << this->nc << "x1" << " - Received: " << Lb.rows() << "x" << Lb.cols() << std::endl;

        return;
    }

    if (Ub.rows() != this->nc_vel || Ub.cols() != 1)
    {
        std::cout << "Wrong Upper Bound dimension" << std::endl;
        std::cout << "Expected: " << this->nc << "x1" << " - Received: " << Ub.rows() << "x" << Ub.cols() << std::endl;
        return;
    }

    this->Lc = Lb;
    this->Uc = Ub;

    for (int i = 0; i < this->N; i++)
    {
        if(this->nc_vel)
        {
            this->Lc_block.block(i * this->nc, 0, this->nc_vel, 1) = this->Lc;
            this->Uc_block.block(i * this->nc, 0, this->nc_vel, 1) = this->Uc;
        }
        this->UpdateObstacleConstraints(obstacles, i);
        this->updateStatesConstraints(i);
        this->Lcv_block.block(i * this->nc, 0, this->nc, this->ncv) = this->C_consV;
        this->Ucv_block.block(i * this->nc, 0, this->nc, this->ncv) = this->C_consV;
    }
    
}

void PredControl::SetReference(int i, const Eigen::MatrixXd &ref)
{
    if (ref.rows() != this->ny || ref.cols() != 1)
    {
        std::cout << "Wrong Ref dimension" << std::endl;
        return;
    }

    // for (int i = 0; i < this->N; i++)
    // {
    // this->Ref_block << ref;
    this->Ref_block.block(i * this->nxa, 0, this->nxa, 1) = ref;
    // }
}

void PredControl::ResetPO()
{
    if (this->solver.isInitialized())
        this->ClearPO();
    if (this->solver.data()->isSet())
        this->ClearData();
    this->last_op = -1;
}

void PredControl::ClearPO()
{
    this->solver.clearSolverVariables();
    this->solver.clearSolver();
}

void PredControl::ClearData()
{
    this->solver.data()->clearLinearConstraintsMatrix();
    this->solver.data()->clearHessianMatrix();
}

void PredControl::UpdateStates(const Eigen::VectorXd &x)
{
    this->x = x;
}

void PredControl::UpdateOptimizationProblem(Eigen::MatrixXd &H,
                                            Eigen::MatrixXd &F,
                                            Eigen::MatrixXd &Ain,
                                            Eigen::VectorXd &lowerBound,
                                            Eigen::VectorXd &upperBound)
{

    // define initial Phi values
    this->DefinePhi();

    // update the prediction model
    this->UpdatePredictionModel();

    // Update the optimization problem

    H = 2 * (this->G.transpose() * this->Q_block * this->G + this->R_block);

    // F = 2 * (((this->Phi * this->x) - this->Ref_block).transpose()) * this->Q_block * this->G;
    // F = 2 * this->G.transpose() * this->Q_block * this-> Phi * (this->x - this->Ref_block.block(0, 0, this->nxa, 1));
    F = 2 * this->G.transpose() * this->Q_block * (this->Phi * this->x - this->Ref_block);
    Ain = this->G_cons;

    lowerBound = this->Lc_block - this->Phi_cons * this->x; //+ this->Lcv_block * this->Lcv_var;

    upperBound = this->Uc_block - this->Phi_cons * this->x; //+ this->Ucv_block * this->Ucv_var;
}

void PredControl::DefinePhi()
{
    this->Phi.block(0, 0, this->ny, this->nxa) = this->C * this->A;
    this->aux_mdl = this->C * this->B;
    
}

void PredControl::DefineConstraintMtxs()
{
    // this->Phi_cons.block(0, 0, this->nc, this->nxa).setZero();
    if(this->nc_vel)
    {
        this->Phi_cons.block(0, 0, this->nc, this->nxa).setZero();
        this->aux_cons.setZero();
        this->aux_cons.block(0, 0, 2, this->nu) = Eigen::MatrixXd::Identity(2, this->nu);
    }
    for(int i=this->nc_vel; i < this->nc; ++i)
    {
        this->Phi_cons.block(i, 0, 1, this->nxa) = this->C_cons.block(i, 0, 1, this->nxa) * this->A;
        this->aux_cons.block(i, 0, 1, this->nu) = this->C_cons.block(i, 0, 1, this->nxa) * this->B;
    }
    for(int i=0; i<this->nc_states; ++i)
    {
        this->Phi_cons.block(this->nc-this->nc_states+i, 0, 1, this->nxa) = this->C_cons.block(this->nc-this->nc_states+i, 0, 1, this->nxa) * this->A;
        this->aux_cons.block(this->nc-this->nc_states+i, 0, 1, this->nu) = this->C_cons.block(this->nc-this->nc_states+i, 0, 1, this->nxa) * this->B;
    }
}

void PredControl::updateStatesConstraints(int iter)
{
    if (!this->nc_states) {
        return;
    }


    for (int i = 0; i<this->nc_states; ++i)
    {
        double d0 = sqrt((this->x(0))*(this->x(0))+(this->x(1)+1.25)*(this->x(1)+1.25));
        // if (d0 > 0.3){
        //     this->C_cons.block(this->nc - (this->nc_states - i), 1, 1, 1) << 1;
        //     this->Lc_block((iter+1)*this->nc-(this->nc_states-i)) = -OsqpEigen::INFTY;
        //     this->Uc_block((iter+1)*this->nc-(this->nc_states-i)) = OsqpEigen::INFTY;
        //     return;
        // }
        this->C_cons.block(this->nc - (this->nc_states - i), 1, 1, 1) << 1;
        this->Lc_block((iter+1)*this->nc-(this->nc_states-i)) = -0.80;
        this->Uc_block((iter+1)*this->nc-(this->nc_states-i)) = OsqpEigen::INFTY; //-0.95;
    }

    // this->C_cons.block(this->nc-this->nc_states - 1, 1, 1, 1) << 1;
    // this->Lc_block((iter+1)*this->nc-1) = -0.05+this->current_iter;
    // this->Uc_block((iter+1)*this->nc-1) = 0.05+this->current_iter;
}
void PredControl::UpdateObstacleConstraints(const std::vector<Obstacle>& obstacles, int iter) 
{
    if (obstacles.empty()) {
        return;
    }

    for (size_t i = 0; i < obstacles.size(); i++) {
        double dx = this->x(0) - obstacles[i].x;
        double dy = this->x(1) - obstacles[i].y;
        double current_dist = sqrt(dx*dx + dy*dy);
        
        if (current_dist < 1e-6) {
            dx = 0.01;
            dy = 0.01;
            current_dist = sqrt(dx*dx + dy*dy);
        }

        double grad_x = dx / current_dist;
        double grad_y = dy / current_dist;

        this->C_cons.block(i+this->nc_vel, 0, 1, 2) << grad_x, grad_y;
        
        double safety_distance = obstacles[i].radius + 0.01 + 0.01;

        this->Lc_block(iter * this->nc + i + this->nc_vel) = safety_distance - current_dist + grad_x*this->x(0) + grad_y*this->x(1);
        this->Uc_block(iter * this->nc + i + this->nc_vel) = OsqpEigen::INFTY;
    }
}

void PredControl::UpdatePredictionModel()
{
    if (this->nc == 0)
    {
        for (int i = 0; i < N; i++)
        {
            int j = 0;

            if (i != 0)
            {
                this->Phi.block(i * this->ny, 0, this->ny, this->nxa) = this->Phi.block((i - 1) * this->ny, 0, this->ny, this->nxa) * this->A;
                this->aux_mdl = this->Phi.block((i - 1) * this->ny, 0, this->ny, this->nxa) * this->B;
            }

            while ((j < this->M) and (i + j < this->N))
            {
                this->G.block((i + j) * this->ny, j * this->nu, this->ny, this->nu) = this->aux_mdl;
                j++;
            }
        }
    }
    else
    {
        this->DefineConstraintMtxs();
        for (int i = 0; i < N; i++)
        {
            int j = 0;

            if (i != 0)
            {
                this->Phi.block(i * this->ny, 0, this->ny, this->nxa) = this->A * this->Phi.block((i - 1) * this->ny, 0, this->ny, this->nxa) ;
                this->aux_mdl = this->Phi.block((i - 1) * this->ny, 0, this->ny, this->nxa) * this->B;

                this->Phi_cons.block(i * this->nc, 0, this->nc, this->nxa) = this->Phi_cons.block((i - 1) * this->nc, 0, this->nc, this->nxa) * this->A;
                this->aux_cons = this->Phi_cons.block((i - 1) * this->nc, 0, this->nc, this->nxa) * this->B;
                this->aux_cons.block(0, 0, 2, this->nu) = Eigen::MatrixXd::Identity(2, this->nu);
            }

            while ((j < this->M) and (i + j < this->N))
            {
                this->G.block((i + j) * this->ny, j * this->nu, this->ny, this->nu) = this->aux_mdl;
                this->G_cons.block((i + j) * this->nc, j * this->nu, this->nc, this->nu) = this->aux_cons;
                j++;
            }
        }
    }
}

void PredControl::updateStatesMatrices(int p)
{
    this->A << 1, 0, -this->v_ref[this->current_iter+p] * sin(this->th_ref[this->current_iter+p]) * this->dt,
                0, 1, this->v_ref[this->current_iter+p] * cos(this->th_ref[this->current_iter+p]) * this->dt,
                0, 0,                    1;

    this->B << cos(this->th_ref[this->current_iter+p]) * this->dt,    -0.0 * this->dt * this->dt * this->v_ref[this->current_iter+p] * sin(this->th_ref[this->current_iter+p]),
            sin(this->th_ref[this->current_iter+p]) * this->dt,    0.0 * this->dt * this->dt * this->v_ref[this->current_iter+p] * sin(this->th_ref[this->current_iter+p]),
                    0                   , this->dt;


}
// void PredControl::UpdatePredictionModel()
// {
//     if (this->nc == 0)
//     {
//         // Eigen::Matrix3d aux_G = Eigen::Matrix3d::Identity();

//         for (int i = 0; i < N; i++)
//         {
//             int j = 0;

//             if (i != 0)
//             {
//                 this->updateStatesMatrices(i);
//                 this->Phi.block(i * this->ny, 0, this->ny, this->nxa) = this->A * this->Phi.block((i - 1) * this->ny, 0, this->ny, this->nxa);
// //                 this->aux_cons = this->Phi_cons.block((i - 1) * this->nc, 0, this->nc, this->nxa) * this->B;
//             }

//             while ((j < this->M) and (i + j < this->N))
//             {
//                 this->updateStatesMatrices(i+j);
//                 if (i != 0){
//                     this->aux_mdl = this->A * G.block((i + j - 1) * this->ny, j * this->nu, this->ny, this->nu);
//                 }
//                 // this->aux_mdl = aux_G * this->B;
//                 this->G.block((i + j) * this->ny, j * this->nu, this->ny, this->nu) = this->aux_mdl;
//                 j++;
//             }
//         }
//     }
//     else
//     {
//         std::cout << "Entró aquí" << std::endl;
//         this->DefineConstraintMtxs();
//         for (int i = 0; i < N; i++)
//         {
//             int j = 0;

//             if (i != 0)
//             {
//                 this->updateStatesMatrices(i);
//                 this->Phi.block(i * this->ny, 0, this->ny, this->nxa) = this->A * this->Phi.block((i - 1) * this->ny, 0, this->ny, this->nxa);
//                 // std::cout << "Entró aquí PHI" << std::endl;
//                 this->Phi_cons.block(i * this->nc, 0, this->nc, this->nxa) = this->Phi_cons.block((i - 1) * this->nc, 0, this->nc, this->nxa) * this->A;
//                 // std::cout << "Entró aquí PHI cons" << std::endl;
//             }

//             while ((j < this->M) and (i + j < this->N))
//             {
//                 this->updateStatesMatrices(i+j);
//                 if (i != 0){
//                     this->aux_mdl = this->A * G.block((i + j - 1) * this->ny, j * this->nu, this->ny, this->nu);
//                     this->aux_cons = G_cons.block((i + j - 1) * this->nc, j * this->nu, this->nc, this->nu) * this->A;
//                 }
//                 // std::cout << "Entró aquí G" << std::endl;
//                 // this->aux_mdl = aux_G * this->B;
//                 this->G.block((i + j) * this->ny, j * this->nu, this->ny, this->nu) = this->aux_mdl;
//                 this->G_cons.block((i + j) * this->nc, j * this->nu, this->nc, this->nu) = this->aux_cons;
//                 j++;
//             }
//         }
//     }
// }
void PredControl::ConfPO()
{
    // then, configure the solver

    this->solver.settings()->setVerbosity(0);

    this->solver.data()->setNumberOfVariables(this->nu * this->M);

    this->hessian_sparse = this->H.sparseView();
    this->solver.data()->clearHessianMatrix();
    this->solver.data()->setHessianMatrix(this->hessian_sparse);

    this->solver.data()->setGradient(F.transpose());

    this->solver.data()->setNumberOfConstraints(this->nc * this->N);
    this->linearMatrix = this->Ain.sparseView();
    this->solver.data()->clearLinearConstraintsMatrix();
    this->solver.data()->setLinearConstraintsMatrix(this->linearMatrix);
    this->solver.data()->setLowerBound(this->Lb);
    this->solver.data()->setUpperBound(this->Ub);

    if (this->nc != 0)
        this->constraints = 1;

    if (!this->first_conf)
    {
        if (!this->solver.initSolver())
            std::cout << "***************** PO  Inicialization Problem ***************** " << std::endl;
        else
            std::cout << "***************** PO OK ***************** " << std::endl;
    }
    else
    {
        if (!this->solver.initSolver())
        {
            std::cout << "Error: " << std::endl;
        }
    }
}

int PredControl::SolvePO()
{

    this->hessian_sparse = this->H.sparseView();
    if (!this->solver.updateHessianMatrix(this->hessian_sparse))
        return -1;

    this->solver.updateGradient(this->F.transpose());

    if (this->constraints != 0)
    {
        this->linearMatrix = this->Ain.sparseView();
        this->solver.updateLinearConstraintsMatrix(this->linearMatrix);
        this->solver.updateBounds(this->Lb, this->Ub);
    }

    if (this->solver.solveProblem() == OsqpEigen::ErrorExitFlag::NoError)
    {
        if (this->solver.getStatus() != OsqpEigen::Status::Solved)
        {
            std::cout << "Not solved - error: " << std::endl;
            std::cout << "Lc_block: " << this->Lc_block << std::endl;
            std::cout << "Lb: " << this->Lb << std::endl;
            return 0;
        }

        this->QPSolution = this->solver.getSolution();
        this->delta_qr = this->QPSolution.block(0, 0, 2, 1);
        return 1;
    }
    else
    {
        if (true)
            std::cout << "Not solved - error" << std::endl;
        return 0;
    }
}