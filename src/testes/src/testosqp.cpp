#include "testes/trajectory.hpp"
#include "testes/pred_control.h"
#include <iostream>

int main(){

    /*  Getting the trajectory  */

    Trajectory path;
    path.setTrajectory(1.25, 5, 0.1, 1, "infinite");
    path.saveTrajectory(path.x_ref, path.y_ref, path.th_ref, path.v_ref, path.w_ref);
    int n_total = path.x_ref.size();
    
    /*  Creating parameters  */

    // float rr = 0.3; // Robot radius 30cm
    PredControl control;
    int N;
    int M;
    Eigen::VectorXd ref(3);
    Eigen::Matrix3d Q;
    Eigen::Matrix2d R;
    Eigen::Vector3d x;
    Eigen::Vector2d delta_qr;
    delta_qr[0] = 0;
    std::vector<double> xrobot;
    std::vector<double> yrobot;
    std::vector<double> throbot;
    std::vector<double> vrobot;
    std::vector<double> wrobot;

    /*  Setting values and shapes  */
    control.x_ref = path.x_ref;
    control.y_ref  = path.y_ref;
    control.th_ref = path.th_ref;
    control.v_ref  = path.v_ref;
    control.dt     = path.dt;
    control.A.resize(3, 3);
    control.B.resize(3, 2);
    control.C.resize(3, 3); 

    N = 2;
    M = 1;

    control.C << 1, 0, 0,
                 0, 1, 0,
                 0, 0, 1;

    Q << 1.0,  0,   0,
          0,  1.0,  0,
          0,   0,  0.5;
          
    R << 1.0,   0,
          0,   1.0;

    x << 0.0, 0.0, 0; // Infinite
    // x << 0.0, 1.0, 0; // Circle
    // x << -1.0, -1.0, 0; // U

    xrobot.push_back(x[0]);
    yrobot.push_back(x[1]);
    throbot.push_back(x[2]);
    vrobot.push_back(0);
    wrobot.push_back(0);

    double mae = 0; // Mean absolute error for position

    /*  Obstacles  */
    std::vector<Obstacle> obstacles;
    // obstacles.push_back({6.0, 6.5, 0.10});
    // obstacles.push_back({0.8, 0.9, 0.05});
    // obstacles.push_back({0.0, -1.25, 0.05});
    int nobs = obstacles.size();
    std::ofstream outFile("obstacles.txt");
    outFile << "x y r\n"; // Header
    for (int i = 0; i<nobs; ++i)
    {
        outFile << obstacles[i].x << " " << obstacles[i].y << " " << obstacles[i].radius << "\n";
    }
    outFile.close();

    /*  States  */
    int nc_states = 1;   // y >= ymin
    control.nc_states = nc_states;
    
    // /*  Vel Constraints  */
    int ncon = 2;
    Eigen::MatrixXd umin(ncon, 1);
    Eigen::MatrixXd umax(ncon, 1);

    if (ncon)
    {
        umin << -1.5, -3.50; // vmin, wmin
        umax << 1.5, 3.50;
    }

    control.C_cons.resize(ncon+nobs+nc_states, 3); // 2 for velocities and 2 for obstacles, # states
    control.C_cons.setZero();
    control.C_consV.resize(ncon+nobs+nc_states, 0); // Center of mass constraints deactivated

    /*  COnfiguring the solver  */

    control.nc_vel = ncon;
    control.solver.settings()->setVerbosity(false);
    control.SetInternalVariables();
    control.SetConstants(path.dt, N, M);
    control.SetWeightMatrices(Q, R);

    for(int i = 0; i < n_total; ++i)
    {
        std::cout << "i: " << i << std::endl;
        control.current_iter = i;
    //     /*
    //         Kinematic Model from paper: https://www.ece.ufrgs.br/~fetter/sbai05_10022.pdf

    //             | 1    0    -v(k)_ref*sin(theta_ref(k))*T |      | cos(theta_ref(k))*T    0 |
    //         A = | 0    1     v(k)_ref*cos(theta_ref(k))*T |  B = | sin(theta_ref(k))*T    0 |
    //             | 0    0                    1             |      |          0             T |
    //     */
        // control.A << 1, 0, -path.v_ref[i] * sin(path.th_ref[i]) * path.dt,
        //              0, 1, path.v_ref[i] * cos(path.th_ref[i]) * path.dt,
        //              0, 0,                       1;

        control.A << 1, 0, -delta_qr[0] * sin(x[2]) * path.dt,
            0, 1, delta_qr[0] * cos(x[2]) * path.dt,
            0, 0,                       1;

        // control.B << cos(path.th_ref[i]) * path.dt,    -0.0 * path.dt * path.dt * path.v_ref[i] * sin(path.th_ref[i]),
        //              sin(path.th_ref[i]) * path.dt,    0.0 * path.dt * path.dt * path.v_ref[i] * sin(path.th_ref[i]),
        //                       0                   , path.dt;

        control.B << cos(x[2]) * path.dt,    -1e-16 * sin(x[2]),//-0.0005 * path.dt * path.dt * delta_qr[0] * sin(x[2]),
                     sin(x[2]) * path.dt,     1e-16 * sin(x[2]), //0.0005 * path.dt * path.dt * delta_qr[0] * sin(x[2]),
                              0                   , path.dt;
    
    //     // control.A = Eigen::MatrixXd::Identity(3, 3) + path.dt * control.A;
    //     // control.B = path.dt * control.B;

        for(int k = 0; k < control.N; ++k)
        {
            if(i+k < n_total)
            {
                ref << path.x_ref[i+k], path.y_ref[i+k], path.th_ref[i+k];
                control.SetReference(k, ref);
            }else{
                ref << path.x_ref[n_total-1], path.y_ref[n_total-1], path.th_ref[n_total-1];
                control.SetReference(k, ref);
            }

        }

        // ref << path.x_ref[i], path.y_ref[i], path.th_ref[i];
        // control.SetReference(1, ref);

        control.ConfPO();
        control.UpdateStates(x);
        control.SetConsBounds(umin, umax, obstacles);
        // control.UpdateObstacleConstraints(obstacles, ncon);
        control.UpdateOptimizationProblem(control.H, control.F, control.Ain, control.Lb, control.Ub);
        control.SolvePO();
        
        mae += sqrt(pow(path.x_ref[i] - x[0], 2) + pow(path.y_ref[i] - x[1], 2));

        for (int j = 0; j < 1; j++)
        {
            x[0] += (control.delta_qr(0)+path.v_ref[i]) * cos(x[2]) * (path.dt);
            x[1] += (control.delta_qr(0)+path.v_ref[i]) * sin(x[2]) * (path.dt);
            x[2] += (control.delta_qr(1)+path.w_ref[i]) * (path.dt);
            // x[0] += (control.delta_qr(0)) * cos(x[2]) * (path.dt);
            // x[1] += (control.delta_qr(0)) * sin(x[2]) * (path.dt);
            // x[2] += (control.delta_qr(1)) * (path.dt);
            // x = control.A * x + control.B * control.delta_qr;
        }

        delta_qr[0] = control.delta_qr(0)+path.v_ref[i];

        if (i<n_total-1){
            xrobot.push_back(x[0]);
            yrobot.push_back(x[1]);
            throbot.push_back(x[2]);
            vrobot.push_back(control.delta_qr(0)+path.v_ref[i]);
            wrobot.push_back(control.delta_qr(1)+path.w_ref[i]);
        }

        control.solver.clearSolver();

        // std::cout << "Phicons: " << control.Phi_cons << std::endl;
        // std::cout << "Gcons: " << control.G_cons << std::endl;
        // std::cout << "aux_cons: " << control.aux_cons << std::endl;
        // std::cout << "Lb: " << control.Lb << std::endl;
        // std::cout << "Ub: " << control.Ub << std::endl;
        // std::cout << "A: " << control.A << std::endl;
        // std::cout << "B: " << control.B << std::endl;

    }

    std::cout << "MAE: " << mae / n_total << std::endl;
    // std::cout << "Phicons: " << control.Phi_cons << std::endl;
    // std::cout << "Gcons: " << control.G_cons << std::endl;
    // std::cout << "aux_cons: " << control.aux_cons << std::endl;
    // std::cout << "Lb: " << control.Lb << std::endl;
    // std::cout << "Ub: " << control.Ub << std::endl;
    // std::cout << "A: " << control.A << std::endl;
    // std::cout << "B: " << control.B << std::endl;

    path.saveTrajectory(xrobot, yrobot, throbot, vrobot, wrobot, "results.txt");

    return 0;
}