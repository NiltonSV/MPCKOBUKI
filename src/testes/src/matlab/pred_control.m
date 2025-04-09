function [solution] = pred_control(x, A, B, C, Ref, Q, R, N, M, lb, ub)

    % #
    nx = size(A, 1);
    nu = size(B, 2);
    ny = size(C, 1);

    Gb      = zeros(ny*N, nu*M);
    Phib    = zeros(ny*N, nx);

    % ************* Q and R block *************
    Qb   = [];
    Rb   = [];

    for i=1:N
        if i-1 < M
            Rb = blkdiag(Rb, R);
        end
        if i < N
            Qb = blkdiag(Qb,power(2, i-1)*Q);
        else
            Qb = blkdiag(Qb,30*power(2, i-1)*Q);
        end
    end

    % ************* Phi and G block *************

    Phib(1:ny, 1:nx)    = C*A;
    aux_mdl = C*B;

    for i=1:N
        j = 0;
        if(i ~= 1)
            Phib((i-1) * ny + (1:ny), 1:nx) = Phib((i - 2) * ny + (1:ny), 1:nx) * A;
            aux_mdl = Phib((i - 2) * ny + (1:ny), 1:nx) * B;
        end

        while (j < M && (i + j - 1) < N)
            Gb((i + j - 1) * ny + (1:ny), j * nu + (1:nu)) = aux_mdl;
            j = j + 1;
        end
    end

    % disp(Gb)

    % Matrices and variables
    H = 2*(Gb'*Qb*Gb+Rb);
    H = (H+H')/2;
    F = 2*Gb'*Qb*Phib*(x-Ref);
    % disp(size(F))
    % F = 2*((Phib*x)-Ref)'*Qb*Gb;

    % disp("********************************")
    % disp(Qb)
    % disp(Gb)
    % disp(Rb)

    % Solver
    options = optimoptions('quadprog', ...
    'Algorithm','interior-point-convex', ...  % Método robusto
    'Display','off', ...                      % Para no imprimir texto
    'MaxIterations',100, ...
    'TolFun',1e-6, ...
    'ConstraintTolerance',1e-6);
    solution = quadprog(H, F, [], [], [], [], lb, ub, [], options);

end
