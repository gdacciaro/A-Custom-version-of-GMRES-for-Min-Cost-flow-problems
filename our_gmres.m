%
% This is our implementation of the GMRES algorithm, exploiting the matrix
% structure:
%       J =   ⎡ D  E' ⎤
%             ⎣ E  0  ⎦
% where D is the vector containing the diagonal elements of the top-left block of J,
%       E is the bottom-left block of J.
%
% The algorithm supports the preconditioning where the preconditioning matrix is defined as:
%       P =   ⎡ D   0 ⎤
%             ⎣ 0  -S ⎦
% where S is the Shur complement matrix and it is defined as : S = E * D^{-1} * E'.
%
% The preconditioning matrix is factorized with the Incomplete Cholesky factorization 
% and the preconditioning is applied to the system as follows:
%                           R^{-T} J R^{-1} Rx = R^{-T} b
% where R is the Incomplete Cholesky factorization of P.
% It is worth noting that the preconditioning matrix is not explicitly computed but block operations are used instead.
%
%
%
% Input: D - the original diagonal vector
%        E - the original E matrix
%        S - the (optional) Shur complement matrix factorized with the Incomplete Cholesky factorization 
%              - set to NaN if preconditioning is not required
%        b - the original b vector
%        starting_point - the starting point of the algorithm
%        threshold - the threshold to stop the algorithm
%        reorth_flag - the flag is used to decide if the reorthogonalization is needed
%        debug - the (optional) flag is used to print the debug information - set to NaN if not used
%
% Output: x - the solution of the system
%         r_rel - the relative residual
%         residuals - the vector of the residuals
%         break_flag - the flag that indicates the reason why the algorithm stopped
%                     0 - the algorithm converged at the threshold
%                     1 - the algorithm converged due to the patient
%                     2 - the algorithm converged due to the lucky breakdown
%                     -1 - the algorithm did not converge
%         k - the number of iterations
%
%
% Authors: 
%   - Gennaro Daniele Acciaro
%   - Alessandro Capurso 
function [x, r_rel, residuals, break_flag, k] = our_gmres(D, E, S, b, starting_point, threshold, reorth_flag, debug)

  % Checks on the input parameters
  assert(size(D,2)==1, "D must be a vector");
  assert(size(starting_point,2)==1, "starting_point must be a vector");
  assert(islogical(reorth_flag), "reorth_flag must be a boolean value");
  if isnan(debug)
      debug = false;
  end
  
  residuals = []; % Vector of the residuals

  dim = size(D,1)+size(E,1);

  % Initialization
  patient_tol = 1e-16;
  patient = 0;
  sn = zeros(dim, 1); cs = zeros(dim, 1); % Instead of saving the whole rotation matrix, we save only the coefficients.
  e1 = zeros(dim, 1);     
  e1(1) = 1;
  last_residual = -1;

  b_norm = norm(b);
  if ~isnan(S)
    S = ichol(-S);
    D_chol = sqrt(sparse(D)); % since D is diagonal, it is equivalent and faster than ichol(sparse(diag(D)));
  else
    D_chol = NaN;
  end

  r = calculate_the_residual_optimized(D, E, S, D_chol, b, starting_point);
  r_norm = norm(r);

  Q(:,1) = r / r_norm;
  e = r_norm * e1;

  for k = 1:dim-1
    % Step 1: Lanczos algorithm
    [H(1:k+1, k), Q(:, k+1), is_breakdown_happened] = lanczos(D, E, S, D_chol, Q, k, reorth_flag);
    if is_breakdown_happened
        break_flag = 2;
        break;
    end

    % Step 2: apply the previous rotations to the newly computed column
    H(:,k) = apply_old_rotations(H(:,k), k, cs, sn);

    % Step 3: apply the current rotation to the newly computed column
    [H(1:k+1, k), e, cs(k), sn(k)] = apply_current_rotation(H(1:k+1,k), e, k);
    
    % Step 4: check the residual
    y = H(1:k, 1:k) \ e(1:k); 
    x = starting_point + Q(:, 1:k) * y;

    r = calculate_the_residual_optimized(D, E, S, D_chol, b, x);
    r_rel = norm(r)/b_norm;

    residuals(end+1) = r_rel;

    if abs(r_rel) < threshold                        % Convergence check
        break_flag = 0;
        break;
    end
    
    if abs(last_residual - r_rel) < patient_tol      % Patient check
        patient = patient+1;
    else
        patient = 0;
    end
    if patient >= 10
        break_flag = 1;
        break;
    end
    
    last_residual = r_rel;

    if debug
        fprintf("Iter: %d Res rel: %e\n", k, last_residual);
    end
  end
    
  if k == dim-1
    break_flag = -1;
  end

  y = H(1:k, 1:k) \ e(1:k);            
  y = starting_point + Q(:, 1:k) * y; 

  % Compute the final solution as x=R^{-1}*y
  if ~isnan(S) 
    m = size(D,1);
    x1 = y(1:m);
    x2 = y(m+1:end);
    x = [D_chol .\ x1 ; S \ x2];
  else
    x = y;
  end

  if debug
    fprintf("The algorithm ends with %d iterations, the break_flag is %d\n", k, break_flag);
    plot_res(residuals,threshold);
  end
end

%
% This function computes the k-th column of the Hessenberg matrix H and the k-th column of the matrix Q
%
% Input: D - the original diagonal vector
%        E - the original E matrix
%        S - the Shur complement matrix
%        D_chol - the Cholesky factorization of the diagonal matrix D
%        Q - the matrix Q of the Lanczos algorithm
%        k - the index of the column to be computed
%        reorth_flag - this flag is used to decide if the reorthogonalization is needed
%
% Output: h - the new column of the Hessenberg matrix H
%         v - the new column of the matrix Q
%         is_breakdown_happened - this flag is used to indicate if the breakdown happened
%
function [h, v, is_breakdown_happened] = lanczos(D, E, S, D_chol, Q, k, reorth_flag)
  q1 = Q(1:size(D,1), k);
  q2 = Q(size(D,1)+1:end, k);
  is_breakdown_happened = false;

  if ~isnan(S)

     B1 = D_chol .\ q1; 
     B2 = S \ q2;
    
     B11 = (D .* B1) + (E' * B2);
     B21 = E * B1;
    
     B12 = D_chol .\ B11;
     B22 = S' \ B21;
    
     v = [B12; B22];
    
  else
     B1 = (D.*q1); 
     B2 = (E'*q2);
     B12 = B1+B2;
     B3 = E*q1;
     v = [B12; B3];
  end

  for i = max(1,k-2):k    
    h(i) = v' * Q(:, i);
    v = v - h(i) * Q(:, i);
  end

  if reorth_flag             % Apply reorthogonalization twice, if necessary
    v = v - Q*(Q'*v);
    v = v - Q*(Q'*v);
  end

  h(k + 1) = norm(v);
  if abs(h(k + 1)) <= 1e-16                                % Lucky breakdown
      is_breakdown_happened = true;
  end
  v = v / h(k + 1);
end

%
% This function applies the previous rotations to the newly computed column 
%
% Input: h - the k-th column of Hn to be rotated
%        k - the index of the column to be computed
%        cs, sn - coefficients of the previous rotations
%
% Output: h - updated vector 
%
function h = apply_old_rotations(h, k, cs, sn)
  % Rotations are computed only on the last items of the vector, since the previous ones are zero

  for i = max(1,k-2):k-1
    temp   =  cs(i) * h(i) + sn(i) * h(i + 1);
    h(i+1) = -sn(i) * h(i) + cs(i) * h(i + 1);
    h(i)   = temp;
  end
end

%
% This function applies the current rotation to the newly computed column
%
% Input: h - the vector to be rotated
%        e - the canonical vector to be rotated
%        k - the index of the column to be computed
% 
% Output: h - updated vector
%         e - updated vector
%         cs_k, sn_k - coefficients of the current rotation
% 
function [h, e, cs_k, sn_k] = apply_current_rotation(h,e,k)
  [cs_k, sn_k] = get_givens_rotation(h,k,k+1);

  tmp       = cs_k * h(k)  + sn_k * h(k + 1);
  h(k + 1)  = -sn_k * h(k) + cs_k * h(k + 1);
  h(k) = tmp;

  temp     =  cs_k * e(k) + sn_k * e(k + 1);
  e(k + 1) = -sn_k * e(k) + cs_k * e(k + 1);
  e(k)     = temp;
end

% 
% This function computes the coefficient s of the Givens rotation matrix
% that zeros the j-th component of the vector v.
%
% Input: v    - vector of size 2
%        i,j  - indices of the vector
% Output: c,s - coefficients of the Givens matrix
% 
function [c,s] = get_givens_rotation(v, i, j) % Coefficients to compute the Givens matrix
    a = v(i);
    b = v(j);

    denominator = sqrt(a^2 + b^2);
    c = a/denominator;
    s = b/denominator;
end

%
% By using the matrix structure, the calculation of the residual can be optimized.
%
% Input: D - the original diagonal vector
%        E - the original E matrix
%        S - the Shur complement matrix
%        D_chol - the Cholesky factorization of the diagonal matrix D
%        b - the original b vector
%        input_vector - the vector on which the residual is computed
%
% Output: r - the residual
%
function r = calculate_the_residual_optimized(D, E, S, D_chol, b, input_vector)
  
  m = size(D,1);

  part_1 = input_vector(1:size(D,1));
  part_2 = input_vector(size(D,1)+1:end);

  if ~isnan(S) 
     
     B1 = D_chol .\ part_1;
     B2 = S \ part_2; 
    
     B11 = (D .* B1) + (E' * B2);
     B21 = E * B1;
    
     B12 = D_chol .\ B11;
     B22 = S' \ B21;

     block1 = D_chol .\ b(1:m);
     block2 = S' \ b(m+1:end);

     r = ([block1;block2] - [B12;B22]);
  else
     r = (b - [(D.*part_1)+(E'*part_2); E*part_1]);
  end
end


%
% This function generate graphs of the residuals in a logarithmic scale
%
%
%
% Input: residuals - the vector containing all the residual values
%        threshold - the threshold used as stopping criteria
%
% Output: the plot representing the evolution of the residual through time
%
function plot_res(residuals,threshold)
    figure;
    p = semilogy(residuals, 'LineWidth',2);
    p.Color=["#0072BD"];
    xlabel('iteration');
    ylabel('residual');
    yline(threshold,'--','Threshold','LineWidth',3);
end 