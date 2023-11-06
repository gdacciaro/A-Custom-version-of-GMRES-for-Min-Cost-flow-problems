## A Custom Version of GMRES Algorithm

<div align="center">
 <p>
    <img style="" src="./logounipi.png" alt="Logo" width="250" >  <br>
    <p>
    Master Degree (Artificial Intelligence)<br>
    Course: Computational Mathematics for Learning and Data Analysis <br> Academic Year: 2022/2023<br> <i>University of Pisa</i>, Italy 🇮🇹
    </p>
  </p>
</div>
<div align="center">
 <p align="center"><h3>👨🏻‍💻 Authors</h3>
    <a href="mailto:g.acciaro@studenti.unipi.it">Gennaro Daniele Acciaro</a> | 
    <a href="mailto:a.capurso1@studenti.unipi.it">Alessandro Capurso</a>
  </p>
    <p align="center">
    <h3><a href="./report.pdf"> 📃 Report</a></h3>
  </p>
</div>

## 🔍 Problem description
###### This project relies on solving the following problems:
(P) is a sparse linear system of the form

```math
\begin{bmatrix}D & E^T\\ E & 0\end{bmatrix}\begin{bmatrix}x \\ y\end{bmatrix} = \begin{bmatrix}b \\ c\end{bmatrix}
```
where $D \in R^{m\times m}$ is a diagonal positive definite matrix (i.e.,  $D=diag(D)>0$ ) and $E \in R^{(n-1)\times m}$ is obtained by removing the last row from the node-arc incidence matrix of a given connected directed **graph**.
These problems arise as the KKT system of the convex quadratic separable Min-Cost Flow Problem, hence you can look, e.g., [here](https://commalab.di.unipi.it/datasets/mcf) for ways to generate meaningful instances of the problem.

(A1) is GMRES, and you must solve the internal problems $\min \; || H_ny-||b||e_1||$ by updating the QR factorization of $H_n$ at each step: given the QR factorization of  $H_{n-1}$  computed at the previous step, apply one more orthogonal transformation to compute that of $H_n$.

(A2) is the same GMRES, but using the so-called *Schur complement preconditioner*

```math
P= \begin{bmatrix}D & 0\\ 0 & -S\end{bmatrix}
```

where $S$ is either $S=-ED^{-1}E^T$ or a sparse approximation of it (to obtain it, for instance, replace the smallest off-diagonal entries of $S$ with zeros). $P$ must be factorized with Incomplete Cholesky factorization.

No off-the-shelf solvers allowed.

## 🔧 Usage
In order to use the algorithm, you need to call the function `our_gmres` in the following way:

    [x, r_rel, residuals, break_flag, k] = our_gmres(D, E, S, b, starting_point, threshold, reorth_flag, debug)

    Input: D - the original diagonal vector
           E - the original E matrix
           S - the (optional) Shur complement matrix factorized with the Incomplete Cholesky factorization 
                - set to NaN if preconditioning is not required
           b - the original b vector
           starting_point - the starting point of the algorithm
           threshold - the threshold to stop the algorithm
           reorth_flag - the flag is used to decide if the reorthogonalization is needed
           debug - the (optional) flag is used to print the debug information - set to NaN if not used
    
    Output: x - the solution of the system
            r_rel - the relative residual
            residuals - the vector of the residuals
            break_flag - the flag that indicates the reason why the algorithm stopped
                0 - the algorithm converged at the threshold
                1 - the algorithm converged due to the patient
                2 - the algorithm converged due to the lucky breakdown
                -1 - the algorithm did not converge
            k - the number of iterations

## 🧪 Tests
To check the correctness and the performance of the algorithm, a suite of tests is provided. In order to run the tests, you need to call the function `run_everything` in the following way:

    test/run_everything.m

The results (.csv and plots) of the tests are stored in their corresponding folders.

To explore a specific case, execute the main.m file and utilize the prompt to select the parameters.

## 🗃️ Main Files

    📦 
     ┣ 📂 graphs       
     ┣ 📂 test
        ┣ 📂 A1_test
        ┣ 📂 A2_test   
        ┗ 📜 run_everything.m   - run all the tests     
     ┣ 📜 our_gmres.m           - our implementation of GMRES      
     ┗ 📜 main.m                
