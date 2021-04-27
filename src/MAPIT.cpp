#include "MAPIT.h"
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]

arma::mat GetLinearKernel(arma::mat X){
    double p = X.n_rows;
    return X.t() * X / p;
}

arma::mat ComputePCs(arma::mat X, int top = 10){
    arma::mat U;
    arma::vec s;
    arma::mat V;
    svd(U, s, V, X);

    arma::mat PCs = U*diagmat(s);
    return PCs.cols(0,top - 1);
}

arma::mat ComputeProjectionMatrix(int n, arma::mat b)
{
    // TODO benchmark this to see if saving the identity to memory is slow
    arma::mat identity = arma::eye<arma::mat>(n, n);
    return identity - b * arma::inv(b.t() * b) * b.t();
}


// Computes Tr(ab)
double ProductTrace(arma::mat a, arma::mat b)
{
    // Computed efficiently using Hadamard (elementwise) product
    // https://en.wikipedia.org/wiki/Trace_(linear_algebra)#Trace_of_a_product
    // https://proofwiki.org/wiki/Trace_of_Matrix_Product
    return arma::as_scalar(accu(a%b));
}

arma::mat ComputeSMatrix(std::vector<arma::mat> matrices)
{
    int num_variance_components = matrices.size();
    arma::mat S = arma::zeros(num_variance_components, num_variance_components); //Create kxk-matrix S to save

    for (int i = 0; i < num_variance_components; i++)
    {
        for (int j = 0; j < num_variance_components; j++)
        {
            if (i <= j) // create upper triangular matrix
            {
                S(i,j) = ProductTrace(matrices[i], matrices[j]);
                S(j,i) = S(i,j);
            }
        }
    }
    return S;
}

arma::vec ComputeqVector(arma::vec yc, std::vector<arma::mat> matrices)
{
    int num_variance_components = matrices.size();
    arma::vec q = arma::zeros(num_variance_components);

    for (int i = 0; i < num_variance_components; i++)
    {
        q(i) = arma::as_scalar(yc.t()*matrices[i]*yc);
    }
    return q;
}

// TODO ComputeHMatrix and ComputeVMatrix can probably be merged
arma::mat ComputeHMatrix(arma::mat Sinv, std::vector<arma::mat> matrices)
{
    arma::mat H = arma::zeros(arma::size(matrices[0]));
    for (int i = 0; i < matrices.size(); i++)
    {
        H = H + Sinv(0, i) * matrices[i];
    }
    return H;
}

arma::mat ComputeVMatrix(arma::vec delta, std::vector<arma::mat> matrices)
{
    arma::mat V = arma::zeros(arma::size(matrices[0]));
    for (int i = 0; i < matrices.size(); i++)
    {
        V = V + delta(i)*matrices[i];
    }
    return V;
}

// TODO should this be called ComputeVarianceDelta?
double ComputeVarianceSigma(arma::vec yc, arma::mat H, arma::mat V)
{
    return arma::as_scalar(2 * yc.t() * H.t() * V * H * yc);
}

double ComputeVarianceSigma(arma::vec yc, arma::mat Sinv, arma::vec delta, std::vector<arma::mat> matrices)
{
    arma::mat H = ComputeHMatrix(Sinv, matrices);
    arma::mat V = ComputeVMatrix(delta, matrices);
    return ComputeVarianceSigma(yc, H, V);
}

arma::vec RemoveFirstElement(arma::vec vector)
{
    arma::vec new_vector = arma::vec(vector.size() - 1);
    for (int i = 1; i < vector.size(); i++)
    {
        new_vector(i - 1) = vector(i);
    }
    return new_vector;
}

////////////////////////////////////////////////////////////////////////////

//Below are functions for MAPIT using two hypothesis testing strategies:
//(1) MAPIT using the Normal or Z-Test
//(2) MAPIT using the Davies Method

//Considered are the following submodels:
//(1) Standard Model ---> y = m+g+e
//(2) Standard + Covariate Model ---> y = Wa+m+g+e
//(3) Standard + Common Environment Model ---> y = m+g+c+e
//(4) Standard + Covariate + Common Environment Model ---> y = Wa+m+g+c+e

//NOTE: delta = {delta(0),delta(1),delta(2)} = {sigma^2,omega^2,tau^2} for models (1) and (2)
//NOTE: delta = {delta(0),delta(1),delta(2),delta(3)} = {sigma^2,omega^2,nu^2,tau^2} for models (3) and (4)

////////////////////////////////////////////////////////////////////////////

// Generalized MAPIT -- should handle all models

// [[Rcpp::export]]
Rcpp::List MAPITCpp(arma::mat X,
                 arma::mat Y,
                 Rcpp::Nullable<Rcpp::NumericMatrix> Z = R_NilValue,
                 Rcpp::Nullable<Rcpp::NumericMatrix> C = R_NilValue,
                 Rcpp::Nullable<Rcpp::NumericVector> variantIndices = R_NilValue,
                 std::string testMethod = "normal",
                 int cores = 1,
                 Rcpp::Nullable<Rcpp::NumericMatrix> GeneticSimilarityMatrix = R_NilValue,
                 std::string phenotypeCovariance = "identity")
{
    std::string logname = "mvMAPIT";
    auto logger = spdlog::get(logname);                         // retrieve existing one
    if (logger == nullptr) logger = spdlog::r_sink_mt(logname);     // or create new one if needed
    spdlog::set_default_logger(logger);                         // and set as default

    spdlog::stopwatch sw;                                   // instantiate a stop watch

    int i;
    const int n = X.n_cols;
    const int p = X.n_rows;
    const int d = Y.n_rows;
    int q = 0;


    logger->info("Number of samples: {}", n);
    logger->info("Number of SNPs: {}", p);
    logger->info("Number of phenotypes: {}", d);
    logger->info("Test method: {}", testMethod);
    logger->info("Phenotype covariance model: {}", phenotypeCovariance);


    if (Z.isNotNull())
    {
        // TODO benchmark this conversion, as we'll do it again below
        q = Rcpp::as<arma::mat>(Z.get()).n_rows; // convert to arma matrix
    }

    Rcpp::NumericVector sigma_est(p);
    Rcpp::NumericVector sigma_se(p);
    Rcpp::NumericVector pve(p);
    arma::mat Lambda(n * d, p);

    arma::mat GSM;
    if (GeneticSimilarityMatrix.isNull())
    {
        GSM = GetLinearKernel(X);
    }
    else
    {
        GSM = Rcpp::as<arma::mat>(GeneticSimilarityMatrix.get());
    }

    arma::vec ind;
    if (variantIndices.isNotNull()) // check if we are provided variants of interest
    {
        ind = Rcpp::as<arma::vec>(variantIndices.get());
    }

    // between phenotype variance
    arma::mat V_M(d, d); V_M.eye(); // effect of error uncorrelated
    arma::mat V_K(d, d);
    arma::mat V_G(d, d);
    if (phenotypeCovariance.compare("covariance") == 0) { // string.compare() returns '0' if equal
        spdlog::info("Covariance of effects proportional to phenotype covariance.");
        V_K = cov(Y.t());
        V_G = V_K;
    } else if (phenotypeCovariance.compare("homogeneous") == 0) {
        spdlog::info("Effect of a variant homogeneous across phenotypes.");
        V_K.ones();
        V_G.ones();
    } else {
        spdlog::info("Effect of a variant uncorrelated across phenotypes.");
        V_K.eye();
        V_G.eye();
    }

#ifdef _OPENMP
    omp_set_num_threads(cores);
#endif
#pragma omp parallel for schedule(dynamic)
    for(i=0; i<p; i++){

        if (variantIndices.isNotNull()) // check if we are provided variants of interest
        {
            // look for i+1 because R uses 1-based indexing
            // if there is no match find == ind.end()
            if (std::find(ind.begin(), ind.end(), i+1) == ind.end())
            {
                continue; // skip to next variant if you don't find a match
            }
        }

        //Compute K and G covariance matrices
        arma::mat K = (GSM * p - arma::trans(X.row(i)) * X.row(i)) / (p - 1); //Create the linear kernel
        arma::mat G = K;
        G.each_row() %= X.row(i);
        G.each_col() %= arma::trans(X.row(i));

        //Transform K and G using projection M
        arma::mat b = arma::zeros(n, q + 2);
        b.col(0) = arma::ones<arma::vec>(n);
        if (q > 0)
        {
            b.cols(1, q) = Rcpp::as<arma::mat>(Z.get()).t();
        }
        b.col(q + 1) = arma::trans(X.row(i));
        arma::mat btb_inv = arma::inv(b.t() * b);
        arma::mat M = ComputeProjectionMatrix(n, b);
        arma::mat Kc = M * K * M;
        arma::mat Gc = M * G * M;
        arma::mat Cc;
        if (C.isNotNull())
        {
            Cc = M * Rcpp::as<arma::mat>(C.get()) * M; // M*C*M, but C is converted to an arma::mat first
        }
        arma::mat Yc = Y * M;
        arma::vec yc = vectorise(Yc); // vectorise the multi-phenotype matrix

        // kronecker products
        Kc = kron(V_K, Kc);
        Gc = kron(V_G, Gc);
        M = kron(V_M, M);
        spdlog::debug("Kronecker product dimensions: {} x {}", M.n_rows, M.n_cols);


        //Compute the quantities q and S
        std::vector<arma::mat> matrices;
        if (C.isNotNull())
        {
            matrices = { Gc, Kc, Cc, M };
        }
        else
        {
            matrices = { Gc, Kc, M };
        }
        arma::vec q = ComputeqVector(yc, matrices);
        arma::mat S = ComputeSMatrix(matrices);

        //Compute delta and Sinv
        arma::mat Sinv = arma::inv(S);
        arma::vec delta = Sinv * q;

        // TODO these blocks should probably be extracted into methods
        if (testMethod == "normal")
        {
            //Compute var(delta(0))
            double V_sigma = ComputeVarianceSigma(yc, Sinv, delta, matrices);

            //Save point estimates and SE of the epistasis component
            sigma_est(i) = delta(0);
            sigma_se(i) = sqrt(V_sigma);
        }
        else if (testMethod == "davies")
        {
            // Find the eigenvalues of the projection matrix
            // TODO this could probably be cleaned up
            arma::vec evals;
            arma::vec eigval;
            arma::mat eigvec;
            if (C.isNull())
            {
                evals = arma::eig_sym((Sinv(0, 0) * Gc + Sinv(0, 1) * M) * q(1) / S(1, 1));
            }
            else
            {
                // Compute P and P^{1/2} matrix
                arma::vec delta_null = arma::inv(S) * q;
                arma::eig_sym(eigval, eigvec, delta_null(0) * Kc + delta_null(1) * Cc + delta_null(2) * M);

                evals = arma::eig_sym((eigvec.cols(find(eigval > 0)) * arma::diagmat(sqrt(eigval(find(eigval > 0)))) * arma::trans(eigvec.cols(find(eigval > 0)))) * (Sinv(0, 0) * Gc + Sinv(0, 1) * Kc + Sinv(0, 2) * Cc + Sinv(0, 3) * M) * (eigvec.cols(find(eigval > 0)) * arma::diagmat(sqrt(eigval(find(eigval > 0)))) * trans(eigvec.cols(find(eigval > 0)))));
            }
            Lambda.col(i) = evals;
        }

        //Compute the PVE
        pve(i) = delta(0) / arma::accu(delta);
    }

    logger->info("Elapsed time: {}", sw);
    //Return a Rcpp::List  of the arguments
    if (testMethod == "davies")
    {
        return Rcpp::List::create(Rcpp::Named("Est") = sigma_est, Rcpp::Named("Eigenvalues") = Lambda,Rcpp::Named("PVE") = pve);
    }
    else //Default to test method "normal"
    {
        //Compute the p-values for each estimate
        Rcpp::NumericVector sigma_pval = 2*Rcpp::pnorm(abs(sigma_est/sigma_se),0.0,1.0,0,0); //H0: sigma = 0 vs. H1: sigma != 0

        return Rcpp::List::create(Rcpp::Named("Est") = sigma_est, Rcpp::Named("SE") = sigma_se, Rcpp::Named("pvalues") = sigma_pval, Rcpp::Named("PVE") = pve);
    }
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

//Here is an alternative version MAPIT with a prespecified genomic regions to analyze

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
Rcpp::List  MAPIT_CisTrans(arma::mat X, arma::vec y, Rcpp::List  regions, bool useCis, int cores = 1){
    int i;
    const int n = X.n_cols;
    const int nsnp = X.n_rows;
    const int p = regions.size();

    //Set up the vectors to save the outputs
    Rcpp::NumericVector sigma_est(p);
    Rcpp::NumericVector pve(p);
    arma::mat Lambda(n, p);

    //Pre-compute the Linear GSM
    arma::mat GSM = GetLinearKernel(X);

#ifdef _OPENMP
    omp_set_num_threads(cores);
#endif
#pragma omp parallel for schedule(dynamic)
    for(i=0; i<p; i++){
        //Pre-compute the Linear GSM
        arma::uvec cis = regions[i];

        //Compute K and G covariance matrices
        arma::mat G;
        if (useCis)
        {
            G = (GetLinearKernel(X.rows(cis - 1)) * cis.n_elem - arma::trans(X.row(i)) * X.row(i)) / (cis.n_elem - 1); //Create the linear kernel
        }
        else
        {
            G = ((GSM * nsnp - GetLinearKernel(X.rows(cis - 1)) * cis.n_elem) - arma::trans(X.row(i)) * X.row(i)) / (nsnp - cis.n_elem - 1); //Create the linear kernel
        }
        G.each_row() %= X.row(i);
        G.each_col() %= arma::trans(X.row(i));

        //Transform K and G using projection M
        arma::mat b = arma::zeros(n, 2);
        b.col(0) = arma::ones<arma::vec>(n);
        b.col(1) = arma::trans(X.row(i));
        arma::mat btb_inv = arma::inv(b.t() * b);
        arma::mat M = ComputeProjectionMatrix(n, b);
        arma::mat Gc = M * G * M;
        arma::vec yc = M * y;

        //Compute the quantities q and S
        std::vector<arma::mat> matrices = { Gc, M };
        arma::vec q = ComputeqVector(yc, matrices);
        arma::mat S = ComputeSMatrix(matrices);

        //Compute delta and Sinv
        arma::mat Sinv = arma::inv(S);
        arma::vec delta = Sinv*q;

        //Save point estimates and SE of the epistasis component
        sigma_est(i) = delta(0);

        //Find the eigenvalues of the projection matrix
        arma::vec evals;
        arma::eig_sym(evals, (Sinv(0,0) * Gc + Sinv(0,1) * M) * q(1) / S(1,1));
        Lambda.col(i) = evals;

        //Compute the PVE
        pve(i) = delta(0) / arma::accu(delta);
    }

    //Return a Rcpp::List  of the arguments
    return Rcpp::List ::create(Rcpp::Named("Est") = sigma_est, Rcpp::Named("Eigenvalues") = Lambda,Rcpp::Named("PVE") = pve);
}