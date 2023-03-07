
#include <RcppArmadillo.h>
#include <algorithm> 
#include <math.h>

using namespace Rcpp;
using namespace std;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::export]]
double fastLikelihood(const NumericVector& covar_params, const double& nu, const arma::mat& Ga, const List& Z, const double& delta, const List& cent, const double& s2, const List& T_i, const List& dd, const NumericVector& ni, const double& N) {
	// Function bess("besselK");
	// Function gam("gamma");
	// const int N = X.size();
	// arma::mat pinv = arma::inv(Psi);
	double dens = 0;
	for (int i = 0; i < N; i++) {
		// double lambda = -(nu + n[i] * p) / 2;
		arma::mat Sigma_gen = T_i[i];
		Sigma_gen.transform([&covar_params](double val) { return(std::pow(covar_params[0], val)); });
		Sigma_gen = arma::pow(Sigma_gen, covar_params[1]);
		//Sigma_gen.diag().fill(covar_params[1]);
		//CS: arma::mat Sigma_gen = arma::ones(n[i], n[i]);
		// CS: Sigma_gen = Sigma_gen * Sigma[1];
		// CS: Sigma_gen.diag().fill(Sigma[0]+Sigma[1]);
		
		// arma::mat sinv = arma::inv(Sigma_gen);

		
		//arma::mat J_n = arma::ones(n[i], 1);
		//arma::mat gen_A = J_n * A.t();
		//arma::mat resid = as<arma::mat>(Y[i]) - as<arma::mat>(X[i]) * B;
		//double Delta = arma::trace(sinv * resid * pinv * resid.t()) + nu;
		//double rho = arma::trace(sinv * gen_A * pinv *gen_A.t());
		// double k = sqrt(Delta*rho);
		// double ratio = (boost::math::cyl_bessel_k(lambda + 1, k) / boost::math::cyl_bessel_k(lambda, k));
		// double ci = log(sqrt(rho / (Delta + nu))) + (log(boost::math::cyl_bessel_k(lambda+0.001, k)) - log(boost::math::cyl_bessel_k(lambda, k))) / 0.001;
		// double bi = sqrt(rho / (Delta + nu))*ratio + (nu + n[i] * p) / (Delta + nu);
		// double ai = sqrt((Delta + nu) / rho)*ratio;
		// dens = dens + log(2) + (nu / 2)*log(nu / 2) - log(2 * 3.1415)* (n[i]*p / 2) - log(arma::det(Sigma_gen)) * (p / 2) - log(arma::det(Psi))*(n[i] / 2) - log(tgamma(nu/2))+ (lambda / 2)*(log(Delta) - log(rho)) + log(bessx) + arma::trace(sinv * resid * pinv * gen_A.t());
		// dens = dens + (nu / 2)*log(nu / 2) - (p / 2)*log(arma::det(Sigma_gen)) - (n[i] / 2)*log(arma::det(Psi)) - log(tgamma(nu / 2)) - (((nu + n[i] * p) / 2) - 1)*ci[i] + arma::trace(sinv * resid * pinv * gen_A.t()) - 0.5*bi[i] * Delta - 0.5*ai[i] * rho;
		// dens + (nu / 2)*log(nu / 2) - (p / 2)*log(arma::det(Sigma_gen)) - (n[i] / 2)*log(arma::det(Psi)) - log(tgamma(nu / 2)) - (nu / 2)*ci + 0.5*arma::trace(sinv * resid * pinv * gen_A.t()) + 0.5*arma::trace(sinv * gen_A * pinv * resid.t()) -0.5*bi*arma::trace(sinv * resid * pinv *resid.t() + nu) - 0.5*ai*arma::trace(sinv * gen_A * pinv * gen_A.t());
	//	if (i == 0) {
//			Rcout << Sigma_gen;
		//}
		arma::mat Zi = as<arma::mat>(Z[i]).t();
		arma::mat ei = as<arma::mat>(cent[i]);
		arma::mat di = as<arma::mat>(dd[i]);

		arma::mat Lambda = Zi.t() * Ga * Zi +Sigma_gen;
		arma::mat OLOinv = arma::inv(Lambda);
		double m = (ei.t() * OLOinv * ei).eval()(0, 0);
		double detL = arma::det(Lambda);
		arma::mat tmp = di.t() * OLOinv;
		double kkk = (1 - tmp * di).eval()(0, 0);
		double A = (tmp * ei).eval()(0, 0) / sqrt(kkk);

		// students_t dist(nu + ni[i]);
		// auto pt = cdf(dist, A*sqrt((nu + ni[i])/(s2*nu + m)));
		double tprob = R::pt(A*sqrt((nu + ni[i]) / (s2*nu + m)), nu + ni[i],true, true);
		dens = dens + log(2) + lgamma(0.5*(nu + ni[i])) - lgamma(0.5*nu) - 0.5*log(detL) - 0.5*ni[i] * log(3.14159265358979323846*nu*s2) - 0.5*(nu + ni[i])*log(1 + (m / s2) / nu) + tprob; //+ log(tprob);
	
	}
	dens = -dens;
	return dens;
}