/* #############################################
 *             This file is part of
 *                    ZERO
 *
 *             Copyright (c) 2020
 *     Released under the Creative Commons
 *        Zero v1.0 Universal License
 *
 *              Find out more at
 *        https://github.com/ds4dm/ZERO
 * #############################################*/


#include "support/utils.h"
#include <armadillo>


arma::sp_mat
Utils::resizePatch(const arma::sp_mat &mat, const unsigned int nR, const unsigned int nC) {
  /**
 @brief Armadillo patch for resizing arma::sp_mat
 @details Armadillo arma::sp_mat::resize() is not robust as it initializes garbage
 values to new columns. This fixes the problem by creating new columns with
 guaranteed zero values. For arma::sp_mat
 */
  arma::sp_mat mMat(nR, nC);
  mMat.zeros();
  if (nR >= mat.n_rows && nC >= mat.n_cols) {
	 if (mat.n_rows >= 1 && mat.n_cols >= 1)
		mMat.submat(0, 0, mat.n_rows - 1, mat.n_cols - 1) = mat;
  } else {
	 if (nR <= mat.n_rows && nC <= mat.n_cols)
		mMat = mat.submat(0, 0, nR, nC);
	 else
		throw ZEROException(ZEROErrorCode::OutOfRange,
								  "Either both dimension should be smaller or larger.");
  }
  return mMat;
}

// For arma::mat
arma::mat Utils::resizePatch(const arma::mat &mat, const unsigned int nR, const unsigned int nC) {
  /**
 @brief Armadillo patch for resizing mat
 @details Armadillo mat::resize() is not robust as it initializes garbage
 values to new columns. This fixes the problem by creating new columns with
 guaranteed zero values. For arma::mat
 */
  arma::mat mMat(nR, nC);
  mMat.zeros();
  if (nR >= mat.n_rows && nC >= mat.n_cols) {
	 if (mat.n_rows >= 1 && mat.n_cols >= 1)
		mMat.submat(0, 0, mat.n_rows - 1, mat.n_cols - 1) = mat;
  } else {
	 if (nR <= mat.n_rows && nC <= mat.n_cols)
		mMat = mat.submat(0, 0, nR, nC);
	 else
		throw ZEROException(ZEROErrorCode::OutOfRange,
								  "Either both dimension should be smaller or larger.");
  }
  return mMat;
}

// For arma::vec
arma::vec Utils::resizePatch(const arma::vec &mat, const unsigned int nR) {
  /**
 @brief Armadillo patch for resizing arma::vec
 @details Armadillo arma::vec::resize() is not robust as it initializes garbage
 values to new columns. This fixes the problem by creating new columns with
 guaranteed zero values. For arma::vec
 */
  arma::vec mMat(nR);
  mMat.zeros();
  if (mat.n_rows > 0) {
	 if (nR > mat.n_rows)
		mMat.subvec(0, mat.n_rows - 1) = mat;
	 else
		mMat = mat.subvec(0, nR - 1);
  }
  return mMat;
}

void Utils::appendSave(const arma::sp_mat &matrix, ///< The arma::sp_mat to be saved
							  const std::string & out,    ///< File name of the output file
							  const std::string & header, ///< A header that might be used to
																	///< check data correctness
							  bool erase                  ///< Should the matrix be appended to the
																	///< current file or overwritten
							  )
/**
 * Utility to append an arma::sp_mat to a data file.
 */
{
  // Using C++ file operations to copy the data into the target given by @out
  unsigned int nR{0}, nC{0}, nnz{0};

  std::ofstream outfile(out, erase ? std::ios::out : std::ios::app);

  nR  = matrix.n_rows;
  nC  = matrix.n_cols;
  nnz = matrix.n_nonzero;

  outfile << header << "\n";
  outfile << nR << "\t" << nC << "\t" << nnz << "\n";
  for (auto it = matrix.begin(); it != matrix.end(); ++it)
	 outfile << it.row() << "\t" << it.col() << "\t" << (*it)
				<< "\n"; // Write the required information of arma::sp_mat
  outfile << "\n";
  outfile.close(); // and close it
}

long int Utils::appendRead(arma::sp_mat &matrix,  ///< Read and store the solution in this matrix.
									const std::string &in, ///< File to read from (could be file very many
																  ///< data is appended one below another)
									long int pos, ///< Position in the long file where reading should start
									const std::string &header ///< Any header to check data sanctity
									)
/**
 * Utility to read an arma::sp_mat from a long file.
 * @returns The end position from which the next data object can be read.
 */
{
  unsigned int nR = 0, nC = 0, nnz = 0;

  std::ifstream infile(in, std::ios::in);
  infile.seekg(pos);

  std::string headerCheckwith;
  infile >> headerCheckwith;

  if (!header.empty() && header != headerCheckwith)
	 throw ZEROException(ZEROErrorCode::InvalidData,
								"Wrong header. Expected " + header + " found " + headerCheckwith);

  infile >> nR >> nC >> nnz;
  if (nR == 0 || nC == 0)
	 matrix.set_size(nR, nC);
  else {
	 arma::umat locations(2, nnz);
	 arma::vec  values(nnz);

	 unsigned int r = 0, c = 0;
	 double       val = 0;

	 for (unsigned int i = 0; i < nnz; ++i) {
		infile >> r >> c >> val;
		locations(0, i) = r;
		locations(1, i) = c;
		values(i)       = val;
	 }
	 matrix = arma::sp_mat(locations, values, nR, nC);
  }

  pos = infile.tellg();
  infile.close();

  return pos;
}

void appendSave(const std::vector<double> &v,
					 const std::string &        out,
					 const std::string &        header,
					 bool                       erase) {
  /**
	* Utility to append an std::vector<double> to a data file.
	*/
  std::ofstream outfile(out, erase ? std::ios::out : std::ios::app);
  outfile << header << "\n" << v.size() << "\n";
  for (const double x : v)
	 outfile << x << "\n";
  outfile.close();
}

long int
appendRead(std::vector<double> &v, const std::string &in, long int pos, const std::string &header) {
  unsigned long int size = 0;
  std::ifstream     infile(in, std::ios::in);
  infile.seekg(pos);
  /**
	* Utility to read an std::vector<double> from a long file.
	* @returns The end position from which the next data object can be read.
	*/

  std::string headerCheckwith;
  infile >> headerCheckwith;

  if (!header.empty() && header != headerCheckwith)
	 throw ZEROException(ZEROErrorCode::InvalidData,
								"Wrong header. Expected " + header + " found " + headerCheckwith);

  infile >> size;

  v.resize(size);
  for (unsigned int i = 0; i < size; ++i)
	 infile >> v[i];
  pos = infile.tellg();
  infile.close();
  return pos;
}

void Utils::appendSave(const arma::vec &  matrix, ///< The arma::vec to be saved
							  const std::string &out,    ///< File name of the output file
							  const std::string &header, ///< A header that might be used to
																  ///< check data correctness
							  bool erase                 ///< Should the arma::vec be appended to the
																  ///< current file or overwritten
) {
  /**
	* Utility to append an arma::vec to a data file.
	*/
  // Using C++ file operations to copy the data into the target given by @out
  unsigned int nR{0};

  std::ofstream outfile(out, erase ? std::ios::out : std::ios::app);

  nR = matrix.n_rows;

  outfile << header << "\n";

  outfile << nR << "\n";
  for (double it : matrix)
	 outfile << it << "\n"; // Write the required information of arma::sp_mat
  outfile << "\n";
  outfile.close(); // and close it
}

long int Utils::appendRead(arma::vec &matrix,     ///< Read and store the solution in this matrix.
									const std::string &in, ///< File to read from (could be file very many
																  ///< data is appended one below another)
									long int pos, ///< Position in the long file where reading should start
									const std::string &header ///< Any header to check data sanctity
) {
  /**
	* Utility to read an arma::vec from a long file.
	* @returns The end position from which the next data object can be read.
	*/
  unsigned int  nR;
  std::string   buffers;
  std::string   checkwith;
  std::ifstream inFile(in, std::ios::in);
  inFile.seekg(pos);

  inFile >> checkwith;
  if (!header.empty() && checkwith != header)
	 throw ZEROException(ZEROErrorCode::InvalidData,
								"Wrong header. Expected " + header + " found " + checkwith);
  inFile >> nR;
  matrix.zeros(nR);
  for (unsigned int i = 0; i < nR; ++i) {
	 double val;
	 inFile >> val;
	 matrix.at(i) = val;
  }

  pos = inFile.tellg();
  inFile.close();

  return pos;
}

void Utils::appendSave(const long int     v,
							  const std::string &out,
							  const std::string &header,
							  bool               erase)
/**
 * Utility to save a long int to file
 */
{
  std::ofstream outfile(out, erase ? std::ios::out : std::ios::app);
  outfile << header << "\n";
  outfile << v << "\n";
  outfile.close();
}

long int
Utils::appendRead(long int &v, const std::string &in, long int pos, const std::string &header) {
  /**
	* Utility to read a long int from a long file.
	* @returns The end position from which the next data object can be read.
	*/
  std::ifstream infile(in, std::ios::in);
  infile.seekg(pos);

  std::string headerCheckwith;
  infile >> headerCheckwith;

  if (!header.empty() && header != headerCheckwith)
	 throw ZEROException(ZEROErrorCode::InvalidData,
								"Wrong header. Expected " + header + " found " + headerCheckwith);

  long int val;
  infile >> val;
  v = val;

  pos = infile.tellg();
  infile.close();

  return pos;
}

void Utils::appendSave(const unsigned int v,
							  const std::string &out,
							  const std::string &header,
							  bool               erase)
/**
 * Utility to save a long int to file
 */
{
  std::ofstream outfile(out, erase ? std::ios::out : std::ios::app);
  outfile << header << "\n";
  outfile << v << "\n";
  outfile.close();
}

long int
Utils::appendRead(unsigned int &v, const std::string &in, long int pos, const std::string &header) {
  std::ifstream infile(in, std::ios::in);
  infile.seekg(pos);

  std::string headerCheckwith;
  infile >> headerCheckwith;

  if (!header.empty() && header != headerCheckwith)
	 throw ZEROException(ZEROErrorCode::InvalidData,
								"Wrong header. Expected " + header + " found " + headerCheckwith);

  unsigned int val;
  infile >> val;
  v = val;

  pos = infile.tellg();
  infile.close();

  return pos;
}

void Utils::appendSave(const std::string &v, const std::string &out, bool erase)
/**
 * Utility to save a long int to file
 */
{
  std::ofstream outfile(out, erase ? std::ios::out : std::ios::app);
  outfile << v << "\n";
  outfile.close();
}

long int Utils::appendRead(std::string &v, const std::string &in, long int pos) {
  /**
	* Utility to read a std::string from a long file.
	* @returns The end position from which the next data object can be read.
	*/
  std::ifstream infile(in, std::ios::in);
  infile.seekg(pos);

  std::string val;
  infile >> val;
  v = val;

  pos = infile.tellg();
  infile.close();

  return pos;
}

bool Utils::containsConstraint(const arma::sp_mat &A,
										 const arma::vec &   b,
										 const arma::vec &   lhs,
										 const double &      rhs,
										 const double        tol) {
  if (lhs.size() != A.n_cols)
	 return false;
  for (int i = 0; i < A.n_rows; ++i) {
	 bool res = true;
	 for (int j = 0; j < A.n_cols; ++j) {
		if (std::abs(lhs.at(j) - A.at(i, j)) > tol) {
		  res = false;
		  break;
		}
	 }
	 if (res && std::abs(b.at(i) - rhs) < tol) {
		return true;
	 }
  }
  return false;
}

bool Utils::containsElement(const arma::vec &b, const double &element, const double tol) {
  for (unsigned int i = 0; i < b.size(); ++i) {
	 if (std::abs(b.at(i) - element) < tol)
		return true;
  }
  return false;
}

bool Utils::containsRow(const arma::sp_mat &A, const arma::vec &row, const double tol) {

  if (row.size() != A.n_cols)
	 return false;
  for (int i = 0; i < A.n_rows; ++i) {
	 bool res = true;
	 for (int j = 0; j < A.n_cols; ++j) {
		if (std::abs(row.at(j) - A.at(i, j)) > tol) {
		  res = false;
		  break;
		}
	 }
	 if (res)
		return true;
  }
  return false;
}
bool Utils::containsConstraint(const arma::sp_mat &A,
										 const arma::vec &   b,
										 const arma::sp_mat &lhs,
										 const double &      rhs,
										 const double        tol) {
  if (lhs.n_rows > 1)
	 return false;
  arma::vec Ai = arma::vec{lhs};
  return Utils::containsConstraint(A, b, Ai, rhs, tol);
}


bool Utils::isZero(const arma::mat &M, double tol) noexcept {
  /**
	* @brief
	* Checking if a given matrix M is a zero matrix
	*
	* @p tol Tolerance, below which a number is treated as 0
	* @warning Tolerance < 0 always returns @p false with no error.
	*
	*/
  return ((abs(M).max() <= tol));
}

bool Utils::isZero(const arma::sp_mat &M, double tol) noexcept {
  /**
	* @brief
	* Checking if a given sparse matrix M is a zero matrix
	*
	* @p tol Tolerance, below which a number is treated as 0
	*
	*/
  if (M.n_nonzero != 0)
	 return false;

  return ((abs(M).max() <= tol));
}
void Utils::sortByKey(perps &set) {
  sort(set.begin(),
		 set.end(),
		 [](std::pair<unsigned int, unsigned int> a, std::pair<unsigned int, unsigned int> b) {
			return a.first < b.first;
		 });
}

VariableBounds Utils::intersectBounds(const VariableBounds &bA, const VariableBounds &bB) {
  VariableBounds bC;
  auto           longest  = bA.size() >= bB.size() ? bA : bB;
  auto           shortest = bA.size() >= bB.size() ? bB : bA;

  // Set the size of the longest
  bC.resize(longest.size());

  for (unsigned int i = 0; i < shortest.size(); ++i) {
	 // Lower bound. The higher, the better
	 bC.at(i).first = bA.at(i).first > bB.at(i).first ? bA.at(i).first : bB.at(i).first;
	 // Upper bound. The lower, the better
	 if (bA.at(i).second < 0)
		bC.at(i).second = bB.at(i).second;
	 else
		bC.at(i).second = (bA.at(i).second < bB.at(i).second || bB.at(i).second < 0)
									 ? bA.at(i).second
									 : bB.at(i).second;
  }

  // Fill remaining element
  for (unsigned int i = shortest.size(); i < longest.size(); ++i) {
	 bC.at(i).first  = longest.at(i).first;
	 bC.at(i).second = longest.at(i).second;
  }

  return bC;
}

std::string Utils::printBounds(const VariableBounds &bounds) {
  std::stringstream r;
  for (unsigned int i = 0; i < bounds.size(); ++i) {
	 r << "var_" << std::to_string(i) << "\t\t\t[" << std::to_string(bounds.at(i).first) << ","
		<< std::to_string(bounds.at(i).second) << "]\n";
  }
  return r.str();
}


/**
 * Round the input @p value to a decimal up to @p to points.
 * @param value The input number
 * @param tol  Decimal precision tolerance
 * @return The rounded value
 */
double Utils::round_nplaces(const double &value, const double &tol)
// from
// https://www.daniweb.com/programming/software-development/code/217332/round-double-to-n-decimal-places
{
  uint32_t to     = round(log10(tol));
  uint32_t places = 1, whole = *(&value);
  for (uint32_t i = 0; i < to; i++)
	 places *= 10;
  double ret = value;

  ret -= whole; // leave decimals

  ret *= places;    // 0.1234 -> 123.4
  ret = round(ret); // 123.4 -> 123
  ret /= places;    // 123 -> .123

  ret += whole; // bring the whole value back
  return ret;
}

/**
 * Normalizes a vector according to the "equilibrium normalization". Namely, we divide for the
 * largest absolute value among the elements of the vector.
 * @p v is the input vector
 * @return The normalized vector
 */
arma::vec Utils::normalizeVec(const arma::vec &v) { return v / arma::max(arma::abs(v)); }

/**
 * Normalizes an inequality according to the "equilibrium normalization". Namely, we divide for
 * the largest absolute value among the elements of the lhs and the rhs.
 * @p rhs is the input and output RHS value
 * @p lhs is the input and output LHS
 */
void Utils::normalizeIneq(arma::vec &lhs, double &rhs, bool force) {
  arma::vec abs  = arma::abs(lhs);
  double    norm = 1;
  for (auto &elem : abs) {
	 if (elem < norm && elem != 0)
		norm = elem;
  }
  if ((abs.max() / norm > 1e3) || force) {
	 LOG_S(5) << "Utils::normalizeIneq:  normalizing inequality.";
	 assert(norm!=0);
	 rhs = rhs / norm;
	 lhs = lhs / norm;
  }
}

/**
 * @brief Given an arma::sp_mat @p A, returns a CoinPackedMatrix.
 * @p A The armadillo sparse matrix
 * @return A CoinPackedMatrix from @p A
 */
CoinPackedMatrix Utils::armaToCoinSparse(const arma::sp_mat &A) {
  CoinPackedMatrix R  = CoinPackedMatrix(true, A.n_cols, 0);
  auto         nnz   = A.n_nonzero;
  std::vector<CoinPackedVector> cols(A.n_cols);
  for (arma::sp_mat::const_iterator it = A.begin(); it != A.end(); ++it)
	 cols.at(it.col()).insert(it.row(),*it);

  for (unsigned int i=0;i<A.n_cols;++i)
	 R.appendCol(cols.at(i));


  assert(A.n_rows == R.getNumRows() && A.n_cols == R.getNumCols());
  return R;
}


/**
 * @brief Given an arma::sp_mat @p A, returns an array of CoinPackedVector(s).
 * @p A The armadillo sparse matrix
 * @return A CoinPackedVector arary from @p A
 */
std::vector<CoinPackedVector> Utils::armaToCoinPackedVector(const arma::sp_mat &A) {
  std::vector<CoinPackedVector> vectors(A.n_rows);

  for (arma::sp_mat::const_iterator it = A.begin(); it != A.end(); ++it)
	 vectors.at(it.row()).insert(it.col(), *it);

  return vectors;
}

/**
 * Create constraints for a given @p model given the matrix @p A, the RHS vector @b, the variables
 * @p x, and an additional RHS of variables @p z. The resulting constraints read:
 *  @f$Ax \quad (sense) \quad b+z@f$
 * @param A The input sparse matrix of LHS
 * @param b The input vector of RHS
 * @param x  The input variables
 * @param basename The basename of these constraints
 * @param model A pointer to the model
 * @param sense As in Gurobi.
 * @param z Additional RHS of variables
 */
void Utils::addSparseConstraints(const arma::sp_mat &A,
											const arma::vec &   b,
											GRBVar *            x,
											const std::string & basename,
											GRBModel *          model,
											int                 sense = GRB_LESS_EQUAL,
											GRBVar *            z     = nullptr) {
  std::vector<GRBLinExpr> Constraints(A.n_rows, 0);
  for (arma::sp_mat::const_iterator it = A.begin(); it != A.end(); ++it) {
	 double coeff = *it;
	 Constraints.at(it.row()).addTerms(&coeff, &x[it.col()], 1);
  }
  if (z != nullptr) {
	 for (unsigned int i = 0; i < A.n_rows; ++i) {
		model->addConstr(Constraints.at(i) - z[i], sense, b(i), basename + "_" + std::to_string(i));
	 }
  } else {
	 for (unsigned int i = 0; i < A.n_rows; ++i) {
		model->addConstr(Constraints.at(i), sense, b(i), basename + "_" + std::to_string(i));
	 }
  }
}
