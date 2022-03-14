/* #############################################
 *             This file is part of
 *                    ZERO
 *
 *             Copyright (c) 2020
 *     Released under the Creative Commons
 *         CC BY-NC-SA 4.0 License
 *
 *              Find out more at
 *        https://github.com/ds4dm/ZERO
 * #############################################*/


#include "support/utils.h"
#include <armadillo>


/**
 * @brief  Armadillo arma::sp_mat::resize() is not robust as it initializes garbage
 values to new columns. This fixes the problem by creating new columns with
 guaranteed zero values. For arma::sp_mat
 * @param mat Input matrix
 * @param nR Number of rows for the output
 * @param nC Number of columns for the output
 * @return The resized @p mat
 */
arma::sp_mat
Utils::resizePatch(const arma::sp_mat &mat, const unsigned int nR, const unsigned int nC) {
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


/**
 * @brief  Armadillo arma::sp_mat::resize() is not robust as it initializes garbage
 values to new columns. This fixes the problem by creating new columns with
 guaranteed zero values. For arma::sp_mat
 * @param mat Input matrix
 * @param nR Number of rows for the output
 * @param nC Number of columns for the output
 * @return The resized @p mat
 */
arma::mat Utils::resizePatch(const arma::mat &mat, const unsigned int nR, const unsigned int nC) {
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

/**
 * @brief  Armadillo arma::mat::resize() is not robust as it initializes garbage
 values to new columns. This fixes the problem by creating new columns with
 guaranteed zero values. For arma::sp_mat
 * @param mat Input vector (row, or column)
 * @param nR Number of rows for the output
 * @return The resized @p mat
 */
arma::vec Utils::resizePatch(const arma::vec &mat, const unsigned int nR) {
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

/**
 * @brief Utility to append an arma::sp_mat to a data file.
 * @param matrix The arma::sp_mat to be saved
 * @param out File name of the output file
 * @param header A header that might be used to check the filetype
 * @param erase Should the matrix be appended to the file or not?
 */
void Utils::appendSave(const arma::sp_mat &matrix,
							  const std::string  &out,
							  const std::string  &header,
							  bool                erase) {
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

/**
 * @brief Utility to read an arma::sp_mat from a long file.
 * @param matrix  Read and store the solution in this matrix.
 * @param in File to read from (could be file very many matrices are appended one below another)
 * @param pos Position in the long file where reading should start
 * @param header Any header to check data sanctity
 * @return The end position from which the next data object can be read.
 */
long int Utils::appendRead(arma::sp_mat      &matrix,
									const std::string &in,
									long int           pos,
									const std::string &header) {
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

/**
 * @brief Utility to append a standard vector to a data file.
 * @param v The vector
 * @param out The output file
 * @param header An optional header
 * @param erase Should the file be over written?
 */
void appendSave(const std::vector<double> &v,
					 const std::string         &out,
					 const std::string         &header,
					 bool                       erase) {
  std::ofstream outfile(out, erase ? std::ios::out : std::ios::app);
  outfile << header << "\n" << v.size() << "\n";
  for (const double x : v)
	 outfile << x << "\n";
  outfile.close();
}

/**
 * @brief Utility to read an std::vector<double> from a long file.
 * @param v The output vector
 * @param in The input file
 * @param pos The position of @p v in the file
 * @param header An optional header check
 * @return The end position from which the next data object can be read.
 */
long int
appendRead(std::vector<double> &v, const std::string &in, long int pos, const std::string &header) {
  unsigned long int size = 0;
  std::ifstream     infile(in, std::ios::in);
  infile.seekg(pos);

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


/**
 * @brief Utility to append an arma::vec to a data file.
 * @param matrix The arma::vec to be saved
 * @param out The output file
 * @param header An optional header
 * @param erase Erase the file?
 */

void Utils::appendSave(const arma::vec   &matrix,
							  const std::string &out,
							  const std::string &header,
							  bool               erase) {
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

/**
 * @brief  Utility to read an arma::vec from a long file.
 * @param matrix The output matrix
 * @param in The input file
 * @param pos The position of @p matrix in the file
 * @param header An optional header
 * @return The end position from which the next data object can be read.
 */
long int Utils::appendRead(arma::vec         &matrix,
									const std::string &in,
									long int           pos,
									const std::string &header) {
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

/**
 * @brief Utility to save a long int to file
 * @param v The int
 * @param out The output file
 * @param header An optional header
 * @param erase Overwrite the file?
 */
void Utils::appendSave(const long int     v,
							  const std::string &out,
							  const std::string &header,
							  bool               erase)

{
  std::ofstream outfile(out, erase ? std::ios::out : std::ios::app);
  outfile << header << "\n";
  outfile << v << "\n";
  outfile.close();
}



/**
 * @brief Utility to read a long int from a file
 * @param v The output number
 * @param in The input file
 * @param pos The position of @p v in the file
 * @param header An optional header check
 * @return The end position from which the next data object can be read.
 */
long int
Utils::appendRead(long int &v, const std::string &in, long int pos, const std::string &header) {
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

/**
 * @brief Utility to save an unsigned int to file
 * @param v The long int to be saved
 * @param out The output file
 * @param header An optional header
 * @param erase Should the file be erased?
 */
void Utils::appendSave(const unsigned int v,
							  const std::string &out,
							  const std::string &header,
							  bool               erase) {
  std::ofstream outfile(out, erase ? std::ios::out : std::ios::app);
  outfile << header << "\n";
  outfile << v << "\n";
  outfile.close();
}

/**
 * @brief An utility to read an unsigned int from a file
 * @param v The output number
 * @param in The input file
 * @param pos The position of the data in the file
 * @param header An optional header
 * @return The end position from which the next data object can be read.
 */
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

/**
 * @brief Utility to write a string to a file
 * @param v The string to be saved
 * @param out The output file
 * @param erase Should the file be erased?
 */
void Utils::appendSave(const std::string &v, const std::string &out, bool erase) {
  std::ofstream outfile(out, erase ? std::ios::out : std::ios::app);
  outfile << v << "\n";
  outfile.close();
}
/**
 * @brief Utility to read a std::string from a long file.
 * @param v The output string
 * @param in The input file
 * @param pos The position of @p v
 * @return The end position from which the next data object can be read.
 */
long int Utils::appendRead(std::string &v, const std::string &in, long int pos) {
  std::ifstream infile(in, std::ios::in);
  infile.seekg(pos);

  std::string val;
  infile >> val;
  v = val;

  pos = infile.tellg();
  infile.close();

  return pos;
}

/**
 * @brief Given a matrix @p A and a vector of LHSs @p b, checks if both @p lhs is in @p A, and @p
 * rhs is in @p b up to a tolerance of @p tol.
 * @param A The matrix containing the constraints
 * @param b The vector containing the RHSs
 * @param lhs The input cut RHS
 * @param rhs The input cut LHS
 * @param tol A numerical tolerance
 * @return True if the constraint was found
 */
bool Utils::containsConstraint(const arma::sp_mat &A,
										 const arma::vec    &b,
										 const arma::vec    &lhs,
										 const double       &rhs,
										 const double        tol) {
  if (lhs.size() != A.n_cols)
	 return false;
  for (int i = 0; i < A.n_rows; ++i) {
	 bool res = true;
	 // The RHS is good
	 if (Utils::isEqual(b.at(i), rhs, tol)) {
		// Okay. Let's check the cut
		bool res = true;
		for (int j = 0; j < A.n_cols; ++j) {
		  if (!Utils::isEqual(lhs.at(j), A.at(i, j), tol)) {
			 // Not equal
			 res = false;
			 break;
		  }
		}
		if (res)
		  return true;
	 }
  }
  return false;
}

/**
 * @brief Check  if the vector @p b contains @p element
 * @param b The vector containing the elements to be checked
 * @param element The input element
 * @param tol A numerical tolerance
 * @return True if the element was found
 */
bool Utils::containsElement(const arma::vec &b, const double &element, const double tol) {
  for (unsigned int i = 0; i < b.size(); ++i) {
	 if (Utils::isEqual(b.at(i), element, tol))
		return true;
  }
  return false;
}


/**
 * @brief Check if @p A contains a row @p row
 * @param A The input matrix
 * @param row The input row
 * @param tol A numerical tolerance
 * @return True if the row was found
 */
bool Utils::containsRow(const arma::sp_mat &A, const arma::vec &row, const double tol) {

  if (row.size() != A.n_cols)
	 return false;
  for (int i = 0; i < A.n_rows; ++i) {
	 bool res = true;
	 for (int j = 0; j < A.n_cols; ++j) {
		if (!Utils::isEqual(row.at(j), A.at(i, j), tol)) {
		  res = false;
		  break;
		}
	 }
	 if (res)
		return true;
  }
  return false;
}
/**
 * @brief Given a matrix @p A and a vector of LHSs @p b, checks if both @p lhs is in @p A, and @p
 * rhs is in @p b up to a tolerance of @p tol.
 * @param A The matrix containing the constraints
 * @param b The vector containing the RHSs
 * @param lhs The input cut RHS (in a matrix form)
 * @param rhs The input cut LHS
 * @param tol A numerical tolerance
 * @return True if the constraint was found
 */
bool Utils::containsConstraint(const arma::sp_mat &A,
										 const arma::vec    &b,
										 const arma::sp_mat &lhs,
										 const double       &rhs,
										 const double        tol) {
  if (lhs.n_rows > 1)
	 return false;
  arma::vec Ai = arma::vec{lhs};
  return Utils::containsConstraint(A, b, Ai, rhs, tol);
}

/**
 * @brief Checks if a matrix is a zero matrix
 * @param M The input matrix
 * @param tol A numerical tolerance
 * @return True if all elements are zero
 */
bool Utils::isZero(const arma::mat &M, double tol) noexcept {

  return (Utils::isEqual(abs(M).max(), 0, tol));
}

/**
 * @brief Checks if a matrix is a zero matrix. Optimized for sparse matrices
 * @param M The input matrix
 * @param tol A numerical tolerance
 * @return True if all elements are zero
 */
bool Utils::isZero(const arma::sp_mat &M, double tol) noexcept {

  if (M.n_nonzero != 0)
	 return false;

  return (Utils::isEqual(abs(M).max(), 0, tol));
}

/**
 * @brief Sort the perps by the first element (key)
 * @param set The input perps
 */
void Utils::sortByKey(perps &set) {
  sort(set.begin(),
		 set.end(),
		 [](std::pair<unsigned int, unsigned int> a, std::pair<unsigned int, unsigned int> b) {
			return a.first < b.first;
		 });
}

/**
 * @brief Given two VariableBounds, it returns the strictest intersection of them
 * @param bA First input bound vector
 * @param bB Second input bound vector
 * @return The intersection
 */
VariableBounds Utils::intersectBounds(const VariableBounds &bA, const VariableBounds &bB) {
  auto longest  = bA.size() >= bB.size() ? bA : bB;
  auto shortest = bA.size() >= bB.size() ? bB : bA;

  // Set the size of the longest
  VariableBounds bC(longest.size());

  for (unsigned int i = 0; i < shortest.size(); ++i) {
	 // Lower bound. The higher, the better
	 if (bA.at(i).first >= 0 || bB.at(i).first >= 0)
		bC.at(i).first = bA.at(i).first > bB.at(i).first ? bA.at(i).first : bB.at(i).first;
	 else
		bC.at(i).first = 0;
	 // Upper bound. The lower, the better
	 if (bA.at(i).second < 0 && bB.at(i).second)
		bC.at(i).second = -1;
	 else {
		if (bA.at(i).second >= 0)
		  bC.at(i).second = bA.at(i).second;
		else
		  bC.at(i).second = bB.at(i).second;
	 }
  }

  // Fill remaining element
  for (unsigned int i = shortest.size(); i < longest.size(); ++i) {
	 bC.at(i).first  = longest.at(i).first;
	 bC.at(i).second = longest.at(i).second;
  }

  return bC;
}

/**
 * @brief Given @p bounds, provides a string with the readable output
 * @param bounds The input bounds
 * @return A readable output for the bounds
 */
std::string Utils::printBounds(const VariableBounds &bounds) {
  std::stringstream r;
  for (unsigned int i = 0; i < bounds.size(); ++i) {
	 r << "var_" << std::to_string(i) << "\t\t\t[" << std::to_string(bounds.at(i).first) << ","
		<< std::to_string(bounds.at(i).second) << "]\n";
  }
  return r.str();
}


/**
 * Round the input @p value to a decimal up to @p numDecimals.
 * @param value The input number
 * @param numDecimals  Decimal precision tolerance
 * @return The rounded value
 */
double Utils::round_nplaces(const double &value, const int &numDecimals) {
  return roundf(value * pow(10, numDecimals)) / pow(10, numDecimals);
}

/**
 * Normalizes a vector according to the "equilibrium normalization". Namely, we divide for the
 * largest absolute value among the elements of the vector.
 * @p v is the input vector
 * @return The normalized vector
 */
arma::vec Utils::normalizeVec(const arma::vec &v) {
  double norm = arma::max(arma::abs(v));
  if (Utils::isEqual(norm, 0, 1e-6))
	 norm = 1;

  return v / arma::max(arma::abs(v));
}

/**
 * Normalizes an inequality according to the "equilibrium normalization". Namely, we divide for
 * the largest absolute value among the elements of the lhs and the rhs. If the ration between the
 * largest non-zero and the smallest non-zero (in abs) is greater than 1e2, we normalize with the
 * former.
 * @p rhs is the input and output RHS value
 * @p lhs is the input and output LHS
 */
void Utils::normalizeIneq(arma::vec &lhs, double &rhs, bool force) {
  // lhs.print("Input with RHS="+std::to_string(rhs));
  arma::vec abs  = arma::abs(lhs);
  double    norm = abs.max();
  for (auto &elem : abs) {
	 if (std::abs(elem) < norm && elem != 0)
		norm = std::abs(elem);
  }
  double ratio = abs.max() / norm;
  if ((ratio > 1e2) || force) {
	 LOG_S(5) << "Utils::normalizeIneq:  normalizing inequality.";
	 // Force this for very bad inequalities with too much of range...
	 if (ratio > 1e2)
		norm = abs.max();

	 assert(norm != 0);
	 rhs = rhs / norm;
	 lhs = lhs / norm;
  }
  // lhs.print("Normalized with RHS="+std::to_string(rhs));
}

/**
 * @brief Given an arma::sp_mat @p A, returns a CoinPackedMatrix.
 * @p A The armadillo sparse matrix
 * @return A CoinPackedMatrix from @p A
 */
CoinPackedMatrix Utils::armaToCoinSparse(const arma::sp_mat &A) {
  CoinPackedMatrix              R   = CoinPackedMatrix(true, A.n_cols, 0);
  auto                          nnz = A.n_nonzero;
  std::vector<CoinPackedVector> cols(A.n_cols);
  for (arma::sp_mat::const_iterator it = A.begin(); it != A.end(); ++it)
	 cols.at(it.col()).insert(it.row(), *it);

  for (unsigned int i = 0; i < A.n_cols; ++i)
	 R.appendCol(cols.at(i));


  assert(A.n_rows == R.getNumRows() && A.n_cols == R.getNumCols());
  return R;
}

/**
 * @brief Given an arma::sp_mat @p A, returns the same matrix where the entries within a
 * relative and absolute tolerance @p tol @p percent are set to zero
 * @p A The armadillo sparse matrix
 * @p tol The absolute tolerance
 * @p percent The relative tolerance
 * @return The new refactored matrix
 */
arma::sp_mat Utils::clearMatrix(const arma::sp_mat &A, double tol, double percent) {
  arma::sp_mat ACopy = A;
  for (arma::sp_mat::const_iterator it = A.begin(); it != A.end(); ++it) {
	 if (Utils::isEqual(*it, 0, tol, percent))
		ACopy.at(it.row(), it.col()) = 0;
  }
  return ACopy;
}

/**
 * @brief Given an arma::vec @p A, returns the same vector where the entries within a
 * relative and absolute tolerance @p tol @p percent are set to zero
 * @p b The armadillo vector
 * @p tol The absolute tolerance
 * @p percent The relative tolerance
 * @return The new refactored matrix
 */
arma::vec Utils::clearVector(const arma::vec &b, double tol, double percent) {
  arma::vec bCopy = b;
  for (unsigned int i = 0; i < bCopy.size(); ++i) {
	 if (Utils::isEqual(b.at(i), 0, tol, percent))
		bCopy.at(i) = 0;
  }
  return bCopy;
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
 * Create constraints for a given @p model given the matrix @p A, the RHS vector @p b, the variables
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
											const arma::vec    &b,
											GRBVar             *x,
											const std::string  &basename,
											GRBModel           *model,
											int                 sense = GRB_LESS_EQUAL,
											GRBVar             *z     = nullptr) {
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


/**
 * @brief Convert a binary vector (0-1 entries) into an unique int
 * @param x The input vector
 * @return The int encoding
 */
int Utils::vecToBin(const arma::vec &x) {
  int output = 0;
  int power  = 1;
  int len    = x.size();

  for (int i = 0; i < len; i++) {
	 output += x.at(len - 1 - i) * power;
	 power *= 2;
  }

  return output;
}
/**
 * @brief Given a number @p num and a bound @p decimalBound, counts the number of non-zero decimals
 * up to @p decimalBound
 * @param num The input number
 * @param decimalBound The maximal bound on decimals to be counted
 * @return The number of non-zero decimals
 */
int Utils::nonzeroDecimals(const double num, const int decimalBound) {

  double integral = 0;
  // Take the fractional
  modf(num, &integral);
  double fractional = num - integral;
  int    count      = 0;

  for (unsigned int i = 0; i < decimalBound; ++i) {
	 fractional *= 10;
	 if (static_cast<long>(fractional) % 10 != 0)
		++count;
  }
  return count;
}
/**
 * @brief Checks if two numbers are equal in absolute tolerance
 * @param a The first number
 * @param b The second number
 * @param tol An absolute tolerance
 * @return True if the two numbers are equal up to the given tolerance
 */
bool Utils::isEqualAbs(const double a, const double b, const double tol) {
  float diff = fabs(a - b);
  if (diff <= tol)
	 return true;
  else
	 return false;
}
/**
 * @brief Checks if two numbers are equal in relative tolerance
 * @param a The first number
 * @param b The second number
 * @param percent Percent of similarity (e.g., 1 and 0.99 are similar at least at 0.99)
 * @return True if the two numbers are equal up to the given tolerance
 */
bool Utils::isEqualRel(const double a, const double b, const double percent) {
  float  diff    = fabs(a - b);
  double A       = fabs(a);
  double B       = fabs(b);
  float  largest = (B > A) ? B : A;

  if (diff <= largest * (1 - percent))
	 return true;
  else
	 return false;
}
/**
 * @brief Checks if two numbers are equal in either absolute or relative tolerance
 * @param a The first number
 * @param b The second number
 * @param tol An absolute tolerance
 * @param percent Percent of similarity (e.g., 1 and 0.99 are similar at least at 0.99)
 * @return
 */
bool Utils::isEqual(const double a, const double b, const double tol, const double percent) {
  if (Utils::isEqualAbs(a, b, tol)) {
	 return true;
  } else if (Utils::isEqualRel(a, b, percent)) {
	 return true;
  }
  return false;
}

/**
 * @brief Gets the sign
 * @param val Input data
 * @return Either -1 or 1
 */
int Utils::getSign(double val) { return val >= 0 ? 1 : -1; }